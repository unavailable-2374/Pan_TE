package RECON::MaskedTrack;

use strict;
use warnings;
use Exporter 'import';
use File::Path qw(make_path);
use Cwd qw(getcwd abs_path);
use RECON::Logger;
use RECON::Utils;
use RECON::Core;

our @EXPORT = qw(
    run_masked_track_independent run_masked_track run_full_genome_track
    extract_unmasked_regions process_sequence_for_unmasked
);

sub run_masked_track_independent {
    my ($genome_file, $output_dir, $threads, $cpu_threads, $genome_size) = @_;
    
    my $current_dir = getcwd();
    chdir $output_dir or die "Cannot change to $output_dir: $!\n";
    
    log_message("INFO", "Starting masked track processing", 
                "genome_size=" . format_size($genome_size));
    
    # Check if the entire process is already completed
    if (-f "masked_track_completed.ok") {
        log_message("INFO", "Masked track already completed", "checkpoint found: masked_track_completed.ok");
        chdir $current_dir;
        return 1;
    }
    
    # Step 1: Extract unmasked regions (with checkpoint)
    my $unmasked_file = "unmasked_regions.fa";
    if (!-f "unmasked_extraction.ok") {
        log_message("INFO", "Extracting unmasked regions", "min_length=50bp");
        
        extract_unmasked_regions($genome_file, $unmasked_file);
        
        if (!-s $unmasked_file) {
            log_message("ERROR", "Failed to extract unmasked regions", "file=$unmasked_file");
            chdir $current_dir;
            return 0;
        }
        
        # Create unmasked extraction checkpoint
        open(my $fh, '>', "unmasked_extraction.ok") or die "Cannot create unmasked extraction checkpoint: $!\n";
        print $fh "Unmasked extraction completed at " . localtime() . "\n";
        my $seq_count = `grep -c '^>' $unmasked_file`;
        chomp $seq_count;
        print $fh "Extracted sequences: $seq_count\n";
        print $fh "File size: " . (-s $unmasked_file) . " bytes\n";
        close $fh;
        log_message("INFO", "Unmasked extraction checkpoint created", "sequences=$seq_count, file=unmasked_extraction.ok");
    } else {
        log_message("INFO", "Unmasked extraction already completed", "checkpoint found: unmasked_extraction.ok");
    }
    
    my $working_genome = $unmasked_file;
    
    # Get sequence count for logging
    my $seq_count = `grep -c '^>' $working_genome`;
    chomp $seq_count;
    log_message("INFO", "Using extracted unmasked regions", "sequences=$seq_count, format=gi|N");
    
    # Step 2: Create sequence name list (with checkpoint)
    if (!-f "seq_naming.ok") {
        log_message("INFO", "Creating sequence name list", "format=gi|N");
        create_seq_name_list($working_genome);
        
        # Create sequence naming checkpoint
        open(my $fh, '>', "seq_naming.ok") or die "Cannot create sequence naming checkpoint: $!\n";
        print $fh "Sequence naming completed at " . localtime() . "\n";
        my $name_count = -s "seq.names" ? `head -1 seq.names` : 0;
        chomp $name_count;
        print $fh "Sequence names created: $name_count\n";
        close $fh;
        log_message("INFO", "Sequence naming checkpoint created", "names=$name_count, file=seq_naming.ok");
    } else {
        log_message("INFO", "Sequence naming already completed", "checkpoint found: seq_naming.ok");
    }
    
    # Step 3: Self-alignment using RMBlastN (with checkpoint)
    my $blast_output = "self_alignment.blast";
    if (!-f "rmblastn.ok") {
        log_message("INFO", "Running self-alignment", "genome=$working_genome");
        
        my $hit_count = run_rmblastn_self_alignment($working_genome, $blast_output, $threads);
        
        if ($hit_count == 0) {
            log_message("WARN", "No BLAST hits found", "skipping_RECON");
            chdir $current_dir;
            return 0;
        }
        
        # Create RMBlastN checkpoint
        open(my $fh, '>', "rmblastn.ok") or die "Cannot create RMBlastN checkpoint: $!\n";
        print $fh "RMBlastN completed at " . localtime() . "\n";
        print $fh "BLAST hits found: $hit_count\n";
        print $fh "Output file size: " . (-s $blast_output || 0) . " bytes\n";
        close $fh;
        log_message("INFO", "RMBlastN checkpoint created", "hits=$hit_count, file=rmblastn.ok");
    } else {
        log_message("INFO", "RMBlastN already completed", "checkpoint found: rmblastn.ok");
    }
    
    # Step 4: Collect MSPs (with checkpoint)
    if (!-f "msp_collection.ok") {
        log_message("INFO", "Collecting MSPs", "input=$blast_output");
        
        # Handle the case where blast output might already be in MSP format
        if (-s $blast_output) {
            # Check if file is already in MSP format by examining first line
            my $first_line = `head -1 $blast_output`;
            chomp $first_line;
            
            if ($first_line =~ /^\d{6}\s+\d+\s+\d{5}\s+\d{5}/) {
                # File is already in MSP format, just copy it
                log_message("INFO", "BLAST output already in MSP format", "copying to msp.out");
                run_cmd("cp $blast_output msp.out");
            } else {
                # File is in BLAST format, need MSPCollect.pl
                log_message("INFO", "Converting BLAST to MSP format", "using MSPCollect.pl");
                run_cmd("MSPCollect.pl $blast_output > msp.out");
            }
        }
        
        if (!-s "msp.out") {
            log_message("WARN", "No MSPs generated", "skipping_RECON");
            chdir $current_dir;
            return 0;
        }
        
        # Create MSP collection checkpoint
        open(my $fh, '>', "msp_collection.ok") or die "Cannot create MSP collection checkpoint: $!\n";
        print $fh "MSP collection completed at " . localtime() . "\n";
        my $msp_lines = `wc -l < msp.out`;
        chomp $msp_lines;
        print $fh "MSP lines: $msp_lines\n";
        print $fh "File size: " . (-s "msp.out") . " bytes\n";
        close $fh;
        log_message("INFO", "MSP collection checkpoint created", "lines=$msp_lines, file=msp_collection.ok");
    } else {
        log_message("INFO", "MSP collection already completed", "checkpoint found: msp_collection.ok");
    }
    
    # Step 5: Run RECON pipeline (with checkpoint)
    if (!-f "recon.ok") {
        my $k_param = determine_k_parameter("msp.out");
        log_message("INFO", "Starting RECON pipeline", "k_parameter=$k_param");
        
        eval {
            run_recon_pipeline($k_param);
            
            # Create RECON checkpoint
            open(my $fh, '>', "recon.ok") or die "Cannot create RECON checkpoint: $!\n";
            print $fh "RECON pipeline completed at " . localtime() . "\n";
            print $fh "K parameter used: $k_param\n";
            if (-d "summary") {
                my $families = -s "summary/families" ? `wc -l < summary/families` : 0;
                chomp $families;
                print $fh "Families found: $families\n";
            }
            close $fh;
            log_message("INFO", "RECON checkpoint created", "k_param=$k_param, file=recon.ok");
        };
        
        if ($@) {
            log_message("ERROR", "RECON pipeline failed", "error=$@");
            chdir $current_dir;
            die "Masked track RECON failure. Cannot proceed.\n";
        }
        
        log_message("INFO", "RECON pipeline completed successfully");
    } else {
        log_message("INFO", "RECON pipeline already completed", "checkpoint found: recon.ok");
    }
    
    # Step 6: Build consensi if families found (with checkpoint)
    if (-d "summary" && -s "summary/families") {
        if (!-f "consensus.ok") {
            log_message("INFO", "Building consensus sequences", "families_found=true");
            
            # Use the working genome for consensus building
            my $local_genome = abs_path($working_genome);
            # Use more threads for build_for_RECON since it's the final step
            my $build_threads = $threads;  # Use all available threads
            run_cmd("build_for_RECON ./ $local_genome $build_threads");
            
            if (-s "consensi.fa") {
                my $count = `grep -c '^>' consensi.fa`;
                chomp $count;
                my $size = -s "consensi.fa";
                
                # Create consensus checkpoint
                open(my $fh, '>', "consensus.ok") or die "Cannot create consensus checkpoint: $!\n";
                print $fh "Consensus building completed at " . localtime() . "\n";
                print $fh "Consensus sequences: $count\n";
                print $fh "File size: $size bytes\n";
                close $fh;
                
                log_message("INFO", "Consensus checkpoint created", "sequences=$count, file=consensus.ok");
                log_message("INFO", "Consensus sequences built", 
                           "count=$count, size=" . format_size($size));
                
                # Now that consensi.fa is successfully generated, cleanup intermediate files
                cleanup_intermediate_files();
            } else {
                log_message("WARN", "No consensus sequences generated despite families found");
            }
        } else {
            log_message("INFO", "Consensus building already completed", "checkpoint found: consensus.ok");
        }
    } else {
        log_message("WARN", "No families found", "skipping_consensus_building");
    }
    
    # Create completion marker
    open(my $fh, '>', "masked_track_completed.ok") or die "Cannot create completion marker: $!\n";
    print $fh "Masked track completed at " . localtime() . "\n";
    
    # Record which checkpoints were successful
    my @checkpoint_files = qw(
        unmasked_extraction.ok seq_naming.ok rmblastn.ok 
        msp_collection.ok recon.ok consensus.ok
    );
    my @completed_checkpoints;
    
    for my $checkpoint (@checkpoint_files) {
        if (-f $checkpoint) {
            push @completed_checkpoints, $checkpoint;
        }
    }
    
    print $fh "Completed checkpoints: " . join(", ", @completed_checkpoints) . "\n";
    
    if (-s "consensi.fa") {
        my $count = `grep -c '^>' consensi.fa`;
        chomp $count;
        print $fh "Final consensus sequences: $count\n";
        print $fh "Final consensus size: " . (-s "consensi.fa") . " bytes\n";
    }
    
    close $fh;
    
    chdir $current_dir;
    
    log_message("INFO", "Masked track completed successfully");
    return 1;
}

sub run_masked_track {
    # Legacy compatibility wrapper
    return run_masked_track_independent(@_);
}

sub run_full_genome_track {
    my ($genome_file, $output_dir, $threads, $cpu_threads, $genome_size) = @_;
    
    log_message("INFO", "Running full genome track", "single_track_mode=true");
    
    # For small genomes, run everything in the main directory
    my $current_dir = getcwd();
    chdir $output_dir or die "Cannot change to $output_dir: $!\n";
    
    # Check if the entire process is already completed
    if (-f "full_genome_completed.ok") {
        log_message("INFO", "Full genome track already completed", "checkpoint found: full_genome_completed.ok");
        chdir $current_dir;
        return 1;
    }
    
    # Step 1: Create sequence name list (with checkpoint)
    if (!-f "seq_naming.ok") {
        log_message("INFO", "Creating sequence name list for full genome");
        create_seq_name_list($genome_file);
        
        # Create sequence naming checkpoint
        open(my $fh, '>', "seq_naming.ok") or die "Cannot create sequence naming checkpoint: $!\n";
        print $fh "Sequence naming completed at " . localtime() . "\n";
        my $name_count = -s "seq.names" ? `head -1 seq.names` : 0;
        chomp $name_count;
        print $fh "Sequence names created: $name_count\n";
        close $fh;
        log_message("INFO", "Sequence naming checkpoint created", "names=$name_count, file=seq_naming.ok");
    } else {
        log_message("INFO", "Sequence naming already completed", "checkpoint found: seq_naming.ok");
    }
    
    # Step 2: Self-alignment (with checkpoint)
    my $blast_output = "self_alignment.blast";
    if (!-f "rmblastn.ok") {
        log_message("INFO", "Running self-alignment on full genome");
        my $hit_count = run_rmblastn_self_alignment($genome_file, $blast_output, $threads);
        
        if ($hit_count == 0) {
            log_message("WARN", "No BLAST hits found in full genome");
            chdir $current_dir;
            return 0;
        }
        
        # Create RMBlastN checkpoint
        open(my $fh, '>', "rmblastn.ok") or die "Cannot create RMBlastN checkpoint: $!\n";
        print $fh "RMBlastN completed at " . localtime() . "\n";
        print $fh "BLAST hits found: $hit_count\n";
        print $fh "Output file size: " . (-s $blast_output || 0) . " bytes\n";
        close $fh;
        log_message("INFO", "RMBlastN checkpoint created", "hits=$hit_count, file=rmblastn.ok");
    } else {
        log_message("INFO", "RMBlastN already completed", "checkpoint found: rmblastn.ok");
    }
    
    # Step 3: MSP collection (with checkpoint)
    if (!-f "msp_collection.ok") {
        log_message("INFO", "Collecting MSPs from full genome");
        
        # Handle the case where blast output might already be in MSP format
        if (-s $blast_output) {
            # Check if file is already in MSP format by examining first line
            my $first_line = `head -1 $blast_output`;
            chomp $first_line;
            
            if ($first_line =~ /^\d{6}\s+\d+\s+\d{5}\s+\d{5}/) {
                # File is already in MSP format, just copy it
                log_message("INFO", "BLAST output already in MSP format", "copying to msp.out");
                run_cmd("cp $blast_output msp.out");
            } else {
                # File is in BLAST format, need MSPCollect.pl
                log_message("INFO", "Converting BLAST to MSP format", "using MSPCollect.pl");
                run_cmd("MSPCollect.pl $blast_output > msp.out");
            }
        }
        
        if (!-s "msp.out") {
            log_message("WARN", "No MSPs generated from full genome");
            chdir $current_dir;
            return 0;
        }
        
        # Create MSP collection checkpoint
        open(my $fh, '>', "msp_collection.ok") or die "Cannot create MSP collection checkpoint: $!\n";
        print $fh "MSP collection completed at " . localtime() . "\n";
        my $msp_lines = `wc -l < msp.out`;
        chomp $msp_lines;
        print $fh "MSP lines: $msp_lines\n";
        print $fh "File size: " . (-s "msp.out") . " bytes\n";
        close $fh;
        log_message("INFO", "MSP collection checkpoint created", "lines=$msp_lines, file=msp_collection.ok");
    } else {
        log_message("INFO", "MSP collection already completed", "checkpoint found: msp_collection.ok");
    }
    
    # Step 4: RECON pipeline (with checkpoint)
    if (!-f "recon.ok") {
        my $k_param = determine_k_parameter("msp.out");
        log_message("INFO", "Starting full genome RECON pipeline", "k_parameter=$k_param");
        
        eval {
            run_recon_pipeline($k_param);
            
            # Create RECON checkpoint
            open(my $fh, '>', "recon.ok") or die "Cannot create RECON checkpoint: $!\n";
            print $fh "RECON pipeline completed at " . localtime() . "\n";
            print $fh "K parameter used: $k_param\n";
            if (-d "summary") {
                my $families = -s "summary/families" ? `wc -l < summary/families` : 0;
                chomp $families;
                print $fh "Families found: $families\n";
            }
            close $fh;
            log_message("INFO", "RECON checkpoint created", "k_param=$k_param, file=recon.ok");
        };
        
        if ($@) {
            log_message("ERROR", "Full genome RECON failed", "error=$@");
            chdir $current_dir;
            die "Full genome RECON failure. Cannot proceed.\n";
        }
    } else {
        log_message("INFO", "RECON pipeline already completed", "checkpoint found: recon.ok");
    }
    
    # Step 5: Build consensus (with checkpoint)
    if (-d "summary" && -s "summary/families") {
        if (!-f "consensus.ok") {
            log_message("INFO", "Building full genome consensus sequences");
            my $local_genome = abs_path($genome_file);
            # Use more threads for build_for_RECON since it's the final step
            my $build_threads = $threads;  # Use all available threads  
            run_cmd("build_for_RECON ./ $local_genome $build_threads");
            
            if (-s "consensi.fa") {
                my $count = `grep -c '^>' consensi.fa`;
                chomp $count;
                
                # Create consensus checkpoint
                open(my $fh, '>', "consensus.ok") or die "Cannot create consensus checkpoint: $!\n";
                print $fh "Consensus building completed at " . localtime() . "\n";
                print $fh "Consensus sequences: $count\n";
                print $fh "File size: " . (-s "consensi.fa") . " bytes\n";
                close $fh;
                
                log_message("INFO", "Full genome consensus built", "sequences=$count");
                log_message("INFO", "Consensus checkpoint created", "sequences=$count, file=consensus.ok");
            } else {
                log_message("WARN", "No consensus sequences generated despite families found");
            }
        } else {
            log_message("INFO", "Consensus building already completed", "checkpoint found: consensus.ok");
        }
    } else {
        log_message("WARN", "No families found in full genome", "skipping_consensus_building");
    }
    
    # Create completion marker
    open(my $fh, '>', "full_genome_completed.ok") or die "Cannot create completion marker: $!\n";
    print $fh "Full genome track completed at " . localtime() . "\n";
    
    # Record which checkpoints were successful
    my @checkpoint_files = qw(seq_naming.ok rmblastn.ok msp_collection.ok recon.ok consensus.ok);
    my @completed_checkpoints;
    
    for my $checkpoint (@checkpoint_files) {
        if (-f $checkpoint) {
            push @completed_checkpoints, $checkpoint;
        }
    }
    
    print $fh "Completed checkpoints: " . join(", ", @completed_checkpoints) . "\n";
    
    if (-s "consensi.fa") {
        my $count = `grep -c '^>' consensi.fa`;
        chomp $count;
        print $fh "Final consensus sequences: $count\n";
        print $fh "Final consensus size: " . (-s "consensi.fa") . " bytes\n";
    }
    
    close $fh;
    
    chdir $current_dir;
    
    log_message("INFO", "Full genome track completed");
    return 1;
}

sub extract_unmasked_regions {
    my ($genome_file, $output_file) = @_;
    
    log_message("INFO", "Extracting unmasked regions", "input=$genome_file, min_length=50bp");
    
    open my $input_fh, '<', $genome_file or die "Cannot open $genome_file: $!";
    open my $output_fh, '>', $output_file or die "Cannot create $output_file: $!";
    
    my ($current_header, $current_seq) = ('', '');
    my ($total_regions, $extracted_regions, $total_bp, $extracted_bp) = (0, 0, 0, 0);
    my $region_counter = 1;
    
    while (<$input_fh>) {
        chomp;
        if (/^>(.*)/) {
            # Process previous sequence
            if ($current_seq) {
                my ($regions, $bp) = process_sequence_for_unmasked($current_header, $current_seq, $output_fh, \$region_counter);
                $extracted_regions += $regions;
                $extracted_bp += $bp;
                $total_bp += length($current_seq);
            }
            
            $current_header = $1;
            $current_seq = '';
            $total_regions++;
        } else {
            $current_seq .= $_;
        }
    }
    
    # Process last sequence
    if ($current_seq) {
        my ($regions, $bp) = process_sequence_for_unmasked($current_header, $current_seq, $output_fh, \$region_counter);
        $extracted_regions += $regions;
        $extracted_bp += $bp;
        $total_bp += length($current_seq);
    }
    
    close $input_fh;
    close $output_fh;
    
    log_message("INFO", "Unmasked region extraction completed", 
                "original_sequences=$total_regions, extracted_regions=$extracted_regions, " .
                "original_bp=" . format_size($total_bp) . ", extracted_bp=" . format_size($extracted_bp));
    
    return $extracted_regions;
}

sub process_sequence_for_unmasked {
    my ($header, $sequence, $output_fh, $counter_ref) = @_;
    
    my $min_length = 50;  # Minimum length for unmasked regions
    my ($regions_count, $total_bp) = (0, 0);
    
    # Find unmasked regions by splitting only on N characters
    # Split on masked regions (N characters only)
    my @unmasked_segments = split /N+/i, $sequence;
    
    foreach my $segment (@unmasked_segments) {
        # Skip empty segments and short segments
        next if length($segment) < $min_length;
        
        # Clean up the segment (remove any remaining Ns)
        $segment =~ s/[Nn]//g;
        next if length($segment) < $min_length;
        
        # Only keep segments with valid DNA bases (allow both upper and lowercase)
        next unless $segment =~ /^[ATGCatgc]+$/;
        
        # Write as gi|N format
        print $output_fh ">gi|" . $$counter_ref . " unmasked_region_" . $$counter_ref . " from_$header\n";
        print $output_fh "$segment\n";
        
        $$counter_ref++;
        $regions_count++;
        $total_bp += length($segment);
    }
    
    return ($regions_count, $total_bp);
}

1;

__END__

=head1 NAME

RECON::MaskedTrack - Masked track implementation for RECON pipeline

=head1 DESCRIPTION

Implements the masked track processing for high-quality consensus building
from masked genome sequences.

=cut