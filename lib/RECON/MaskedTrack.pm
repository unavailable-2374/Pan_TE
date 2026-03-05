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
    run_masked_track_independent
    extract_unmasked_regions extract_unmasked_regions_extended
);

sub run_masked_track_independent {
    my @args = @_;
    my $current_dir = getcwd();
    my $success;
    eval { $success = _run_masked_track_impl(@args); };
    my $err = $@;
    chdir $current_dir;  # ALWAYS restore working directory
    die $err if $err;
    return $success;
}

sub _run_masked_track_impl {
    my ($genome_file, $output_dir, $threads, $cpu_threads, $genome_size, $original_genome, $masking_type) = @_;
    $masking_type ||= 'hard';

    chdir $output_dir or die "Cannot change to $output_dir: $!\n";

    log_message("INFO", "Starting masked track processing",
                "genome_size=" . format_size($genome_size) . ", masking_type=$masking_type");

    # Check if the entire process is already completed
    if (-f "masked_track_completed.ok") {
        log_message("INFO", "Masked track already completed", "checkpoint found: masked_track_completed.ok");
        return 1;
    }

    # Step 1: Prepare working genome
    my $working_genome;

    if ($masking_type eq 'soft') {
        # Soft masking: use the entire soft-masked genome directly (no extraction)
        # rmblastn -lcase_masking will skip seeding in lowercase regions
        # but can extend alignments through them — preserving full context
        my $renamed_file = "genome_renamed.fa";
        if (!-f "genome_prep.ok") {
            log_message("INFO", "Soft masking mode: rename to gi|N and fragment >20MB sequences",
                        "genome=$genome_file, max_fragment=20MB");
            # Rename to gi|N and split any sequence >20MB into multiple gi|N entries
            # This ensures all downstream files (seq.names, database, query chunks) have consistent IDs
            # and every chunk file stays ≤20MB
            my $total_gi = rename_and_fragment_fasta($genome_file, $renamed_file, 20);

            if (!-s $renamed_file) {
                log_message("ERROR", "Failed to prepare genome sequences", "file=$renamed_file");
                return 0;
            }

            open(my $fh, '>', "genome_prep.ok") or die "Cannot create genome prep checkpoint: $!\n";
            print $fh "Genome preparation completed at " . localtime() . "\n";
            print $fh "Mode: soft masking (rename + fragment)\n";
            print $fh "Sequences (gi|N): $total_gi\n";
            print $fh "File size: " . (-s $renamed_file) . " bytes\n";
            close $fh;
            log_message("INFO", "Genome prepared for RECON",
                        "gi_sequences=$total_gi, mode=soft_masking");
        } else {
            log_message("INFO", "Genome preparation already completed", "checkpoint found: genome_prep.ok");
        }
        $working_genome = $renamed_file;
    } else {
        # Hard masking: extract unmasked (non-N) regions
        my $unmasked_file = "unmasked_regions.fa";
        if (!-f "unmasked_extraction.ok") {
            if ($original_genome && -f $original_genome) {
                log_message("INFO", "Extracting unmasked regions with flanking extension",
                            "min_length=50bp, flank=500bp, original_genome=$original_genome");
                extract_unmasked_regions_extended($genome_file, $unmasked_file, $original_genome, 'hard');
            } else {
                log_message("INFO", "Extracting unmasked regions", "min_length=50bp");
                extract_unmasked_regions($genome_file, $unmasked_file, 'hard');
            }

            if (!-s $unmasked_file) {
                log_message("ERROR", "Failed to extract unmasked regions", "file=$unmasked_file");
                return 0;
            }

            open(my $fh, '>', "unmasked_extraction.ok") or die "Cannot create unmasked extraction checkpoint: $!\n";
            print $fh "Unmasked extraction completed at " . localtime() . "\n";
            my $seq_count = `grep -c '^>' $unmasked_file`;
            chomp $seq_count;
            print $fh "Extracted sequences: $seq_count\n";
            print $fh "File size: " . (-s $unmasked_file) . " bytes\n";
            close $fh;
            log_message("INFO", "Unmasked extraction checkpoint created", "sequences=$seq_count");
        } else {
            log_message("INFO", "Unmasked extraction already completed", "checkpoint found: unmasked_extraction.ok");
        }
        $working_genome = $unmasked_file;
    }

    my $seq_count = `grep -c '^>' $working_genome`;
    chomp $seq_count;
    log_message("INFO", "Working genome ready", "sequences=$seq_count, file=$working_genome");

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
            log_message("INFO", "No BLAST hits found - masked genome has no new repeats", "creating empty consensi.fa");
            open(my $empty_fh, '>', "consensi.fa") or die "Cannot create empty consensi.fa: $!\n";
            close($empty_fh);
            open(my $ok_fh, '>', "masked_track_completed.ok") or die "Cannot create completion marker: $!\n";
            print $ok_fh "Completed at " . localtime() . "\nNo BLAST hits - no new repeats found\n";
            close($ok_fh);
            return 1;
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

            if ($first_line =~ /^\d+\s+\d+\s+\d+\s+\d+/) {
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
        genome_prep.ok unmasked_extraction.ok
        seq_naming.ok rmblastn.ok msp_collection.ok recon.ok consensus.ok
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

    log_message("INFO", "Masked track completed successfully");
    return 1;
}

sub extract_unmasked_regions_extended {
    my ($masked_genome_file, $output_file, $original_genome_file, $masking_type) = @_;
    $masking_type ||= 'hard';

    my $flank_size = 500;       # Extend each side by 500bp
    my $merge_gap  = 200;       # Merge regions within 200bp of each other
    my $min_length = 50;        # Minimum unmasked region length

    # Soft masking: extract from the masked genome itself (preserves lowercase for rmblastn)
    # Hard masking: extract from original genome (N regions have no sequence info)
    my $extraction_source = ($masking_type eq 'soft') ? $masked_genome_file : $original_genome_file;

    log_message("INFO", "Extracting unmasked regions with flanking extension (streaming)",
                "masked=$masked_genome_file, extraction_source=$extraction_source, " .
                "masking_type=$masking_type, flank=${flank_size}bp, merge_gap=${merge_gap}bp");

    # Phase 1: Stream masked genome line-by-line to find unmasked region coordinates
    # Hard mask: N/n = masked, everything else = unmasked
    # Soft mask: lowercase (a-z) = masked, uppercase (A-Z) = unmasked
    open my $mask_fh, '<', $masked_genome_file
        or die "Cannot open $masked_genome_file: $!";

    my @all_regions;    # [{chr => name, start => N, end => N}, ...]
    my %chr_lengths;    # chr_name => length
    my $current_chr = '';
    my $position = 0;
    my $in_unmasked = 0;
    my $unmasked_start = 0;

    while (my $line = <$mask_fh>) {
        chomp $line;
        if ($line =~ /^>(\S+)/) {
            # Flush previous chromosome's last unmasked region
            if ($in_unmasked && $current_chr) {
                if (($position - $unmasked_start) >= $min_length) {
                    push @all_regions, { chr => $current_chr, start => $unmasked_start, end => $position };
                }
                $in_unmasked = 0;
            }
            $chr_lengths{$current_chr} = $position if $current_chr;
            $current_chr = $1;
            $position = 0;
        } else {
            for my $i (0 .. length($line) - 1) {
                my $c = substr($line, $i, 1);
                my $is_masked;
                if ($masking_type eq 'soft') {
                    $is_masked = ($c ge 'a' && $c le 'z');  # lowercase = soft-masked
                } else {
                    $is_masked = ($c eq 'N' || $c eq 'n');   # N = hard-masked
                }
                if ($is_masked) {
                    if ($in_unmasked) {
                        if (($position - $unmasked_start) >= $min_length) {
                            push @all_regions, { chr => $current_chr, start => $unmasked_start, end => $position };
                        }
                        $in_unmasked = 0;
                    }
                } else {
                    if (!$in_unmasked) {
                        $unmasked_start = $position;
                        $in_unmasked = 1;
                    }
                }
                $position++;
            }
        }
    }
    # Flush last chromosome
    if ($in_unmasked && $current_chr) {
        if (($position - $unmasked_start) >= $min_length) {
            push @all_regions, { chr => $current_chr, start => $unmasked_start, end => $position };
        }
    }
    $chr_lengths{$current_chr} = $position if $current_chr;
    close $mask_fh;

    log_message("INFO", "Found unmasked regions before extension",
                "regions=" . scalar(@all_regions));

    # Phase 2: Extend and merge regions per chromosome (coordinate-only, no sequence data)
    my %chr_regions;
    for my $r (@all_regions) {
        push @{$chr_regions{$r->{chr}}}, $r;
    }

    my @final_regions;
    for my $chr (sort keys %chr_regions) {
        my $chr_len = $chr_lengths{$chr} || 0;
        my @regions = sort { $a->{start} <=> $b->{start} } @{$chr_regions{$chr}};

        # Extend each region by flank_size, clamped to chromosome bounds
        for my $r (@regions) {
            $r->{start} = ($r->{start} - $flank_size > 0) ? $r->{start} - $flank_size : 0;
            $r->{end}   = ($r->{end} + $flank_size < $chr_len) ? $r->{end} + $flank_size : $chr_len;
        }

        # Merge regions within merge_gap of each other
        my @merged;
        my $prev = $regions[0];
        for my $i (1 .. $#regions) {
            if ($regions[$i]->{start} - $prev->{end} <= $merge_gap) {
                $prev->{end} = $regions[$i]->{end} if $regions[$i]->{end} > $prev->{end};
            } else {
                push @merged, $prev;
                $prev = $regions[$i];
            }
        }
        push @merged, $prev;
        push @final_regions, @merged;
    }

    log_message("INFO", "Regions after extension and merging",
                "merged_regions=" . scalar(@final_regions));

    # Phase 3: Write BED file and extract from appropriate source
    # Soft mask: extract from masked genome (preserves lowercase for rmblastn soft_masking)
    # Hard mask: extract from original genome (recovers real DNA at N-masked positions)
    my $bed_file = "${output_file}.regions.bed";

    if (! -f "${extraction_source}.fai") {
        log_message("INFO", "Creating samtools index",
                    "file=$extraction_source");
        run_cmd("samtools faidx $extraction_source");
    }

    open my $bed_fh, '>', $bed_file or die "Cannot create $bed_file: $!";
    for my $r (@final_regions) {
        print $bed_fh "$r->{chr}\t$r->{start}\t$r->{end}\n";
    }
    close $bed_fh;

    my $raw_fasta = "${output_file}.raw";
    run_cmd("bedtools getfasta -fi $extraction_source -bed $bed_file > $raw_fasta");

    # Phase 4: Reformat bedtools output headers to gi|N format and apply quality filters
    open my $raw_fh, '<', $raw_fasta or die "Cannot open $raw_fasta: $!";
    open my $out_fh, '>', $output_file or die "Cannot create $output_file: $!";

    my $region_counter = 1;
    my ($total_bp, $extracted_count) = (0, 0);
    my ($current_header, $current_seq_data) = ('', '');

    while (my $fasta_line = <$raw_fh>) {
        chomp $fasta_line;
        if ($fasta_line =~ /^>(.+)/) {
            # Process previous sequence
            if ($current_header && $current_seq_data) {
                my $len = length($current_seq_data);
                if ($len >= $min_length
                    && $current_seq_data !~ /^N+$/i
                    && $current_seq_data =~ /^[ATGCNatgcn]+$/) {
                    print $out_fh ">gi|$region_counter $current_header extended_unmasked\n";
                    print $out_fh "$current_seq_data\n";
                    $region_counter++;
                    $extracted_count++;
                    $total_bp += $len;
                }
            }
            $current_header = $1;
            $current_seq_data = '';
        } else {
            $current_seq_data .= $fasta_line;
        }
    }
    # Process last sequence
    if ($current_header && $current_seq_data) {
        my $len = length($current_seq_data);
        if ($len >= $min_length
            && $current_seq_data !~ /^N+$/i
            && $current_seq_data =~ /^[ATGCNatgcn]+$/) {
            print $out_fh ">gi|$region_counter $current_header extended_unmasked\n";
            print $out_fh "$current_seq_data\n";
            $region_counter++;
            $extracted_count++;
            $total_bp += $len;
        }
    }

    close $raw_fh;
    close $out_fh;

    # Clean up temporary files
    unlink($bed_file, $raw_fasta);

    log_message("INFO", "Extended unmasked region extraction completed (streaming)",
                "masking_type=$masking_type, extracted_regions=$extracted_count, " .
                "total_bp=" . format_size($total_bp) .
                ", avg_length=" . ($extracted_count > 0 ? int($total_bp / $extracted_count) : 0));

    return $extracted_count;
}

sub extract_unmasked_regions {
    my ($genome_file, $output_file, $masking_type) = @_;
    $masking_type ||= 'hard';

    log_message("INFO", "Extracting unmasked regions",
                "input=$genome_file, masking_type=$masking_type, min_length=50bp");

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
                my ($regions, $bp) = process_sequence_for_unmasked($current_header, $current_seq, $output_fh, \$region_counter, $masking_type);
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
        my ($regions, $bp) = process_sequence_for_unmasked($current_header, $current_seq, $output_fh, \$region_counter, $masking_type);
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
    my ($header, $sequence, $output_fh, $counter_ref, $masking_type) = @_;
    $masking_type ||= 'hard';

    my $min_length = 50;  # Minimum length for unmasked regions
    my ($regions_count, $total_bp) = (0, 0);

    # Split on masked regions:
    #   Hard mask: N runs
    #   Soft mask: lowercase runs
    my @unmasked_segments;
    if ($masking_type eq 'soft') {
        @unmasked_segments = split /[a-z]+/, $sequence;
    } else {
        @unmasked_segments = split /N+/i, $sequence;
    }

    foreach my $segment (@unmasked_segments) {
        # Skip empty segments and short segments
        next if length($segment) < $min_length;

        # Clean up the segment
        if ($masking_type eq 'soft') {
            # Soft mask: segments are already pure uppercase from split
            next unless $segment =~ /^[ATGCN]+$/;
        } else {
            # Hard mask: remove any remaining Ns
            $segment =~ s/[Nn]//g;
            next if length($segment) < $min_length;
            next unless $segment =~ /^[ATGCatgc]+$/;
        }
        
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