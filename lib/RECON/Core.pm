package RECON::Core;

use strict;
use warnings;
use Exporter 'import';
use File::Path qw(make_path);
use FindBin qw($Bin);
use RECON::Logger;
use RECON::Utils;

# RepeatModeler2-compatible scoring matrix for rmblastn
my $MATRIX_DIR = "$Bin/Matrices/ncbi/nt";
$ENV{BLASTMAT} = $MATRIX_DIR if -d $MATRIX_DIR;

our @EXPORT = qw(
    determine_k_parameter run_recon_pipeline
    run_rmblastn_self_alignment merge_msp_files
    check_checkpoint create_checkpoint remove_checkpoint remove_all_checkpoints
    cleanup_intermediate_files
);

sub check_checkpoint {
    my ($checkpoint_file) = @_;
    return -f $checkpoint_file;
}

sub create_checkpoint {
    my ($checkpoint_file, $step_name) = @_;
    
    open(my $fh, '>', $checkpoint_file) or die "Cannot create checkpoint $checkpoint_file: $!\n";
    print $fh "Completed at " . localtime() . "\n";
    print $fh "Step: $step_name\n";
    close $fh;
    
    log_message("INFO", "Checkpoint created", "file=$checkpoint_file, step=$step_name");
}

sub remove_checkpoint {
    my ($checkpoint_file) = @_;
    
    if (-f $checkpoint_file) {
        unlink($checkpoint_file);
        log_message("INFO", "Checkpoint removed", "file=$checkpoint_file");
        return 1;
    }
    return 0;
}

sub remove_all_checkpoints {
    my @checkpoint_files = ("imagespread.ok", "eledef.ok", "eleredef.ok", "edgeredef.ok", "famdef.ok");
    my $removed_count = 0;
    
    for my $checkpoint (@checkpoint_files) {
        if (remove_checkpoint($checkpoint)) {
            $removed_count++;
        }
    }
    
    log_message("INFO", "All RECON checkpoints removed", "count=$removed_count");
    return $removed_count;
}

sub cleanup_intermediate_files {
    my ($force) = @_;

    # Only cleanup if consensi.fa was successfully generated or force flag is set
    if (!$force && !(-s "consensi.fa")) {
        log_message("INFO", "Skipping cleanup", "consensi.fa not found or empty");
        return 0;
    }

    # Check debug flag
    my $debug = $ENV{'RECON_DEBUG'} || 0;
    if ($debug) {
        log_message("INFO", "Skipping cleanup", "RECON_DEBUG=$debug");
        return 0;
    }

    my @cleanup_dirs = ("ele_def_res", "ele_redef_res", "edge_redef_res");
    my $cleaned_count = 0;
    my $cleaned_files = 0;

    # Clean up RECON intermediate directories
    for my $dir (@cleanup_dirs) {
        if (-d $dir) {
            run_cmd("rm -rf $dir");
            $cleaned_count++;
            log_message("INFO", "Cleaned up intermediate directory", "dir=$dir");
        }
    }

    # Clean up RMBlastN chunk files (residual .blast and .msp files)
    my @rmblast_blast_files = glob("rmblast_chunk_*.blast");
    my @rmblast_msp_files = glob("rmblast_chunk_*.msp");
    my @rmblast_chunk_fa_files = glob("rmblast_chunk_*.fa");

    for my $file (@rmblast_blast_files, @rmblast_msp_files, @rmblast_chunk_fa_files) {
        if (-f $file) {
            unlink($file);
            $cleaned_files++;
        }
    }
    if (@rmblast_blast_files || @rmblast_msp_files || @rmblast_chunk_fa_files) {
        log_message("INFO", "Cleaned up RMBlastN chunk files",
                    "blast=" . scalar(@rmblast_blast_files) .
                    ", msp=" . scalar(@rmblast_msp_files) .
                    ", fa=" . scalar(@rmblast_chunk_fa_files));
    }

    # Clean up BLAST database files
    my @blast_db_files = glob("sample_db.n*");
    for my $file (@blast_db_files) {
        if (-f $file) {
            unlink($file);
            $cleaned_files++;
        }
    }
    if (@blast_db_files) {
        log_message("INFO", "Cleaned up BLAST database files", "count=" . scalar(@blast_db_files));
    }

    # Clean up RepeatMasker intermediate files
    my @rm_cat_files = glob("*.cat *.cat.gz");
    my @rm_tbl_files = glob("*.tbl");
    my @rm_out_files = glob("sample_masked.fa.out sample.fa.out");
    my @rm_dirs = glob("RM_*");

    for my $file (@rm_cat_files, @rm_tbl_files, @rm_out_files) {
        if (-f $file) {
            unlink($file);
            $cleaned_files++;
        }
    }
    for my $dir (@rm_dirs) {
        if (-d $dir) {
            run_cmd("rm -rf $dir");
            $cleaned_count++;
        }
    }
    if (@rm_cat_files || @rm_tbl_files || @rm_out_files || @rm_dirs) {
        log_message("INFO", "Cleaned up RepeatMasker intermediate files",
                    "cat=" . scalar(@rm_cat_files) .
                    ", tbl=" . scalar(@rm_tbl_files) .
                    ", out=" . scalar(@rm_out_files) .
                    ", dirs=" . scalar(@rm_dirs));
    }

    # Clean up merged consensus temporary file
    if (-f "merged_consensus.fa") {
        unlink("merged_consensus.fa");
        $cleaned_files++;
        log_message("INFO", "Cleaned up merged consensus file", "file=merged_consensus.fa");
    }

    if ($cleaned_count > 0 || $cleaned_files > 0) {
        log_message("INFO", "Intermediate file cleanup completed",
                    "directories_removed=$cleaned_count, files_removed=$cleaned_files, " .
                    "consensi_size=" . (-s "consensi.fa" || 0));
    }

    return $cleaned_count + $cleaned_files;
}

sub determine_k_parameter {
    my ($msp_file) = @_;
    
    return 1 unless -s $msp_file; # Default fallback
    
    # Get MSP file size for K parameter determination (RECON documentation recommendation)
    my $msp_size_bytes = -s $msp_file;
    my $msp_size_gb = $msp_size_bytes / (1024 * 1024 * 1024);
    
    # Count MSP entries for additional information
    my $msp_count = `wc -l < $msp_file`;
    chomp $msp_count;
    
    # RECON documentation recommends K values based on MSP file size for stability
    my $k_param;
    if ($msp_size_gb < 1) {
        $k_param = 1;
    } elsif ($msp_size_gb < 5) {
        $k_param = 4;
    } elsif ($msp_size_gb < 20) {
        $k_param = 16;
    } else {
        $k_param = 64;
    }
    
    log_message("INFO", "Determined K parameter based on MSP size", 
                sprintf("msp_size=%.2fGB, msp_lines=%d, k=%d", 
                        $msp_size_gb, $msp_count, $k_param));
    
    return $k_param;
}

sub run_recon_pipeline {
    my ($k_param) = @_;
    $k_param ||= 1;  # Use default K=1 for small MSP files
    
    log_message("INFO", "Starting RECON pipeline", "k_parameter=$k_param");
    
    # Validate input files and create RECON-compatible names
    if (!-f "seq.names") {
        die "ERROR: seq.names not found. Cannot proceed with RECON.\n";
    }
    if (!-f "msp.out") {
        die "ERROR: msp.out not found. Cannot proceed with RECON.\n";
    }
    
    # Create RECON-compatible symlinks
    unlink("seqnames", "msps.out");
    symlink("seq.names", "seqnames") or die "Cannot create seqnames symlink: $!\n";
    symlink("msp.out", "msps.out") or die "Cannot create msps.out symlink: $!\n";
    
    # Create necessary directories (following RepeatModeler pattern)
    run_cmd("mkdir images") unless -d "images";
    run_cmd("mkdir summary") unless -d "summary";
    
    my @recon_steps = (
        { cmd => "imagespread seqnames msps.out $k_param", 
          check => "images/spread1", name => "Image spreading", checkpoint => "imagespread.ok", post => \&process_imagespread },
        { cmd => "eledef seqnames msps.out single", 
          check => "ele_def_res", name => "Element definition", checkpoint => "eledef.ok", pre => \&setup_eledef, post => \&process_eledef },
        { cmd => "eleredef seqnames", 
          check => "ele_redef_res", name => "Element redefinition", checkpoint => "eleredef.ok", pre => \&setup_eleredef, post => \&cleanup_eleredef },
        { cmd => "edgeredef seqnames", 
          check => "edge_redef_res", name => "Edge redefinition", checkpoint => "edgeredef.ok", pre => \&setup_edgeredef, post => \&cleanup_edgeredef },
        { cmd => "famdef seqnames", 
          check => "families", name => "Family definition", checkpoint => "famdef.ok", pre => \&setup_famdef, post => \&cleanup_famdef }
    );
    
    # Log existing checkpoints for resume information
    my @existing_checkpoints;
    for my $step (@recon_steps) {
        if (check_checkpoint($step->{checkpoint})) {
            push @existing_checkpoints, $step->{name};
        }
    }
    
    if (@existing_checkpoints) {
        log_message("INFO", "Found existing checkpoints", 
                    "completed_steps=" . join(", ", @existing_checkpoints));
    } else {
        log_message("INFO", "No checkpoints found", "starting_from_beginning");
    }
    
    for my $step (@recon_steps) {
        # Check if step already completed (checkpoint exists)
        if (check_checkpoint($step->{checkpoint})) {
            log_message("INFO", "Skipping $step->{name}", "checkpoint_found=$step->{checkpoint}");
            next;
        }
        
        log_message("INFO", "Running: $step->{name}");
        
        # Pre-processing if defined
        if ($step->{pre}) {
            $step->{pre}->();
        }
        
        eval {
            # Handle gmon.out cleanup like RepeatModeler
            my $gmon_backup = lc($step->{name});
            $gmon_backup =~ s/\s+//g;
            $gmon_backup .= "-gmon.out";
            
            run_cmd($step->{cmd});
            
            # Move gmon.out if it exists (matching RepeatModeler pattern)
            if (-f "gmon.out") {
                run_cmd("mv gmon.out $gmon_backup");
            }
        };
        
        if ($@ || !check_step_success($step)) {
            log_message("ERROR", "$step->{name} failed", "error=$@");
            die "RECON step failed: $step->{name}\n";
        }
        
        # Post-processing if defined
        if ($step->{post}) {
            if ($step->{name} eq "Image spreading") {
                $step->{post}->($k_param);
            } else {
                $step->{post}->();
            }
        }
        
        # Create checkpoint after successful completion
        create_checkpoint($step->{checkpoint}, $step->{name});
        
        log_message("INFO", "$step->{name} completed", "output=$step->{check}");
    }
    
    # Move families file to summary if it exists
    if (-f "families") {
        run_cmd("mv families summary/") if -d "summary";
    }
    
    # Only cleanup symlinks, preserve intermediate directories for consensus building
    unlink("seqnames", "msps.out");
    
    log_message("INFO", "RECON pipeline completed", "results in summary/");
    
    # Note: Intermediate directories (ele_def_res, ele_redef_res, edge_redef_res) 
    # are preserved until final consensus is built by build_for_RECON
}

sub process_imagespread {
    # Process images following RepeatModeler pattern exactly
    # $sect should equal $k_param from imagespread step
    my ($k_param) = @_;
    my $sect = $k_param || 1;
    
    for (my $i = 1; $i <= $sect; $i++) {
        my $spread = "images/spread$i";
        if (-f $spread) {
            run_cmd("sort -T . -k 3,3 -k 4n,4n -k 5nr,5nr $spread >> images/images_sorted");
        }
    }
    
    # Remove spread files like RepeatModeler
    unlink glob "images/spread*";
    
    if (-f "images/images_sorted") {
        log_message("INFO", "Imagespread results processed", "sorted=images/images_sorted, sections=$sect");
    } else {
        log_message("WARN", "Imagespread did not produce results", "file=images/spread*");
    }
}

sub setup_eledef {
    # Remove and recreate ele_def_res directory (following RepeatModeler pattern)
    run_cmd("rm -rf ele_def_res") if -d "ele_def_res";
    run_cmd("mkdir ele_def_res");
    log_message("INFO", "Created ele_def_res directory");
}

sub setup_eleredef {
    # Remove and recreate ele_redef_res directory, then create symlinks
    run_cmd("rm -rf ele_redef_res") if -d "ele_redef_res";
    run_cmd("mkdir ele_redef_res");
    
    # Create symlinks for eleredef (required by RECON workflow)
    # tmp -> ele_def_res, tmp2 -> ele_redef_res
    unlink("tmp", "tmp2");
    symlink("ele_def_res", "tmp") or die "Cannot create symlink tmp -> ele_def_res: $!\n";
    symlink("ele_redef_res", "tmp2") or die "Cannot create symlink tmp2 -> ele_redef_res: $!\n";
    log_message("INFO", "Created eleredef symlinks", "tmp->ele_def_res, tmp2->ele_redef_res");
}

sub cleanup_eleredef {
    # Remove symlinks after eleredef
    unlink("tmp", "tmp2");
}

sub setup_edgeredef {
    # Remove and recreate edge_redef_res directory, then create symlinks
    run_cmd("rm -rf edge_redef_res") if -d "edge_redef_res";  
    run_cmd("mkdir edge_redef_res");
    
    # Create symlinks for edgeredef (required by RECON workflow)  
    # tmp -> ele_redef_res, tmp2 -> edge_redef_res
    unlink("tmp", "tmp2");
    symlink("ele_redef_res", "tmp") or die "Cannot create symlink tmp -> ele_redef_res: $!\n";
    symlink("edge_redef_res", "tmp2") or die "Cannot create symlink tmp2 -> edge_redef_res: $!\n";
    log_message("INFO", "Created edgeredef symlinks", "tmp->ele_redef_res, tmp2->edge_redef_res");
}

sub cleanup_edgeredef {
    # Remove symlinks after edgeredef
    unlink("tmp", "tmp2");
}

sub setup_famdef {
    # Create symlink for famdef (required by RECON workflow)
    # tmp -> edge_redef_res
    unlink("tmp");
    symlink("edge_redef_res", "tmp") or die "Cannot create symlink tmp -> edge_redef_res: $!\n";
    log_message("INFO", "Created famdef symlink", "tmp->edge_redef_res");
}

sub cleanup_famdef {
    # Remove symlink after famdef
    unlink("tmp");
}


sub check_step_success {
    my ($step) = @_;
    
    if ($step->{check} eq "ele_def_res" || $step->{check} eq "ele_redef_res" || $step->{check} eq "edge_redef_res") {
        return -d $step->{check} && glob("$step->{check}/*");
    } elsif ($step->{check} eq "images/spread1") {
        return -f "images/spread1";
    } elsif ($step->{check} eq "families") {
        # famdef creates families file in summary directory
        return -f "summary/families" || -f "families";
    } else {
        return -f $step->{check};
    }
}

sub process_eledef {
    # eledef creates files in ele_def_res directory (following RepeatModeler pattern)
    if (-d "ele_def_res" && glob("ele_def_res/*")) {
        my @ele_files = glob("ele_def_res/*");
        log_message("INFO", "Element definition completed", "files_created=" . scalar(@ele_files));
    } else {
        log_message("WARN", "Element definition did not produce results", "dir=ele_def_res");
    }
}

sub run_rmblastn_self_alignment {
    my ($query_file, $output_file, $threads, $database) = @_;
    $threads ||= 4;
    
    log_message("INFO", "Starting parallel RMBlastN self-alignment", 
                "query=$query_file, threads=$threads");
    
    # Create database if not provided
    unless ($database) {
        $database = "${query_file}_db";
        create_blast_database($query_file, $database);
    }
    
    # Check query file size to decide strategy
    my $query_size = -s $query_file;
    my $size_mb = $query_size / (1024 * 1024);

    if ($size_mb <= 20) {
        # ≤20MB: single process, no chunking needed
        return run_single_rmblastn_process($query_file, $output_file, $threads, $database);
    }

    # >20MB: chunk into ≤20MB pieces to control memory
    return run_parallel_rmblastn_chunked($query_file, $output_file, $threads, $database);
}

sub run_single_rmblastn_process {
    my ($query_file, $output_file, $threads, $database) = @_;

    log_message("INFO", "Using single-process RMBlastN", "small_file_optimization");

    # Run RMBlastN with RepeatModeler2-compatible parameters
    # comparison.matrix + complexity_adjust: proven RM2 scoring scheme
    # word_size 9: balances sensitivity vs memory for genome-scale self-alignment
    # lcase_masking: honor lowercase soft masking in query (skip seeding, allow extension)
    my $blast_tmp = "${output_file}.blast_tmp";
    my $cmd = "rmblastn -query $query_file -db $database " .
              "-out $blast_tmp " .
              "-num_threads $threads " .
              "-word_size 9 -matrix comparison.matrix " .
              "-gapopen 20 -gapextend 5 " .
              "-complexity_adjust -lcase_masking -evalue 1e-5 " .
              "-max_target_seqs 50000";

    run_cmd($cmd);

    # Convert BLAST output to MSP format (consistent with parallel/sequential paths)
    if (-s $blast_tmp) {
        run_cmd("MSPCollect.pl $blast_tmp > $output_file");
        unlink($blast_tmp);
    } else {
        # No hits: create empty output
        open(my $fh, '>', $output_file);
        close($fh);
        unlink($blast_tmp);
    }

    my $hit_count = `wc -l < $output_file`;
    chomp $hit_count;
    log_message("INFO", "Single RMBlastN completed", "msp_lines=$hit_count, output=$output_file");

    return $hit_count;
}

sub run_parallel_rmblastn_chunked {
    my ($query_file, $output_file, $threads, $database) = @_;

    # Split query into ≤20MB chunks to cap per-process memory
    my @chunk_files = split_fasta_by_size($query_file, "rmblast_chunk", 20);
    my $actual_chunks = @chunk_files;

    if ($actual_chunks <= 1) {
        return run_single_rmblastn_process($query_file, $output_file, $threads, $database);
    }

    # Process-level parallelism: 1 thread per process, maximize concurrency
    # rmblastn internal threading (-num_threads) is unreliable,
    # process-level parallelism via ForkManager is more effective
    my $num_processes = $threads;
    $num_processes = $actual_chunks if $num_processes > $actual_chunks;
    my $threads_per_process = 1;

    log_message("INFO", "Using parallel RMBlastN strategy",
                "total_threads=$threads, processes=$num_processes, threads_per_process=$threads_per_process, chunks=$actual_chunks");
    
    # Use Parallel::ForkManager for parallel processing
    eval { require Parallel::ForkManager; };
    if ($@) {
        log_message("WARN", "Parallel::ForkManager not available, using sequential processing", "error=$@");
        return run_sequential_chunked($query_file, $output_file, $threads, $database, \@chunk_files);
    }
    
    my $pm = Parallel::ForkManager->new($num_processes);
    my @msp_files;
    
    # Process chunks in parallel
    for my $i (0..$#chunk_files) {
        my $chunk_file = $chunk_files[$i];
        my $chunk_blast = "rmblast_chunk_${i}.blast";
        my $chunk_msp = "rmblast_chunk_${i}.msp";
        
        push @msp_files, $chunk_msp;
        
        $pm->start and next;  # Fork process
        
        # Child process
        eval {
            # Run rmblastn with RepeatModeler2-compatible parameters
            my $cmd = "rmblastn -query $chunk_file -db $database " .
                     "-out $chunk_blast " .
                     "-num_threads $threads_per_process " .
                     "-word_size 9 -matrix comparison.matrix " .
                     "-gapopen 20 -gapextend 5 " .
                     "-complexity_adjust -lcase_masking -evalue 1e-5 " .
                     "-max_target_seqs 50000";

            run_cmd($cmd);

            # Convert to MSP immediately and cleanup
            if (-s $chunk_blast) {
                run_cmd("MSPCollect.pl $chunk_blast > $chunk_msp");
                unlink($chunk_blast);  # Delete blast output to save space
            } else {
                # Create empty MSP file
                open(my $fh, '>', $chunk_msp);
                close($fh);
            }
        };
        
        if ($@) {
            log_message("ERROR", "Chunk processing failed", "chunk=$i, error=$@");
            # Create empty MSP file to avoid breaking the merge
            open(my $fh, '>', $chunk_msp);
            close($fh);
        }
        
        $pm->finish;  # Exit child process
    }
    
    # Wait for all processes to complete
    $pm->wait_all_children;
    
    log_message("INFO", "All parallel RMBlastN processes completed", "merging_results");
    
    # Merge all MSP files directly - no need to convert back to BLAST format
    # The output_file parameter is expected to be the final MSP file (msp.out)
    merge_msp_files(\@msp_files, $output_file);
    
    # Cleanup chunk files and intermediate MSPs
    for my $chunk_file (@chunk_files) {
        unlink($chunk_file);
    }
    for my $msp_file (@msp_files) {
        unlink($msp_file);
    }
    
    # Count final hits
    my $hit_count = `wc -l < $output_file`;
    chomp $hit_count;
    log_message("INFO", "Parallel RMBlastN completed", 
                "chunks_processed=$actual_chunks, total_hits=$hit_count");
    
    return $hit_count;
}

sub run_sequential_chunked {
    my ($query_file, $output_file, $threads, $database, $chunk_files_ref) = @_;
    
    log_message("INFO", "Running sequential chunked processing", "fallback_mode");
    
    my @chunk_files = @$chunk_files_ref;
    my @msp_files;
    
    for my $i (0..$#chunk_files) {
        my $chunk_file = $chunk_files[$i];
        my $chunk_blast = "rmblast_chunk_${i}.blast";
        my $chunk_msp = "rmblast_chunk_${i}.msp";
        
        push @msp_files, $chunk_msp;
        
        log_message("INFO", "Processing chunk $i", "file=$chunk_file");
        
        eval {
            # Run rmblastn with RepeatModeler2-compatible parameters
            my $cmd = "rmblastn -query $chunk_file -db $database " .
                     "-out $chunk_blast " .
                     "-num_threads $threads " .
                     "-word_size 9 -matrix comparison.matrix " .
                     "-gapopen 20 -gapextend 5 " .
                     "-complexity_adjust -lcase_masking -evalue 1e-5 " .
                     "-max_target_seqs 50000";

            run_cmd($cmd);

            # Convert to MSP and cleanup
            if (-s $chunk_blast) {
                run_cmd("MSPCollect.pl $chunk_blast > $chunk_msp");
                unlink($chunk_blast);
            } else {
                open(my $fh, '>', $chunk_msp);
                close($fh);
            }
        };
        
        if ($@) {
            log_message("WARN", "Chunk $i processing failed", "error=$@");
            open(my $fh, '>', $chunk_msp);
            close($fh);
        }
    }
    
    # Merge MSP files directly to output (same as parallel path)
    merge_msp_files(\@msp_files, $output_file);

    # Cleanup
    for my $chunk_file (@chunk_files) {
        unlink($chunk_file);
    }
    for my $msp_file (@msp_files) {
        unlink($msp_file);
    }

    my $hit_count = `wc -l < $output_file`;
    chomp $hit_count;
    return $hit_count;
}

sub merge_msp_files {
    my ($msp_files_ref, $output_file) = @_;
    my @msp_files = @$msp_files_ref;
    
    log_message("INFO", "Merging MSP files", "count=" . scalar(@msp_files) . ", output=$output_file");
    
    open(my $out_fh, '>', $output_file) or die "Cannot create $output_file: $!\n";
    
    for my $msp_file (@msp_files) {
        next unless -f $msp_file && -s $msp_file;
        
        open(my $in_fh, '<', $msp_file) or next;
        while (my $line = <$in_fh>) {
            print $out_fh $line;
        }
        close($in_fh);
    }
    
    close($out_fh);
}

1;

__END__

=head1 NAME

RECON::Core - Core RECON algorithm implementation

=head1 DESCRIPTION

Provides the core RECON pipeline functionality including element definition,
refinement, family definition, and edge refinement steps.

=cut