package RECON::SamplingTrack;

use strict;
use warnings;
use Exporter 'import';
use File::Path qw(make_path);
use Cwd qw(getcwd abs_path);
use RECON::Logger;
use RECON::Utils;
use RECON::Core;
use RECON::MaskedTrack;

our @EXPORT = qw(
    run_sampling_track_independent run_sampling_track_progressive run_sampling_track
    analyze_round_metrics evaluate_stopping_criteria update_accumulated_mask
    cleanup_failed_round log_debugging_info extract_unmasked_regions_with_gi_format
);

sub run_sampling_track_independent {
    my ($genome_file, $bed_files_ref, $masked_consensi, $output_dir, 
        $threads, $cpu_threads, $genome_size) = @_;
    
    my @bed_files = @$bed_files_ref;
    my $current_dir = getcwd();
    chdir $output_dir or die "Cannot change to $output_dir: $!\n";
    
    log_message("INFO", "Starting independent sampling track", "progressive_masking=enabled");
    
    # Check if sampling track is already completed
    if (-f "sampling_track_completed.ok") {
        log_message("INFO", "Sampling track already completed", "checkpoint found: sampling_track_completed.ok");
        chdir $current_dir;
        return 1;
    }
    
    # Initialize cumulative exclusion mask (will be created in first round)
    my $exclusion_mask;
    
    # Initialize cumulative consensi list for progressive masking
    my @cumulative_consensi;
    push @cumulative_consensi, $masked_consensi if -s $masked_consensi;
    
    # Also include RepeatScout high-quality consensus if available
    my $repeatscout_consensus = find_repeatscout_consensus($output_dir);
    push @cumulative_consensi, $repeatscout_consensus if -s $repeatscout_consensus;
    
    # Enhanced sampling parameters - starting from 90MB with intelligent stopping
    my $round = 1;
    my @sample_sizes = (90, 270, 540, 1200);  # MB - progressive scaling
    my $accumulated_mask;
    my %previous_metrics;
    my %cumulative_families;  # Track families across rounds
    
    for my $sample_size_mb (@sample_sizes) {
        last if $round > 4;  # Max 4 rounds
        
        my $round_dir = "round_$round";
        
        # Check if this round is already completed
        if (-f "$round_dir/round_completed.ok" && -s "$round_dir/consensi.fa") {
            log_message("INFO", "Round $round already completed", "checkpoint found: $round_dir/round_completed.ok");
            
            # Add this round's consensi to cumulative list
            push @cumulative_consensi, abs_path("$round_dir/consensi.fa");
            
            # Load metrics for stopping criteria evaluation
            my %current_metrics = analyze_round_metrics($round, $round_dir, \%previous_metrics, 
                                                       $sample_size_mb, 'masked');
            %previous_metrics = %current_metrics;
            
            # Update accumulated mask if sampled regions exist
            if (-s "$round_dir/sampled_regions.bed") {
                update_accumulated_mask(\$accumulated_mask, "$round_dir/sampled_regions.bed");
            }
            
            $round++;
            next;
        }
        
        log_message("INFO", "Sampling Round $round", "target_size=${sample_size_mb}MB");
        make_path($round_dir);
        
        # Change to round directory for all work
        chdir $round_dir or die "Cannot change to $round_dir: $!\n";
        
        # Create initial exclusion mask from BED files (in round directory)
        if ($round == 1) {
            $exclusion_mask = "initial_exclusion.bed";
            if (@bed_files) {
                my $bed_list = join(" ", @bed_files);
                run_cmd("cat $bed_list | bedtools sort | bedtools merge > $exclusion_mask");
            } else {
                open(my $fh, '>', $exclusion_mask);
                close($fh);
            }
            $accumulated_mask = $exclusion_mask;
        } else {
            # For subsequent rounds, copy accumulated mask to this round
            $exclusion_mask = "accumulated_mask.bed";
            if (defined $accumulated_mask && -f "../$accumulated_mask") {
                run_cmd("cp ../$accumulated_mask $exclusion_mask");
            } else {
                # Fallback to empty mask
                open(my $fh, '>', $exclusion_mask);
                close($fh);
            }
        }
        
        # Check available genome space
        my $available_size = calculate_available_size($genome_file, $exclusion_mask);
        my $required_size = $sample_size_mb * 1000000;
        
        if ($available_size < $required_size * 1.5) {
            log_message("WARN", "Insufficient genome space for sampling", 
                       "available=" . format_size($available_size) . ", required=" . format_size($required_size));
            log_message("INFO", "Falling back to full masked track processing", "switching to complete genome analysis");
            
            # Run full masked track processing as fallback
            my $success = run_masked_track_independent(
                $genome_file,
                ".",  # Current round directory
                $threads,
                $cpu_threads,
                calculate_genome_size($genome_file)
            );
            
            if ($success) {
                # Mark this round as completed using masked track approach
                open(my $fh, '>', "round_completed.ok") or die "Cannot create round completion marker: $!\n";
                print $fh "Round $round completed at " . localtime() . "\n";
                print $fh "Method: Full masked track (sampling skipped due to insufficient space)\n";
                print $fh "Available space: " . format_size($available_size) . "\n";
                print $fh "Required space: " . format_size($required_size) . "\n";
                close $fh;
                
                log_message("INFO", "Round $round completed with masked track fallback", 
                           "consensi generated successfully");
            } else {
                log_message("ERROR", "Masked track fallback also failed", 
                           "round $round could not be completed");
            }
            
            chdir "..";
            last;
        }
        
        # Step 1: Sample genome (with checkpoint)
        my $sample_file = "sample.fa";
        if (!-f "sampling.ok") {
            log_message("INFO", "Starting genome sampling", "target_size=${sample_size_mb}MB");
            perform_adaptive_sampling($genome_file, $exclusion_mask, $accumulated_mask, 
                                     $sample_file, $sample_size_mb, RECON::Utils::FRAGMENT_SIZE);
            
            # Create sampling checkpoint
            open(my $fh, '>', "sampling.ok") or die "Cannot create sampling checkpoint: $!\n";
            print $fh "Sampling completed at " . localtime() . "\n";
            print $fh "Sample size: " . (-s $sample_file || 0) . " bytes\n";
            close $fh;
            log_message("INFO", "Sampling checkpoint created", "file=sampling.ok");
        } else {
            log_message("INFO", "Sampling already completed", "checkpoint found: sampling.ok");
        }
        
        # Step 2: Enhanced progressive masking with comprehensive library integration (with checkpoint)
        my $masked_sample = "sample_masked.fa";
        if (!-f "pre-mask.ok") {
            log_message("INFO", "Starting enhanced progressive masking", 
                       "total_libraries=" . scalar(@cumulative_consensi));
            
            # Categorize and log all consensus libraries being used
            my @library_categories;
            for my $i (0..$#cumulative_consensi) {
                my $consensi = $cumulative_consensi[$i];
                my $size = -s $consensi || 0;
                my $seq_count = `grep -c '^>' '$consensi' 2>/dev/null || echo 0`;
                chomp $seq_count;
                
                my ($type, $priority);
                if ($consensi =~ /masked_track/) {
                    $type = "RECON_masked_track";
                    $priority = 1;  # High quality from masked genome
                } elsif ($consensi =~ /high_quality_consensus\.fasta/) {
                    $type = "RepeatScout_enhanced_refiner";
                    $priority = 2;  # Enhanced Refiner.py results
                } elsif ($consensi =~ /RepeatScout/) {
                    $type = "RepeatScout_standard";
                    $priority = 3;  # Standard RepeatScout
                } elsif ($consensi =~ /round_\d+/) {
                    $type = "RECON_sampling_round";
                    $priority = 4;  # Previous sampling rounds
                } else {
                    $type = "unknown";
                    $priority = 9;  # Lowest priority
                }
                
                push @library_categories, {
                    file => $consensi,
                    type => $type,
                    priority => $priority,
                    size => $size,
                    sequences => $seq_count,
                    index => $i
                };
                
                log_message("INFO", "Consensus library [$i]", 
                           "type=$type, priority=$priority, file=$consensi, " .
                           "size=$size bytes, sequences=$seq_count");
            }
            
            # Sort libraries by priority for optimal masking order
            @library_categories = sort { $a->{priority} <=> $b->{priority} } @library_categories;
            
            my $rm_threads = int($threads / 2) || 1;  # Use threads/2 for RepeatMasker
            
            # Initial TRF masking for tandem repeats
            run_trf_masking($sample_file, $masked_sample, $rm_threads);
            log_message("INFO", "TRF masking completed", "file=$masked_sample");
            
            # Merge all consensus libraries into a single file for RepeatMasker
            my $merged_library = "merged_consensus.fa";
            my $total_sequences = 0;
            my $total_masked_regions = 0;
            
            if (@library_categories > 0) {
                log_message("INFO", "Merging consensus libraries", 
                           "total_libraries=" . scalar(@library_categories) . 
                           ", output=$merged_library");
                
                open(my $merged_fh, '>', $merged_library) or die "Cannot create merged library: $!\n";
                
                my $lib_counter = 0;
                for my $lib_info (@library_categories) {
                    my $consensi = $lib_info->{file};
                    next unless -s $consensi;
                    
                    $lib_counter++;
                    log_message("INFO", "Adding library $lib_counter to merged file", 
                               "type=$lib_info->{type}, sequences=$lib_info->{sequences}, " .
                               "source=$consensi");
                    
                    # Add sequences with modified headers to track source
                    open(my $lib_fh, '<', $consensi) or die "Cannot read $consensi: $!\n";
                    while (my $line = <$lib_fh>) {
                        if ($line =~ /^>/) {
                            chomp $line;
                            # Add source information to header
                            $line .= "_src_$lib_info->{type}_lib${lib_counter}";
                            print $merged_fh "$line\n";
                            $total_sequences++;
                        } else {
                            print $merged_fh $line;
                        }
                    }
                    close($lib_fh);
                }
                close($merged_fh);
                
                log_message("INFO", "Library merging completed", 
                           "merged_library=$merged_library, total_sequences=$total_sequences, " .
                           "file_size=" . (-s $merged_library) . " bytes");
                
                # Run RepeatMasker once with the merged library
                if (-s $merged_library) {
                    log_message("INFO", "Running RepeatMasker with merged library", 
                               "sequences=$total_sequences, threads=$rm_threads");
                    
                    my $temp_out = "${masked_sample}.merged_masked";
                    my $masking_success = run_repeatmasker($masked_sample, $merged_library, $temp_out, $rm_threads);
                    
                    if ($masking_success && -s $temp_out) {
                        # Count newly masked regions
                        my $prev_n_count = `grep -o 'N' '$masked_sample' | wc -l` || 0;
                        chomp $prev_n_count;
                        my $new_n_count = `grep -o 'N' '$temp_out' | wc -l` || 0;
                        chomp $new_n_count;
                        $total_masked_regions = $new_n_count - $prev_n_count;
                        
                        run_cmd("mv $temp_out $masked_sample");
                        
                        log_message("INFO", "RepeatMasker with merged library completed", 
                                   "newly_masked=${total_masked_regions}bp, " .
                                   "libraries_applied=" . scalar(@library_categories));
                    } else {
                        log_message("WARN", "RepeatMasker with merged library failed", 
                                   "library=$merged_library, temp_output=$temp_out");
                        unlink($temp_out) if -f $temp_out;
                    }
                } else {
                    log_message("WARN", "Merged library is empty or invalid", "file=$merged_library");
                }
            } else {
                log_message("INFO", "No consensus libraries available for masking", "skipping RepeatMasker");
            }
            
            # Create comprehensive pre-masking checkpoint
            open(my $fh, '>', "pre-mask.ok") or die "Cannot create pre-mask checkpoint: $!\n";
            print $fh "Enhanced pre-masking with merged libraries completed at " . localtime() . "\n";
            print $fh "Masked file size: " . (-s $masked_sample || 0) . " bytes\n";
            print $fh "Total libraries merged: " . scalar(@cumulative_consensi) . "\n";
            print $fh "Total sequences in merged library: $total_sequences\n";
            print $fh "Merged library file: $merged_library\n";
            print $fh "Total masked regions: ${total_masked_regions} bp\n";
            print $fh "RepeatMasker runs: 1 (merged library approach)\n";
            
            # Record detailed library information with priority order
            print $fh "\n=== Library Merge Order (by priority) ===\n";
            for my $lib_info (@library_categories) {
                print $fh sprintf("Priority %d: %s\n", $lib_info->{priority}, $lib_info->{type});
                print $fh sprintf("  File: %s (%d bytes, %d sequences)\n", 
                                 $lib_info->{file}, $lib_info->{size}, $lib_info->{sequences});
            }
            
            print $fh "\n=== Integration Summary ===\n";
            my $has_repeatscout = scalar(grep { $_->{type} eq "RepeatScout_enhanced_refiner" } @library_categories);
            my $has_recon = scalar(grep { $_->{type} eq "RECON_masked_track" } @library_categories);
            my $previous_rounds = scalar(grep { $_->{type} eq "RECON_sampling_round" } @library_categories);
            
            print $fh "RepeatScout enhanced consensus: " . ($has_repeatscout ? "YES" : "NO") . "\n";
            print $fh "RECON masked track: " . ($has_recon ? "YES" : "NO") . "\n";
            print $fh "Previous rounds: $previous_rounds\n";
            print $fh "Optimization: Single RepeatMasker run with merged library (improved efficiency)\n";
            
            close $fh;
            log_message("INFO", "Enhanced pre-mask checkpoint created", 
                       "file=pre-mask.ok, merged_sequences=$total_sequences, total_masked=${total_masked_regions}bp");
        } else {
            log_message("INFO", "Pre-masking already completed", "checkpoint found: pre-mask.ok");
        }
        
        # Step 2.5: Extract unmasked regions and filter short fragments (with checkpoint)
        my $unmasked_sample = "sample_unmasked_filtered.fa";
        if (!-f "unmasked_extract.ok") {
            log_message("INFO", "Starting unmasked region extraction and filtering", 
                       "min_length=50bp, input=$masked_sample");
            
            my ($extracted_count, $filtered_count) = extract_unmasked_regions_with_gi_format(
                $masked_sample, $unmasked_sample
            );
            
            # Create unmasked extraction checkpoint
            open(my $fh, '>', "unmasked_extract.ok") or die "Cannot create unmasked extraction checkpoint: $!\n";
            print $fh "Unmasked region extraction completed at " . localtime() . "\n";
            print $fh "Input file: $masked_sample\n";
            print $fh "Output file: $unmasked_sample\n";
            print $fh "Total extracted fragments: $extracted_count\n";
            print $fh "Fragments after filtering (≥50bp): $filtered_count\n";
            print $fh "Filtered out (short fragments): " . ($extracted_count - $filtered_count) . "\n";
            print $fh "Final file size: " . (-s $unmasked_sample || 0) . " bytes\n";
            close $fh;
            
            log_message("INFO", "Unmasked extraction checkpoint created", 
                       "extracted=$extracted_count, filtered=$filtered_count, file=unmasked_extract.ok");
        } else {
            log_message("INFO", "Unmasked region extraction already completed", "checkpoint found: unmasked_extract.ok");
        }
        
        # Update the sample file to use unmasked, filtered sequences for downstream processing
        $masked_sample = $unmasked_sample;
        
        # Step 3: Sequence renaming (with checkpoint)
        my $renamed_sample = "sample_renamed.fa";
        if (!-f "renaming.ok") {
            log_message("INFO", "Starting sequence renaming", "format=gi|N");
            my $seq_count = rename_fasta_ids($masked_sample, $renamed_sample, "gi");
            
            # Create renaming checkpoint
            open(my $fh, '>', "renaming.ok") or die "Cannot create renaming checkpoint: $!\n";
            print $fh "Renaming completed at " . localtime() . "\n";
            print $fh "Sequences renamed: $seq_count\n";
            close $fh;
            log_message("INFO", "Renaming checkpoint created", "sequences=$seq_count, file=renaming.ok");
        } else {
            log_message("INFO", "Sequence renaming already completed", "checkpoint found: renaming.ok");
        }
        $masked_sample = $renamed_sample;  # Use renamed file for downstream processing
        
        # Step 4: RMBlastN self-alignment with direct MSP output (with checkpoint)
        if (!-f "rmblastn.ok") {
            log_message("INFO", "Starting RMBlastN self-alignment", "database=sample_db, direct_msp=true");
            create_blast_database($masked_sample, "sample_db");
            
            # Call rmblastn with direct MSP output (msp.out)
            run_rmblastn_self_alignment($masked_sample, "msp.out", $threads, "sample_db");
            
            my $msp_count = -s "msp.out" ? `wc -l < msp.out` : 0;
            chomp $msp_count;
            
            # Create RMBlastN checkpoint (now includes MSP generation)
            open(my $fh, '>', "rmblastn.ok") or die "Cannot create RMBlastN checkpoint: $!\n";
            print $fh "RMBlastN with MSP generation completed at " . localtime() . "\n";
            print $fh "MSP output size: " . (-s "msp.out" || 0) . " bytes\n";
            print $fh "MSP lines: $msp_count\n";
            close $fh;
            log_message("INFO", "RMBlastN checkpoint created", "msp_lines=$msp_count, file=rmblastn.ok");
        } else {
            log_message("INFO", "RMBlastN with MSP already completed", "checkpoint found: rmblastn.ok");
        }
        
        if (!-s "msp.out") {
            log_message("WARN", "No repeats found in round $round");
            chdir "..";
            $round++;
            next;
        }
        
        # Step 6: RECON pipeline preparation (with checkpoint)
        if (!-f "recon_prep.ok") {
            log_message("INFO", "Starting RECON preparation", "creating sequence list and K-parameter");
            create_seq_name_list($masked_sample);
            run_cmd("cp msp.out msp.out.backup");  # Keep original
            my $k_param = determine_k_parameter("msp.out");
            
            # Create preparation checkpoint
            open(my $fh, '>', "recon_prep.ok") or die "Cannot create RECON prep checkpoint: $!\n";
            print $fh "RECON preparation completed at " . localtime() . "\n";
            print $fh "K parameter: $k_param\n";
            print $fh "Sequence list: seq.names\n";
            close $fh;
            log_message("INFO", "RECON prep checkpoint created", "k_param=$k_param, file=recon_prep.ok");
        } else {
            log_message("INFO", "RECON preparation already completed", "checkpoint found: recon_prep.ok");
        }
        
        # Step 7: RECON pipeline execution (with checkpoint)
        if (!-f "recon.ok") {
            log_message("INFO", "Starting RECON pipeline execution");
            my $k_param = determine_k_parameter("msp.out");  # Re-determine in case of restart
            
            eval {
                run_recon_pipeline($k_param);
                
                # Create RECON execution checkpoint
                open(my $fh, '>', "recon.ok") or die "Cannot create RECON checkpoint: $!\n";
                print $fh "RECON pipeline completed at " . localtime() . "\n";
                print $fh "K parameter used: $k_param\n";
                if (-d "summary") {
                    my $families = -s "summary/families" ? `wc -l < summary/families` : 0;
                    chomp $families;
                    print $fh "Families found: $families\n";
                }
                close $fh;
                log_message("INFO", "RECON pipeline checkpoint created", "file=recon.ok");
            };
            
            if ($@) {
                log_message("WARN", "RECON failed in round $round", "error=$@");
                
                # Log debugging information
                log_debugging_info($round, $round_dir);
                log_message("INFO", "Preserving failed round data for debugging", "round_dir=$round_dir");
                
                # DO NOT clean up failed data - preserve for debugging
                # Just remove completion checkpoints to prevent false positive detection
                unlink("recon.ok") if -f "recon.ok";
                unlink("consensus.ok") if -f "consensus.ok"; 
                unlink("round_completed.ok") if -f "round_completed.ok";
                
                chdir "..";
                $round++;
                next;
            }
        } else {
            log_message("INFO", "RECON pipeline already completed", "checkpoint found: recon.ok");
        }
        
        # Step 8: Consensus building (with checkpoint)
        if (-d "summary" && -s "summary/families") {
            if (!-f "consensus.ok") {
                log_message("INFO", "Starting consensus building", "families_found=true");
                my $local_genome = abs_path($masked_sample);
                run_cmd("build_for_RECON ./ $local_genome $cpu_threads");
                
                if (-s "consensi.fa") {
                    my $consensi_count = `grep -c '^>' consensi.fa`;
                    chomp $consensi_count;
                    
                    # Create consensus checkpoint
                    open(my $fh, '>', "consensus.ok") or die "Cannot create consensus checkpoint: $!\n";
                    print $fh "Consensus building completed at " . localtime() . "\n";
                    print $fh "Consensus sequences: $consensi_count\n";
                    print $fh "File size: " . (-s "consensi.fa") . " bytes\n";
                    close $fh;
                    
                    push @cumulative_consensi, abs_path("consensi.fa");
                    log_message("INFO", "Consensus checkpoint created", "sequences=$consensi_count, file=consensus.ok");
                } else {
                    log_message("WARN", "No consensus sequences generated despite families found");
                }
            } else {
                log_message("INFO", "Consensus building already completed", "checkpoint found: consensus.ok");
                # Still need to add to cumulative list if not already there
                if (-s "consensi.fa") {
                    my $consensi_path = abs_path("consensi.fa");
                    push @cumulative_consensi, $consensi_path unless grep { $_ eq $consensi_path } @cumulative_consensi;
                }
            }
            
            # Create final round completion marker only if consensus exists
            if (-s "consensi.fa" && !-f "round_completed.ok") {
                open(my $round_fh, '>', "round_completed.ok") or die "Cannot create round completion marker: $!\n";
                print $round_fh "Round $round completed at " . localtime() . "\n";
                my $consensi_count = `grep -c '^>' consensi.fa`;
                chomp $consensi_count;
                print $round_fh "Consensus sequences: $consensi_count\n";
                print $round_fh "All checkpoints: sampling.ok, pre-mask.ok, renaming.ok, rmblastn.ok, recon_prep.ok, recon.ok, consensus.ok\n";
                close $round_fh;
                log_message("INFO", "Round $round completion marker created", "file=round_completed.ok");
            }
        } else {
            log_message("INFO", "No families found, skipping consensus building", "summary_missing=true");
        }
        
        # Generate sampled regions BED file for this round (in round directory)
        # This will be used to update the accumulated mask for next round
        my $sampled_regions_file = "sampled_regions.bed";
        # TODO: This should extract sampled regions from the sampling step
        # For now, create empty file as placeholder
        unless (-f $sampled_regions_file) {
            open(my $fh, '>', $sampled_regions_file);
            close($fh);
        }
        
        # Return to sampling_track directory
        chdir "..";
        
        # Calculate comprehensive metrics for stopping decision
        my %current_metrics = analyze_round_metrics($round, $round_dir, \%previous_metrics, 
                                                   $sample_size_mb, 'masked');
        
        # Evaluate stopping criteria using new intelligent system
        my ($should_continue, $reason) = evaluate_stopping_criteria($round, \%current_metrics, 
                                                                   \%previous_metrics, undef);
        
        log_message("INFO", "Round $round sampling decision", 
                    "decision=" . ($should_continue ? "CONTINUE" : "STOP") . ", reason=$reason");
        
        # Update metrics for next round comparison
        %previous_metrics = %current_metrics;
        
        # Stop if criteria met (intelligent stopping system)
        if (!$should_continue) {
            log_message("INFO", "Intelligent stopping criteria triggered", 
                       "round=$round, sample_size=${sample_size_mb}MB, reason=$reason");
            last;
        }
        if (-s "$round_dir/$sampled_regions_file") {
            update_accumulated_mask(\$accumulated_mask, "$round_dir/$sampled_regions_file");
        }
        $round++;
    }
    
    # Create completion marker for checkpoint system
    open(my $fh, '>', "sampling_track_completed.ok") or die "Cannot create completion marker: $!\n";
    print $fh "Completed at " . localtime() . "\n";
    print $fh "Rounds completed: " . ($round - 1) . "\n";
    close $fh;
    
    chdir $current_dir;
    
    log_message("INFO", "Sampling track completed", "rounds=" . ($round - 1));
}

sub run_sampling_track_progressive {
    # Legacy compatibility - redirect to independent track
    return run_sampling_track_independent(@_);
}

sub run_sampling_track {
    # Legacy compatibility - redirect to independent track
    return run_sampling_track_independent(@_);
}

sub analyze_round_metrics {
    my ($round, $round_dir, $prev_metrics_ref, $sample_mb, $track_type) = @_;
    my %metrics;
    
    $metrics{sample_mb} = $sample_mb;
    $metrics{track_type} = $track_type || 'masked';  # 'masked' or 'raw'
    
    # 1. Count families (recon_families_total)
    if (-f "$round_dir/summary/families") {
        $metrics{recon_families_total} = `wc -l < $round_dir/summary/families`;
        chomp $metrics{recon_families_total};
        
        # 2. Count singletons and calculate family size distribution
        if (-f "$round_dir/summary/eles") {
            open(my $fh, '<', "$round_dir/summary/eles") or return %metrics;
            my %family_counts;
            while (<$fh>) {
                next if /^#/;
                my @fields = split;
                next if @fields < 1;
                $family_counts{$fields[0]}++;
            }
            close($fh);
            
            # 5. Singletons count (copy_number = 1)
            $metrics{singletons} = scalar(grep { $family_counts{$_} == 1 } keys %family_counts);
            
            # 6. P95 family size (95th percentile of copy numbers)
            my @sizes = sort { $b <=> $a } values %family_counts;
            if (@sizes > 0) {
                my $p95_index = int(@sizes * 0.05);
                $metrics{p95_family_size} = $sizes[$p95_index] || $sizes[-1];
            } else {
                $metrics{p95_family_size} = 0;
            }
        }
    }
    
    # 2. MSP lines (HSP/edge count)
    if (-f "$round_dir/msp.out") {
        $metrics{msp_lines} = `wc -l < $round_dir/msp.out`;
        chomp $metrics{msp_lines};
    } else {
        $metrics{msp_lines} = 0;
    }
    
    # 3. MSP density (lines per Gb of sampled sequence)
    if ($sample_mb > 0) {
        $metrics{msp_density_lines_per_GbSeq} = ($metrics{msp_lines} || 0) / ($sample_mb / 1024);
    } else {
        $metrics{msp_density_lines_per_GbSeq} = 0;
    }
    
    # Calculate incremental metrics vs previous rounds
    if ($prev_metrics_ref && %$prev_metrics_ref) {
        # 7. New families (真正新增的家族数)
        $metrics{new_families} = ($metrics{recon_families_total} || 0) - 
                                 ($prev_metrics_ref->{recon_families_total} || 0);
        
        # 8. New families per 100Mb sampling efficiency
        $metrics{new_families_per_100Mb} = $sample_mb > 0 ? 
            100 * ($metrics{new_families} || 0) / $sample_mb : 0;
            
        # Singleton ratio (for stopping criteria)
        $metrics{singleton_ratio} = $metrics{recon_families_total} > 0 ?
            ($metrics{singletons} || 0) / $metrics{recon_families_total} : 0;
            
        # P95 growth rate (for large family detection)
        if (($prev_metrics_ref->{p95_family_size} || 0) > 0) {
            $metrics{p95_growth_rate} = (($metrics{p95_family_size} || 0) - 
                                       ($prev_metrics_ref->{p95_family_size} || 0)) / 
                                       $prev_metrics_ref->{p95_family_size};
        } else {
            $metrics{p95_growth_rate} = 0;
        }
        
        # MSP density growth rate (for graph enrichment detection)
        if (($prev_metrics_ref->{msp_density_lines_per_GbSeq} || 0) > 0) {
            $metrics{msp_density_growth_rate} = (($metrics{msp_density_lines_per_GbSeq} || 0) - 
                                               ($prev_metrics_ref->{msp_density_lines_per_GbSeq} || 0)) / 
                                               $prev_metrics_ref->{msp_density_lines_per_GbSeq};
        } else {
            $metrics{msp_density_growth_rate} = 0;
        }
    } else {
        # First round - no comparison available
        $metrics{new_families} = $metrics{recon_families_total} || 0;
        $metrics{new_families_per_100Mb} = $sample_mb > 0 ? 
            100 * ($metrics{new_families} || 0) / $sample_mb : 0;
        $metrics{singleton_ratio} = $metrics{recon_families_total} > 0 ?
            ($metrics{singletons} || 0) / $metrics{recon_families_total} : 0;
        $metrics{p95_growth_rate} = 0;
        $metrics{msp_density_growth_rate} = 0;
    }
    
    # Log comprehensive metrics
    log_message("INFO", "Round $round comprehensive metrics", 
                sprintf("sample_mb=%.1f, families=%d, singletons=%d (%.1f%%), " .
                        "new_families=%d (%.1f/100Mb), p95_size=%d, msp_lines=%d, msp_density=%.1f",
                        $metrics{sample_mb} || 0,
                        $metrics{recon_families_total} || 0,
                        $metrics{singletons} || 0,
                        ($metrics{singleton_ratio} || 0) * 100,
                        $metrics{new_families} || 0,
                        $metrics{new_families_per_100Mb} || 0,
                        $metrics{p95_family_size} || 0,
                        $metrics{msp_lines} || 0,
                        $metrics{msp_density_lines_per_GbSeq} || 0));
    
    return %metrics;
}

sub evaluate_stopping_criteria {
    my ($round, $current_metrics, $prev_metrics, $round_history) = @_;
    
    my $track_type = $current_metrics->{track_type} || 'masked';
    my $continue = 0;
    my $stop_reason = "";
    my @reasons;
    
    # Thresholds based on track type
    my $continue_threshold = $track_type eq 'raw' ? 20 : 30;
    my $stop_threshold = $track_type eq 'raw' ? 10 : 15;
    
    log_message("INFO", "Evaluating stopping criteria for round $round",
                sprintf("track_type=%s, sample_mb=%.1f, families=%d, new_families=%d",
                        $track_type, 
                        $current_metrics->{sample_mb} || 0,
                        $current_metrics->{recon_families_total} || 0,
                        $current_metrics->{new_families} || 0));
    
    # === CONTINUE CONDITIONS (any one satisfied) ===
    
    # 1. High new family discovery rate
    my $new_families_rate = $current_metrics->{new_families_per_100Mb} || 0;
    if ($new_families_rate >= $continue_threshold) {
        $continue = 1;
        push @reasons, sprintf("High discovery rate: %.1f families/100Mb (≥%d)", 
                              $new_families_rate, $continue_threshold);
    }
    
    # 2. P95 family size still growing significantly (≥20% increase)
    my $p95_growth = $current_metrics->{p95_growth_rate} || 0;
    if ($p95_growth >= 0.20) {
        $continue = 1;
        push @reasons, sprintf("P95 family size growing: +%.1f%% (≥+20%%)",
                              $p95_growth * 100);
    }
    
    # 3. MSP density still increasing significantly (≥15% increase) 
    my $msp_growth = $current_metrics->{msp_density_growth_rate} || 0;
    if ($msp_growth >= 0.15) {
        $continue = 1;
        push @reasons, sprintf("MSP density growing: +%.1f%% (≥+15%%)",
                              $msp_growth * 100);
    }
    
    # === STOP CONDITIONS (both must be satisfied for confident stopping) ===
    
    my $should_stop = 0;
    my @stop_conditions;
    
    # Stop condition 1: Low new family rate for current AND previous round
    my $current_low_rate = $new_families_rate < $stop_threshold;
    my $prev_low_rate = 0;
    
    if ($prev_metrics && exists $prev_metrics->{new_families_per_100Mb}) {
        $prev_low_rate = ($prev_metrics->{new_families_per_100Mb} || 0) < $stop_threshold;
    }
    
    my $consecutive_low_rate = $current_low_rate && $prev_low_rate;
    
    if ($consecutive_low_rate) {
        push @stop_conditions, sprintf("Low discovery rate for 2 rounds: %.1f < %d",
                                     $new_families_rate, $stop_threshold);
    }
    
    # Stop condition 2: High singleton ratio (≥60%) AND stable P95 (±5%)
    my $singleton_ratio = $current_metrics->{singleton_ratio} || 0;
    my $high_singleton_ratio = $singleton_ratio >= 0.60;
    
    my $stable_p95 = abs($p95_growth) <= 0.05;  # ±5%
    
    if ($high_singleton_ratio && $stable_p95) {
        push @stop_conditions, sprintf("High singletons (%.1f%% ≥60%%) + stable P95 (%.1f%% ±5%%)",
                                     $singleton_ratio * 100, $p95_growth * 100);
    }
    
    # Decide: STOP only if both stop conditions met AND no continue conditions
    if (@stop_conditions >= 2 && !$continue) {
        $should_stop = 1;
        $stop_reason = "Both stop conditions met: " . join("; ", @stop_conditions);
    }
    
    # Override: Always continue if this is round 1 (90MB)
    if ($round == 1) {
        $continue = 1;
        $should_stop = 0;
        $stop_reason = "First round - always continue to 270MB";
    }
    
    # Log detailed metrics
    log_message("INFO", "Round $round detailed metrics",
                sprintf("msp_lines=%d, msp_density=%.1f lines/Gb, " .
                        "new_families=%d (%.1f/100Mb), singletons=%d (%.1f%%), " .
                        "p95_size=%d (+%.1f%%), msp_growth=+%.1f%%",
                        $current_metrics->{msp_lines} || 0,
                        $current_metrics->{msp_density_lines_per_GbSeq} || 0,
                        $current_metrics->{new_families} || 0,
                        $new_families_rate,
                        $current_metrics->{singletons} || 0,
                        $singleton_ratio * 100,
                        $current_metrics->{p95_family_size} || 0,
                        $p95_growth * 100,
                        $msp_growth * 100));
    
    # Final decision
    my $decision = $should_stop ? "STOP" : "CONTINUE";
    my $reason = $should_stop ? $stop_reason : 
                 (@reasons ? join("; ", @reasons) : "Default continue");
    
    log_message("INFO", "Sampling decision for round $round", 
                "$decision - $reason");
    
    return (!$should_stop, $reason);
}

sub update_accumulated_mask {
    my ($mask_ref, $new_regions_file) = @_;
    
    if (!$$mask_ref) {
        $$mask_ref = $new_regions_file;
        return;
    }
    
    # Merge BED files (simple concatenation for now)
    my $temp_file = "temp_mask.bed";
    run_cmd("cat $$mask_ref $new_regions_file > $temp_file");
    run_cmd("bedtools sort -i $temp_file | bedtools merge > merged_mask.bed");
    run_cmd("mv merged_mask.bed $$mask_ref");
    unlink($temp_file);
}

sub cleanup_failed_round {
    my ($round_dir) = @_;
    
    # Remove all checkpoint files for fresh restart
    my @checkpoint_files = qw(
        sampling.ok pre-mask.ok renaming.ok rmblastn.ok 
        recon_prep.ok recon.ok consensus.ok 
        round_completed.ok
    );
    
    for my $checkpoint (@checkpoint_files) {
        if (-f $checkpoint) {
            unlink($checkpoint);
            log_message("DEBUG", "Removed checkpoint", "file=$checkpoint");
        }
    }
    
    # Remove potentially corrupted output files but keep logs for debugging
    my @cleanup_files = qw(
        consensi.fa msp.out msp.out.backup
        sample_renamed.fa sample_masked.fa sample.fa
        seq.names
        summary/families summary/eles
    );
    
    for my $file (@cleanup_files) {
        if (-f $file) {
            unlink($file);
            log_message("DEBUG", "Removed output file", "file=$file");
        }
    }
    
    # Remove database files
    for my $db_file (glob("sample_db*")) {
        unlink($db_file);
    }
    
    # Remove summary directory if it exists and is empty/corrupted
    if (-d "summary") {
        # Remove all files in summary directory
        for my $summary_file (glob("summary/*")) {
            unlink($summary_file);
        }
        rmdir("summary");
    }
    
    log_message("INFO", "Cleaned up failed round", "round_dir=$round_dir, checkpoints_removed=" . scalar(@checkpoint_files));
}

sub log_debugging_info {
    my ($round, $round_dir) = @_;
    
    log_message("DEBUG", "Failed round debugging info", "round=$round, dir=$round_dir");
    
    # Log file sizes and existence
    my @debug_files = qw(
        seq.names msp.out sample_renamed.fa
        eledef.log edgeredef.log famdef.log
    );
    
    for my $file (@debug_files) {
        if (-f $file) {
            my $size = -s $file;
            log_message("DEBUG", "File exists", "file=$file, size=$size bytes");
        } else {
            log_message("DEBUG", "File missing", "file=$file");
        }
    }
    
    # Check if seq.names has correct format
    if (-f "seq.names") {
        my $first_line = `head -1 seq.names`;
        chomp $first_line;
        if ($first_line =~ /^\d+$/) {
            log_message("DEBUG", "seq.names format correct", "count=$first_line");
        } else {
            log_message("DEBUG", "seq.names format incorrect", "first_line=$first_line");
        }
    }
    
    # Log MSP file statistics
    if (-f "msp.out") {
        my $msp_lines = `wc -l < msp.out`;
        chomp $msp_lines;
        log_message("DEBUG", "MSP statistics", "lines=$msp_lines");
        
        # Show first few MSP lines
        if ($msp_lines > 0) {
            my $sample_msp = `head -3 msp.out`;
            chomp $sample_msp;
            $sample_msp =~ s/\n/ | /g;
            log_message("DEBUG", "MSP sample", "first_3_lines=$sample_msp");
        }
    }
}

sub find_repeatscout_consensus {
    my ($base_dir) = @_;
    
    # Get current working directory to understand context
    my $current_wd = getcwd();
    
    # Look for RepeatScout high-quality consensus file in various possible locations
    # The challenge is that we might be called from different contexts:
    # 1. From sampling_track directory: need ../RepeatScout/...
    # 2. From sampling_track/round_1: need ../../../RepeatScout/...
    
    my @possible_paths = (
        # From current working directory (most reliable approach)
        "../../../RepeatScout/refiner_output/high_quality_consensus.fasta",  # From round_1
        "../../RepeatScout/refiner_output/high_quality_consensus.fasta",     # From sampling_track
        "../RepeatScout/refiner_output/high_quality_consensus.fasta",        # From RECON dir
        "../../../../RepeatScout/refiner_output/high_quality_consensus.fasta", # Extra level up
        
        # Using base_dir parameter (fallback)
        "$base_dir/../RepeatScout/refiner_output/high_quality_consensus.fasta",
        "$base_dir/../../RepeatScout/refiner_output/high_quality_consensus.fasta",
        "$base_dir/../../../RepeatScout/refiner_output/high_quality_consensus.fasta",
        
        # Enhanced Refiner output paths (with improved consensus building)
        "../../../RepeatScout/refiner_output/consensus_masking.fa",
        "../../RepeatScout/refiner_output/consensus_masking.fa",
        "../RepeatScout/refiner_output/consensus_masking.fa",
        "$base_dir/../RepeatScout/refiner_output/consensus_masking.fa",
        "$base_dir/../../RepeatScout/refiner_output/consensus_masking.fa",
        "$base_dir/../../../RepeatScout/refiner_output/consensus_masking.fa", 
        
        # Fallback to standard RepeatScout output
        "../../../RepeatScout/consensi.fa",
        "../../RepeatScout/consensi.fa", 
        "../RepeatScout/consensi.fa",
        "$base_dir/../RepeatScout/consensi.fa",
        "$base_dir/../../RepeatScout/consensi.fa",
        "$base_dir/../../../RepeatScout/consensi.fa",
    );
    
    # Debug: log current context
    log_message("DEBUG", "RepeatScout consensus search context", 
               "current_wd=$current_wd, base_dir=$base_dir");
    
    for my $path (@possible_paths) {
        log_message("DEBUG", "Checking RepeatScout path", "path=$path, exists=" . (-f $path ? "yes" : "no"));
        
        if (-f $path && -s $path) {
            # Validate that this is a properly formatted consensus file
            my $seq_count = `grep -c '^>' '$path' 2>/dev/null || echo 0`;
            chomp $seq_count;
            
            if ($seq_count > 0) {
                log_message("INFO", "Found RepeatScout consensus", 
                           "file=$path, sequences=$seq_count, absolute=" . abs_path($path));
                return abs_path($path);
            } else {
                log_message("WARN", "Empty RepeatScout consensus file found", "file=$path");
            }
        }
    }
    
    log_message("WARN", "RepeatScout consensus not found", 
               "current_wd=$current_wd, searched paths: " . join(", ", @possible_paths));
    return "";
}

sub extract_unmasked_regions_with_gi_format {
    my ($input_file, $output_file) = @_;
    
    log_message("INFO", "Extracting unmasked regions with gi|N formatting", 
               "input=$input_file, output=$output_file, min_length=50bp");
    
    # Create temporary file for MaskedTrack extraction
    my $temp_file = "${output_file}.tmp";
    
    # Use MaskedTrack's extract_unmasked_regions function (now 50bp minimum, N-only splitting)
    my $extracted_regions = extract_unmasked_regions($input_file, $temp_file);
    
    # Reformat the output to use simple gi|N format
    open(my $temp_fh, '<', $temp_file) or die "Cannot open temp file $temp_file: $!\n";
    open(my $out_fh, '>', $output_file) or die "Cannot create output file $output_file: $!\n";
    
    my $sequence_id = 0;
    my $total_sequences = 0;
    my $current_sequence = "";
    
    while (my $line = <$temp_fh>) {
        chomp $line;
        
        if ($line =~ /^>/) {
            # Process previous sequence if exists
            if ($current_sequence) {
                $total_sequences++;
                $sequence_id++;
                
                # Output with simple gi|N format as requested
                print $out_fh ">gi|$sequence_id\n";
                
                # Format sequence with 80 characters per line
                for (my $i = 0; $i < length($current_sequence); $i += 80) {
                    print $out_fh substr($current_sequence, $i, 80) . "\n";
                }
            }
            $current_sequence = "";
        } else {
            $current_sequence .= $line;
        }
    }
    
    # Process the last sequence
    if ($current_sequence) {
        $total_sequences++;
        $sequence_id++;
        
        print $out_fh ">gi|$sequence_id\n";
        for (my $i = 0; $i < length($current_sequence); $i += 80) {
            print $out_fh substr($current_sequence, $i, 80) . "\n";
        }
    }
    
    close($temp_fh);
    close($out_fh);
    
    # Clean up temporary file
    unlink($temp_file);
    
    log_message("INFO", "Unmasked region extraction with gi|N formatting completed", 
               "extracted_regions=$extracted_regions, formatted_sequences=$total_sequences");
    
    # Return the same count for both since MaskedTrack already does 50bp filtering
    return ($extracted_regions, $total_sequences);
}

1;

__END__

=head1 NAME

RECON::SamplingTrack - Enhanced sampling track implementation for RECON pipeline

=head1 DESCRIPTION

Implements an intelligent sampling track with:
- Progressive multi-round sampling (90MB → 270MB → 540MB → 1.2GB)
- Adaptive stopping criteria based on discovery metrics
- Comprehensive family discovery analysis
- MSP density and P95 family size monitoring
- Cumulative masking and progressive enrichment

=head1 SAMPLING STRATEGY

Starting from 90MB with intelligent stopping based on:
- New families per 100MB discovery rate
- P95 family size growth (large family detection)
- MSP density changes (graph enrichment)
- Singleton ratio and stability metrics

=cut
