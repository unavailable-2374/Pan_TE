#!/usr/bin/env perl
use strict;
use warnings;
use File::Path qw(make_path);
use File::Spec;
use File::Basename;
use Cwd qw(getcwd abs_path);
use List::Util qw(max);

# Sampling track with progressive masking and intelligent stopping criteria
sub run_sampling_track_independent {
    my ($genome_file, $bed_files_ref, $masked_consensi, $output_dir, 
        $threads, $cpu_threads, $genome_size) = @_;
    
    my @bed_files = @$bed_files_ref;
    my $current_dir = getcwd();
    chdir $output_dir or die "Cannot change to $output_dir: $!\n";
    
    print "Starting independent sampling track with progressive masking...\n";
    
    # Create initial exclusion mask from BED files
    my $exclusion_mask = "initial_exclusion.bed";
    if (@bed_files) {
        my $bed_list = join(" ", @bed_files);
        system("cat $bed_list | bedtools sort | bedtools merge > $exclusion_mask");
    } else {
        open(my $fh, '>', $exclusion_mask);
        close($fh);
    }
    
    # Initialize progressive masking
    my @cumulative_consensi;
    push @cumulative_consensi, $masked_consensi if -s $masked_consensi;
    
    # Sampling parameters
    my $round = 1;
    my @sample_sizes = (30, 90, 270, 540);  # MB
    my $accumulated_mask = $exclusion_mask;
    
    # Metrics tracking for stopping criteria
    my %previous_metrics;
    
    for my $sample_size_mb (@sample_sizes) {
        last if $round > 4;  # Max 4 rounds
        
        print "\n=== Sampling Round $round (${sample_size_mb}MB) ===\n";
        
        my $round_dir = "round_$round";
        make_path($round_dir);
        
        # Check available genome space
        my $available_size = calculate_available_genome_size($genome_file, $accumulated_mask);
        my $required_size = $sample_size_mb * 1000000;
        
        if ($available_size < $required_size * 1.5) {
            print "Insufficient genome space for ${sample_size_mb}MB sampling\n";
            last;
        }
        
        # Sample genome
        my $sample_file = "$round_dir/sample.fa";
        perform_sampling($genome_file, $accumulated_mask, $sample_file, $sample_size_mb);
        
        # Apply progressive masking
        my $masked_sample = "$round_dir/sample_masked.fa";
        apply_progressive_masking($sample_file, \@cumulative_consensi, $masked_sample, $cpu_threads);
        
        # Run BLAST and RECON
        my $blast_output = "$round_dir/self_alignment.blast";
        run_blast_and_recon($masked_sample, $blast_output, $round_dir, $threads, $cpu_threads);
        
        # Build consensi for this round
        if (-d "$round_dir/summary" && -s "$round_dir/summary/families") {
            chdir $round_dir;
            my $main_dir = File::Spec->rel2abs("../../../..");
            my $genome_path = File::Spec->catfile($main_dir, "genome", "genome.fa");
            system("build_for_RECON ./ $genome_path $cpu_threads");
            chdir "..";
            
            push @cumulative_consensi, "$round_dir/consensi.fa" if -s "$round_dir/consensi.fa";
        }
        
        # Calculate metrics for stopping decision
        my %current_metrics = calculate_round_metrics($round_dir, $sample_size_mb);
        
        # Stopping criteria evaluation
        if ($round >= 2 && should_stop_sampling(\%previous_metrics, \%current_metrics, $round)) {
            print "Stopping criteria met. Ending sampling at round $round\n";
            last;
        }
        
        # Update for next round
        %previous_metrics = %current_metrics;
        update_accumulated_mask(\$accumulated_mask, "$round_dir/sampled_regions.bed");
        $round++;
    }
    
    chdir $current_dir;
}

sub calculate_available_genome_size {
    my ($genome_file, $mask_bed) = @_;
    
    # Get total genome size
    my $total_size = `seqkit stats -T $genome_file | tail -1 | cut -f5`;
    chomp $total_size;
    $total_size =~ s/,//g;
    
    # Get masked size
    my $masked_size = 0;
    if (-s $mask_bed) {
        $masked_size = `bedtools merge -i $mask_bed | awk '{sum += \$3-\$2} END {print sum}'`;
        chomp $masked_size;
        $masked_size ||= 0;
    }
    
    return $total_size - $masked_size;
}

sub perform_sampling {
    my ($genome_file, $mask_bed, $output_file, $size_mb) = @_;
    
    print "Sampling ${size_mb}MB from genome...\n";
    
    # Use bedtools to sample regions not in mask
    my $target_size = $size_mb * 1000000;
    my $fragment_size = 40000;
    
    # Create genome file for bedtools
    my $genome_sizes = "genome.sizes";
    system("seqkit seq -n $genome_file | xargs -I {} sh -c 'echo -e \"{}\t\$(seqkit stats -T $genome_file | grep {} | cut -f5 | tr -d ,)\"' > $genome_sizes");
    
    # Generate random windows avoiding masked regions
    system("bedtools random -l $fragment_size -n " . int($target_size/$fragment_size) . 
           " -g $genome_sizes | bedtools subtract -a - -b $mask_bed > sampled_regions.bed");
    
    # Extract sequences
    system("bedtools getfasta -fi $genome_file -bed sampled_regions.bed -fo $output_file");
}

sub apply_progressive_masking {
    my ($input_fa, $consensi_ref, $output_fa, $threads) = @_;
    
    # First apply TRF masking
    system("trf $input_fa 2 5 7 80 10 50 2000 -m -h > /dev/null 2>&1");
    system("mv ${input_fa}.2.5.7.80.10.50.2000.mask $output_fa");
    
    # Apply each consensi file sequentially
    for my $consensi (@$consensi_ref) {
        next unless -s $consensi;
        print "Applying masking with: $consensi\n";
        my $temp_out = "${output_fa}.tmp";
        system("RepeatMasker -pa $threads -lib $consensi -xsmall -nolow -norna $output_fa > /dev/null 2>&1");
        system("mv ${output_fa}.masked $temp_out 2>/dev/null && mv $temp_out $output_fa");
    }
}

sub run_blast_and_recon {
    my ($genome_file, $blast_output, $work_dir, $threads, $cpu_threads) = @_;
    
    chdir $work_dir;
    
    # Create BLAST database
    system("makeblastdb -in $genome_file -dbtype nucl -out genome_db");
    
    # Run BLAST self-alignment
    system("blastn -query $genome_file -db genome_db -out $blast_output " .
           "-evalue 1e-5 -num_threads $threads -outfmt 6");
    
    # Convert to MSP
    system("MSPCollect.pl $blast_output > msp.out");
    
    return unless -s "msp.out";
    
    # Create sequence name list
    system("grep '^>' $genome_file | sed 's/>//' > seq_name.list");
    
    # Determine K parameter
    my $k_param = 10;  # Default
    my $msp_lines = `wc -l < msp.out`;
    chomp $msp_lines;
    $k_param = 15 if $msp_lines > 100000;
    $k_param = 20 if $msp_lines > 500000;
    
    # Run RECON
    system("imagespread seq_name.list msp.out");
    system("eledef msp.out.spread");
    system("eleredef seq_name.list ele_def_res -K $k_param");
    system("edgeredef seq_name.list");
    system("famdef edge_redef_res seq_name.list");
    
    chdir "..";
}

sub calculate_round_metrics {
    my ($round_dir, $sample_size_mb) = @_;
    my %metrics;
    
    # Count families
    if (-f "$round_dir/summary/families") {
        $metrics{total_families} = `wc -l < $round_dir/summary/families`;
        chomp $metrics{total_families};
        
        # Count singletons
        $metrics{singletons} = `grep -c '^1\t' $round_dir/summary/families`;
        chomp $metrics{singletons};
    }
    
    # MSP density
    if (-f "$round_dir/msp.out") {
        $metrics{msp_lines} = `wc -l < $round_dir/msp.out`;
        chomp $metrics{msp_lines};
        $metrics{msp_density} = $metrics{msp_lines} / ($sample_size_mb * 1000);  # per kb
    }
    
    # Family size distribution
    if (-f "$round_dir/summary/eles") {
        my @family_sizes;
        open(my $fh, '<', "$round_dir/summary/eles");
        my $current_family = "";
        my $family_size = 0;
        
        while (<$fh>) {
            next if /^#/;
            my @fields = split;
            if ($fields[0] ne $current_family) {
                push @family_sizes, $family_size if $family_size > 0;
                $current_family = $fields[0];
                $family_size = 1;
            } else {
                $family_size++;
            }
        }
        push @family_sizes, $family_size if $family_size > 0;
        close($fh);
        
        @family_sizes = sort {$b <=> $a} @family_sizes;
        $metrics{p95_family_size} = $family_sizes[int(@family_sizes * 0.05)] || 0;
    }
    
    return %metrics;
}

sub should_stop_sampling {
    my ($prev_metrics, $curr_metrics, $round) = @_;
    
    return 0 if $round < 2;  # Always do at least 2 rounds
    
    # Calculate derived metrics
    my $new_families = ($curr_metrics->{total_families} || 0) - ($prev_metrics->{total_families} || 0);
    my $new_families_per_100mb = $new_families * 100 / (270 - 90);  # Assuming 90->270 transition
    
    my $singleton_ratio = 0;
    if ($curr_metrics->{total_families} && $curr_metrics->{total_families} > 0) {
        $singleton_ratio = $curr_metrics->{singletons} / $curr_metrics->{total_families};
    }
    
    my $p95_change = 0;
    if ($prev_metrics->{p95_family_size} && $prev_metrics->{p95_family_size} > 0) {
        $p95_change = ($curr_metrics->{p95_family_size} - $prev_metrics->{p95_family_size}) 
                     / $prev_metrics->{p95_family_size};
    }
    
    my $msp_density_change = 0;
    if ($prev_metrics->{msp_density} && $prev_metrics->{msp_density} > 0) {
        $msp_density_change = ($curr_metrics->{msp_density} - $prev_metrics->{msp_density})
                             / $prev_metrics->{msp_density};
    }
    
    print "Round $round metrics:\n";
    print "  New families per 100MB: $new_families_per_100mb\n";
    print "  Singleton ratio: " . sprintf("%.2f", $singleton_ratio) . "\n";
    print "  P95 family size change: " . sprintf("%.2f%%", $p95_change * 100) . "\n";
    print "  MSP density change: " . sprintf("%.2f%%", $msp_density_change * 100) . "\n";
    
    # Stopping criteria
    # Stop if: low new families AND high singleton ratio
    if ($new_families_per_100mb < 15 && $singleton_ratio > 0.6) {
        print "Stopping: Low new families and high singleton ratio\n";
        return 1;
    }
    
    # Stop if: P95 not growing AND MSP density not increasing much
    if (abs($p95_change) < 0.05 && $msp_density_change < 0.15) {
        print "Stopping: Family sizes stabilized and MSP density plateau\n";
        return 1;
    }
    
    # Continue conditions
    # Continue if: good new family discovery rate
    if ($new_families_per_100mb >= 30) {
        print "Continuing: High new family discovery rate\n";
        return 0;
    }
    
    # Continue if: P95 still growing significantly
    if ($p95_change >= 0.20) {
        print "Continuing: Large families still expanding\n";
        return 0;
    }
    
    # Continue if: MSP density increasing significantly
    if ($msp_density_change >= 0.15) {
        print "Continuing: MSP density still increasing\n";
        return 0;
    }
    
    # Default: stop if none of the continue conditions met
    return 1;
}

sub update_accumulated_mask {
    my ($mask_ref, $new_regions) = @_;
    
    return unless -s $new_regions;
    
    my $temp_mask = "temp_mask.bed";
    system("cat $$mask_ref $new_regions | bedtools sort | bedtools merge > $temp_mask");
    $$mask_ref = $temp_mask;
}

1;