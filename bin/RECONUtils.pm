package RECONUtils;

################################################################################
# RECONUtils - Utility functions for Advanced RECON Pipeline
# Version: 3.0
################################################################################

use strict;
use warnings;
use 5.010;

use base 'Exporter';
use File::Temp qw(tempfile tempdir);
use List::Util qw(min max sum);
use POSIX qw(floor ceil);

our @EXPORT = qw(
    calculate_gc_content
    calculate_gc_content_string
    estimate_repeat_content
    calculate_kmer_complexity
    generate_sampling_layer
    generate_coarse_sample
    identify_te_hotspots
    generate_focused_samples
    check_for_tirs
    check_for_ltrs
    check_for_tsds
    assess_coding_potential
    perform_msa
    analyze_alignment_edges
    trim_consensus
    calculate_sequence_distances
    hierarchical_clustering
    build_consensus_sequence
    calculate_average_divergence
    split_chimeric_family
    families_are_redundant
    merge_family_info
    classify_ltr_family
    classify_dna_family
    calculate_std_dev
);

# ===================================================================
# SEQUENCE ANALYSIS FUNCTIONS
# ===================================================================

sub calculate_gc_content {
    my ($file) = @_;
    
    my $total_bases = 0;
    my $gc_bases = 0;
    
    open(my $fh, '<', $file) or return 0.5;
    while (my $line = <$fh>) {
        next if $line =~ /^>/;
        chomp $line;
        $total_bases += length($line);
        $gc_bases += ($line =~ tr/GCgc//);
    }
    close $fh;
    
    return $total_bases > 0 ? $gc_bases / $total_bases : 0.5;
}

sub calculate_gc_content_string {
    my ($seq) = @_;
    
    my $total = length($seq);
    return 0.5 unless $total > 0;
    
    my $gc = ($seq =~ tr/GCgc//);
    return $gc / $total;
}

sub estimate_repeat_content {
    my ($file) = @_;
    
    # Quick repeat estimation using k-mer frequency
    my %kmers;
    my $k = 15;
    my $total_kmers = 0;
    
    open(my $fh, '<', $file) or return 0;
    while (my $line = <$fh>) {
        next if $line =~ /^>/;
        chomp $line;
        
        for (my $i = 0; $i <= length($line) - $k; $i++) {
            my $kmer = substr($line, $i, $k);
            $kmers{$kmer}++;
            $total_kmers++;
        }
    }
    close $fh;
    
    # Calculate repeat score based on k-mer redundancy
    my $unique_kmers = scalar(keys %kmers);
    my $repeat_score = 1 - ($unique_kmers / $total_kmers);
    
    return $repeat_score;
}

sub calculate_kmer_complexity {
    my ($file, $k) = @_;
    $k ||= 15;
    
    my %kmers;
    my $total = 0;
    
    open(my $fh, '<', $file) or return 0.5;
    while (my $line = <$fh>) {
        next if $line =~ /^>/;
        chomp $line;
        
        for (my $i = 0; $i <= length($line) - $k; $i++) {
            my $kmer = substr($line, $i, $k);
            $kmers{$kmer}++;
            $total++;
        }
    }
    close $fh;
    
    # Shannon entropy of k-mer distribution
    my $entropy = 0;
    foreach my $count (values %kmers) {
        my $p = $count / $total;
        $entropy -= $p * log($p) / log(2);
    }
    
    # Normalize by maximum possible entropy
    my $max_entropy = log(scalar(keys %kmers)) / log(2);
    
    return $max_entropy > 0 ? $entropy / $max_entropy : 0;
}

# ===================================================================
# SAMPLING FUNCTIONS
# ===================================================================

sub generate_sampling_layer {
    my ($genome_file, $layer) = @_;
    
    if ($layer->{method} && $layer->{method} eq 'random') {
        # Random sampling
        my $num_samples = int($layer->{sample_size} / 10000);  # Assume 10kb per sample
        system("seqkit sample -s 42 -n $num_samples $genome_file > $layer->{file} 2>/dev/null");
    }
    else {
        # Window-based sampling
        my $window_cmd = "seqkit sliding " .
                        "-s $layer->{overlap} " .
                        "-W $layer->{window_size} " .
                        "$genome_file | " .
                        "seqkit head -n " . int($layer->{sample_size} / $layer->{window_size}) .
                        " > $layer->{file} 2>/dev/null";
        system($window_cmd);
    }
}

sub generate_coarse_sample {
    my ($genome_file, $genome_size) = @_;
    
    my $sample_size = min(1_000_000_000, $genome_size * 0.1);  # 10% or 1GB max
    my $sample_file = "$ENV{PWD}/tmp/coarse_sample.fa";
    
    # Use stratified sampling
    my $num_windows = 1000;
    my $window_size = int($genome_size / $num_windows);
    
    system("seqkit sliding -s $window_size -W $window_size $genome_file | " .
           "seqkit sample -s 42 -p 0.1 > $sample_file 2>/dev/null");
    
    return $sample_file;
}

sub identify_te_hotspots {
    my ($sample_file) = @_;
    
    # Run quick RepeatMasker to identify TE-rich regions
    my $rm_output = "$ENV{PWD}/tmp/hotspot_rm.out";
    system("RepeatMasker -pa 4 -q -no_is -norna -nolow -div 40 $sample_file 2>/dev/null");
    
    my %hotspots;
    
    if (-f "$sample_file.out") {
        open(my $fh, '<', "$sample_file.out");
        while (my $line = <$fh>) {
            next if $line =~ /^(SW|score|\s*$)/;
            
            my @fields = split(/\s+/, $line);
            next unless @fields >= 15;
            
            my $seq_id = $fields[4];
            my $start = $fields[5];
            my $end = $fields[6];
            
            push @{$hotspots{$seq_id}}, {
                start => $start,
                end => $end,
                repeat_class => $fields[10],
            };
        }
        close $fh;
    }
    
    return \%hotspots;
}

sub generate_focused_samples {
    my ($genome_file, $hotspots) = @_;
    
    my @focused_samples;
    my $sample_idx = 0;
    
    foreach my $seq_id (keys %$hotspots) {
        my $regions = $hotspots->{$seq_id};
        next unless scalar(@$regions) > 5;  # Focus on TE-rich sequences
        
        # Merge overlapping regions
        my @merged = merge_overlapping_regions($regions);
        
        foreach my $region (@merged) {
            $sample_idx++;
            my $sample_file = "$ENV{PWD}/tmp/focused_sample_$sample_idx.fa";
            
            # Extract region with flanking sequences
            my $start = max(1, $region->{start} - 5000);
            my $end = $region->{end} + 5000;
            
            system("seqkit subseq -r ${start}:${end} $genome_file > $sample_file 2>/dev/null");
            
            push @focused_samples, $sample_file if -s $sample_file;
        }
    }
    
    return \@focused_samples;
}

sub merge_overlapping_regions {
    my ($regions) = @_;
    
    # Sort by start position
    my @sorted = sort { $a->{start} <=> $b->{start} } @$regions;
    
    my @merged;
    my $current = $sorted[0];
    
    for my $i (1..$#sorted) {
        my $next = $sorted[$i];
        
        if ($next->{start} <= $current->{end} + 1000) {  # Allow 1kb gap
            # Merge
            $current->{end} = max($current->{end}, $next->{end});
        }
        else {
            push @merged, $current;
            $current = $next;
        }
    }
    
    push @merged, $current;
    
    return @merged;
}

# ===================================================================
# STRUCTURAL FEATURE DETECTION  
# ===================================================================

sub check_for_tirs {
    my ($sequence) = @_;
    
    return { found => 0 } unless $sequence && length($sequence) > 20;
    
    # Check for terminal inverted repeats
    my $max_tir_length = min(100, int(length($sequence) / 2));
    
    for (my $tir_len = 10; $tir_len <= $max_tir_length; $tir_len++) {
        my $left_term = substr($sequence, 0, $tir_len);
        my $right_term = substr($sequence, -$tir_len);
        
        # Reverse complement of right terminal
        my $right_rc = reverse_complement($right_term);
        
        # Calculate identity
        my $matches = 0;
        for (my $i = 0; $i < $tir_len; $i++) {
            $matches++ if substr($left_term, $i, 1) eq substr($right_rc, $i, 1);
        }
        
        my $identity = $matches / $tir_len;
        
        if ($identity >= 0.8) {  # 80% identity threshold
            return {
                found => 1,
                length => $tir_len,
                identity => $identity,
                left => $left_term,
                right => $right_term,
            };
        }
    }
    
    return { found => 0 };
}

sub check_for_ltrs {
    my ($sequence) = @_;
    
    return { found => 0 } unless $sequence && length($sequence) > 500;
    
    # LTRs are typically 100-5000bp
    my $seq_len = length($sequence);
    my $min_ltr = 100;
    my $max_ltr = min(5000, int($seq_len * 0.4));
    
    for (my $ltr_len = $max_ltr; $ltr_len >= $min_ltr; $ltr_len -= 50) {
        my $left_ltr = substr($sequence, 0, $ltr_len);
        my $right_ltr = substr($sequence, -$ltr_len);
        
        # Calculate similarity
        my $matches = 0;
        for (my $i = 0; $i < $ltr_len; $i++) {
            $matches++ if substr($left_ltr, $i, 1) eq substr($right_ltr, $i, 1);
        }
        
        my $identity = $matches / $ltr_len;
        
        if ($identity >= 0.8) {  # 80% identity threshold
            return {
                found => 1,
                length => $ltr_len,
                identity => $identity,
            };
        }
    }
    
    return { found => 0 };
}

sub check_for_tsds {
    my ($family) = @_;
    
    # TSDs are typically 2-20bp
    my @tsd_lengths = (4, 5, 9, 2, 3, 6, 7, 8, 10, 11, 12);
    
    foreach my $tsd_len (@tsd_lengths) {
        my $tsd_counts = {};
        
        # Check each family member for TSDs
        foreach my $member (@{$family->{members}}) {
            next unless $member->{flanking_left} && $member->{flanking_right};
            
            my $left_tsd = substr($member->{flanking_left}, -$tsd_len);
            my $right_tsd = substr($member->{flanking_right}, 0, $tsd_len);
            
            if ($left_tsd eq $right_tsd) {
                $tsd_counts->{$left_tsd}++;
            }
        }
        
        # Check if we found consistent TSDs
        foreach my $tsd (keys %$tsd_counts) {
            if ($tsd_counts->{$tsd} >= scalar(@{$family->{members}}) * 0.5) {
                return {
                    found => 1,
                    length => $tsd_len,
                    sequence => $tsd,
                    frequency => $tsd_counts->{$tsd} / scalar(@{$family->{members}}),
                };
            }
        }
    }
    
    return { found => 0 };
}

sub assess_coding_potential {
    my ($sequence) = @_;
    
    return { has_orfs => 0 } unless $sequence;
    
    my $longest_orf = 0;
    my @orfs;
    
    # Check all three reading frames
    for my $frame (0, 1, 2) {
        my $orf_start = -1;
        
        for (my $i = $frame; $i < length($sequence) - 2; $i += 3) {
            my $codon = substr($sequence, $i, 3);
            
            # Start codon
            if ($codon eq 'ATG' && $orf_start < 0) {
                $orf_start = $i;
            }
            # Stop codon
            elsif (($codon eq 'TAA' || $codon eq 'TAG' || $codon eq 'TGA') && $orf_start >= 0) {
                my $orf_length = $i - $orf_start + 3;
                push @orfs, {
                    start => $orf_start,
                    end => $i + 2,
                    length => $orf_length,
                    frame => $frame,
                };
                $longest_orf = max($longest_orf, $orf_length);
                $orf_start = -1;
            }
        }
    }
    
    return {
        has_orfs => scalar(@orfs) > 0,
        num_orfs => scalar(@orfs),
        longest_orf => $longest_orf,
        orfs => \@orfs,
    };
}

# ===================================================================
# ALIGNMENT AND CLUSTERING FUNCTIONS
# ===================================================================

sub perform_msa {
    my ($sequences, $output_file) = @_;
    
    # Write sequences to temp file
    my ($temp_fh, $temp_file) = tempfile();
    for (my $i = 0; $i < @$sequences; $i++) {
        print $temp_fh ">seq_$i\n$sequences->[$i]\n";
    }
    close $temp_fh;
    
    # Run MAFFT
    system("mafft --quiet --auto $temp_file > $output_file 2>/dev/null");
    
    unlink $temp_file;
}

sub analyze_alignment_edges {
    my ($alignment_file) = @_;
    
    # Read alignment
    my @aligned_seqs;
    my $current_seq = '';
    
    open(my $fh, '<', $alignment_file);
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ /^>/) {
            push @aligned_seqs, $current_seq if $current_seq;
            $current_seq = '';
        }
        else {
            $current_seq .= $line;
        }
    }
    push @aligned_seqs, $current_seq if $current_seq;
    close $fh;
    
    return { needs_trimming => 0 } unless @aligned_seqs;
    
    my $aln_length = length($aligned_seqs[0]);
    my $num_seqs = scalar(@aligned_seqs);
    
    # Calculate conservation at each position
    my @conservation;
    for (my $pos = 0; $pos < $aln_length; $pos++) {
        my %bases;
        for my $seq (@aligned_seqs) {
            my $base = substr($seq, $pos, 1);
            $bases{$base}++ unless $base eq '-';
        }
        
        my $max_count = 0;
        foreach my $count (values %bases) {
            $max_count = max($max_count, $count);
        }
        
        push @conservation, $max_count / $num_seqs;
    }
    
    # Find trim points (conservation < 50%)
    my $trim_start = 0;
    my $trim_end = $aln_length;
    
    # Find start trim point
    for (my $i = 0; $i < @conservation; $i++) {
        if ($conservation[$i] >= 0.5) {
            $trim_start = $i;
            last;
        }
    }
    
    # Find end trim point
    for (my $i = $#conservation; $i >= 0; $i--) {
        if ($conservation[$i] >= 0.5) {
            $trim_end = $i + 1;
            last;
        }
    }
    
    return {
        needs_trimming => ($trim_start > 0 || $trim_end < $aln_length),
        trim_start => $trim_start,
        trim_end => $trim_end,
    };
}

sub calculate_sequence_distances {
    my ($members) = @_;
    
    my $n = scalar(@$members);
    my @distance_matrix;
    
    for (my $i = 0; $i < $n; $i++) {
        for (my $j = 0; $j < $n; $j++) {
            if ($i == $j) {
                $distance_matrix[$i][$j] = 0;
            }
            elsif ($j > $i) {
                my $dist = calculate_pairwise_distance(
                    $members->[$i]{sequence},
                    $members->[$j]{sequence}
                );
                $distance_matrix[$i][$j] = $dist;
                $distance_matrix[$j][$i] = $dist;
            }
        }
    }
    
    return \@distance_matrix;
}

sub calculate_pairwise_distance {
    my ($seq1, $seq2) = @_;
    
    return 1.0 unless $seq1 && $seq2;
    
    # Use simple Hamming distance for speed
    my $len = min(length($seq1), length($seq2));
    my $mismatches = 0;
    
    for (my $i = 0; $i < $len; $i++) {
        $mismatches++ if substr($seq1, $i, 1) ne substr($seq2, $i, 1);
    }
    
    return $mismatches / $len;
}

sub hierarchical_clustering {
    my ($distance_matrix, $threshold) = @_;
    
    my $n = scalar(@$distance_matrix);
    
    # Initialize clusters (each sequence is its own cluster)
    my @clusters;
    for (my $i = 0; $i < $n; $i++) {
        push @clusters, [$i];
    }
    
    # Single-linkage clustering
    while (scalar(@clusters) > 1) {
        # Find minimum distance between clusters
        my ($min_dist, $merge_i, $merge_j) = (1.0, -1, -1);
        
        for (my $i = 0; $i < @clusters; $i++) {
            for (my $j = $i + 1; $j < @clusters; $j++) {
                my $dist = cluster_distance($clusters[$i], $clusters[$j], $distance_matrix);
                if ($dist < $min_dist) {
                    $min_dist = $dist;
                    $merge_i = $i;
                    $merge_j = $j;
                }
            }
        }
        
        # Stop if minimum distance exceeds threshold
        last if $min_dist > $threshold;
        
        # Merge clusters
        if ($merge_i >= 0 && $merge_j >= 0) {
            push @{$clusters[$merge_i]}, @{$clusters[$merge_j]};
            splice(@clusters, $merge_j, 1);
        }
        else {
            last;
        }
    }
    
    return \@clusters;
}

sub cluster_distance {
    my ($cluster1, $cluster2, $distance_matrix) = @_;
    
    # Single linkage: minimum distance between any two members
    my $min_dist = 1.0;
    
    foreach my $i (@$cluster1) {
        foreach my $j (@$cluster2) {
            $min_dist = min($min_dist, $distance_matrix->[$i][$j]);
        }
    }
    
    return $min_dist;
}

sub build_consensus_sequence {
    my ($members) = @_;
    
    return '' unless @$members;
    return $members->[0]{sequence} if @$members == 1;
    
    # Get sequences
    my @sequences = map { $_->{sequence} } @$members;
    
    # Find minimum length
    my $min_len = min(map { length($_) } @sequences);
    
    # Build consensus
    my $consensus = '';
    for (my $pos = 0; $pos < $min_len; $pos++) {
        my %bases;
        for my $seq (@sequences) {
            my $base = substr($seq, $pos, 1);
            $bases{$base}++;
        }
        
        # Get most frequent base
        my ($best_base) = sort { $bases{$b} <=> $bases{$a} } keys %bases;
        $consensus .= $best_base;
    }
    
    return $consensus;
}

sub calculate_average_divergence {
    my ($members) = @_;
    
    return 0 unless @$members > 1;
    
    my $total_dist = 0;
    my $count = 0;
    
    for (my $i = 0; $i < @$members; $i++) {
        for (my $j = $i + 1; $j < @$members; $j++) {
            $total_dist += calculate_pairwise_distance(
                $members->[$i]{sequence},
                $members->[$j]{sequence}
            );
            $count++;
        }
    }
    
    return $count > 0 ? $total_dist / $count : 0;
}

# ===================================================================
# UTILITY FUNCTIONS
# ===================================================================

sub reverse_complement {
    my ($seq) = @_;
    
    $seq = reverse($seq);
    $seq =~ tr/ATGCatgc/TACGtacg/;
    
    return $seq;
}

sub trim_consensus {
    my ($sequence, $trim_start, $trim_end) = @_;
    
    return substr($sequence, $trim_start, $trim_end - $trim_start);
}

sub split_chimeric_family {
    my ($family) = @_;
    
    # Simple splitting based on sequence similarity clustering
    my $distance_matrix = calculate_sequence_distances($family->{members});
    my $clusters = hierarchical_clustering($distance_matrix, 0.25);  # 25% divergence
    
    return undef unless @$clusters > 1;
    
    my @split_families;
    for (my $i = 0; $i < @$clusters; $i++) {
        my $cluster = $clusters->[$i];
        next unless @$cluster >= 2;  # Minimum family size
        
        my @cluster_members = map { $family->{members}[$_] } @$cluster;
        
        push @split_families, {
            id => $i + 1,
            members => \@cluster_members,
            consensus => build_consensus_sequence(\@cluster_members),
            parent_id => $family->{id},
        };
    }
    
    return \@split_families;
}

sub families_are_redundant {
    my ($family1, $family2) = @_;
    
    return 0 unless $family1->{consensus} && $family2->{consensus};
    
    my $len1 = length($family1->{consensus});
    my $len2 = length($family2->{consensus});
    
    # Length difference check
    my $len_ratio = min($len1, $len2) / max($len1, $len2);
    return 0 if $len_ratio < 0.8;
    
    # Sequence similarity check
    my $distance = calculate_pairwise_distance(
        $family1->{consensus},
        $family2->{consensus}
    );
    
    return $distance < 0.05;  # 95% similarity
}

sub merge_family_info {
    my ($target_family, $source_family) = @_;
    
    # Merge members
    push @{$target_family->{members}}, @{$source_family->{members}};
    
    # Update consensus if we have more members now
    if (scalar(@{$target_family->{members}}) > 10) {
        $target_family->{consensus} = build_consensus_sequence($target_family->{members});
    }
    
    # Merge structural annotations
    foreach my $key (qw(structure_type tir_length ltr_length tsd_length)) {
        $target_family->{$key} ||= $source_family->{$key} if $source_family->{$key};
    }
}

sub classify_ltr_family {
    my ($family) = @_;
    
    # Simple classification based on length and features
    my $length = length($family->{consensus} || '');
    
    if ($length > 4000) {
        if ($family->{has_coding} && $family->{longest_orf} > 1000) {
            return 'Gypsy';  # Likely Gypsy element
        }
        return 'Copia';  # Likely Copia element
    }
    
    return 'LTR_unknown';
}

sub classify_dna_family {
    my ($family) = @_;
    
    # Classification based on TIR and TSD features
    if ($family->{tsd_length}) {
        if ($family->{tsd_length} == 2) {
            return 'hAT';
        }
        elsif ($family->{tsd_length} == 9) {
            return 'MuDR';
        }
        elsif ($family->{tsd_length} == 5) {
            return 'PIF_Harbinger';
        }
    }
    
    return 'DNA_unknown';
}

sub calculate_std_dev {
    my ($values) = @_;
    
    return 0 unless @$values;
    
    my $mean = sum(@$values) / @$values;
    my $variance = sum(map { ($_ - $mean) ** 2 } @$values) / @$values;
    
    return sqrt($variance);
}

1;