package RECON::Utils;

use strict;
use warnings;
use Exporter 'import';
use File::Path qw(make_path);
use Cwd qw(getcwd abs_path);
use RECON::Logger;

our @EXPORT = qw(
    run_cmd create_seq_name_list split_fasta_by_size rename_fasta_ids
    perform_adaptive_sampling calculate_bed_size perform_sampling
    calculate_available_size filter_self_hits
    run_trf_masking run_repeatmasker create_blast_database
    format_size classify_genome_size check_required_tools
);

# Genome size thresholds
our %GENOME_THRESHOLDS = (
    TINY   => 500_000_000,    # 500MB
    SMALL  => 1_000_000_000,  # 1GB
    MEDIUM => 5_000_000_000,  # 5GB
    LARGE  => 10_000_000_000, # 10GB
);

# Constants
use constant {
    FRAGMENT_SIZE => 40000,
    MAX_RECON_ITERATIONS => 50000,
};

sub run_cmd {
    my ($cmd) = @_;
    log_message("DEBUG", "Executing command", $cmd);
    
    my $exit_code = system($cmd);
    if ($exit_code != 0) {
        die "Command failed with exit code $exit_code: $cmd\n";
    }
    
    return $exit_code;
}

sub create_seq_name_list {
    my ($fasta_file) = @_;
    
    # First pass: count sequences
    my $seq_count = 0;
    my @seq_names;
    
    open(my $fasta_fh, '<', $fasta_file) or die "Cannot open $fasta_file: $!\n";
    while (my $line = <$fasta_fh>) {
        if ($line =~ /^>(\S+)/) {
            push @seq_names, $1;
            $seq_count++;
        }
    }
    close($fasta_fh);
    
    # Sort sequence names lexicographically for imagespread compatibility
    # This produces the order: gi|1, gi|10, gi|100, gi|1000, etc.
    @seq_names = sort @seq_names;
    
    # Write seq.names file with count on first line (required by RECON tools)
    open(my $list_fh, '>', 'seq.names') or die "Cannot create seq.names: $!\n";
    print $list_fh "$seq_count\n";  # First line must be sequence count
    
    for my $name (@seq_names) {
        print $list_fh "$name\n";
    }
    
    close($list_fh);
    
    log_message("INFO", "Created sequence name list", "file=seq.names, count=$seq_count");
}

sub split_fasta_by_size {
    my ($input_file, $output_prefix, $max_size_mb) = @_;
    $max_size_mb ||= 100; # Default 100MB chunks
    
    my $max_size_bytes = $max_size_mb * 1024 * 1024;
    my $current_size = 0;
    my $file_count = 1;
    my $current_output;
    my @output_files;
    my $current_record = "";
    my $in_sequence = 0;
    
    open(my $input_fh, '<', $input_file) or die "Cannot open $input_file: $!\n";
    
    while (my $line = <$input_fh>) {
        if ($line =~ /^>/) {
            # Header line - start of new sequence
            
            # If we have a previous record and size limit reached, create new chunk
            if ($current_record && $current_size >= $max_size_bytes) {
                # Write current record to current chunk
                print $current_output $current_record if $current_output;
                close($current_output) if $current_output;
                
                # Start new chunk
                my $output_file = "${output_prefix}_part${file_count}.fa";
                open($current_output, '>', $output_file) or die "Cannot create $output_file: $!\n";
                push @output_files, $output_file;
                $file_count++;
                $current_size = 0;
                $current_record = "";
            }
            
            # If no current output file, create first one
            if (!$current_output) {
                my $output_file = "${output_prefix}_part${file_count}.fa";
                open($current_output, '>', $output_file) or die "Cannot create $output_file: $!\n";
                push @output_files, $output_file;
                $file_count++;
            }
            
            # Write previous record if exists
            if ($current_record) {
                print $current_output $current_record;
                $current_size += length($current_record);
            }
            
            # Start new record
            $current_record = $line;
            $in_sequence = 1;
        } elsif ($in_sequence) {
            # Sequence line
            $current_record .= $line;
        } else {
            # Non-FASTA line (shouldn't happen in proper FASTA)
            $current_record .= $line;
        }
    }
    
    # Write final record
    if ($current_record && $current_output) {
        print $current_output $current_record;
    }
    
    close($input_fh);
    close($current_output) if $current_output;
    
    log_message("INFO", "FASTA file split by complete records", 
                "input=$input_file, chunks=" . scalar(@output_files) . 
                ", target_size_mb=$max_size_mb");
    
    return @output_files;
}

sub rename_fasta_ids {
    my ($input_file, $output_file, $prefix) = @_;
    $prefix ||= "seq";
    
    open(my $input_fh, '<', $input_file) or die "Cannot open $input_file: $!\n";
    open(my $output_fh, '>', $output_file) or die "Cannot create $output_file: $!\n";
    
    my $seq_count = 1;
    while (my $line = <$input_fh>) {
        if ($line =~ /^>/) {
            if ($prefix eq "gi") {
                # Special formatting for gi format: >gi|1, >gi|2, etc.
                print $output_fh ">${prefix}|${seq_count}\n";
            } else {
                # Default formatting: >prefix_1, >prefix_2, etc.
                print $output_fh ">${prefix}_${seq_count}\n";
            }
            $seq_count++;
        } else {
            print $output_fh $line;
        }
    }
    
    close($input_fh);
    close($output_fh);
    
    log_message("INFO", "Renamed FASTA IDs", "sequences=$seq_count, output=$output_file");
    return $seq_count - 1;
}

sub perform_adaptive_sampling {
    my ($genome_file, $exclusion_bed, $accumulated_mask, $output_file, $sample_size_mb, $fragment_size) = @_;
    
    log_message("INFO", "Starting adaptive sampling", 
                "target_size=${sample_size_mb}MB, fragment_size=$fragment_size");
    
    # Use bedtools to exclude regions and sample
    my $cmd = "bedtools random -l $fragment_size -n " . int($sample_size_mb * 1000000 / $fragment_size) . 
              " -g ${genome_file}.fai";
    
    if ($exclusion_bed && -s $exclusion_bed) {
        $cmd .= " | bedtools subtract -a stdin -b $exclusion_bed";
    }
    
    $cmd .= " | bedtools getfasta -fi $genome_file -bed stdin > $output_file";
    
    run_cmd($cmd);
    
    my $actual_size = -s $output_file || 0;
    log_message("INFO", "Adaptive sampling completed", 
                "output_size=" . format_size($actual_size) . ", file=$output_file");
}

sub calculate_available_size {
    my ($genome_file, $exclusion_bed) = @_;
    
    # Calculate total genome size
    my $total_size = 0;
    if (-f "${genome_file}.fai") {
        $total_size = `awk '{sum += \$2} END {print sum}' ${genome_file}.fai`;
        chomp $total_size;
    } else {
        $total_size = -s $genome_file || 0;
    }
    
    # Calculate excluded size
    my $excluded_size = 0;
    if ($exclusion_bed && -s $exclusion_bed) {
        $excluded_size = calculate_bed_size($exclusion_bed);
    }
    
    return $total_size - $excluded_size;
}

sub calculate_bed_size {
    my ($bed_file) = @_;
    return 0 unless -s $bed_file;
    
    my $total_size = `awk '{sum += (\$3 - \$2)} END {print sum+0}' $bed_file`;
    chomp $total_size;
    return $total_size || 0;
}

sub perform_sampling {
    my ($genome_file, $bed_file, $output_file) = @_;
    
    my $cmd = "bedtools getfasta -fi $genome_file -bed $bed_file > $output_file";
    run_cmd($cmd);
    
    log_message("INFO", "Sampling completed", "output=$output_file");
}

sub filter_self_hits {
    my ($blast_file, $output_file, $min_length, $max_evalue) = @_;
    $min_length ||= 50;
    $max_evalue ||= 1e-5;
    
    open(my $input_fh, '<', $blast_file) or die "Cannot open $blast_file: $!\n";
    open(my $output_fh, '>', $output_file) or die "Cannot create $output_file: $!\n";
    
    my $filtered_count = 0;
    while (my $line = <$input_fh>) {
        chomp $line;
        my @fields = split /\t/, $line;
        next if @fields < 12;
        
        my ($query, $subject, $identity, $length, $mismatches, $gaps,
            $qstart, $qend, $sstart, $send, $evalue, $bitscore) = @fields;
        
        # Filter self-hits and low quality hits
        next if $query eq $subject;
        next if $length < $min_length;
        next if $evalue > $max_evalue;
        
        print $output_fh "$line\n";
        $filtered_count++;
    }
    
    close($input_fh);
    close($output_fh);
    
    log_message("INFO", "Filtered BLAST hits", "kept=$filtered_count hits");
    return $filtered_count;
}

sub run_trf_masking {
    my ($input_file, $output_file, $threads) = @_;
    $threads ||= 4;
    
    # Try RepeatMasker first, then TRF, then dustmasker
    my @masking_tools = (
        "RepeatMasker -species simple -nolow -no_is -norna -parallel $threads -dir . $input_file && mv ${input_file}.masked $output_file",
        "trf $input_file 2 7 7 80 10 50 500 -m -h && mv ${input_file}.2.7.7.80.10.50.500.mask $output_file",
        "dustmasker -in $input_file -out $output_file -outfmt fasta"
    );
    
    for my $cmd (@masking_tools) {
        eval {
            run_cmd($cmd);
            return; # Success
        };
    }
    
    # If all masking tools fail, just copy the input
    run_cmd("cp $input_file $output_file");
    log_message("WARN", "TRF masking failed, using unmasked sequence", "file=$output_file");
}

sub run_repeatmasker {
    my ($input_file, $library_file, $output_file, $threads) = @_;
    $threads ||= 4;
    
    my $success = 0;
    my $cmd = "RepeatMasker -lib $library_file -nolow -no_is -norna -parallel $threads -dir . $input_file";
    
    eval {
        run_cmd($cmd);
        if (-f "${input_file}.masked") {
            run_cmd("mv ${input_file}.masked $output_file");
            $success = 1;
        }
    };
    
    if ($@ || !-s $output_file) {
        # Fallback: copy original if RepeatMasker fails
        run_cmd("cp $input_file $output_file");
        log_message("WARN", "RepeatMasker failed, using unmasked sequence", "file=$output_file, error=" . ($@ || "output file empty"));
        $success = 0;
    }
    
    return $success;
}

sub create_blast_database {
    my ($fasta_file, $db_name) = @_;
    
    my $cmd = "makeblastdb -in $fasta_file -dbtype nucl -out $db_name";
    run_cmd($cmd);
    
    log_message("INFO", "Created BLAST database", "database=$db_name");
}

sub format_size {
    my ($bytes) = @_;
    return "0B" unless $bytes;
    
    my @units = qw(B KB MB GB TB);
    my $unit_index = 0;
    
    while ($bytes >= 1024 && $unit_index < @units - 1) {
        $bytes /= 1024;
        $unit_index++;
    }
    
    return sprintf("%.1f%s", $bytes, $units[$unit_index]);
}

sub classify_genome_size {
    my ($size) = @_;
    
    return 'TINY' if $size < $GENOME_THRESHOLDS{TINY};
    return 'SMALL' if $size < $GENOME_THRESHOLDS{SMALL};
    return 'MEDIUM' if $size < $GENOME_THRESHOLDS{MEDIUM};
    return 'LARGE' if $size < $GENOME_THRESHOLDS{LARGE};
    return 'HUGE';
}

sub check_required_tools {
    my @required_tools = qw(
        edgeredef eleredef eledef famdef
        blastn makeblastdb seqkit mafft bedtools
    );
    
    my @missing_tools;
    for my $tool (@required_tools) {
        unless (system("which $tool > /dev/null 2>&1") == 0) {
            push @missing_tools, $tool;
        }
    }
    
    if (@missing_tools) {
        log_message("ERROR", "Missing required tools", join(", ", @missing_tools));
        die "Required tools not found. Please ensure RECON and dependencies are installed.\n";
    }
    
    log_message("INFO", "All required tools found");
}

1;

__END__

=head1 NAME

RECON::Utils - Utility functions for RECON pipeline

=head1 DESCRIPTION

Provides common utility functions for genome processing, sampling, masking,
and BLAST operations used throughout the RECON pipeline.

=cut