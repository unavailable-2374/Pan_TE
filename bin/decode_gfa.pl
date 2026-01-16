#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use autodie;
use File::Path qw(make_path remove_tree);
use Cwd qw(abs_path);
use POSIX qw(strftime);
use File::Spec;

# Constants
use constant {
    MIN_SEQ_LENGTH => 50,
    DEFAULT_THREADS => 4,
    MAX_RETRIES => 3,
};

our $VERSION = '1.0.0';

# Configuration
my %config = (
    vcf_file => undef,
    threads => DEFAULT_THREADS,
    out_dir => undef,
    debug => 0,
);

# Setup
init_environment();

# Main execution
process_data();
cleanup_tmp_files($config{out_dir});

exit(0);

###############################################################################
# Core Functions
###############################################################################

sub init_environment {
    parse_arguments();
    validate_environment();
    setup_directories();
    setup_logging();
    $config{out_dir} = abs_path($config{out_dir});
    $config{vcf_file} = abs_path($config{vcf_file});
}

sub process_data {
    my $base_name = basename($config{vcf_file}, ".vcf");
    my $output_file = File::Spec->catfile($config{out_dir}, "$base_name.fa");

    log_message("Starting parallel processing with $config{threads} workers");

    # Step 1: Count data lines (single pass)
    my $total_lines = count_data_lines($config{vcf_file});
    log_message("VCF has $total_lines data lines");

    # Step 2: Calculate line ranges for each worker
    my $num_workers = $config{threads};
    $num_workers = $total_lines if $num_workers > $total_lines;  # Don't create more workers than lines

    my $lines_per_worker = int(($total_lines + $num_workers - 1) / $num_workers);

    # Step 3: Fork workers to process line ranges
    my @output_files;
    my @pids;

    for my $i (0 .. $num_workers - 1) {
        my $start_line = $i * $lines_per_worker;
        my $end_line = ($i + 1) * $lines_per_worker - 1;
        $end_line = $total_lines - 1 if $end_line >= $total_lines;

        # Skip if this worker has no lines to process
        last if $start_line >= $total_lines;

        my $worker_output = File::Spec->catfile($config{out_dir}, "tmp_${base_name}_worker_${i}.fa");
        push @output_files, $worker_output;

        my $pid = fork();
        if (!defined $pid) {
            die "Fork failed: $!\n";
        } elsif ($pid == 0) {
            # Child process: read original file, process assigned line range
            process_line_range($config{vcf_file}, $worker_output, $start_line, $end_line, $i);
            exit(0);
        } else {
            push @pids, $pid;
        }
    }

    # Step 4: Wait for all workers
    for my $pid (@pids) {
        waitpid($pid, 0);
        if ($? != 0) {
            warn "Worker process $pid exited with status $?\n";
        }
    }

    # Step 5: Merge outputs
    merge_outputs(\@output_files, $output_file);
    log_message("Merged outputs to $output_file");

    # Step 6: Cleanup worker outputs
    for my $f (@output_files) {
        unlink $f if -e $f;
    }

    print STDERR "PROGRESS:100\n";
}

sub count_data_lines {
    my $file = shift;
    my $count = 0;

    open my $fh, '<', $file;
    while (<$fh>) {
        $count++ unless /^#/;
    }
    close $fh;

    return $count;
}

###############################################################################
# Setup Functions
###############################################################################

sub parse_arguments {
    GetOptions(
        'vcf=s'     => \$config{vcf_file},
        'threads=i' => \$config{threads},
        'out=s'     => \$config{out_dir},
        'debug'     => \$config{debug},
        'help|h'    => \my $help,
        'version|v' => \my $version,
    ) or pod2usage(2);

    pod2usage(1) if $help;
    
    if ($version) {
        print "decode_gfa.pl version $VERSION\n";
        exit(0);
    }

    # Validate required arguments
    pod2usage("Error: VCF file is required") unless $config{vcf_file};
    pod2usage("Error: Number of threads is required") unless $config{threads};
    
    # Set default output directory
    $config{out_dir} ||= dirname($config{vcf_file});
}

sub validate_environment {
    # Check input file
    die "Error: Input file '$config{vcf_file}' does not exist\n" 
        unless -e $config{vcf_file};
    die "Error: Input file '$config{vcf_file}' is not readable\n" 
        unless -r $config{vcf_file};
    
    # Validate thread count
    die "Error: Number of threads must be positive\n" 
        unless $config{threads} > 0;
        
    # Check output directory permissions
    die "Error: Cannot write to output directory '$config{out_dir}'\n"
        unless -w $config{out_dir};
}

sub setup_directories {
    make_path($config{out_dir}) unless -d $config{out_dir};
}

sub setup_logging {
    my $log_file = "$config{out_dir}/decode_gfa.log";
    open my $log_fh, '>>', $log_file;
    print $log_fh strftime("%Y-%m-%d %H:%M:%S", localtime), 
          " Starting decode_gfa.pl v$VERSION\n";
    close $log_fh;
}

###############################################################################
# Worker Processing Functions
###############################################################################

sub process_line_range {
    my ($vcf_file, $output_file, $start_line, $end_line, $worker_id) = @_;
    my $base_name = basename($vcf_file, ".vcf");
    my $temp_file = File::Spec->catfile($config{out_dir}, "tmp_${base_name}_allele_${worker_id}.fa");

    # Initialize worker output file
    open my $out_fh, '>', $output_file;

    # Open VCF and skip to our assigned range
    open my $in_fh, '<', $vcf_file;

    my $current_line = 0;  # Data line counter (excluding headers)
    my $processed = 0;

    while (my $line = <$in_fh>) {
        # Skip headers
        next if $line =~ /^#/;

        # Skip lines before our range
        if ($current_line < $start_line) {
            $current_line++;
            next;
        }

        # Stop after our range
        last if $current_line > $end_line;

        chomp $line;
        my @fields = split(/\t/, $line);
        $current_line++;

        next unless @fields >= 5;

        my $alt_allele = $fields[4];
        $processed++;

        # Report progress every 1000 entries
        if ($processed % 1000 == 0) {
            print STDERR "Worker $worker_id: processed $processed entries\n";
        }

        next unless length($alt_allele) > MIN_SEQ_LENGTH;

        if ($alt_allele =~ /,/) {
            # Multi-allelic: write to temp, refine, append to output
            process_multiple_alleles_worker($alt_allele, $out_fh, $temp_file);
        } else {
            print_sequence($out_fh, "SV", $alt_allele);
        }
    }

    close $in_fh;
    close $out_fh;

    # Cleanup temp file if exists
    unlink $temp_file if -e $temp_file;

    print STDERR "Worker $worker_id: completed $processed entries (lines $start_line-$end_line)\n";
}

sub process_multiple_alleles_worker {
    my ($alt_allele, $out_fh, $temp_file) = @_;
    my @alleles = split(/,/, $alt_allele);
    my $long_alleles = 0;

    open my $tmp_fh, '>', $temp_file;

    for my $i (0 .. $#alleles) {
        if (length($alleles[$i]) > MIN_SEQ_LENGTH) {
            print_sequence($tmp_fh, $i, $alleles[$i]);
            $long_alleles++;
        }
    }

    close $tmp_fh;

    if ($long_alleles > 1) {
        refine_alleles_worker($temp_file, $out_fh);
    } elsif ($long_alleles == 1) {
        # Append single allele to output
        open my $in_fh, '<', $temp_file;
        while (<$in_fh>) {
            print $out_fh $_;
        }
        close $in_fh;
    }

    unlink $temp_file;
}

sub refine_alleles_worker {
    my ($temp_file, $out_fh) = @_;

    my $refiner_output = "${temp_file}.out";

    my $cmd = join(" ",
        "Refiner_for_Graph",
        $temp_file,
        $refiner_output,
        "--distance-threshold", "0.8",
        "-t", "1",
        "-v"
    );

    system($cmd) == 0
        or die "Feature Extraction Failure: $?\n";

    # Append refiner output to worker output
    if (-e $refiner_output) {
        open my $in_fh, '<', $refiner_output;
        while (<$in_fh>) {
            print $out_fh $_;
        }
        close $in_fh;
        unlink $refiner_output;
    }
}

###############################################################################
# Output Merge Function
###############################################################################

sub merge_outputs {
    my ($output_files, $final_output) = @_;

    open my $out_fh, '>', $final_output;

    for my $file (@$output_files) {
        next unless -e $file && -s $file;  # Skip empty or non-existent files
        open my $in_fh, '<', $file;
        while (<$in_fh>) {
            print $out_fh $_;
        }
        close $in_fh;
    }

    close $out_fh;
}

###############################################################################
# Utility Functions
###############################################################################

sub print_sequence {
    my ($fh, $id, $sequence) = @_;
    print $fh ">$id\n$sequence\n";
}

sub log_message {
    my $message = shift;
    my $timestamp = strftime("%Y-%m-%d %H:%M:%S", localtime);
    print STDERR "[$timestamp] $message\n" if $config{debug};
    
    open my $log_fh, '>>', "$config{out_dir}/decode_gfa.log";
    print $log_fh "[$timestamp] $message\n";
    close $log_fh;
}

sub cleanup_tmp_files {
    my $dir = shift;
    opendir(my $dh, $dir) or die "Cannot open directory $dir: $!\n";
    while (my $file = readdir($dh)) {
        next if ($file eq '.' or $file eq '..');
        if ($file =~ /^tmp_/) {
            my $path = File::Spec->catfile($dir, $file);
            if (-d $path) {
                remove_tree($path) or warn "Could not remove directory $path: $!";
            } else {
                unlink($path) or warn "Could not remove file $path: $!";
            }
            print "Removed $path\n";
        }
    }
    closedir($dh);
}

__END__

=head1 NAME

decode_gfa.pl - Process VCF files for genomic feature analysis

=head1 SYNOPSIS

decode_gfa.pl --vcf <vcf_file> --threads <num_threads> [options]

Options:
  --vcf        Input VCF file (required)
  --threads    Number of threads for parallel processing (required)
  --out        Output directory (default: VCF file directory)
  --debug      Enable debug logging
  --help       Show this help message
  --version    Show version number

=head1 DESCRIPTION

This script processes VCF files to extract and analyze genomic features.
It handles both single and multiple alleles, with support for parallel
processing and refined analysis of complex variants.

=head1 AUTHOR

AGIS Zhou-Lab

=cut
