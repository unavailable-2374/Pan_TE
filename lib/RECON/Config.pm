package RECON::Config;

use strict;
use warnings;
use Exporter 'import';
use Getopt::Long;
use File::Spec;
use Cwd 'abs_path';
use RECON::Logger;
use RECON::Utils;

our @EXPORT = qw(
    parse_arguments validate_inputs find_input_files
    get_pipeline_config
);

# Default configuration
our %DEFAULT_CONFIG = (
    threads => 4,
    cpu_threads => undef,  # Will be calculated
    fragment_size => 40000,
    use_full_genome => 0,
    memory_limit => undef,
    
    # Genome size thresholds
    dual_track_threshold => 1_000_000_000,  # 1GB
    
    # Sampling parameters
    sample_sizes => [30, 90, 270, 540],  # MB
    max_sampling_rounds => 4,
    
    # RECON parameters
    default_k => 18,
    max_recon_iterations => 50000,
    
    # Resource limits
    max_blast_chunk_size => 100,  # MB
    max_memory_per_process => 8000,  # MB
);


sub parse_arguments {
    my %config = %DEFAULT_CONFIG;
    
    GetOptions(
        'threads=i' => \$config{threads},
        'memory|M=i' => \$config{memory_limit},
        'fragment-size=i' => \$config{fragment_size},
        'full-genome' => \$config{use_full_genome},
        'help|h' => \$config{help},
        'version|v' => \$config{version},
    );
    
    if ($config{help}) {
        print_usage();
        exit 0;
    }
    
    if ($config{version}) {
        print_version();
        exit 0;
    }
    
    # Parse positional arguments
    if (@ARGV < 3) {
        print_usage();
        exit 1;
    }
    
    ($config{threads}, $config{genome_file}, $config{genome_size}) = @ARGV;
    
    # Validate and process inputs
    validate_inputs(\%config);
    
    return \%config;
}

sub validate_inputs {
    my ($config) = @_;
    
    # Validate threads
    if ($config->{threads} < 1 || $config->{threads} > 256) {
        die "Error: threads must be between 1 and 256\n";
    }
    
    # Calculate CPU threads if not specified
    $config->{cpu_threads} = int($config->{threads} / 4) || 1;
    $config->{cpu_threads} = 20 if $config->{cpu_threads} > 20;  # Cap at 20
    
    # Validate genome file
    unless (-f $config->{genome_file}) {
        die "Error: Genome file does not exist: $config->{genome_file}\n";
    }
    
    # Convert to absolute path
    $config->{genome_file} = abs_path($config->{genome_file});
    
    # Validate genome size
    unless ($config->{genome_size} =~ /^\d+$/ && $config->{genome_size} > 0) {
        die "Error: Genome size must be a positive integer\n";
    }
    
    # Check memory constraints
    if ($config->{memory_limit}) {
        if ($config->{memory_limit} < 1000) {
            log_message("WARN", "Low memory limit specified", "memory=$config->{memory_limit}MB");
        }
        
        # Adjust thread count based on memory
        my $memory_per_thread = $config->{memory_limit} / $config->{threads};
        if ($memory_per_thread < 500) {
            my $recommended_threads = int($config->{memory_limit} / 500);
            $recommended_threads = 1 if $recommended_threads < 1;
            
            log_message("WARN", "Insufficient memory per thread", 
                       "current=${memory_per_thread}MB/thread, reducing to $recommended_threads threads");
            $config->{threads} = $recommended_threads;
            $config->{cpu_threads} = int($config->{threads} / 4) || 1;
        }
    }
    
    # Determine processing mode
    my $genome_class = classify_genome_size($config->{genome_size});
    $config->{genome_class} = $genome_class;
    
    if ($config->{use_full_genome} || $config->{genome_size} < $config->{dual_track_threshold}) {
        $config->{processing_mode} = 'single_track';
    } else {
        $config->{processing_mode} = 'dual_track';
    }
    
    log_message("INFO", "Configuration validated", 
                "mode=$config->{processing_mode}, genome_class=$genome_class, " .
                "threads=$config->{threads}, cpu_threads=$config->{cpu_threads}");
}

sub find_input_files {
    my ($config) = @_;
    
    # Find the original genome file for size calculation
    # Based on Pan_TE pipeline structure: out_dir/genome/genome.fa
    # Current working directory is: out_dir/RECON/
    my $main_output_dir = File::Spec->rel2abs("..");
    my $original_genome_file = File::Spec->catfile($main_output_dir, "genome", "genome.fa");
    
    # Use original genome size for classification, but keep masked genome for processing
    if (-f $original_genome_file) {
        my $original_genome_size = -s $original_genome_file;
        $config->{original_genome_file} = $original_genome_file;
        $config->{original_genome_size} = $original_genome_size;
        
        log_message("INFO", "Found original genome file", 
                    "file=$original_genome_file, size=$original_genome_size bytes");
        
        # Re-classify based on original genome size
        my $genome_class = classify_genome_size($original_genome_size);
        $config->{genome_class} = $genome_class;
        
        # Update processing mode based on original size
        if ($config->{use_full_genome} || $original_genome_size < $config->{dual_track_threshold}) {
            $config->{processing_mode} = 'single_track';
        } else {
            $config->{processing_mode} = 'dual_track';
        }
        
        log_message("INFO", "Genome classification updated", 
                    "original_size=$original_genome_size, genome_class=$genome_class, mode=$config->{processing_mode}");
    } else {
        log_message("WARN", "Original genome file not found", 
                    "expected=$original_genome_file, using masked genome size for classification");
    }
    
    # Find masked genome and BED files
    my $genome_file = $config->{genome_file};
    my $genome_dir = File::Spec->rel2abs(File::Spec->updir($genome_file));
    
    # Look for BED exclusion files
    my @bed_files;
    
    # Common BED file patterns
    my @bed_patterns = (
        "*.bed",
        "tmp/*.bed", 
        "../tmp/*.bed",
        "exclusion*.bed",
        "mask*.bed"
    );
    
    for my $pattern (@bed_patterns) {
        my @files = glob($pattern);
        push @bed_files, grep { -s $_ } @files;
    }
    
    # Remove duplicates and get absolute paths
    my %seen;
    @bed_files = grep { !$seen{$_}++ } map { abs_path($_) } @bed_files;
    
    $config->{bed_files} = \@bed_files;
    
    log_message("INFO", "Input files identified", 
                "genome=$genome_file, bed_files=" . scalar(@bed_files));
    
    return $config;
}

sub get_pipeline_config {
    my ($genome_size, $threads) = @_;
    
    my %pipeline_config = %DEFAULT_CONFIG;
    
    # Adjust configuration based on genome size
    my $genome_class = classify_genome_size($genome_size);
    
    if ($genome_class eq 'HUGE') {
        $pipeline_config{max_blast_chunk_size} = 50;  # Smaller chunks
        $pipeline_config{sample_sizes} = [20, 60, 180, 360];  # Smaller samples
    } elsif ($genome_class eq 'LARGE') {
        $pipeline_config{max_blast_chunk_size} = 75;
        $pipeline_config{sample_sizes} = [25, 75, 225, 450];
    } elsif ($genome_class eq 'SMALL' || $genome_class eq 'TINY') {
        $pipeline_config{sample_sizes} = [50, 150, 400];  # Fewer rounds
        $pipeline_config{max_sampling_rounds} = 3;
    }
    
    # Adjust for thread count
    if ($threads <= 4) {
        $pipeline_config{max_blast_chunk_size} = 200;  # Larger chunks for fewer threads
    } elsif ($threads >= 32) {
        $pipeline_config{max_blast_chunk_size} = 25;   # Smaller chunks for many threads
    }
    
    return \%pipeline_config;
}

sub print_usage {
    print <<EOF;
Usage: run_RECON_advanced threads genome_file genome_size

Advanced RECON pipeline with dual-track processing for transposable element discovery.

Positional Arguments:
    threads        Number of processing threads (1-256)
    genome_file    Input genome file in FASTA format (required)
    genome_size    Genome size in bytes (required)

Options:
    --memory, -M SIZE      Memory limit in MB
    --fragment-size SIZE   Genome fragment size (default: 40000)
    --full-genome         Force single-track full genome processing
    --help, -h            Show this help message
    --version, -v         Show version information

Examples:
    # Standard dual-track processing
    run_RECON_advanced 80 genome.fa 2500000000
    
    # Memory-limited processing
    run_RECON_advanced 16 genome.fa 1000000000 --memory 8000
    
    # Force full-genome mode for small genomes
    run_RECON_advanced 8 small_genome.fa 500000000 --full-genome

Output:
    The pipeline creates consensus TE libraries in the current directory:
    - masked_track/consensi.fa    (high-quality consensus)
    - sampling_track/round_N/     (sampling results by round)
    - RECON.ok                    (completion marker)

For more information, see the Pan_TE documentation.
EOF
}

sub print_version {
    print <<EOF;
Pan_TE Advanced RECON Pipeline v2.0.0

A comprehensive pipeline for transposable element detection and annotation
supporting both linear and graph genomes with adaptive dual-track processing.

Copyright (C) 2025 Pan_TE Development Team
Licensed under GPL v3.0
EOF
}

1;

__END__

=head1 NAME

RECON::Config - Configuration management for RECON pipeline

=head1 DESCRIPTION

Handles command-line argument parsing, input validation, and pipeline
configuration for different genome sizes and processing modes.

=cut