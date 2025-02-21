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

# Constants
use constant {
    MIN_SEQ_LENGTH => 50,
    DEFAULT_THREADS => 4,
    MAX_RETRIES => 3,
    TMP_DIR => '.',
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
    chdir($config{out_dir}) or die "Cannot change directory to $config{out_dir}: $!\n";
    print "DEBUG: Changed working directory to $config{out_dir}\n";
}

sub process_data {
    my $vcf_entries = read_vcf_file($config{vcf_file});
    my $output_file = File::Spec->catfile($config{out_dir}, "processed.fa");
    process_vcf_entries($vcf_entries, $output_file);
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
    make_path(TMP_DIR) unless -d TMP_DIR;
}

sub setup_logging {
    my $log_file = "$config{out_dir}/decode_gfa.log";
    open my $log_fh, '>>', $log_file;
    print $log_fh strftime("%Y-%m-%d %H:%M:%S", localtime), 
          " Starting decode_gfa.pl v$VERSION\n";
    close $log_fh;
}

###############################################################################
# VCF Processing Functions
###############################################################################

sub read_vcf_file {
    my $file = shift;
    my @entries;
    
    log_message("Reading VCF file: $file");
    open my $fh, '<', $file;
    
    my $total_lines = 0;
    while (<$fh>) {
        $total_lines++ unless /^#/;
    }
    
    seek($fh, 0, 0);
    
    select((select(STDERR), $| = 1)[0]);
    
    my $processed = 0;
    while (my $line = <$fh>) {
        next if $line =~ /^#/;  # Skip headers
        chomp $line;
        
        $processed++;
        if ($processed % 100 == 0) {
            my $percent = int(($processed / $total_lines) * 100);
            print STDERR "PROGRESS:$percent\n";
        }
        
        my @fields = split(/\t/, $line);
        next unless @fields >= 5;  # Ensure minimum fields
        
        push @entries, {
            chrom => $fields[0],
            pos   => $fields[1],
            ref   => $fields[3],
            alt   => $fields[4],
        };
    }
    
    print STDERR "PROGRESS:100\n";
    
    close $fh;
    log_message("Read " . scalar(@entries) . " VCF entries");
    return \@entries;
}

sub process_vcf_entries {
    my $entries = shift;
    my $base_name = basename($config{vcf_file}, ".vcf");
    my $output_file = TMP_DIR . "/$base_name.fa";
    my $temp_file = TMP_DIR . "/$base_name.tmp.fa";
    
    # Force immediate output flush
    $| = 1;
    select((select(STDERR), $| = 1)[0]);
    
    open my $out_fh, '>>', $output_file;
    
    my $total_entries = scalar(@$entries);
    my $processed = 0;
    my $last_percent = -1;
    
    foreach my $entry (@$entries) {
        my $alt_allele = $entry->{alt};
        $processed++;
        
        # Calculate and display progress
        my $percent = int(($processed / $total_entries) * 100);
        if ($percent != $last_percent) {
            print STDERR "PROGRESS:$percent\n";  # 特殊格式的进度信息
            $last_percent = $percent;
        }
        
        next unless length($alt_allele) > MIN_SEQ_LENGTH;
        
        if ($alt_allele =~ /,/) {
            process_multiple_alleles($alt_allele, $base_name, $temp_file);
        } else {
            print_sequence($out_fh, "SV", $alt_allele);
        }
    }
    
    close $out_fh;
}

sub process_multiple_alleles {
    my ($alt_allele, $base_name, $temp_file) = @_;
    my @alleles = split(/,/, $alt_allele);
    my $long_alleles = 0;

    my $tmp_dir = TMP_DIR;
    make_path($tmp_dir) unless -d $tmp_dir;
    
    open my $tmp_fh, '>', $temp_file;
    
    for my $i (0 .. $#alleles) {
        if (length($alleles[$i]) > MIN_SEQ_LENGTH) {
            print_sequence($tmp_fh, $i, $alleles[$i]);
            $long_alleles++;
        }
    }
    
    close $tmp_fh;
    
    if ($long_alleles > 1) {
        refine_alleles($temp_file, $base_name);
    } elsif ($long_alleles == 1) {
        system("cat $temp_file >> " . TMP_DIR . "/$base_name.fa") == 0
            or die "Failed to append single sequence: $!\n";
    }
    
    unlink $temp_file;
}

sub refine_alleles {
    my ($temp_file, $base_name) = @_;

    $config{out_dir} = abs_path($config{out_dir});

    chdir($config{out_dir}) or die "Cannot change directory to $config{out_dir}: $!\n";
    
    my $new_temp_file = File::Spec->catfile($config{out_dir}, "$base_name.tmp.fa");
    my $output_fasta  = File::Spec->catfile($config{out_dir}, "$base_name.tmp.fa.out");

    my $cmd = join(" ",
        "Refiner_for_Graph",
        $new_temp_file,
        $output_fasta,
        "--distance-threshold", "0.8",
        "-t", $config{threads},
	"-v"
    );
    
    system($cmd) == 0
        or die "Feature Extraction Failure: $?\n";

    my $final_fa = File::Spec->catfile($config{out_dir}, "$base_name.fa");
    system("cat $output_fasta >> $final_fa")
        == 0 or die "Failed to append Refine results: $!\n";
}

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
