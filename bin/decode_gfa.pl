#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use autodie;  # Automatically die if file operations fail
use File::Path qw(make_path);
use constant {
    MIN_SEQ_LENGTH => 50,
    CDHIT_IDENTITY => 0.8,
    CDHIT_COVERAGE => 0.8,
};

our $VERSION = '1.0.0';

my %opts = parse_arguments();

make_path('tmp') unless -d 'tmp';

process_vcf($opts{vcf}, $opts{threads});

sub parse_arguments {
    my %opts;
    GetOptions(
        'vcf=s'     => \$opts{vcf},
        'threads=i' => \$opts{threads},
        'help|h'    => \$opts{help},
        'version|v' => \$opts{version},
    ) or pod2usage(2);

    pod2usage(1) if $opts{help};
    if ($opts{version}) {
        print "decode_gfa.pl version $VERSION\n";
        exit(0);
    }

    pod2usage("Error: VCF file is required") unless $opts{vcf};
    pod2usage("Error: Number of threads is required") unless $opts{threads};
    
    die "Error: Input file '$opts{vcf}' does not exist\n" unless -e $opts{vcf};
    die "Error: Input file '$opts{vcf}' is not readable\n" unless -r $opts{vcf};
    
    die "Error: Number of threads must be positive\n" unless $opts{threads} > 0;

    return %opts;
}

sub process_vcf {
    my ($vcf_file, $threads) = @_;
    
    my $base_name = basename($vcf_file, ".vcf");
    my $output_file = "tmp/$base_name.fa";
    my $temp_file = "tmp/$base_name.tmp.fa";
    
    open(my $vcf_fh, '<', $vcf_file) or die "Cannot open $vcf_file: $!\n";
    open(my $out_fh, '>>', $output_file) or die "Cannot open $output_file: $!\n";

    while (my $line = <$vcf_fh>) {
        chomp $line;
        next if $line =~ /^#/;  # Skip header lines
        
        my @fields = split(/\t/, $line);
        next unless @fields >= 5;  # Ensure minimum fields exist
        
        my $alt_allele = $fields[4];
        next unless length($alt_allele) > MIN_SEQ_LENGTH;
        
        if ($alt_allele =~ /,/) {
            process_multiple_alleles($alt_allele, $base_name, $threads, $temp_file);
        } else {
            print_sequence($out_fh, "SV", $alt_allele);
        }
    }
    
    close $vcf_fh;
    close $out_fh;
}

sub process_multiple_alleles {
    my ($alt_allele, $base_name, $threads, $temp_file) = @_;
    
    my @alleles = split(/,/, $alt_allele);
    my $long_alleles = 0;
    
    open(my $tmp_fh, '>', $temp_file) or die "Cannot open $temp_file: $!\n";
    for my $i (0 .. $#alleles) {
        if (length($alleles[$i]) > MIN_SEQ_LENGTH) {
            print_sequence($tmp_fh, $i, $alleles[$i]);
            $long_alleles++;
        }
    }
    close $tmp_fh;
    
    if ($long_alleles > 1) {
        my $output_fasta = "$temp_file.out";
        my $cmd = join(" ",
            "Refiner_for_Graph",
             $temp_file,
             $output_fasta,
            "--distance-threshold", 0.7,
	    "-t", $threads
	    #"2>>", "tmp/refine.log"
        );
        
        system($cmd) == 0
            or die " Feature Extraction Failure: $?\n";
        
        system("cat $output_fasta >> tmp/$base_name.fa") == 0
            or die "Failed to append Refine results: $?\n";
            
        unlink $output_fasta;
	unlink glob ('tmp/*.fa.n*');
    } elsif ($long_alleles == 1) {
        system("cat $temp_file >> tmp/$base_name.fa") == 0
            or die "Failed to append single sequence: $?\n";
    }
    
    unlink $temp_file;
}

sub print_sequence {
    my ($fh, $id, $sequence) = @_;
    print $fh ">$id\n$sequence\n";
}

__END__

=head1 NAME

decode_gfa.pl - Process VCF files and extract sequences

=head1 SYNOPSIS

decode_gfa.pl --vcf <vcf_file> --threads <num_threads> [options]

Options:
  --vcf        Input VCF file (required)
  --threads    Number of threads for cd-hit-est (required)
  --help       Show this help message
  --version    Show version number

=head1 DESCRIPTION

This script processes VCF files to extract sequences and perform clustering using cd-hit-est
for sequences longer than 50 bases.

=head1 AUTHOR

Lab Name

=cut
