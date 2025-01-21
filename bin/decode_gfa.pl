#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use autodie;
use File::Path qw(make_path);
use File::Which;
use File::Temp qw(tempdir);

use constant {
    MIN_SEQ_LENGTH => 50,
    CLUSTER_IDENTITY => 0.8,
};

our $VERSION = '1.0.0';

my $mmseqs_path = which('mmseqs') 
    or die "Error: mmseqs is not installed or not in PATH\n";

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
    
    open(my $vcf_fh, '<', $vcf_file);
    open(my $out_fh, '>>', $output_file);

    while (my $line = <$vcf_fh>) {
        chomp $line;
        next if $line =~ /^#/;
        
        my @fields = split(/\t/, $line);
        next unless @fields >= 5;
        
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
    
    open(my $tmp_fh, '>', $temp_file);
    for my $i (0 .. $#alleles) {
        if (length($alleles[$i]) > MIN_SEQ_LENGTH) {
            print_sequence($tmp_fh, "seq_$i", $alleles[$i]);
            $long_alleles++;
        }
    }
    close $tmp_fh;
    
    if ($long_alleles > 1) {
        my $tmp_dir = tempdir(CLEANUP => 1);
        my $output_fasta = "$temp_file.out";
        my $log_file = "tmp/mmseqs.log";
        
        my $db_prefix = "$tmp_dir/seqdb";
        system("mmseqs createdb $temp_file $db_prefix >>$log_file 2>&1") == 0
            or die "MMseqs2 createdb failed: $?\n";
    
        my $clu_prefix = "$tmp_dir/clu";
        my $cmd = join(" ",
            "mmseqs cluster",
            $db_prefix,
            $clu_prefix,
            $tmp_dir,
            "--min-seq-id", "0.8",          
            "--threads", $threads,           
            "--cluster-mode", 0,            
            
            "--cov-mode", 2,                
            "--min-coverage", "0.8",        
            
            "--alignment-mode", 3,          
            "--max-seqs", 300,              
            "--min-aln-len", 30,            
            
            "--k-score", "seq:5",           

            "--max-seq-len", 32768,         
            "--split-memory-limit", "0",    
            
            "--remove-tmp-files", 1,        
            "--filter-hits", 1,             
            
            ">>$log_file 2>&1"
        );
        
        system($cmd) == 0
            or die "MMseqs2 clustering failed: $?\n";
        
        system("mmseqs createseqfiledb $db_prefix $clu_prefix $db_prefix"."_rep >>$log_file 2>&1") == 0
            or die "MMseqs2 createseqfiledb failed: $?\n";
        
        system("mmseqs result2flat $db_prefix $db_prefix $db_prefix"."_rep $output_fasta --use-fasta-header >>$log_file 2>&1") == 0
            or die "MMseqs2 result2flat failed: $?\n";
        
        system("cat $output_fasta >> tmp/$base_name.fa") == 0
            or die "Failed to append MMseqs2 results: $?\n";
            
        unlink $output_fasta;
        
    } elsif ($long_alleles == 1) {
        system("cat $temp_file >> tmp/$base_name.fa") == 0
            or die "Failed to append single sequence: $!\n";
    }
    
    unlink $temp_file;
}

sub print_sequence {
    my ($fh, $id, $sequence) = @_;
    print $fh ">$id\n$sequence\n";
}

__END__

=head1 NAME

decode_gfa.pl - Process VCF files and extract sequences using MMseqs2 for clustering

=head1 SYNOPSIS

decode_gfa.pl --vcf <vcf_file> --threads <num_threads> [options]

Options:
  --vcf        Input VCF file (required)
  --threads    Number of threads for MMseqs2 (required)
  --help       Show this help message
  --version    Show version number

=head1 DESCRIPTION

This script processes VCF files to extract sequences and perform clustering using MMseqs2
for sequences longer than 50 bases. MMseqs2 is used for its superior speed, sensitivity
and memory efficiency compared to other clustering tools.

Key features of MMseqs2 clustering:
- Very fast: 10-100x faster than CD-HIT
- Memory efficient: Uses compressed index structures
- Sensitive: Uses profile-based clustering
- Scalable: Can handle millions of sequences

=head1 REQUIREMENTS

- MMseqs2 must be installed and available in PATH

=head1 AUTHOR

Lab Name

=cut


sub process_multiple_alleles {
    my ($alt_allele, $base_name, $threads, $temp_file) = @_;
    
    my @alleles = split(/,/, $alt_allele);
    my $long_alleles = 0;
    
    open(my $tmp_fh, '>', $temp_file);
    for my $i (0 .. $#alleles) {
        if (length($alleles[$i]) > MIN_SEQ_LENGTH) {
            print_sequence($tmp_fh, "seq_$i", $alleles[$i]);
            $long_alleles++;
        }
    }
    close $tmp_fh;
    
    if ($long_alleles > 1) {
        my $tmp_dir = "tmp/mmseqs_tmp";
        make_path($tmp_dir) unless -d $tmp_dir;
        
        my $output_fasta = "$temp_file.out";
        my $log_file = "tmp/mmseqs.log";
        
        my $db_prefix = "$tmp_dir/seqdb";
        print STDERR "Creating MMseqs2 database...\n";
        my $createdb_cmd = "mmseqs createdb $temp_file $db_prefix";
        print STDERR "Running: $createdb_cmd\n";
        system("$createdb_cmd >>$log_file 2>&1");
        
        if ($? != 0) {
            die "MMseqs2 createdb failed with status $?: $!\n" .
                "Check log file: $log_file\n";
        }
   
        print STDERR "Running MMseqs2 clustering...\n";
        my $clu_prefix = "$tmp_dir/clu";
        my $cmd = join(" ",
            "mmseqs cluster",
            $db_prefix,
            $clu_prefix,
            $tmp_dir,
            "--min-seq-id", "0.8",
            "--threads", $threads,
            "--cluster-mode", 0,
            "--cov-mode", 2,
            "-c", "0.8"
        );
        
        print STDERR "Running: $cmd\n";
        my $cluster_exit = system("$cmd >>$log_file 2>&1");
        
        if ($cluster_exit != 0) {
            print STDERR "Contents of log file:\n";
            system("cat $log_file");
            die "MMseqs2 clustering failed with status " . ($cluster_exit >> 8) . 
                "\nCommand was: $cmd\n";
        }
        
        print STDERR "Extracting cluster representatives...\n";
        my $createseqdb_cmd = "mmseqs createseqfiledb $db_prefix $clu_prefix $db_prefix"."_rep";
        print STDERR "Running: $createseqdb_cmd\n";
        system("$createseqdb_cmd >>$log_file 2>&1");
        
        if ($? != 0) {
            die "MMseqs2 createseqfiledb failed: $?\nCheck log file: $log_file\n";
        }
        
        print STDERR "Converting to FASTA format...\n";
        my $result2flat_cmd = "mmseqs result2flat $db_prefix $db_prefix $db_prefix"."_rep $output_fasta --use-fasta-header";
        print STDERR "Running: $result2flat_cmd\n";
        system("$result2flat_cmd >>$log_file 2>&1");
        
        if ($? != 0) {
            die "MMseqs2 result2flat failed: $?\nCheck log file: $log_file\n";
        }

        if (-e $output_fasta && -s $output_fasta) {
            system("cat $output_fasta >> tmp/$base_name.fa") == 0
                or die "Failed to append MMseqs2 results: $!\n";
        } else {
            die "Output FASTA file is missing or empty: $output_fasta\n";
        }
        
        unlink $output_fasta;
        
    } elsif ($long_alleles == 1) {
        system("cat $temp_file >> tmp/$base_name.fa") == 0
            or die "Failed to append single sequence: $!\n";
    }
    
    unlink $temp_file;
}