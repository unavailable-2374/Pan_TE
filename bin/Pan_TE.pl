#!/usr/bin/perl
use strict;
use Cwd;
use Getopt::Long;
use threads;
use POSIX qw(:sys_wait_h ceil floor);
use Cwd qw/abs_path getcwd cwd/;
use File::Basename;

my $command_line_geta = join " ", @ARGV;
my $dirname = dirname($0);
$dirname =~ s/\/bin$//;

my $usage = <<USAGE;
Usage:
    perl $0 [options]

For example:
    perl $0 --genome genome.fasta --cpu 80

Parameters::
[General]

    --ref <string>     Required
    genome file in Fasta format.

[other]
    --GFA <string> default:NA
    Pan-genome file in GFA format.

    --out <string>    default: .
    the work dir.

    -M <int>    
    memory limit (in MB) for the program, default 0; 0 for unlimitted;

    --model <string>
    P or M or F or O. P:Plants, M:Metazoans, F:Fungi, and O: Others.

    --model_dir <string>
    Provide model_dir that could be downloaded from website (optional requirements). 

    --hmmscan <int>
    path to hmmscan

    --cpu <int>    default: 4
    the number of threads, preferably in multiples of 4.

    --style <string>
    PGGB or Minigraph-cactus

    --sensitive defaut:normal
    Sensitivity represents a parameter when merging genomes; 
    the more sensitive is more cautious for merging, 
    the longer the procedure takes and the more TEs may be obtained.
    super,hight,normal,low

    --help|-h Display this help info

Version: 1.0.0
USAGE

if (@ARGV==0){die $usage}

my ( $hmm_dir, $genome, $lmerSize, $model, $model_dir, $mem, $sensitive, $gn_ls, $CPU, $cmdString, $out, $style, $help);
GetOptions(
    "hmmscan:s" =>\$hmm_dir,
    "help|h" =>\$help,
    "genome:s" => \$genome,
    "cpu:s" => \$CPU,
    "out:s" => \$out,
    "GFA:s" => \$gn_ls,
    "sensitive:s" => \$sensitive,
    "M:s" => \$mem,
    "model_dir:s" => \$model_dir,
    "model:s" => \$model,
    "style:s" => \$styleq
);

my $Luck = int(rand(9999999))+1;
$out ||= "panTE_".$Luck;
$CPU ||= 4;
$model ||= "P";
$mem ||=0;
$sensitive ||= "normal";

if($genome && $gn_ls){print "please input genome only or list only\n"}

if($genome){
    print STDERR "\n============================================\n";
    print STDERR "RUN Genome Model $out" . "(" . (localtime) . ")" . "\n";
	$genome = abs_path($genome);
	my @filestat = stat ($genome);
	$lmerSize = ceil( ( log( $filestat[7] ) / log( 4 ) + 1 ) );
	mkdir "$out" unless -e "$out";
	chdir $out;
    print STDERR "Step 1: Prepare For Genome $genome" . "(" . (localtime) . ")" . "\n";
	unless (-e "genome.ok") {
        mkdir "genome" unless -e "genome";
        $cmdString="awk '{print \$1}' $genome > genome/genome.fa";
        if (-e "genome/genome.fa"){
            print STDERR "Skip genome.fa existed " . "(" . (localtime) . ")" . "\n";
        }else{
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            `samtools faidx genome/genome.fa`;
			$cmdString="index genome/genome";
			system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        }
        open OUT, ">", "genome.ok" or die $!; close OUT;
    }else{
        print STDERR "Skip Step 1 for the file genome.ok exists\n";
    }
	
	my $dir = getcwd;
	$genome = $dir."/genome";

    print STDERR "Step 2: Looking For LTRs " . "(" . (localtime) . ")" . "\n";

    unless (-e "LTR.ok") {
		mkdir "look4ltrs" unless -e "look4ltrs";
		chdir "look4ltrs";
        my $timeout = 2 * 60 * 60;
        $cmdString="timeout $timeout look4ltrs --fasta ../RECON/repeat/ --out ./ --parallel $CPU > ltr.stat";
        my $thr1 = threads->create(\&run_cmd1, $cmdString);
    }else{
        print STDERR "Skip Step 2 for the file LTR.ok exists\n";
    }

    print STDERR "Step 3: Running RepeatScout " . "(" . (localtime) . ")" . "\n";

    unless (-e "RepeatScout.ok") {
        mkdir "RepeatScout" unless -e "RepeatScout";
        chdir "RepeatScout";
        mkdir "tmp" unless -e "tmp";
        if($filestat[7] > 1000000000){
            $cmdString="build_RS tmp $genome/genome $CPU 600 200000> build.log";
            my $thr2 = threads->create(\&run_cmd2, $cmdString);
        }else{
            $cmdString="build_RS tmp $genome/genome $CPU all > build.log";
            my $thr2 = threads->create(\&run_cmd2, $cmdString);
        }
        chdir "../";
        open OUT, ">", "RepeatScout.ok" or die $!; close OUT;
    }else{
        print STDERR "Skip Step 3 for the file RepeatScout.ok exists\n";
    }
    
    print STDERR "Step 4: Running RECON " . "(" . (localtime) . ")" . "\n";
	my $simple_size=0;
    unless (-e "RECON.ok") {
        mkdir "RECON" unless -e "RECON";
        chdir "RECON";
        my $num = ceil($CPU/4);
        `cat ../RepeatScout/re.fa `
        $cmdString="RepeatMasker -lib ../RepeatScout/consensi.fa -e rmblast -dir repeat -pa $num ../genome/genome.fa -lcambig -gff -s -div 20 > rmod.log";
		unless(-e "repeat/genome.fa.tbl"){
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
		}
		if($filestat[7] > "10000000000"){
            print STDERR "Use Huge genome model, it will spend a lot of time\n";
            print STDERR "Huge genome model will run 6 rounds\n";
	        my @sample_size = (0,240,480,720,1400,2800,5600);
	        for(my $i=1;$i<=6;$i++){
                if(-e "round-$i.ok"){next;}
		        mkdir "round-$i" unless -e "round-$i";
		        chdir "round-$i";
                $cmdString="run_RECON $i $CPU ../repeat/genome.fa.masked $sample_size[$i] 200000";
                system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
                `cat consensi.fa >> ../consensi.fa`;
                chdir "../";
	        }
        }elsif($filestat[7] > "5000000000"){
            print STDERR "Use Big genome model, it will spend a lot of time\n";
            print STDERR "Big genome model will run 5 rounds\n";
	        my @sample_size = (0,120,240,480,860,2800);
	        for(my $i=1;$i<=5;$i++){
                if(-e "round-$i.ok"){next;}
		        mkdir "round-$i" unless -e "round-$i";
		        chdir "round-$i";
                $cmdString="run_RECON $i $CPU ../repeat/genome.fa.masked $sample_size[$i] 200000";
                system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
                `cat consensi.fa >> ../consensi.fa`;
                chdir "../";
	        }            
        }elsif($filestat[7] > "1000000000"){
            print STDERR "Use middle genome model, it will spend a lot of time\n";
            print STDERR "middle genome model will run 5 rounds\n";
	        my @sample_size = (0,80,120,240,480,860);
	        for(my $i=1;$i<=5;$i++){
                if(-e "round-$i.ok"){next;}
		        mkdir "round-$i" unless -e "round-$i";
		        chdir "round-$i";
                $cmdString="run_RECON $i $CPU ../repeat/genome.fa.masked $sample_size[$i] 100000";
                system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
                `cat consensi.fa >> ../consensi.fa`;
                chdir "../";
	        }            
        }elsif($filestat[7] > "500000000"){
            print STDERR "Use middle genome model, it will spend a lot of time\n";
            print STDERR "middle genome model will run 4 rounds\n";
	        my @sample_size = (0,80,120,240,360);
	        for(my $i=1;$i<=4;$i++){
                if(-e "round-$i.ok"){next;}
		        mkdir "round-$i" unless -e "round-$i";
		        chdir "round-$i";
                $cmdString="run_RECON $i $CPU ../repeat/genome.fa.masked $sample_size[$i] 80000";
                system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
                `cat consensi.fa >> ../consensi.fa`;
                chdir "../";
	        }            
        }else{
            print STDERR "Use small genome model, it will spend a lot of time\n";
            if($filestat[7] > "300000000"){
                print STDERR "small genome model will run 3 rounds\n";
                my @sample_size = (0,80,120,240);
                for(my $i=1;$i<=3;$i++){
                    if(-e "round-$i.ok"){next;}
                    mkdir "round-$i" unless -e "round-$i";
                    chdir "round-$i";
                    $cmdString="run_RECON $i $CPU ../repeat/genome.fa.masked $sample_size[$i] 40000";
                    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
                    `cat consensi.fa >> ../consensi.fa`;
                    chdir "../";
                }  
            }else{
                print STDERR "small genome model will run 1 rounds\n";
                my $size = length($filestat[7]) - 1;
                my $len = 0 x $size;
                $size = "1".$len;
                unless(-e "round-1.ok"){
                mkdir "round-1" unless -e "round-1";
                chdir "round-1";
                $cmdString="run_RECON 1 $CPU ../repeat/genome.fa.masked $size";
                system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
                `cat consensi.fa >> ../consensi.fa`;
                chdir "../"; }
            }          
        }
		chdir "../"; 
		open OUT, ">", "RECON.ok" or die $!; close OUT;
    }else{
		print STDERR "Skip Step 4 for the file RECON.ok exists\n";
	}
    
    print STDERR "Step 4: Combining For results " . "(" . (localtime) . ")" . "\n";

    unless (-e "Combine.ok") {
    mkdir "Combine" unless -e "Cmbine";
    chdir "Combine";

    $cmdString="cat ../look4ltrs/consensi.fa ../RepeatScout/consensi.fa ../RECON/consensi.fa > raw.fa";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    $cmdString="cd-hit-est -i raw.fa -o raw_TEs.fa -c 0.8 -aS 0.8 -M 0 -d 0 -T 0 -n 5 -g 1 -b 500 -G 0 -A 80 -l 30";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    $cmdString="sh run_Classifier";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    }
    print STDERR "\n============================================\n";
    print STDERR "Finish run Pan_TE" . "(" . (localtime) . ")" . "\n";
}

if($gn_ls){
    print STDERR "\n============================================\n";
    print STDERR "RUN Graph-Genome Model $out" . "(" . (localtime) . ")" . "\n";
	$gn_ls = abs_path($gn_ls);
	mkdir "$out" unless -e "$out";
	chdir $out;
    my $paths = abs_path($out);
    print STDERR "Step 1: Process For Graph-Genome $gn_ls" . "(" . (localtime) . ")" . "\n";
	unless (-e "genome.ok") {
        mkdir "genome" unless -e "genome";
        if($style eq "PGGB"){
            $cmdString="Process_For_PGGB $gn_ls $paths $min_length $max_length $core_min_length $core_max_length $core_uniq_length";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            $genome = $gn_l."/pggb.fa";
        }else{
            $cmdString="Process_For_MC $gn_ls $paths $min_length $max_length $core_min_length $core_max_length $core_uniq_length";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            $genome = $gn_l."/mc.fa";
        }
        my @filestat = stat ($genome);
        $lmerSize = ceil( ( log( $filestat[7] ) / log( 4 ) + 1 ) );
        open OUT, ">", "genome.ok" or die $!; close OUT;
    }else{
        print STDERR "Skip Step 1 for the file genome.ok exists\n";
    }
}


sub run_ltr {
    my ($cmd) = @_;
    system($cmd) == 0 or die "failed to execute: $cmd\n";

    $cmdString="grep -v \"chrom\" Rtr/genome.rtr | awk '\$16>0.9&&\$6!=\"NA\"{print \$1\"\t\"\$3\"\t\"\$6}' | sort | uniq >  LTR.bed";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n"; 

    $cmdString="grep -v \"chrom\" Rtr/genome.rtr | awk '\$16>0.9&&\$6==\"NA\"{print \$1\"\t\"\$3\"\t\"\$4}' | sort | uniq >> LTR.bed";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    $cmdString="seqtk subseq $genome/genome.fa LTR.bed > LTR.fa";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    mkdir "tmp" unless -e "tmp";
    $cmdString="build_LL tmp consensi.fa $CPU > build.log";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    `renameTE LTR consensi.fa re.fa`;
    chdir "../";
    open OUT, ">", "LTR.ok" or die $!; close OUT;
    print "LTR completed\n";
}

sub run_rs {
    my ($cmd) = @_;
    system($cmd) == 0 or die "failed to execute: $cmd\n";
    `renameTE RS consensi.fa re.fa`;
    print "RepeatScout completed\n";
}
