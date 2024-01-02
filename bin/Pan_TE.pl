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
    perl $0 --genome genome.fasta --cpu 80 --ClassifyTE_dir Your_Path_To_ClassifyTE

Parameters::
[General]

    --ref <string>     Required
    genome file in fasta format.

    --ClassifyTE_dir <string>
    Provide ClassifyTE's dir for run classification.

[other]
    --list <string> default:NA
    path file for GFA, see gfa.list for formate.

    --out <string>    default: .
    the work dir.

    -M <int>    
    memory limit (in MB) for the program, default 0; 0 for unlimitted;

    --flag <string>
    PGGB or Minigraph-Cactus.

    --cpu <int>    default: 4
    the number of threads, preferably in multiples of 4.

    --fragment_size <int>    default: 40000
    the length for fragment.

    --help|-h Display this help info

Version: 1.0.0
USAGE

if (@ARGV==0){die $usage}

my ($genome, $style, $model_dir, $mem, $gn_ls, $CPU, $cmdString, $out, $fragment_size, $lmerSize, $help);
GetOptions(
    "help|h" =>\$help,
    "genome:s" => \$genome,
    "cpu:s" => \$CPU,
    "out:s" => \$out,
    "list:s" => \$gn_ls,
    "M:s" => \$mem,
    "ClassifyTE_dir:s" => \$model_dir,
    "flag:s" => \$style,
	"fragment_size" => \$fragment_size
);

my $Luck = int(rand(9999999))+1;
$out ||= "panTE_".$Luck;
$CPU ||= 4;
$mem ||=0;
$fragment_size ||="40000";
my @filestat;
my $thr1;
my $thr2;
my $dir;

if($genome && $gn_ls){print "please input genome only or list only\n"}

if($genome){ 
    print STDERR "\n============================================\n";
    print STDERR "RUN Genome Model $out" . "(" . (localtime) . ")" . "\n";
	$genome = abs_path($genome);
	@filestat = stat ($genome);
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
}if($gn_ls){
    print STDERR "\n============================================\n";
    print STDERR "RUN Graph-Genome Model $out" . "(" . (localtime) . ")" . "\n";
	$gn_ls = abs_path($gn_ls);
	mkdir "$out" unless -e "$out";
	chdir $out;
    my ($min_length,$max_length);
    my $core_min_length = 30;
    my $core_uniq_length = 0.4;
    print STDERR "Step 1: Process For Graph-Genome $gn_ls" . "(" . (localtime) . ")" . "\n";
	unless (-e "genome.ok") {
        if($filestat[7] > "10000000000"){
            $min_length = 100000;
            $max_length = 200000;
        }elsif($filestat[7] > "5000000000"){
            $min_length = 80000;
            $max_length = 200000;
        }elsif($filestat[7] > "1000000000"){
            $min_length = 50000;
            $max_length = 100000;
        }else{
            $min_length = 40000;
            $max_length = 80000;
        }
        mkdir "genome" unless -e "genome";
        chdir "genome";
        if($style eq "PGGB"){
            $cmdString="Process_For_PGGB $gn_ls ./ $min_length $max_length $core_min_length $min_length $core_uniq_length $CPU";
            print $cmdString."\n";  
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        }else{
            $cmdString="Process_For_MC $gn_ls ./ $min_length $max_length $core_min_length $min_length $core_uniq_length $CPU";
            print $cmdString."\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        }
        $genome = $gn_ls."/genome.fa";
        @filestat = stat ($genome);
        $lmerSize = ceil( ( log( $filestat[7] ) / log( 4 ) + 1 ) );
        open OUT, ">", "genome.ok" or die $!; close OUT;
    }else{
        print STDERR "Skip Step 1 for the file genome.ok exists\n";
    }
}
	
$dir = getcwd;
$genome = $dir."/genome";

print STDERR "Step 2: Looking For LTRs " . "(" . (localtime) . ")" . "\n";

unless (-e "LTR.ok") {
    chdir $dir;
    mkdir "look4ltrs" unless -e "look4ltrs";
    my $timeout = 2 * 60 * 60;
    $cmdString="timeout $timeout look4ltrs --fasta $genome --out ./ --parallel $CPU > ltr.stat";
    $thr1 = threads->create(\&run_ltr, $cmdString);
}else{
    print STDERR "Skip Step 2 for the file LTR.ok exists\n";
}

print STDERR "Step 3: Running RepeatScout " . "(" . (localtime) . ")" . "\n";

unless (-e "RepeatScout.ok") {
    chdir $dir;
    mkdir "RepeatScout" unless -e "RepeatScout";
    chdir "RepeatScout";
    mkdir "tmp" unless -e "tmp";
    if($filestat[7] > 1000000000){
        $cmdString="build_RS $dir/RepeatScout/tmp $genome/genome $CPU 600 > build.log";
        $thr2 = threads->create(\&run_rs, $cmdString);
    }else{
        $cmdString="build_RS $dir/RepeatScout/tmp $genome/genome $CPU all > build.log";
        $thr2 = threads->create(\&run_rs, $cmdString);
    }
    if($thr1){$thr1->join();$thr2->join();}
    else{$thr2->join();}
}else{
    print STDERR "Skip Step 3 for the file RepeatScout.ok exists\n";
}

print STDERR "Step 4: Running RECON " . "(" . (localtime) . ")" . "\n";
my $simple_size=0;
unless (-e "RECON.ok") {
    chdir $dir;
    mkdir "RECON" unless -e "RECON";
    chdir "RECON";
    my $num = ceil($CPU/4);
    my $sample_time;
    $cmdString="cat ../RepeatScout/consensi.fa ../look4ltrs/consensi.fa > lib.fa";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            `renameTE RS lib.fa re.fa`;
    $cmdString="RepeatMasker -lib re.fa -e rmblast -dir repeat -pa $num ../genome/genome.fa -lcambig -gff -s -div 20 > rmod.log";
    unless(-e "repeat/genome.fa.tbl"){
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        $cmdString="cat lib.fa | filter-stage-2.prl --cat repeat/genome.fa.out > stag2.fa";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        $cmdString="filter_by_stage2 stag2.fa repeat/genome.fa.out";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        $cmdString="bedtools merge -i mask.bed > out.bed";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        $cmdString="bedtools maskfasta -fi ../genome/genome.fa -fo repeat/genome.masked -bed out.bed";
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
            $cmdString="run_RECON $i $CPU ../repeat/genome.masked $i 200000";
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
            my $sz = substr($filestat[7],0,1);
            $size = int($sz.$len)/1000000;
            unless(-e "round-1.ok"){
            mkdir "round-1" unless -e "round-1";
            chdir "round-1";
            $cmdString="run_RECON 1 $CPU ../repeat/genome.masked $size 40000";
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

    $cmdString="sh run_Classifier $model_dir";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
}
print STDERR "\n============================================\n";
print STDERR "Finish run Pan_TE" . "(" . (localtime) . ")" . "\n";


sub run_ltr {
    chdir $dir;
    chdir "look4ltrs";
    my ($cmd) = @_;
    system($cmd) == 0 or die "failed to execute: $cmd\n";
    
    chdir $dir;
    chdir "look4ltrs";
    `uniq $dir/look4ltrs/Bed/genome.bed > $dir/look4ltrs/Bed/LTR.bed`;
    $cmdString="seqtk subseq $genome/genome.fa $dir/look4ltrs/Bed/LTR.bed > $dir/look4ltrs/LTR.fa";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    mkdir "tmp" unless -e "tmp";
    $cmdString="build_LL $dir/look4ltrs/tmp consensi.fa $CPU > build.log";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    open OUT, ">", "$dir/LTR.ok" or die $!; close OUT;
    print "LTR completed\n";
}

sub run_rs {
    my ($cmd) = @_;
    chdir $dir;
    chdir "RepeatScout";
    system($cmd) == 0 or die "failed to execute: $cmd\n";
    print "RS completed\n";
    open OUT, ">", "$dir/RepeatScout.ok" or die $!; close OUT;
}

