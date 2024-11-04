#!/usr/bin/perl
use strict;
use warnings;

open(LIST,$ARGV[0]);
open(OUT,">>SV.fa");
while(<LIST>){
    chomp;
    if(/^#/){next;}
    my @arr = split(/\t/);
    if(length($arr[4]) > 50){
        if($arr[4] =~ /,/){
            my @array = split(/,/,$arr[4]);
            my $num =0;
            open(TMP,">tmp.fa");
            for(my $i = 0;$i<=$#array;$i++){
                if(length ($array[$i])>50){
                    print TMP ">".$i."\n".$array[$i]."\n";
                    $num++;
                }
            }if($num > 1){
                system("cd-hit-est -i tmp.fa -o tmp.fasta -c 0.8 -aS 0.8 -M 0 -d 0 -T 0 -n 5 -g 1 -b 500 -G 0 -A 80 -l 30 &>> tmp.log");
                `cat tmp.fasta >> SVs.fa`;
            }elsif($num==1){`cat tmp.fa >> SVs.fa`;}
            unlink "tmp.fa";
        }else{print OUT ">SV\n".$arr[4]."\n";}
    }
}
