#!/usr/bin/perl

use warnings;
use strict;

open(IN,$ARGV[0]);
open(OUT,">TEs.fa");

my $num = 0;
while (<IN>) {
    chomp;
    if(/>/){
        $num = sprintf("%08d", $num);
        if(/ClassI_/ || /ClassII_/ || /ClassIII_/){
            s/.*ClassIII_/>TE$num#/;
            s/.*ClassII_/>TE$num#/;
            s/.*ClassI_/>TE$num#/;
            s/_/\//;
        }else{
            $_=">TE".$num."#unknow";
        }s/TE/TE_/;
        print OUT $_."\n";
        $num++;
    }else{print OUT $_."\n";}
}
