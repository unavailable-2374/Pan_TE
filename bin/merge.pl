#!/usr/bin/perl

use warnings;
use strict;

open(IN,$ARGV[0]);
open(OUT,">TE_lib.fa");

my %hash;
my $flag = 0;
my $key;
my @arr;
while (<IN>) {
    chomp;
    if(/>/ && /\//){
        $key = $_;
        $flag = 1;
    }elsif (/>/ && $_ !~ /\//) {
        $flag = 0;
    }elsif($flag == 1){
        $hash{$_} = $key;
    }
}
open(FA,$ARGV[1]);

while(<FA>){
    chomp;
    if(/>/){
		$key = $_;
    }elsif($hash{$_}){print OUT $hash{$_}."\n".$_."\n";}
	else{print OUT $key."\n".$_."\n";}
}
