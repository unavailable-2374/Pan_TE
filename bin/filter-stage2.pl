#!/usr/bin/perl
use strict;
use warnings;

open(IN,$ARGV[0]);

my $min;
my %hash;

while(<IN>){
    chomp;
    if(/RS/){
        my @arr = split(/RS/);
        my @array = split(/\s+/,$arr[1]);
        if(exists $hash{$array[0]}){
            $hash{$array[0]}++;
        }else{
            $hash{$array[0]}=1;
        }
    }
}

open(FA,$ARGV[1]);

while(<FA>){
    chomp;
    if(/>/){
        s/>RS//;
        my @arr = split(/#/);
        if(exists $hash{$arr[0]} && $hash{$arr[0]}>10){
            print ">RS".$_."\n";
            $min = 1; 
        }else{$min = 0;}
    }elsif($min == 1 && /A/){
        print $_."\n";
    }
}
