#!/usr/bin/perl
use strict;
use warnings;

my %hash;
my $key;

open(LIST,$ARGV[0]);
while(<LIST>){
    chomp;
    if(/>/){$key = $_;$hash{$key}="";}
    else{$hash{$key}=$_;}
}

foreach my $keys (keys %hash){
    print $keys."\n".$hash{$keys}."\n";
}
