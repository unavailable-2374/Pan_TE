use strict;
use warnings;

open(IN,$ARGV[0]);
open(FILE,$ARGV[1]);
open(OUT,">$ARGV[2]");

my %hash;
my $i=0;
my $num = 0;
my $used = 0;

while(<IN>){
    chomp;
    my @arr = split(/\t/);
    if($arr[0] eq "S"){
        $hash{$arr[1]} = $arr[2];
		$num+=length($arr[2]);
    }
}

while(<FILE>){
    chomp;
    my @arr = split(/\t/);
    print OUT ">segment_".$i."\n";
    for(my $j=0;$j<=$#arr;$j++){
        print OUT $hash{$arr[$j]};
		$used += length($hash{$arr[$j]});
    }print OUT "\n";
	$i++;
}
print $used/$num."\n";
