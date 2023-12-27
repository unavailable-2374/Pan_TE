#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(shuffle);

my $input_file = $ARGV[0];
my $fixed_length =$ARGV[2];  # 提取序列片段的固定长度
my $MAX_SIZE = $ARGV[1]; # sample size
$MAX_SIZE = $MAX_SIZE * 1024 * 1024;  # 500MB

unless ($input_file && $fixed_length) {
    die "Usage: $0 <fasta file> <fixed sequence length>\n";
}

# 读取fasta文件
open my $fh, '<', $input_file or die "Could not open $input_file: $!";
my @sequences;
my ($header, $sequence);

while (<$fh>) {
    chomp;
    if (/^>(.*)/) {
        if ($header) {
            push @sequences, { header => $header, sequence => $sequence };
        }
        $header = $1;
        $sequence = "";
    } else {
        $sequence .= $_;
    }
}
push @sequences, { header => $header, sequence => $sequence } if $header;
close $fh;

# 检查每个序列的长度是否大于固定长度
@sequences = grep { length($_->{sequence}) >= $fixed_length } @sequences;
if (scalar(@sequences) == 0) {
    die "No sequences found with length greater than or equal to $fixed_length\n";
}

# 随机抽取序列
my $current_size = 0;
while ($current_size < $MAX_SIZE) {
    my $random_seq = $sequences[rand @sequences];
    my $start_pos = rand(length($random_seq->{sequence}) - $fixed_length + 1);
    my $extracted_seq = substr($random_seq->{sequence}, $start_pos, $fixed_length);

    last if $current_size + length($extracted_seq) > $MAX_SIZE;
    
    print ">".$random_seq->{header}."_pos_$start_pos\n";
    print $extracted_seq."\n";
    $current_size += length($extracted_seq);
}

