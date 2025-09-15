#!/usr/bin/perl

if (!@ARGV) {
    die "usage: MSPCollect BLAST_output_file\n";
}

open (BLST, "$ARGV[0]") || die "usage: MSPCollect BLAST_output_file\nCan not open the BLAST_output_file $ARGV[0]\n";

$score_cutoff = 80;
$iden_cutoff = 70;
$score = -1;

# Initialize variables
$head_query = 0;
$tail_query = 0;
$head_sbjct = 0;
$tail_sbjct = 0;
$query = "";
$sbjct = "";

while (<BLST>) {
    chomp;
    if (/^Query=\s+(.+?)(?:\s|$)/) {
        $query = $1;
        $query =~ s/\s+$//;  # Remove trailing whitespace
    }
    elsif (/^>(.+?)(?:\s|$)/) {
        $sbjct_current = $1;
        $sbjct_current =~ s/\s+$//;  # Remove trailing whitespace
    }
    elsif (/^\s+Score =\s*(\S+)/) {
        if (&filter) {
            &report;
        }
        # Reset for new alignment
        $head_q_boolean = 1;
        $head_s_boolean = 1;
        $head_query = 0;
        $tail_query = 0;
        $head_sbjct = 0;
        $tail_sbjct = 0;
        
        $sbjct = $sbjct_current;
        $score_string = $1;
        if ($score_string =~ /(\S+)(e|E)\+(\d+)/) {
            $score = $1 * exp($3*log(10));
        }
        else {$score = $score_string;}
    }
    elsif (/^\s+Identities =\s+(\d+)\/(\d+)/) {
        $iden = $1/$2*100;
    }
    elsif (/^Query\s+(\d+)\s+\S+\s+(\d+)/) {
        if ($head_q_boolean) {
            $head_query = $1;
            $head_q_boolean = 0;
        }
        $tail_query = $2;
    }
    elsif (/^Sbjct\s+(\d+)\s+\S+\s+(\d+)/) {
        if ($head_s_boolean) {
            $head_sbjct = $1;
            $head_s_boolean = 0;
        }
        $tail_sbjct = $2;
    }
    elsif (/^Parameters:/) {
        if (&filter) {&report;}
        $score = -1;
    }
}

sub filter {
    # Basic score and identity filtering
    if ($score <= 0 || $score <= $score_cutoff || $iden <= $iden_cutoff) {
        return 0;
    }
    
    # Calculate alignment length and filter by minimum length (80bp)
    my $query_length = abs($tail_query - $head_query) + 1;
    my $sbjct_length = abs($tail_sbjct - $head_sbjct) + 1;
    my $min_length = ($query_length < $sbjct_length) ? $query_length : $sbjct_length;
    
    if ($min_length < 80) {
        return 0;
    }
    
    # Filter out perfect self-matches (same sequence, same coordinates)
    # This removes trivial matches like: gi|1 1-1000 vs gi|1 1-1000
    # But keeps useful matches like: gi|1 1-500 vs gi|1 600-1100 (internal repeats)
    if ($query eq $sbjct && 
        $head_query == $head_sbjct && 
        $tail_query == $tail_sbjct) {
        return 0;
    }
    
    return 1;
}

sub report {
    # Convert lcl| format to gi| format for both query and subject IDs
    my $formatted_query = convert_id_format($query);
    my $formatted_sbjct = convert_id_format($sbjct);
    
    # Ensure coordinates are in forward order (start < end)
    my ($query_start, $query_end) = ($head_query < $tail_query) ? 
                                   ($head_query, $tail_query) : 
                                   ($tail_query, $head_query);
    my ($sbjct_start, $sbjct_end) = ($head_sbjct < $tail_sbjct) ? 
                                   ($head_sbjct, $tail_sbjct) : 
                                   ($tail_sbjct, $head_sbjct);
    
    printf("%06d %03d %05d %05d %s %05d %05d %s\n", $score, $iden, $query_start, $query_end, $formatted_query, $sbjct_start, $sbjct_end, $formatted_sbjct);
}

sub convert_id_format {
    my ($id) = @_;
    
    # Convert lcl|sequence_name to gi|sequence_name
    if ($id =~ /^lcl\|(.+)$/) {
        return "gi|$1";
    }
    
    # If already in gi| format or other format, keep as is
    if ($id =~ /^gi\|/) {
        return $id;
    }
    
    # For any other format, add gi| prefix
    return "gi|$id";
}