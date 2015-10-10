#!/usr/bin/perl
use strict;
use warnings;

# declaring some variables
my ($dna_seq,$AA,$C,$GG,$TT,$OTHER);
my %dinucleotides = ();
# initializing variables
$dna_seq  = "ATGCATTGGGGAACCCTGTGCGGATTCTTGTGGCTTTGGCCCTATCTTTTCTATGTCCAAGCTG".
            "TGCCCATCCAAAAAGTCCAAGATGACACCAAAACCCTCATCAAGACAATTGTCACCAGGATCAA";

$AA = $C = $GG = $TT = $OTHER = 0;
# Looping through the sequence
my $seqlen = length($dna_seq);

for (my $i = 0; $i < $seqlen; $i++) {
    my $char;
    $char = uc(substr($dna_seq,$i,2));
    $dinucleotides{$char}++;
    if ($char eq "TT") {
    	print "TT found at position $i!\n";
    }
};

# Printing results
my @sorted_di = sort {
	$dinucleotides{$b} <=> $dinucleotides{$a}
} keys %dinucleotides;

foreach my $di (@sorted_di) {
	if ($di eq "TT") {
		print "------\n",
		      "$di\t$dinucleotides{$di}\n",
		      "------\n";
	} else {
		print "$di\t$dinucleotides{$di}\n";
	}
}