#!/usr/bin/perl

use warnings;
use strict;

my $input = shift @ARGV;

my @numbers    = split /[,]/, $input;
@numbers       = sort {$a <=> $b} @numbers;

die "Give me more than one number separated by commas.\n" 
	if (@numbers == 1);

my $median     =  0;

if (scalar(@numbers) % 2 == 0) {
	# EVEN NUMBER
	my $first_idx  = (scalar(@numbers) / 2) - 1;
	my $second_idx = $first_idx + 1;
	$median        = ($numbers[$first_idx] + $numbers[$second_idx]) / 2;
} else {
	# ODD NUMBER
	my $median_idx = int(scalar(@numbers) / 2);
	$median        = $numbers[$median_idx];
}

print "These are your sorted numbers:\n",
      join(",", @numbers), "\n",
      "This is the median: $median\n";
