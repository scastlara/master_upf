#!/usr/bin/perl

use warnings;
use strict;

print "Gimme a number: ";
my $number = <STDIN>;
chomp $number;

my $int = 1;

for (my $i = 1; $i <= $number; $i++) {
	$int *= $i;
}

print "Factorial of $number = $int\n";
