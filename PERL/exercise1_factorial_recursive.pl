#!/usr/bin/perl

use warnings;
use strict;

print "Gimme a number: ";
my $number = <STDIN>;
chomp $number;

print "Factorial of $number = ", factorial($number), "\n";


sub factorial {
	my $num = shift;

	if ($num == 1) {
		return 1;
	} else {
		return $num * factorial($num - 1);
	}
}