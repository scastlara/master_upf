#!/usr/bin/perl

use warnings;
use strict;

for (my $i = 10; $i > 0; $i--) {
	print "$i\n";
	sleep(1);
}

print "BOOM!\n";