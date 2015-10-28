#!/usr/bin/perl

use warnings;
use strict;

for (my $i = 10; $i > 0; $i--) {
	print STDERR "$i \r";    # STDERR channel avoids the use of buffers,
				 # so you can print out directly and avoid
				 # sleeping issues.
	sleep(1);
}

print "BOOM!\n";
