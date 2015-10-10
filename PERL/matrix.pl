#!/usr/bin/perl

use warnings;
use strict;


#-- VARIABLES --#
my @matrix    = ([1,2,3], [5,7,12], [10,12,13]);
my @trasposed = ();


#-- MAIN PROGRAM --#
print "ORIGINAL MATRIX:\n";
my $trasposed = matrix_trasposomatic(\@matrix, 1);

print "\n\nTRASPOSED:\n";
matrix_trasposomatic($trasposed);


#-- FUNCTIONS --#
sub matrix_trasposomatic {
	my $matrix           = shift;
	my $traspose_flag    = shift;
	my $trasposed_matrix = ();
	my @diagonal         = ();

	for (my $i = 0; $i <= $#{$matrix}; $i++) {
		for (my $j = 0; $j <= $#{$matrix[$i]}; $j++) {
			print "$matrix->[$i]->[$j]\t";
			$trasposed_matrix->[$j]->[$i] = $matrix->[$i]->[$j];
			push @diagonal, $matrix->[$i]->[$j] if $i == $j;
		}
		print "\n";
	}
	
	print "DIAGONAL: ", join(",", @diagonal), "\n";
	return $traspose_flag ? $trasposed_matrix : 0;
}


