#!/usr/bin/perl
use warnings;
use strict;
use GD::Simple;

#-- VARIABLES --#
my $nmax  = shift @ARGV;
my $l     = 50;
my $s_dir = [[0,1], [1,0], [-1, 0], [0,-1]];
my $o_dir = [[1,1], [-1,-1], [-1, 1], [1,-1]];
my $init  = [[200,125]];
my $img   = GD::Simple->new(400,250);
my $n = 0;

#-- MAIN --#
cross($n, $nmax, $init, $s_dir, $o_dir, $l, $img);
print $img->png;


#-- FUNCTIONS --#
sub cross {
	my ($n, $nmax, $init, $s_dir, $o_dir, $l, $img) = @_;
	my $dirs;
	my $force = 0;
	
	my @newstart = ();
  
	if ($n % 2 == 0) {
	$dirs  = $s_dir;
	$force = $l;
	} else {
	$dirs  = $o_dir;
	$force = $l / 1.41; # Pythagorean theorem
	}

	foreach my $point (@$init) {
	foreach my $dir (@$dirs) {
		my $endpoint = [
			$force*$dir->[0] + $point->[0], 
			$force*$dir->[1] + $point->[1]
		];
		$img->moveTo($point->[0],$point->[1]);
		$img->lineTo($endpoint->[0],$endpoint->[1]);
		push @newstart, $endpoint;
	}
	}
	
	return if $n == $nmax;
	$n++;
	cross($n, $nmax, \@newstart, $s_dir, $o_dir, $l/2, $img);
}
