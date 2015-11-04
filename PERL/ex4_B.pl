#!/usr/bin/perl
#
use strict;
use warnings;

my %mixture = ();

my ($line, @row, $i);

open (IN, "<fileA.tbl");

$line = <IN>;

while ($line) {

	chomp ($line);

	@row = split (" ", $line);
	
    	$mixture{$row[0]}[0] = $row[1];

	$line = <IN>;
};
close (IN);

open (IN, "<fileB.tbl");

$line = <IN>;

while ($line) {

	chomp ($line);

	@row = split (" ", $line);
	
    	$mixture{$row[0]}[1] = $row[1]; 

	$line = <IN>;
};
close (IN);

# Printing results

print ("Identifier", "\t", "File A", "\t", "File B", "\n");

foreach my $identifier (sort {$a <=> $b} keys %mixture){
	for ($i = 0; $i < 3; $i++){
		# Prints the identifier number
		if ($i == 0) {print $identifier, "\t"x2;}
		else {
			# Prints a defined value
			if (defined($mixture{$identifier}[$i-1])) {
				print ($mixture{$identifier}[$i-1]);
			} else {
				# Prints NA when the value is undefined
				print ("NA");
			}
			# Prints a tab char when value is from File A
			if ($i == 1) {print "\t";}
			# Prints return when value is from File B
			else {print "\n";}
		}
	}
}
