use warnings;
use strict;

local $/ = ">"; # This magic variable changes the record separator.
                # It will make the diamond operator read from > to >
                # instead of from \n to \n (line by line)
my $motif = "C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H";
<>; # Skip first blank "line"

while(<>){
	chomp; 
	my ($id, @seq) = split /\n/; # split by lines, 
                                 # the first is the ID and the rest is the sequence
	my $seq = join("", @seq); 
	print ">$id\n$seq\n" 
		if ($seq =~ s/($motif)/\-\[$1\]\-/gi);
}