use warnings;
use strict;

local $/ = ">"; 
<>; 

while(<>){
	chomp; 
	my ($id, @seq) = split /\n/; 
	my $seq = join("", @seq); 
	print ">$id\n$seq\n" 
		if ($seq =~ s/(C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H)/\-\[$1\]\-/gi);
}