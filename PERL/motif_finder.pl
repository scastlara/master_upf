use warnings;
use strict;

local $/ = ">"; 
<>; 
my $motif = "C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H";

while(<>){
	chomp; 
	my ($id, @seq) = split /\n/; 
	my $seq = join("", @seq); 
	print ">$id\n$seq\n" 
		if ($seq =~ s/($motif)/\-\[$1\]\-/gi);
}