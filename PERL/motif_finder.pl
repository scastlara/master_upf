use warnings;
use strict;

local $/ = ">"; # This magic variable changes the record separator.
                # It will make the diamond operator read from > to >
                # By default it is set to "\n", that's why <> reads "line" by "line".
                
my $motif = "C[.{2}.{4}]C.{3}[LIVMFYWC].{8}H.{3,5}H";
<>; # Skip first blank "line"

while(<>){
	chomp; 
	my ($id, @seq) = split /\n/; # split by lines, 
                                 # the first is the ID and the rest is the sequence
	my $seq = join("", @seq); 
	print ">$id\n$seq\n" 
		if ($seq =~ s/($motif)/\-\[$1\]\-/gi); # the /g modifier means "match more than one time if possible"
                                               # the /i modifier makes the match case insensitive
}
