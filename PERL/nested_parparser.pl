use strict;
use warnings;

my %substring = (); 
my @saved;
my $iop; 
my $whole;

while(<>){
    chomp;  
    $whole .=$_;
}

for my $i (0..length($whole)) {
    my $substr = substr($whole, $i, 1);

    if ($substr eq "("){
        $iop++; 
        $substring{$iop}->{start} = $i;
    } elsif ($substr eq ")"){
        $substring{$iop}->{end} = $i+1; 
        my $length = $substring{$iop}->{end} - $substring{$iop}->{start};
        my $str = substr($whole, $substring{$iop}->{start}, $length);
        push @saved, $str;
        $iop--;
    }
    
}

print join("\n", @saved), "\n";