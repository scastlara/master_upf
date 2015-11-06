#!/usr/bin/perl

# This program doesn't do anything useful.
# It's only a way to parse a text like this:
# (aaa(bbb(ccc)bbb)aaa)
# It can be extended to do interesting things (maybe some Perl-RNA project :P)

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