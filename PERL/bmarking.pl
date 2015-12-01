#! perl -slw
use strict;
use Benchmark qw[ cmpthese ];

our $string = 'x'x5000;

cmpthese -1, {
    substr      => q[
        for ( 0..length( $string)-1 ) {
            my $c = substr $string, $_, 1;
        }
    ],
    splitArray  => q[
        my @c=split'',$string;
        for( 0 .. $#c ){
            my $c = $c[$_];
        }
    ],
    splitFor    => q[
        for( split'', $string ){
            my $c = $_;
        }
    ],
    unpack      => q[
        for( unpack 'C*', $string ) {
            my $c = chr;
        }
    ],
    reverseChop => q[
        my $s = reverse $string;
        my $c;
        $c = $_ while chop $s;
    ],
    chop        => q[
        my $s = $string;
        my $c;
        $c = $_ while chop $s;
    ],
    ramfile        => q[
        open my $ram, '<', \$string;
        my $c;
        $c = $_ while $_ = getc( $ram );
    ],
};