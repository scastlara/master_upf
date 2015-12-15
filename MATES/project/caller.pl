
use warnings;
use strict;
use lib '/home/sergio/code/lib';
use Gauss qw(solve_system);
my $A = [

    [1, 1, 2],
    [0, 1, -1],


];


my $variables = solve_system($A);

for my $i (0.. $#{$variables}) {
    print "var$i = $variables->[$i]\n";
}
