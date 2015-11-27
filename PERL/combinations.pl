use warnings;
use strict;

my @list = @ARGV;
my @all_combinations = ();

for my $i (1..$#list)  {
    my @combinations = combinations(\@list, $i);
    push @all_combinations, join(" ", @$_) for @combinations;
}

print join("\n", @all_combinations), "\n";

sub combinations {
    my $list = shift;
    my $n    = shift;

    error("Something went wrong when getting the combinations of your files", 1) 
        if $n > @$list;

    return map [$_], @$list if $n <= 1;

    my @comb;

    for (my $i = 0; $i+$n <= @$list; ++$i) {
        my $val  = $list->[$i];
        my @rest = @$list[$i+1..$#$list];
        push @comb, [$val, @$_] for combinations(\@rest, $n-1);
    }

    return @comb;
}