#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use lib '/home/sergio/code/lib';
use Gauss "solve_system";

my $data_file = shift @ARGV;
my $data      = read_data($data_file);
my $n         = scalar(@{ $data });
my @x_list    = map {$_->[0]} @{$data};
my @y_list    = map {$_->[1]} @{$data};;

my $x_sum = summatory(\@x_list);
my $y_sum = summatory(\@y_list);

my @x_2 = map { expo($_, 2) } @x_list;
my @x_3 = map { expo($_, 3) } @x_list;
my @x_4 = map { expo($_, 4) } @x_list;

my $sum_x2 = summatory(\@x_2);
my $sum_x3 = summatory(\@x_3);
my $sum_x4 = summatory(\@x_4);

my @x_times_y = ();
my $y_mean = $y_sum / $n;

foreach my $i (0..$#{$data}) {
    push @x_times_y, $data->[$i]->[0] * $data->[$i]->[1];
}

my $x_times_y_sum = summatory(\@x_times_y);

my @x2_times_y = ();

foreach my $i (0..$#{$data}) {
    push @x2_times_y, ($data->[$i]->[0]**2) * $data->[$i]->[1];
}
my $x2_times_y_sum = summatory(\@x2_times_y);

my $SS_y   = 0;
foreach my $y_i (@y_list) {
    $SS_y += ($y_i - $y_mean)**2;
}





# LINEAR REGRESSION

my $simple_matrix = [
    [$n , $x_sum,  $y_sum],
    [$x_sum,  $sum_x2    , $x_times_y_sum ]
];

my $linear_results = solve_system($simple_matrix);
my $SS_linear_reg = 0;
foreach my $x_i (@x_list) {
    my $y_pred = $linear_results->[0] + $linear_results->[1] * $x_i;
    $SS_linear_reg += ($y_pred - $y_mean)**2;
}

my $R_linear_reg = $SS_linear_reg / $SS_y;
print "R square linear = $R_linear_reg con ss linear $SS_linear_reg\n";
print "a = $linear_results->[0] b = $linear_results->[1]\n";

# MULTIPLE REGRESSION

my $matrix = [
    [$n,      $x_sum,  $sum_x2, $y_sum          ],
    [$x_sum,  $sum_x2, $sum_x3, $x_times_y_sum  ],
    [$sum_x2, $sum_x3, $sum_x4, $x2_times_y_sum ]
];

my $results = solve_system($matrix);


# TOTAL SUM OF SQUARES

my $SS_quad_reg = 0;
foreach my $x_i (@x_list) {
    my $y_pred = $results->[0] + $results->[1] * $x_i + $results->[2] * $x_i ** 2;
    $SS_quad_reg += ($y_pred - $y_mean)**2;
}

my $R_quad_square = $SS_quad_reg / $SS_y;
print "R square quadratic = $R_quad_square\n";





foreach my $i (0.. $#{$results}) {
    print "var", ++$i, " = $results->[$i -1]\n";
}

# FUNCTIONS
# -------------------
sub read_data {
    my $file = shift;
    my $data;

    open my $fh, "<", $file
        or die "Can't open data file $file : $!\n";

    <$fh>; # skip first line
    while (<$fh>) {
        chomp;
        my ($x, $y) = split /\s+/;
        push @{$data}, [$x, $y];
    }
    return $data;
}

# -------------------
sub summatory {
    my $list = shift;
    my $sum  = 0;

    foreach my $number (@{ $list }) {
        $sum += $number;
    }

    return $sum;
}

# -------------------
sub expo {
    my $number = shift;
    my $n      = shift;
    my $result = 1;

    foreach my $i (1..$n) {
        $result *= $number;
    }

    return $result;
}
