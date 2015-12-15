package Gauss;
use strict;
use warnings;
use Exporter qw(import);

#===============================================================================
# VARIABLES AND OPTIONS
#===============================================================================
our @ISA         = qw(Exporter);
our @EXPORT_OK   = qw(solve_system);
our %EXPORT_TAGS = ( DEFAULT => [qw(solve_system)]);

#===============================================================================
# MAIN
#===============================================================================
#-----
sub solve_system {
    my $matrix = shift;
    foreach my $row (0..$#{$matrix}) {

        foreach my $variable ($row + 1 .. $#{$matrix}){

            if ($matrix->[$variable]->[$row] != 0){

                my $diagonal_coeff = $matrix->[$row]->[$row];
                my $factor = - $diagonal_coeff / $matrix->[$variable]->[$row];

                foreach my $col (0..@{$matrix}) {
                    $matrix->[$variable][$col] = $matrix->[$variable][$col]*$factor+$matrix->[$row][$col];
                }

            } # if

        } # foreach variable

    } # foreach main

    my $vars = get_variables($matrix);
    return $vars;
}

#-----
sub get_variables {
    my $matrix    = shift;
    my @variables = ();
    my $last_row_idx = $#{$matrix};

    my $last_constant = $matrix->[$last_row_idx]->[$last_row_idx +1];
    my $last_coeff    = $matrix->[$last_row_idx]->[$last_row_idx];

    # Last variable
    push @variables, $last_constant / $last_coeff;

    # Iterator helps us identify the last non-zero variable
    my $iterator = 1;

    # From last row -1 to the first one
    for (my $row = $#{$matrix} -1 ; $row >= 0 ; $row--) {
        my $constant   = $matrix->[$row]->[-1];
        my $last_idx   = -2 - $iterator; # first variable coeff that is not zero

        my $summ = 0;
        # From the last variable to the first one that is not zero
        for (my $col = -2 ; $col > $last_idx ; $col--) {
            $summ += $matrix->[$row]->[$col] * $variables[$col +1];
        }

        my $newvar = ($constant - $summ) / $matrix->[$row]->[$last_idx];

        # Add variable to list
        unshift @variables, $newvar;
        $iterator++;
    }

    return(\@variables);
}
