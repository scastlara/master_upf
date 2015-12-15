
use warnings;
use strict;
use Gauss qw(solve_system);


my $x_system = [
#    a         b       c    v1   b1
    [3.023, 15.761, 22.868, 1, 3.023  ],
    [2.012, 15.419, 21.825, 1, 2.012  ],
    [2.239, 16.218, 20.516, 1, 2.239  ],
    [2.398, 15.62, 19.457, 1,  2.398  ],
];

my $y_system = [
#        d       e        f      v2    b2
    [  3.023,  15.761,  22.868,  1, -12.59  ],
    [  2.012,  15.419,  21.825,  1, -12.09 ],
    [  2.239,  16.218,  20.516,  1, -10.60 ],
    [  2.398,  15.62,   19.457,  1, -10.28 ],
];

my $z_system = [
#       g      h       i       v3   b3
    [ 3.023, 15.761, 22.868,   1, 12.627 ],
    [ 2.012, 15.419, 21.825,   1, 12.240 ],
    [ 2.239, 16.218, 20.516,   1, 12.369 ],
    [ 2.398, 15.62,  19.457,   1, 11.851 ],

];

my $variables = ();


 foreach my $point_coord ($x_system, $y_system, $z_system) {
     my $results = solve_system($point_coord);
     push @{$variables}, $results;
 }


my ($LM, $TV) = build_affinity($variables);
my $determinant = det3x3($LM);

if ($determinant == 1) {
    print "This could be a rotation. Det: $determinant\n";
} else {
    print "This is not a rotation. Det: $determinant\n";
}

# CALCULAR COORDENADAS DEL LIGANDO (APLICAR AFFINE MAP)

my $coords_file = shift @ARGV;

my $init_coords = read_coords($coords_file);

my %final_coords = ();
foreach my $atom (keys %{$init_coords}) {
    affine_map(
        $atom,
        $init_coords->{$atom},
        \%final_coords,
        $LM,
        $TV
    );
}




# --------
sub build_affinity {
    my $results_array = shift;
    my $LM = ();
    my $TV = ();

    foreach my $row (@{$results_array}) {
        my $newrow = ();
        for my $i (0..3) {
            if ($i != 3) {
                push @{$newrow}, $row->[$i];
            } else {
                push @{$TV}, $row->[$i];
            }
        }
        push @{$LM}, $newrow;
    }
    return($LM, $TV);
}

# --------
sub det3x3 {
    my $matrix = shift;
    my ($P1, $P2, $P3, $N1, $N2, $N3, $det);

    $P1 = $matrix->[0][0] * $matrix->[1][1] * $matrix->[2][2];
    $P2 = $matrix->[1][0] * $matrix->[2][1] * $matrix->[0][2];
    $P3 = $matrix->[0][1] * $matrix->[1][2] * $matrix->[2][0];
    $N1 = $matrix->[0][2] * $matrix->[1][1] * $matrix->[2][0];
    $N2 = $matrix->[0][0] * $matrix->[2][1] * $matrix->[1][2];
    $N3 = $matrix->[1][0] * $matrix->[0][1] * $matrix->[2][2];

    $det = $P1 + $P2 + $P3 - $N1 - $N2 - $N3;

    return($det);
}


sub read_coords {
    my $file = shift;
    my %atom_coords = ();

    open my $fh, "<", $file
        or die "Can't open $file : $!\n";

    <$fh> for (1..4); # skip 4 lines
    while (<$fh>) {
        chomp;
        my ($atom, $x, $y, $z) = split /\s+/;
        $atom_coords{$atom} = [$x, $y, $z];
    }

    return (\%atom_coords);
}

sub affine_map {
    my $atom      = shift;
    my $i_coords  = shift;
    my $f_coords  = shift;
    my $lm_matrix = shift;
    my $tr_vector = shift;

    
}


use Data::Dumper;
