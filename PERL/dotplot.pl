#!/usr/bin/perl
use strict;
#use warnings;

# Variables Declaration and Initialization
my %SEQS = ();
my ($name1, $name2, $file);

# Get sequence names and input filename
# from command-line arguments
if (scalar @ARGV > 2){ # Separete each argument
    ($name1, $name2, $file) = @ARGV;
} else{
    # Split the first argument to the two names of the sequences
    ($name1, $name2) = split(":", $ARGV[0]);
    $file = $ARGV[1]
}


# Load sequences from input file into a simple hash
#  (sequence names as keys,
#   sequence strings as the corresponding values)
open(FASTA, $file); #"short_seqs.fa.txt"
local $/ = ">";
while (<FASTA>) {
    chomp;
    my ($id, @seq) = split /\n/; # split by lines,
                                 # the first is the ID and the rest is the sequence
    my $seq = join("", @seq);
    $SEQS{$id} = $seq;
}; # while
close(FASTA);


# Create the nxm matrix,
# load sequences on the left column and top row,
# set 0/1 values for match/mismatch respectively
#Choosing sequences to work with
my ($seqA, $seqB) = ($SEQS{$name1}, $SEQS{$name2});

$seqA = uc($seqA); # Ensuring all characters
$seqB = uc($seqB); # from the sequences are in upper case

my $n = length($seqA);
my $m = length($seqB);

my @MATRIX = (); # Initialize the matrix holder

$MATRIX[0][0] = q{*}; # The upper and leftmost cell
                      # does not contain data at all

# Now load sequences into the "zeroes" cells
for (my $j = 1; $j < $n; $j++) {
    # Loops through columns and loads sequence A in row 0
    $MATRIX[0][$j] = substr($seqA, $j, 1);
}; # for $j

for (my $i = 1; $i < $m; $i++) {
    # Loops through rows and loads sequence B in column 0
    $MATRIX[$i][0] = substr($seqB, $i, 1);
}; # for $j

# Check data
#for (my $j = 0; $j < $n; $j++){
    #for (my $i = 0; $i < $m; $i++){
        #print $MATRIX[$j][$i];
    #}
    #print "\n";
#}

# And finally, calculate identities/mismatches for each cell
$n++; # Increase length to include the sequence
$m++; # cells in the matrix, see below

for (my $i = 1; $i < $m; $i++) {
    # loop through rows (skip row 0 -> sequence holder)
    for (my $j = 1; $j < $n; $j++) {
    # loop through columns (skip col 0 -> sequence holder)
        if ($MATRIX[$i][0] eq $MATRIX[0][$j]) {
            $MATRIX[$i][$j] = 1; # match
        } else {
            $MATRIX[$i][$j] = 0; # mismatch
        };
     }; # for $j
}; # for $i


#Print to the terminal the Dot-Plot Matrix
for (my $i = 0; $i < $m; $i++) {
    #loop through rows
    for (my $j = 0; $j < $n; $j++) {
        #loop through columns
        unless ($i == 0 || $j == 0) {
            #print dot-plot cells
            print STDOUT ($MATRIX[$i][$j] ? q{*} : q{ });
        } else {
            #print nucleotides
            print STDOUT $MATRIX[$i][$j];
        };
        #print horizontal spacer
        print STDOUT q{ };
    }; # for $j
    #print rows
    print STDOUT "\n";
}; # for $i

exit(0);
