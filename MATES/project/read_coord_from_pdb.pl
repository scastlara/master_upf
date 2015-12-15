#!/usr/bin/perl
#This script read the coordinates of a protein and ligand from a PDB file, and copies their coordinates to a new coordinate file.

use strict;
use warnings;

my %prot;
my %ligand;

open(I, "<2WY4.pdb");
while(<I>){
    next if not (/^ATOM/ or /^HETATM/);
    if(/^ATOM/){                        #We save the coordinates for each atom in the globin protein
        my $pos = substr($_,6,5)+0;
        my $x = substr($_,30,8)+0;
        my $y = substr($_,38,8)+0;
        my $z = substr($_,46,8)+0;
        $prot{$pos} = "$x;$y;$z";
    }

    if(/^HETATM/ and /HEM/){            #We save the coordinates of each atom in the ligand HEM
        my $pos = substr($_,6,5)+0;
        my $x = substr($_,30,8)+0;
        my $y = substr($_,38,8)+0;
        my $z = substr($_,46,8)+0;
        $ligand{$pos} = "$x;$y;$z";
    }
}
close(I);


#Here we do magic!
my %new_prot = %prot;
my %new_ligand = %ligand;


#We save the new coordinates in a new files
open(O, ">protein_2WY4.coord");
print O "###The format of this file is (coordinates only have 6 positions):\n###ATOM_number\tX_coord\tY_coord\tZ_coord\n\n\n";
my @positions = sort { $a <=> $b } keys %new_prot;          #In order to write the atoms in the right order
foreach my $pos (@positions){
    my ($x,$y,$z) = split(';',$new_prot{$pos});             #The coordinates are splitted in x, y, z components before saving them to the file
    print O "$pos\t$x\t$y\t$z\n";
}

open(O, ">ligand_2WY4.coord");
print O "###The format of this file is (coordinates only have 6 positions):\n###ATOM_number\tX_coord\tY_coord\tZ_coord\n\n\n";
my @positions = sort { $a <=> $b } keys %new_ligand;          #In order to write the atoms in the right order
foreach my $pos (@positions){
    my ($x,$y,$z) = split(';',$new_ligand{$pos});             #The coordinates are splitted in x, y, z components before saving them to the file
    print O "$pos\t$x\t$y\t$z\n";
}
