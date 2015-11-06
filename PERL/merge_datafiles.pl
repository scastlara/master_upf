#!/usr/bin/perl
use strict;
use warnings;

my %input;
my @values;

my $tblA = shift or die "Usage: $0 FILE\n";
open my $tblA_fh, '<', $tblA or die "Could not open '$tblA' $!";

my $tblB = shift or die "Usage: $0 FILE\n";
open my $tblB_fh, '<', $tblB or die "Could not open '$tblB' $!";

while (<$tblA_fh>) {
	chomp $_;
		(my $id,@values) = split /\s+/, $_;
		$input{$id}= [@values];
		
};
close $tblA_fh;

while (<$tblB_fh>){
	chomp $_;
	my ($id,$value) = split /\s+/, $_;
	
	if (exists $input{$id}){
		push (@{$input{$id}}, $value);
	}elsif (not exists $input{$id}){
		@values = ('NA', split (/\b/,$value))	;
		$input{$id}= [@values]; 	
	};		
};
close $tblB_fh;


open (my $outfh, '>', 'fileOut.tbl') or die "Could not create 'fileOut.tbl' $!";
foreach my $c (sort {$a <=> $b} keys %input){
	print $outfh $c."\t";
	if (@{ $input{$c} } == 2) {
		print $outfh $input{$c}->[0], "\t", $input{$c}->[1], "\n";
	}else{
		print $outfh $input{$c}->[0], "\t", "NA", "\n";
	};
	
};
