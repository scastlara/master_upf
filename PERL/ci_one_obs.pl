#!/usr/bin/perl

use Statistics::R;

# This module allows the script to interact with R.
# Installation: invoke the CPAN module in bash 'perl -MCPAN -e shell'
# Then type the command 'install Statistics::R'

use strict;
use warnings;

#### This script computes a confidence interval (CI) for a one 
#### observation sample produced by a normal distribution with 
#### unknown mean and variance. Then it tests the method with a 
#### randomly generated resamplings of a normal distribution with
#### mean = $truemean and sd = $truesd. The number of resamplings
#### is submitted by the user.


my $R = Statistics::R->new();

# Creates a communication bridge with R and starts R

my $confidence;    # Confidence level for your CI
my $error;         
my $observation;   # The single observation of the sample 

my $tup = 100;     # This is to include high confidence values
my $tdown = 1.1;   # This is to include low confidence values
my $t;

print "Type a confidence: ";
$confidence = <STDIN>;

$t=($tup+$tdown)/2;
$R->set('t', $t);
$R->run(q`arrel <- sqrt((1/(2*t))*log((t+1)/(t-1)))`);
$R->run(q`extrem_b <- (t+1)*arrel`);
$R->run(q`extrem_a <- (t-1)*arrel`);
$R->run(q`alpha <- pnorm(extrem_b)-pnorm(extrem_a)`);
my $alpha = $R->get('alpha');

$error = $alpha - (1-$confidence);
$error = abs($error);

# The following loop computes a constant t such that [X-t*|X-X0|, X+t*|X+X0|] 
# is a CI with confidence level = $confidence. It makes use of interaction 
# with R via the module "Statistics::R".

while ($error >= 0.0001) {

	if ($alpha <= 1-$confidence){
		$tup=$t;
		$t=($tup+$tdown)/2;
	}
	else {
	$tdown=$t;	
	$t=($tup+$tdown)/2;
	}
	
$R->set('t', $t);
$R->run(q`arrel <- sqrt((1/(2*t))*log((t+1)/(t-1)))`);
$R->run(q`extrem_b <- (t+1)*arrel`);
$R->run(q`extrem_a <- (t-1)*arrel`);
$R->run(q`alpha <- pnorm(extrem_b)-pnorm(extrem_a)`);
$alpha = $R->get('alpha');
$error = abs($alpha-(1-$confidence));
}

print "t = $t\n";

# Yes! t only depends on the confidence sought.

# Now we carry out a random test showing the suitability of the method.
# We generate a random sample of a normal distribution, we use the first 
# element to produce a CI with confidence level = $confidence, then test
# what is the proportion of times the resampling meets the CI.

my $truemean=19;
my $truesd=1.2;

print "Number of 1-observation resamplings for the test: ";
my $resamplings = <STDIN>;
chomp($resamplings);

$R->set('truemean', $truemean);
$R->set('truesd', $truesd);
$R->set('samplesize', $resamplings);
$R->run(q`sample <- rnorm(samplesize,truemean,truesd)`);

my $current_sample; 
my $hits = 0;
my $length = 0;
$observation = $R->get('sample[1]');

# CI = [X-t*|X-x0|, X+t*|X-x0|]
# where x0 = $observation = first 1-observation sampling

for (my $i=1; $i <= $resamplings; $i++){
$R->set('i', $i);
$current_sample = $R->get('sample[i]');
$length = $length + 2*$t*abs($current_sample-$observation)/$resamplings;
	if ($truemean > ($current_sample - $t*abs($current_sample-$observation))){
		if ($truemean < ($current_sample + $t*abs($current_sample-$observation))){
		$hits++;
		}
	}
}

# $hits counts for how many 1-observation resamplings
# the observation belongs to the CI it produces

print "(hits:resamplings) ratio = ", $hits/$resamplings, "\n"; 

($hits/$resamplings >= $confidence) && print "Passed!\n";
($hits/$resamplings < $confidence) && print "Not passed!\n";

# Verification is attained as the ratio $hits / $resamplings is
# higher than $confidence as the $resamplings gets larger.

print "Average CI length = ", $length;
$R->stop();
