#!/usr/bin/perl
use warnings;
use strict;
use Benchmark;

my $file = shift @ARGV;

my $dot_data = _slurp($file);

 timethese(5,
     {"substring" => "submethod",
     "chop" => "chopmethod",
     "regexp" => "regexpmethod"});

#--------------------------------------------------------------------------------
sub _slurp {
    my $file   = shift;

    local $/ = undef;

    open my $fh, "<", $file;

    my $dot = <$fh>;

    $dot =~ s/\s+=\s+/=/g;
    $dot =~ s/(\->|\-\-)/ $1 /g;
    $dot =~ s/\n\s*#.*?\n/\n/g;

    return $dot;
}


sub submethod {

    for (my $i = 0; $i < length($dot_data); $i++) {
        my $substring = substr($dot_data, $i, 1);
    }

    return;
}


sub chopmethod {
    my $s = reverse($dot_data);
    my $char;
    while (length($s)) {
        $char = chop($s);
    }
}

sub regexpmethod {
    
    while ($dot_data =~ m/\G(.)/gsc) {
    }
}