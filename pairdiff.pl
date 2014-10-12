#! /usr/bin/perl -w

=head DESCRIPTION

## this is analogous to Tajima's nucleotide diversity calculation.
## pi=sigma(x_i*x_j*pi_ij) where xi and xj are the respective frequencies of the ith and jth sequences, pi_ij is the number of nucleotide differences per nucleotide site between the ith and jth sequences, and n is the number of sequences in the sample.

=head SYNOPSIS
Usage: perl pairdiff.pl oltype

=cut

use strict;

## data input and 2d matrix construction specifying the dimension, since we have 12 species xi and xj is a constant which equals to 1/12. the only dimension that matters is the pair numbers to calculate πij later.

my $file = "sankoffin-" . $ARGV[0] . ".txt";
my @in;
my @pair_pi;
open IN, "$file" or die "file not here";
while (<IN>) {
    chomp;
    push @in, [ split /\s/ ];
}
close IN;

my $pairnumber = $#{ $in[0] };

## for loop going through all combinations and caculate the pairwise difference
for my $i ( 0 .. 10 ) {
    my $diff = 0;
    for my $j ( 0 .. $pairnumber ) {
        if ( $in[ $i + 1 ][$j] - $in[$i][$j] != 0 ) { $diff++; }
    }
    push @pair_pi, $diff / $pairnumber;
}

## output the segregating site proportion and pi.
my $pi;
for my $i ( 0 .. $#pair_pi ) {
    $pi += $pair_pi[$i];
}
$pi = 2 * $pi / (12*(11));
my $err=sqrt($pi*(1-$pi)/(66*$pairnumber));

my $invariable_total = 0;
for my $i ( 0 .. $pairnumber ) {
    my @single_site = map $_->[$i], @in;
    if ( @single_site == grep { $_ == 1 } @single_site ) {

        # all equal
        $invariable_total++;
    }
}
my $segregating_prop = $invariable_total / $pairnumber;

print "$segregating_prop\n$pi\n$err\n";
