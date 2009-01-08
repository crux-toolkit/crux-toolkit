#!/usr/bin/perl

# FILE: compare_scores.pl
# AUTHOR: Barbara Frewen
# DATE: December 23, 2008
# DESCRIPTION:  Compare search scores for crux and sequest. 
#               Assumes searches have already been run.  Crux results
#               are in an sqt file and sequest results are in
#               individual files, one for each spectrum.  Seqeust
#               result files must end in <scan #>.<charge>.scores


use strict;

my $usage = 
"USAGE: compare_scores <sqt file> <sequest-score directory> <results name>";

if( @ARGV != 3 ){
    die "$usage\n";
}

# get arguments
my $sqtfile = $ARGV[0];
my $sequestDir = $ARGV[1];
$sequestDir =~ s/([^\/]$)/$1\//; # put a / on the end if not there
                                # already
my $results = $ARGV[2];


# parse sqt file with crux scores
print "Parsing sqt file $sqtfile\n";

# store as a hash of lists. Key=sequence, data=[scan:charge:xcorr]
my %cruxScores;
my $line;
my $scan = 0;
my $charge = 0;

open SQT, $sqtfile;
while( $line = <SQT> ){
    chomp $line;

    if( $line =~ /^S/ ){
        my($s,$n1,$n2,$z) = split /\t/, $line;
        $scan = $n1;
        $charge = $z;
        next;
    }elsif( $line =~ /^M/ ){

        my($m,$xrank,$srank,$mass,$dcn,$xcorr,$sp,$mions,$tions,$seq) = 
            split /\t/, $line;
        
        $seq =~ s/^..([\*\@\#A-Z]*)..$/$1/;

        # sequest changes negative scores to 0.
        $xcorr = 0 if $xcorr < 0;

        #print "scan: $scan, charge: $charge, xcorr: $xcorr, seq: $seq\n";
        my $data = "$scan:$charge:$xcorr";

        if( $cruxScores{$seq} ){
            push @{ $cruxScores{$seq} }, $data ;   
        }else{
            $cruxScores{$seq} = [ $data ];
        }

    }# else skip
}

close SQT;

# for each sequest result file

my @sequestfiles = <$sequestDir/*scores>;

my $numfiles = @sequestfiles;
print "There are $numfiles *scores files in the directory $sequestDir\n";

# store pairs of scores in three lists, one for each charge
my @scorePairs;
my @scorePairsZ1;
my @scorePairsZ2;
my @scorePairsZ3;

for my $filename (@sequestfiles){

    next if (! -e $filename);
    open FILE, $filename or warn "Couldn't open sequest file $filename\n";
    (my $seqScan = $filename) =~ s/^.*\.([0-9]*)\.[0-3]\.scores/$1/;
    $seqScan =~ s/^0*//;
    (my $seqCharge = $filename) =~ s/^.*\.[0-9]*\.([0-3])\.scores/$1/;

    #print "$seqScan $seqCharge\n";

    # skip ahead to data
    while(my $line = <FILE> ){
        if( $line =~ / *---/ ){
            last;
        }
    }

    while(my $line = <FILE> ){
        chomp $line;
        $line =~ s,\/, ,;      # remove the / in Rank/Sp column
        $line =~ s/\+[0-9]*//; # remove multiple copies column
        $line =~ s/ +/\t/g;    # fix variable spacing
        
        my($blank,$num,$r1, $r1,$mass,$dcn,$xcorr,$sp,$ions,$prot,$seq) = 
            split /\t/, $line;
        $seq =~ s/^..([\*\@\#A-Z]*)..$/$1/;
        
        # look up sequence crux list
        next if ! $cruxScores{$seq};
        
        my @entryList = @{ $cruxScores{$seq} };
        
        # find matching scan/charge
        my $cruxEntry;
        for( @entryList ){
            $cruxEntry = $_ if( $_ =~ /$seqScan:$seqCharge/ );
        }
        next if ! defined $cruxEntry;
        
        # store score pairs
        my($cscan,$cz,$cruxScore) = split /:/, $cruxEntry;
        #print "seq $seq $cruxScore\t$xcorr\n";

        push @scorePairs, "$cruxScore\t$xcorr\t$cscan\t$cz\t$seq";

        if( $seqCharge == 1 ){
            push @scorePairsZ1, "$cruxScore\t$xcorr";
        }elsif( $seqCharge == 2 ){
            push @scorePairsZ2, "$cruxScore\t$xcorr";
        }if( $seqCharge == 3 ){
            push @scorePairsZ3, "$cruxScore\t$xcorr";
        }
    }
    close FILE;
}

# return results for plotting externally
# create a file with a list of the differences of the scores
#     and one with the score pairs
#open DIFFS, ">$results.diffs" or 
#                  die "Can't write results to $results.diffs";
open PAIRS, ">$results.score-pairs" or 
                  die "Can't write results to $results.score-pairs";
open PLOTME, ">plotme.$$" or die "Can't open plotme.$$\n";

for( @scorePairs ){
    print PAIRS  "$_\n";

    my($s1,$s2) = split /\t/, $_;
    print PLOTME "$s1\t$s2\n" if ($s1 != 0 and $s2 != 0); #plot non-zero scores
#    my $diff = abs($s1 - $s2);
#    print DIFFS  "$diff\n";
}
close PAIRS;
#close DIFFS;

=pod
open DIFFS, ">$results.diffs.z1" or 
                  die "Can't write results to $results.diffs.z1";
open PAIRS, ">$results.score-pairs.z1" or 
                  die "Can't write results to $results.score-pairs.z1";

for( @scorePairsZ1 ){
    print PAIRS  "$_\n";

    my($s1,$s2) = split /\t/, $_;
    my $diff = abs($s1 - $s2);
    print DIFFS  "$diff\n";
}
close PAIRS;
close DIFFS;


open DIFFS, ">$results.diffs.z2" or 
                  die "Can't write results to $results.diffs.z2";
open PAIRS, ">$results.score-pairs.z2" or 
                  die "Can't write results to $results.score-pairs.z2";

for( @scorePairsZ2 ){
    print PAIRS  "$_\n";

    my($s1,$s2) = split /\t/, $_;
    my $diff = abs($s1 - $s2);
    print DIFFS  "$diff\n";
}
close PAIRS;
close DIFFS;

open DIFFS, ">$results.diffs.z3" or 
                  die "Can't write results to $results.diffs.z3";
open PAIRS, ">$results.score-pairs.z3" or 
                  die "Can't write results to $results.score-pairs.z3";

for( @scorePairsZ3 ){
    print PAIRS  "$_\n";

    my($s1,$s2) = split /\t/, $_;
    my $diff = abs($s1 - $s2);
    print DIFFS  "$diff\n";
}
close PAIRS;
close DIFFS;

=cut

# plot scores, crux vs sequest

my $plotname = "$results.png";
open PLOT, "|gnuplot" or die "Can't open gnuplot\n";

print PLOT "set term png\n";
print PLOT "set output \"$plotname\"\n";
print PLOT "set xlabel \"Crux xcorr\"\n";
print PLOT "set ylabel \"SEQUEST xcorr\"\n";
print PLOT "plot 'plotme.$$' using 1:2 with points notitle, x lt 3\n";
print PLOT "\n";


close PLOT;

`rm -f plotme.$$`;

open EXACT, ">$results.exact";

for my $entry ( @scorePairs ){

    my($s1,$s2) = split /\t/, $entry;

    my $dif = abs($s1 - $s2);
    print EXACT "$entry\n" if ($dif < 0.005 and $s1 != 0); 
}

close EXACT;









