#!/usr/bin/perl

# FILE: compare_best_match.pl
# AUTHOR: Barbara Frewen
# DATE: December 23, 2008
# DESCRIPTION:  Compare search results for crux and sequest to see if best
#               match for each is same sequence. 
#               Assumes searches have already been run.  Crux results
#               are in an sqt file and sequest results are in
#               individual files, one for each spectrum.  Seqeust
#               result files must end in <scan #>.<charge>.scores


use strict;

my $usage = 
"USAGE: compare_best_match <sqt file> <sequest-score dir/prefix> <results name>";

if( @ARGV != 3 ){
    die "$usage\n";
}

# get arguments
my $sqtfile = $ARGV[0];
my $sequestDir = $ARGV[1];
#$sequestDir =~ s/([^\/]$)/$1\//; # put a / on the end if not there
                                # already
my $results = $ARGV[2];


# parse sqt file with crux scores
print "Parsing sqt file $sqtfile\n";

# store crux results in hash. Key=scan:charge, data=sequence:xcorr
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

        # only keep rank=1 matches
        next if $xrank != 1;
        
        $seq =~ s/^..([\*\@\#A-Z]*)..$/$1/;
        #print "scan: $scan, charge: $charge, xcorr: $xcorr, seq: $seq\n";
        my $key = "$scan:$charge";
        my $data = "$seq:$xcorr";

        if( $cruxScores{$key} ){
            #die "There are two results for $key.\n" if 
            #    $cruxScores{$key} ne $data ;
            push @{ $cruxScores{$key} }, $data;
        }else{
            $cruxScores{$key} = [ $data ];
        }

    }# else skip
}

close SQT;

# store pairs of scores in three lists, one for each charge
my @scorePairs;

# for each spec in sqt file
for my $scan_charge ( keys(%cruxScores) ){
    my ($scan, $charge) = split /:/, $scan_charge;

    # pad the scan number with appropriate number of zeros
    $scan = "000000" . $scan;
    $scan =~ s/^0*(......)$/$1/;

  #my $filename = $sequestDir."spectra_random_1000.$scan.$scan.$charge.scores";
    my $filename = $sequestDir.".$scan.$scan.$charge.scores";
    open FILE, $filename or warn "Couldn't open sequest file $filename\n";


    # skip ahead to data
    while(my $line = <FILE> ){
        if( $line =~ / *---/ ){
            last;
        }
    }

    # get line with best match
    my $line = <FILE>;
    chomp $line;
    $line =~ s,\/, ,;      # remove the / in Rank/Sp column
    $line =~ s/\+[0-9]*//; # remove multiple copies column
    $line =~ s/ +/\t/g;    # fix variable spacing
        
    my($blank,$num,$r1, $r1,$mass,$dcn,$xcorr,$sp,$ions,$prot,$seq) = 
        split /\t/, $line;
    $seq =~ s/^..([\*\@\#A-Z]*)..$/$1/;
    close FILE;
        
    # look up sequence crux list
    my @seq_score_list = @{ $cruxScores{$scan_charge} };
    #my $seq_score = $cruxScores{$scan_charge};

    for my $seq_score (@seq_score_list){


        my ($crux_seq, $crux_score) = split /:/, $seq_score;        
        #print "crux seq: $crux_seq sequest seq: $seq\n";
        if( $crux_seq eq $seq ){
            push @scorePairs, "$crux_score\t$xcorr\t$scan\t$charge\t$seq";
            last;
        }else{
            #print "$scan_charge $crux_seq $seq\n";
        }
    }
}# next psm

# return results for plotting externally
# create a file with score pairs
open PAIRS, ">$results-firstrank.score-pairs" or 
                  die "Can't write results to $results-firstrank.score-pairs";

for( @scorePairs ){
    print PAIRS  "$_\n";
}
close PAIRS;

# plot the results
my $plotname = "$results-firstrank.png";
open PLOT, "|gnuplot" or die "Can't open gnuplot\n";

print PLOT "set term png\n";
print PLOT "set output \"$plotname\"\n";
print PLOT "set xlabel \"Crux xcorr\"\n";
print PLOT "set ylabel \"SEQUEST xcorr\"\n";
print PLOT "set xrange [0:*]\n";
print PLOT "set yrange [0:*]\n";
print PLOT "plot '$results-firstrank.score-pairs' using 1:2 with points notitle, x lt 3\n";
print PLOT "\n";


close PLOT;


open EXACT, ">$results.exact";

for my $entry ( @scorePairs ){

    my($s1,$s2) = split /\t/, $entry;

    my $dif = abs($s1 - $s2);
    print EXACT "$entry\n" if ($dif < 0.005 and $s1 != 0); 
}

close EXACT;















