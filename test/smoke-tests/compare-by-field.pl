#!/usr/bin/perl

# A script to compare to files line by line
# on lines that differ, compare field by field

# TODO: first compare col 1 (scan num) and if not the same, try sorting
# if still not the same, just return with different psms error
# or make that a different script, find-same-psms which returns
# true (all the same) or returns two files with those that are 
# in common (scan,charge,seq)

use strict;
use Getopt::Long;

my $usage = "USAGE: compare-by-field.pl [-full] <file 1> <file 2>";

die "$usage\n" if @ARGV < 2;

my $print_full_details;

GetOptions("full"=>\$print_full_details) or die "$usage\n";

my $filename1 = shift @ARGV;
my $filename2 = shift @ARGV;

open FILE1, "$filename1" or die "Can't open file $filename1\n";
open FILE2, "$filename2" or die "Can't open file $filename1\n";

# store how many rows have differences by column number
my %cols;
my $line_count = 0;
my $success = 0;

# confirm two headers are the same,
# store for reporting differences
my $header = <FILE1>;
my $header2 = <FILE2>;
if( $header ne $header2 ){
    print "The two headers are not the same. Not doing file comparison.\n";
    exit(1);
}

while(my $line1 = <FILE1>){

    my $line2 = <FILE2>;
    $line_count++;

    if( ! $line2 ){
        print "File $filename1 contains more lines than $filename2\n";
        last;
    }

    next if $line1 eq $line2 ;

    # else look at each field
    chomp $line1;
    chomp $line2;
    my @parts1 = split /\s/, $line1;
    my @parts2 = split /\s/, $line2;
    my $col_count = 0;
    foreach my $token1 (@parts1){
        my $token2 = shift @parts2;
        $col_count++;

        next if $token1 == $token2;

        # else report it
        my $dif = abs($token1 - $token2);
        print "line $line_count, col $col_count $token1\t $token2\t$dif\n"
            if $print_full_details;
        $cols{$col_count} += 1;
        $success++;
    }

}

if( <FILE2> ){
    print "File $filename2 is longer than  $filename1\n";
    $success++;
}

close FILE1;
close FILE2;

# report differences
# get column names from the header
my @col_names = split /\t/, $header;
my @sorted_keys = sort { $a <=> $b } keys(%cols);
for my $c ( @sorted_keys){
    print "column: $col_names[$c-1] ($c), rows differ: $cols{$c}\n";
}

exit $success;































