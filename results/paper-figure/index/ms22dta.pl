#!/usr/bin/perl -w
die "Usage: $0 <ms2-file>\n\n\tCreates dta files from an ms2 file.\n" 
	unless @ARGV==1;
$/="S";
my $ms2_file  = shift;
open(MS2, $ms2_file)||die "open:$!\n";
<MS2>;
my $dta_prefix = $ms2_file;
$dta_prefix=~s/\.ms2//;
my $dta_count=0;
my $spec=0;

while(<MS2>){
	chomp;
	if (/\r/){
		print "Please remove carriage returns from file.\n";
		exit;
	}
	my @lines = split /\n/, $_;
	my $first = shift @lines;
	# $first has format S scan1 scan2 precursor_mass
	my ($s, $scan1, $scan2, $prec_mass) = split /\s+/, $first;

  my @zlines = grep {/^Z/} @lines;
  my @nozlines = grep {!/^Z/ and !/^I/} @lines; # for now strip out ^I as well
	next unless @zlines;

	my ($d, $charge, $mass_plus_H) =  split /\s+/, $zlines[0];

	# $mass_plus_H from ms2 is M+H+
	
	# dta name has name,scan1,scan2,charge
	# first line has M+H+ mass,charge
	open(ONE, ">${dta_prefix}.$scan1.$scan2.$charge.dta")||die "open:$!";
	$dta_count++;
	print ONE map{"$_\n"} "$mass_plus_H\t$charge", @nozlines; 
	close ONE;

	if(@zlines==2){
		my ($d, $charge, $mass_plus_H) =  split /\s+/, $zlines[1];
		open(TWO, ">${dta_prefix}.$scan1.$scan2.$charge.dta")||die "open:$!";
		$dta_count++;
		print TWO map{"$_\n"} "$mass_plus_H\t$charge",  @nozlines;
		close TWO;
	}
        #only create 10 dta files
        #if ($dta_count > 100){
        #    close MS2;
        #    exit;
        #}

}
close MS2;
