#!/usr/bin/perl
#FILE: spectacle
#AUTHOR: Barbara Frewen
#CREATE DATE: 01/10/06
#DESCRIPTION: Read in an annotation file and 
#    produce plots of spectra therein

use strict;
use Getopt::Long;
use FileHandle;

# !!!!! see TODO items !!!!!!!!
#TODO label size bigger for png
#TODO not printing the last peak?
#TODO mark the precursor m/z

my $usage = "Usage: spectacle [options] <annotations>

Options: 
  -format <eps|fig|png> - format of output file. Default is eps. 
  -html <file>    - write an html file with all of the figures. 
  -names <file>   - specifiy names of output files. Only enabled for png.
  -pre-mz         - mark the precursor m/z on the x-axis.
  -title          - create a title above the figure of the form 
                     <identifier> m/z=<m/z> <comment>.    
  -labeltop       - labels at the top of figure
  -nobox          - do not draw a border at top and left of plot
  -script <file>  - print commands to file instead of running
                    gnuplot directly\n";


#read the command line
my $format = "eps";
my $names_file;
my $label_up;
my $html;
my $premz;
my $title;
my $gnuplot_file;
my $annot_file;
my $no_box;
GetOptions("names=s"=>\$names_file,
           "format=s"=>\$format,
           "html=s"=>\$html,
           "pre-mz"=>\$premz,
           "labeltop"=>\$label_up,
           "script=s"=>\$gnuplot_file,
           "title"=>\$title);   
#TODO add verbosity?
if(@ARGV != 1){
  print "Missing annotation file or option arguments\n";
  print $usage;
  exit;
}else{
  $annot_file = $ARGV[0];
}

#check for option errors...

#format: must be eps png fig; write terminal string
my $terminal;
if($format eq "eps"){
  $terminal = "set terminal postscript eps color solid 24\n";
}elsif($format eq "png"){
  #$terminal = "set terminal png medium color\n";
  $terminal = "set terminal png medium\n";
}elsif($format eq "fig"){
  $terminal = "set terminal fig color small 24\n";
}else{
  print "Format $format not recognized. Use png, eps, or fig.";
  exit;
}
#print "Format is $format\n";

#check for gnuplot and exit if not present
if( `which gnuplot` eq "" ){
    die "Could not find required gnuplot.\n";
}

#html: change to false and write warning if format!=png
if($html && $format ne "png" ){
  print STDERR "Warning: HTML file cannot be created with $format format\n";
  undef $html;
}

#names: open and put in a hash (id as key, file name as entry)
#I haven't checked to see if any names are repeated
(open NAMES, $names_file or 
 die "Could not open names file $names_file\n") if $names_file;; 
my %spec_names;
if($names_file){
  my $line;
  while($line = <NAMES>){
    my ($id, $filename) = split /\s/, $line;
    $spec_names{$id} = $filename;
  }
}
close NAMES;

open ANNOT, $annot_file or die "Couldn't open annotation file $annot_file\n";
my $line = <ANNOT>;

#open gnuplot and initialize settings for all plots
($gnuplot_file) ? (open GNUPLOT, ">$gnuplot_file") : (open GNUPLOT, "|gnuplot");
print GNUPLOT $terminal;
#print GNUPLOT "set noytics\n"
print GNUPLOT "set xtics nomirror\n";
print GNUPLOT "set xlabel \"m/z\"\nset ylabel \"Intensity\"\n";
print GNUPLOT "set border 3\n" if $no_box;
#save fig names for html page, temp file names to delete
my @figNames;
my @tmpFileNames;

while( $line =~ /^>/){  #each loop is one spectrum
  #remove > and whitespace
  chomp $line;
  $line =~ s/>\s*//;  
  my($id, $mz) = split /\s/, $line;
  #print "id is $id, mz is $mz\n";
  my $comment =~ s/$id\s+$mz\s+//;

  #check names
  my $fig_file_name = "$id.$format";
  if( $names_file and !$spec_names{$id} ){
    #print "Don't make fig for $id\n";
    next;  
  }elsif( $names_file ){
    $fig_file_name = "$spec_names{$id}.$format";
  }
  push @figNames, $fig_file_name;

  #title string
  if($title){
    $title = $id." mz: ".$mz." ".$comment;
    #print "title is $title\n"; 
  }

  #list of series, store text to be printed to file
  my %series;  #indexed by color
  open BLACK, ">black.$id.$$.tmp" 
    or die "Can't create temporary data file\n";
  push @tmpFileNames, "black.$id.$$.tmp";

  #list of labels
  my @labels;

  #list of warnings
  my %warnings;

  #for each peak
  my $peak_line;
  #my $max_y = 0;
  my $max_y = -10;
  while( ($peak_line = <ANNOT>) !~ /^>/ and (not eof ANNOT)){  
    chomp $peak_line;
    my($mz, $in, $label, $color) = split /\s+/, $peak_line;
    
    #print "mz: $mz int: $in label: $label color: $color\n";
    #keep track of largest intensity for adjusting y-range
    $max_y = $in if ($in > $max_y);

    #add label to list
    (push @labels, "$mz $in, $label") if $label;

    #check series, create new if necessary
    if(!$color){
      print BLACK "$mz $in\n";
    }elsif( $series{$color} ){
      $series{$color} = $series{$color}."$mz $in\n";
    }else{
      #check that color is legal
      #add to list
      if($color eq "red" or
       $color eq "green" or
       $color eq "blue" or
       $color eq "magenta"){
      $series{$color} = $series{$color}."$mz $in\n";
      }else{
      $warnings{$color} = "Warning: unrecognized color ($color)\n";
      print BLACK "$mz $in\n";
      }
    }
  }#last peak read

  $line = $peak_line;
  #print "end of spec, line is $line";
  
  #close series files
  close BLACK;
  #print colors to files
  foreach my $color (keys %series){
    open FILE, ">$color.$id.$$.tmp";
    print FILE $series{$color};
    close FILE;
    push @tmpFileNames, "$color.$id.$$.tmp";
  }

  #print warnings  
  foreach my $color (keys %warnings){
    print STDERR $warnings{$color};
  }

  #adjust the labels to be above the peak
  my $offset = $max_y * 0.04;
  my $abs_y = $max_y if $label_up;
  #print "max_y is $max_y, offset is $offset, and abs_y is $abs_y\n";
  my $label_string;
  for (@labels){
    my ($mz, $in, $label) = split /\s/;
    ($label_up) ? $in = $abs_y : $in +=$offset;
    $label_string .= "set label \"$label\" at $mz, $in center\n";
  }

  #clear labels from last plot
  print GNUPLOT "set nolabel\n";
  #write filename
  print GNUPLOT "set output \"$fig_file_name\"\n";

  #set yrange as greater than highest peak 
  $max_y *= 1.1; # if $label_string;
  #print GNUPLOT "set yrange [-10001:$max_y]\n";
  #print GNUPLOT "set yrange [-100:10]\n";
  print GNUPLOT "set yrange [0:$max_y]\n";

  #write title
  print GNUPLOT "set title \"$title\"\n" if $title;

  #write labels
  #print GNUPLOT $label_string;
  #print GNUPLOT "set nolabel\n";
  #write plot command
  my $plot = "plot \"black.$id.$$.tmp\" notitle with impulses lt -1, ";
  foreach my $color (keys %series){
    my $color_number = get_color_number($color, $format);
    $plot = $plot.
      " \"$color.$id.$$.tmp\" notitle with impulses lt $color_number,";
  }
  $plot =~ s/,\s*$//;  #remove last comma
  print GNUPLOT $plot."\n";

  #skip over any blank lines between spec
  while( ($line !~ /^>/) and (not eof ANNOT) ){
    $line = <ANNOT>;
  }
#end while, next spec
}continue{  #here if a spec was skipped b/c not in names file
  while( ($line !~ /^>/) and (not eof ANNOT) ){
    $line = <ANNOT>;
  }
} 

close GNUPLOT;
for (@tmpFileNames){
  #print "temp filename is $_\n";
  `rm $_`;
}

#write html file
if($html){
  open HTML, ">$html" or die "Couldn't open $html\n";
  print HTML "<html><body>\n";
  for (@figNames){
    print HTML "<p>$_</p><p><img src=\"$_\"></p>\n";
  }
  print HTML "</body></html>";
  close HTML;
}

########## END OF MAIN #####################


#number    fig    png    eps      
#-1      black          
#1    black    red    red      
#2    blue    green    green      
#3    green    blue    blue      
#4      magenta    magenta
#5    red    
#6    magenta    
#7             black


sub get_color_number{ #color, format
  my $color = shift;
  my $format = shift;
  my $code = 1;

  if($color eq "black"){
    if($format eq "png"){
      $code = -1;
    }elsif($format eq "eps"){
      $code = 7;
    }elsif($format eq "fig"){
      $code = 1;
    }
  }elsif($color eq "red"){
    if($format eq "png"){
      $code = 1;
    }elsif($format eq "eps"){
      $code = 1;
    }elsif($format eq "fig"){
      $code = 5;
    }
  }elsif($color eq "green"){
    if($format eq "png"){
      $code = 2;
    }elsif($format eq "eps"){
      $code = 2;
    }elsif($format eq "fig"){
      $code = 3;
    }
  }elsif($color eq "blue"){
    if($format eq "png"){
      $code = 3;
    }elsif($format eq "eps"){
      #$code = 3;
      $code = 21;
    }elsif($format eq "fig"){
      $code = 2;
    }
  }elsif($color eq "magenta"){
    if($format eq "png"){
      $code = 4;
    }elsif($format eq "eps"){
      $code = 4;
    }elsif($format eq "fig"){
      $code = 6;
    }
  }
  return $code;
}
