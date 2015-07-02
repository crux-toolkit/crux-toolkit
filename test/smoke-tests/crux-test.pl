#!/usr/bin/perl -w

#
# Note: copied from Charles Grant's meme-test.pl

=head1 crux-test.pl

Usage: crux-test.pl [-up] <test list file>

The crux-test.pl script provides a framework for testing the overall 
functionality of programs in the crux distribution. It read a list of
tests from a '='-delimited file, one test per line. Each test is described
by four fields:

=over 4

=item 0
a Boolean value (0 or 1) indicating whether the test should be run under a
Darwin OS

=item o
the name of the test

=item o
the path to the file containing the known good output

=item o
the command line implementing the test

=back

The command line for each test is executed with the STDOUT captured to
a file. STDERR is written as is.  The test output is compared to the
known good output contained in the file specified in the second
field. If the test output matches the known good output, then the test
succeeds; otherwise, it fails. The results are printed to the standard
output.  When a test fails, a copy of the observed output (with the 
filename extension ".observed") is placed alongside the known good output
file.

If the '-u' option is specified, then the files containing the
known good output will be replaced by the output of the test.

It is assumed that the command to be tested sends its output to
standard out.

=cut

use strict;
use Getopt::Std;
use File::Temp qw/ tempfile tempdir /;

# Get the name of the architecture.
my $arch = `uname`;
print(STDERR "Running Crux smoke tests under $arch.\n");

# Handle the command line arguments.
my $usage = "Usage: crux-test.pl [-up] <file name>\n";
my (%options, $update);
if (!getopts('up:', \%options)) {
  die($usage);
};
if (scalar @ARGV != 1) {
  die("Incorrect number of arguments.\n", $usage);
}
if (defined $options{'u'}) {
  $update = 1;
} else {
  $update = 0;
}

if( defined $options{'p'}){
    my $path = $options{'p'};
    $ENV{'PATH'} = $path . ":" . $ENV{'PATH'};
    print "Prepended to PATH $path\n";
}

# Set an interrupt handler
my $caught_interrupt = 0;
$SIG{'INT'} = 'sigint_handler';

# Read each test from the file and execute it.
my $ignore_string = "";
my $num_tests = 0;
my $num_successful_tests = 0;
while (my $line = <ARGV>) {

  if ($caught_interrupt) {
    die("Testing was interrupted.");
  }  

  # Skip comments and blank lines
  next if $line =~ /^\#/;
  next if $line =~ /^\s*$/;

  # Parse the test parameters
  chomp $line;
  my @fields = split "=", $line;
  my $do_darwin = $fields[0];
  my $test_name = $fields[1];
  my $standard_filename = $fields[2];
  my $cmd = $fields[3];
  $cmd =~ s/^ //;
  $ignore_string = $fields[4];

  # Remove whitespace from the filename.
  $standard_filename =~ s/ //g;

  # Skip this test if we're on a Darwin OS and it doesn't pass there.
  if (($arch eq "Darwin\n") && ($do_darwin == "0")) {
    print "\n----- Skipping test $test_name \n";
    print STDERR "\n----- Skipping test $test_name \n";
    next;
  }

  # Execute the test
  print "\n----- Running test $test_name \n";
  print STDERR "\n----- Running test $test_name \n";
  my ($output_fh, $output_filename) = tempfile();
  if (!$output_filename) {
      die("Unable to create output file.\n");
  }

  my $result = &test_cmd($cmd, $standard_filename, 
			 $output_filename, $ignore_string);
  $num_tests++;

  if ($result == 0) {
      print "----- SUCCESS - test $test_name\n";
      $num_successful_tests++;
  } else {
      print "----- FAILURE - test $test_name\n";
  }    


  # Clean up
  close $output_fh;
  if ($update) {
    # If the update flag was given replace the existing standard
    # with the current test ouput.
    system("cp $output_filename $standard_filename");
  }
  unlink $output_filename;

}

# Print the summary of the tests.
print "\n-------------------------------\n";
print "----- Successful Tests: $num_successful_tests\n";
my $num_failed_tests = $num_tests - $num_successful_tests;
print "----- Failed Tests: $num_failed_tests\n";
print "----- Total Tests: $num_tests\n";
print "-------------------------------\n\n";
exit($num_failed_tests);

=head2 test_cmd($cmd, $standard_filename, $output_filename, $ignore_string)

Arguments:

$cmd                A string containing a shell command to execute.
$standard_filename  The name of a file containing the expected ouput
                    for the shell command.
$output_filename    The name of the file used to store the output of the
                    command.
$ignore_string      A string of regexes to be ignored by diff, each regex
                    in single quotes separated by a space. eg 'H' 'Time ex' 

The test_cmd function executes the command line passed in the $cmd variable
storing the output in $output_filename. If the return status indicates the
command was successful the output of the command is compared to the expected
output as stored in $standard_filename using diff. If there are differences
the output of diff is sent to standard out.

Returns 0 if the command executed successfully and the output exactly matched
the expected output, 1 otherwise.

=cut

sub test_cmd() {
  my ($cmd, $standard_filename, $output_filename, $ignore_string) = @_;
  my ($diff_fh, $diff_filename) = tempfile();
  if (!$diff_filename) {
    die("unable to create diff file.\n");
  }
  $cmd .= " > $output_filename";
  my $result = system($cmd);
  if ($result == 0) {
    # The command was successful, now vet the output.
      if( $ignore_string ){
	  $ignore_string =~ s/ \'/ -I \'/g; #add -I to each regex
      }else{
	  $ignore_string = " ";
      }

    my $diff_cmd = "diff $ignore_string $standard_filename $output_filename" .
	                                                  "> $diff_filename";
    $result = system($diff_cmd);
    if ($result == 0) {
      # The output of the command matches the expected output.
    } else {
      # The output of the command doesn't match the expected output.

      # Store a copy of the observed output.
      my $copy_command = "cp $output_filename $standard_filename.observed";
      print("$copy_command\n");
      system($copy_command);

      # Print the diff output.
      open(DIFF, $diff_filename) || die("Unable to read diff file.\n");
      my $diff_line;
      while ($diff_line = <DIFF>) {
        print $diff_line;
      }
      close DIFF;

      # Also check to see which columns changed.
      $diff_cmd = "./compare-by-field.pl $standard_filename $output_filename" .
          "> $diff_filename" ;
      system($diff_cmd);
      print "Compare by field results\n";
      open(DIFF, $diff_filename) || die("Unable to read diff file.\n");
      while ($diff_line = <DIFF>) {
        print $diff_line;
      }
      close DIFF;

      unlink $diff_filename;

    }
  } else {
    # The command was not succesful
    die("Testing failed.");
  }
  return $result;
}

sub sigint_handler {
  $caught_interrupt = 1;
}

