#!/usr/bin/env python
# AUTHOR: William Stafford Noble
# CREATE DATE: May Day, 2013
import sys

USAGE = """USAGE: convertCometPHP.py <php file>

  This is a goofy little helper script that takes as input one of the
  parameter PHP files distributed with Jimmy Eng's Comet search
  engine, strips off some of the unnecessary text, adds HTML headers
  and prints the result to standard output.  The idea is to allow
  automatic incorporation of any changes Jimmy makes to his
  documentation into Crux.  This is brittle, because it's based on
  hand-designed parsing code.  A better solution would be for Jimmy to
  store his parameter docs in, say, XML, and then use an existing
  parser.

"""

# Parse the command line.
if (len(sys.argv) != 2):
  sys.stderr.write(USAGE)
  exit(1)
phpFileName = sys.argv[1]

phpFile = open(phpFileName, "r")

# Find the start of the interesting section.
for line in phpFile:
  if ("<h2>" in line):

    # Get the name of the parameter.
    words = line.rstrip().split()
    if ((len(words) != 3) or 
        (words[0] != "<h2>Comet") or 
        (words[2][-5:] != "</h2>")):
      sys.stderr.write("Can't parse %s." % line)
      sys.exit(1)
    parameter = words[2][:-5]

    # Print the HTML header.
    sys.stdout.write("<html><head><title>comet parameter: %s</title>\n" % parameter)
    sys.stdout.write("<link href=\"../crux.css\" rel=\"styleSheet\" type=\"text/css\">\n")
    sys.stdout.write("</title><body>\n")
    sys.stdout.write("<h2>Comet parameter: %s</h2>\n" % parameter)
    break

# Print the body of the file.
for line in phpFile:
  if "</div" in line:
    break
  sys.stdout.write(line)
phpFile.close()

sys.stdout.write("</body></html>\n")
