#!/usr/bin/env python
# AUTHOR: William Stafford Noble
# CREATE DATE: 13 Sep 2015
import sys

USAGE = """USAGE: tsv2html.py [options] <tab-delimited file>

  Convert a tab-delimited text file an HTML table.  Use "-" to read
  from stdin.

  Options:
    --title <string>   
    --table-only
"""

title = ""
tableOnly = False

# Parse the command line.
sys.argv = sys.argv[1:]
while (len(sys.argv) > 1):
  nextArg = sys.argv[0]
  sys.argv = sys.argv[1:]
  if (nextArg == "--title"):
    title = sys.argv[1]
    sys.argv = sys.argv[1:]
  elif (nextArg == "--table-only"):
    tableOnly = True
  else:
    sys.stderr.write("Invalid option (%s).\n" % nextArg)
    sys.exit(1)
if (len(sys.argv) != 1):
  sys.stderr.write(USAGE)
  sys.exit(1)
tsvFileName = sys.argv[0]

# Open the input file.
if (tsvFileName == "-"):
  tsvFile = sys.stdin
else:
  tsvFile = open(tsvFileName, "r")

# Print the header.
if (not tableOnly):
  sys.stdout.write("<html>\n")
  sys.stdout.write("<head>\n")
  if (title != ""):
    sys.stdout.write("<title>%s</title>\n" % title)
  sys.stdout.write("</head>\n")
  sys.stdout.write("<body bgcolor=\"#ffffff\">\n")
  if (title != ""):
    sys.stdout.write("<h1>%s</h1>" % title)
sys.stdout.write("<table border cellspacing=0>\n")

# Read and display header row
fields = tsvFile.readline().rstrip().split("\t")
for field in fields:
  sys.stdout.write("<td><b>%s</b></td>" % field)
sys.stdout.write("</tr>\n")

# Process lines.
for line in tsvFile:

  fields = line.rstrip().split("\t")
  sys.stdout.write("<tr>\n")
  for field in fields:
    if (field == ""):
      field = "&nbsp"
    sys.stdout.write("<td>%s</td>" % field)
  sys.stdout.write("</tr>\n")
tsvFile.close()

# Print tail info.
sys.stdout.write("</table>\n")
if not tableOnly:
  sys.stdout.write("</body>\n")
  sys.stdout.write("</html>\n")
