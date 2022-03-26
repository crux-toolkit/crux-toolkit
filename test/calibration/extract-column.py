#!/usr/bin/env python
# AUTHOR: William Stafford Noble
# CREATE DATE: 4 Feb 2022
import sys

USAGE = """USAGE: extract-column.py <file> <column>

Extract from a tab-delimited file all values from the column
with the specified header.
"""

# Parse the command line.
if (len(sys.argv) != 3):
    sys.stderr.write(USAGE)
    sys.exit(1)
my_filename = sys.argv[1]
column = sys.argv[2]

with open(my_filename, "r") as my_file:

    # Get the column index.
    header_line = my_file.readline()
    column_index = header_line.rstrip().split("\t").index(column)
    if (column_index == -1):
        sys.stderr.write(f"Cannot find column {column}.\n{header_line}")
        sys.exit(1)
    sys.stderr.write(f"Extracting {column} from column {column_index}.\n")

    # Print that column.
    for line in my_file:
        words = line.rstrip().split("\t")
        sys.stdout.write(f"{words[column_index]}\n")
