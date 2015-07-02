#!/bin/env python
# AUTHOR: Paul Draghicesu
# EDITED BY: Sean McIlwain
# CREATE DATE: 19 January 2010

from optparse import OptionParser
import os
import sys

usage = """usage: dta2ms2 <dta dir> <ms2 file>

Given a directory of dta files, create an
ms2 file.  Parses scan number from the
dta file.  Each dta file has a scan number and
charge information within the filename itself.
dta filename format expected is file.scan#.scan#.charge#.dta

"""


def main():
        parser = OptionParser(usage)
        (option, args) = parser.parse_args()
        if len(args) != 2:
                parser.error("incorrect number of arguments")
        datadir=args[0]
        filenames = os.listdir(datadir)
        ms2_filename = args[1]
        mass_p = 1.00727646688
        dta_numbers = {}
        for filename in filenames:
                tokens = filename.split(".")
                key = int(tokens[1])
                if not dta_numbers.has_key(key):
                        dta_numbers[key] = []

                dta_numbers[key].append(filename)
        ms2_file = open(ms2_filename, "w")
        ms2_file.write("H\tCreationDate\t8/03/09\n")
        ms2_file.write("H\tExtractor\tPaul Draghicescu\n")
        ms2_file.write("H\tExtractorVersion\t1\n")
        ms2_file.write("H\tExtractorOptions\t1\n")
        keys = dta_numbers.keys()
        keys.sort()
        for key in keys:
                dta_file = open(datadir + dta_numbers[key][0], "r")
                first_line = dta_file.readline().split()
                charge1 = int(first_line[1])
                mph1 = float(first_line[0]) #mass line is M+H
		mass1 = mph1 - mass_p
                if (len(dta_numbers[key]) == 1):
                        mz = (mass1 + (mass_p * charge1)) / charge1
                else:
                        mz1 = (mass1 + (mass_p * charge1)) / charge1
                        dta_file2 = open(datadir + dta_numbers[key][1], "r")
                        first_line2 = dta_file2.readline().split()
                        charge2 = int(first_line2[1])
			mph2 = float(first_line2[0])
                        mass2 = mph2 - mass_p
                        mz2 = (mass2 + (mass_p * charge2)) / charge2
                        mz = (mz1 + mz2) / 2
                        dta_file2.close()
                ms2_file.write("S\t" + str(key) + "\t" + str(key) + "\t" + str(mz) + "\n")
                ms2_file.write("Z\t" + str(charge1) + "\t" + str(mph1) + "\n")
                if (len(dta_numbers[key]) == 2):
                        ms2_file.write("Z\t" + str(charge2) + "\t" + str(mph2) + "\n")
                for line in dta_file:
                        ms2_file.write(line)
                dta_file.close()
        ms2_file.close()	
if __name__ == "__main__":
        main()
