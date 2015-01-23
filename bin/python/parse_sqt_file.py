#!PYTHON
# FILE: parse_sqt_file.py
# AUTHOR: CHRIS PARK
# CREATE DATE: 20 OCTOBER 2006

"""
This script parses a .sqt file and returns objects that contain the necessary information
"""

import os
import sys

#-------------------
class CruxObject:
    """an object that contains the crux result fields"""

    def __init__(self, file, mode="sp", fromFile=True):
        self.file = file
        self.mode = mode
        self.spectrums = []

        if fromFile:
            if not self.__parse():
                print "%s\n" % "failed to create sqt object"
        else:
            if not self.__parse_string():
                print "%s\n" % "failed to create sqt object"
                
    def __parse(self):
        """ starts from filename, parse the crux result file and creates spectrum and peptide objects"""
        try:
            crux_file = open(self.file, "r")
            
        except IOError:
            print "%s\n" % "failed to open file"
            return False

        if not self.__parse_crux_object(crux_file):
            crux_file.close()
            return False

        crux_file.close()
        return True
    
    def __parse_string(self):
        """ starts from string, parse the crux result file and creates spectrum and peptide objects"""
        return self.__parse_crux_object_2(self.file)

    def __parse_crux_object(self, crux_file):
        count = 0
        for line in crux_file:
            # new spectrum, create a spectrum object
            if line.startswith('#'):
                fields = line.rstrip('\n').split()
                if fields[2] == "SCAN":
                    scan_num = int(fields[4])
                elif fields[2] == "PRECURSOR":
                    mz = float(fields[4])
                elif fields[2] == "CHARGE:":
                    charge = int(fields[3])
                    #create spectrum object
                    self.spectrums.append(Spectrum(scan_num, charge, mz))
                    # peptide objects of the spectrum
                    count += 1
                    print "%s: %d" % ("spectrum added", count)
            elif line.startswith('P'):
                fields = line.rstrip('\n').split()
                # xcorr output format
                if self.mode == "xcorr":
                    peptide =  Peptide(int(fields[1]), float(fields[4]),float(fields[5]), \
                                       float(fields[3]), fields[7], int(fields[2]))
                    # sp output format
                elif self.mode == "sp":
                    peptide =  Peptide(None, None, float(fields[3]), float(fields[2]), \
                                       fields[5], int(fields[1]))
                else:
                    print "%s\n" % "In correct mode, must be sp or xcorr"
                    return False
                # add peptide to spectrum
                self.spectrums[len(self.spectrums)-1].addPeptide(peptide)
        return True

    def __parse_crux_object_2(self, crux_file):
        count = 0
        lines = crux_file.split('\n')
        
        for line in lines:
            # new spectrum, create a spectrum object
            if line.startswith('#'):
                fields = line.split()
                if fields[2] == "SCAN":
                    scan_num = int(fields[4])
                elif fields[2] == "PRECURSOR":
                    mz = float(fields[4])
                elif fields[2] == "CHARGE:":
                    charge = int(fields[3])
                    #create spectrum object
                    self.spectrums.append(Spectrum(scan_num, charge, mz))
                    # peptide objects of the spectrum
                    count += 1
                    print "%s: %d" % ("spectrum added", count)
            elif line.startswith('P'):
                fields = line.split()
                # Peptide object arguments(xcorr_rank, xcore, sp, mass, sequence, sp_rank=None)
                if self.mode == "xcorr":
                    peptide =  Peptide(int(fields[1]), float(fields[4]),float(fields[5]), \
                                       float(fields[3]), fields[6], int(fields[2]))
                    # sp output format
                elif self.mode == "sp":
                    peptide =  Peptide(None, None, float(fields[5]), float(fields[3]), \
                                       fields[6], int(fields[2]))
                elif self.mode == "logp_exp_sp":
                    peptide =  Peptide(None, None, float(fields[5]), float(fields[3]), \
                                       fields[6], int(fields[2]))
                else:
                    print "%s\n" % "In correct mode, must be sp or xcorr"
                    return False
                # add peptide to spectrum
                self.spectrums[len(self.spectrums)-1].addPeptide(peptide)
        return True

#-------------------

class SqtObject:
    """an object that contains the spectrum object"""

    def __init__(self, file):
        self.file = file
        self.spectrums = []

        if not self.__parse():
            print "%s\n" % "failed to create sqt object"
            self = None
            
    def __parse(self):
        """ parse the squest file and creates spectrum and peptide objects"""
        
        try:
            sqt_file = open(self.file, "r")

        except IOError:
            print "%s\n" % "failed to open file"
            return False

        # might have to split lines
        # iterate over all line of sqt file
        for line in sqt_file:
            if line.startswith('H') or line.startswith('L') : continue
            elif line.startswith('S'):
                #parse spectrum
                fields = line.rstrip('\n').split()
                id = fields[1].lstrip('0')
                self.spectrums.append(Spectrum(int(id), int(fields[3]), float(fields[6])))
                
            elif line.startswith('M'):
                #parse peptide
                fields = line.rstrip('\n').split()
                peptide = Peptide(int(fields[1]), float(fields[5]),float(fields[6]), \
                                   float(fields[3]), fields[9])
                self.spectrums[len(self.spectrums)-1].addPeptide(peptide)
        sqt_file.close()
        return True

#-------------------

class Spectrum:
    """an object that stores the spectrum"""
    #the spectrums scan num, charge, observed mass
    #array of peptide & its score
    def __init__(self, scan, charge, observed_mass):
        """ initialize the spectrum object """
        self.fields = {}
        self.peptides = []
        self.fields["scan"] = scan
        self.fields["charge"] = charge
        self.fields["observed_mass"] = observed_mass
        
    def addPeptide(self, peptide):
        """add a peptide to the spectrum"""
        self.peptides.append(peptide)

    def numPeptides(self):
        """how many peptides are there in the spectrum? """
        return len(self.peptides)

#-------------------
        
class Peptide:
    """an object that contains info for the peptide"""
    # xcore, Sp, mass, sequence,  
    def __init__(self, xcorr_rank, xcore, sp, mass, sequence, sp_rank=None):
        """  fill in the peptide object with all its components """
        self.components = {}
        self.components["xcorr_rank"] = xcorr_rank
        self.components["xcorr"] = xcore
        self.components["sp"] = sp
        self.components["mass"] = mass
        self.components["sequence"] = self.__extractSequence(sequence)
        self.components["sp_rank"] = sp_rank
        
    def __extractSequence(self, sequence):
        """ extracts the sequence from the two dots in between """
        real_sequence = sequence.split('.')
        if len(real_sequence) == 3:
            return real_sequence[1]
        return real_sequence[0]

