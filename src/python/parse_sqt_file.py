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

class SqtObject:
    """an object that contains the spectrum object"""

    def __init__(self, file):
        self.file = file
        self.spectrums = []

        if not self.__parse():
            print "%s\n" % "failed to create sqt object"
            self = None;
            
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
    #the spectrums ID num, charge, observed mass
    #array of peptide & its score
    def __init__(self, id, charge, observed_mass):
        """ initialize the spectrum object """
        self.fields = {}
        self.peptides = []
        self.fields["id"] = id
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
    def __init__(self, xcore_rank, xcore, sp, mass, sequence, sp_rank=None):
        """  fill in the peptide object with all its components """
        self.components = {}
        self.components["xcore_rank"] = xcore_rank
        self.components["xcore"] = xcore
        self.components["sp"] = sp
        self.components["mass"] = mass
        self.components["sequence"] = self.__extractSequence(sequence)
        self.components["sp_rank"] = sp_rank
        
    def __extractSequence(self, sequence):
        """ extracts the sequence from the two dots in between """
        real_sequence = sequence.split('.')
        return real_sequence[1]


