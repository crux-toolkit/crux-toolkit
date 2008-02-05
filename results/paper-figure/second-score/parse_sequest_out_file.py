#!PYTHON
# FILE: parse_sequest_out_file.py
# AUTHOR: CHRIS PARK
# CREATE DATE: NOVEMBER 29 2006

"""
This script parses the sequest output file and returns objects that contain the necessary information
"""

import os
import sys
from parse_sqt_file import Spectrum
from parse_sqt_file import Peptide

#-------------------

class SequestObject:
    """an object that contains the sequest object"""

    def __init__(self, result):
        self.result = result
        self.peptides = []

        if not self.__parse():
            print "%s\n" % "failed to create sequest object"
            self = None;
            
    def __parse(self):
        """ parse the sequest file and creates peptide objects"""
        rank_area = False
        
        # iterate over all line of sqt file
        result_line = self.result.split('\n')
        
        for element in result_line:
            
            #debug
            #print fields

            #skip if nothing is here
            if len(element) == 0:
                continue
            
            #are we in the ranking area?
            #if not continue and go on...
            if not rank_area:
                temp_fields = element.split()
                
                if temp_fields[0] == "---" :
                    rank_area = True                    
                continue


            temp_fields = element.split('/')
            
            if len(temp_fields) != 3 :
                continue

            fields = []
            fields.extend(temp_fields[0].split())
            fields.extend(temp_fields[1].split())
            fields.extend(temp_fields[2].split())
                        
            #test this...
            temp = fields[len(fields)-1]
            if len(temp) < 2 or temp[len(temp)-2] != '.':
                continue

            #print fields
            
            #ok now we are in ranking land        
            #parse peptide
            peptide = Peptide(int(fields[1]), float(fields[5]),float(fields[6]), \
                              float(fields[3]), fields[len(fields)-1], int(fields[2]))

            #add peptide to list
            self.peptides.append(peptide)  
        return True
