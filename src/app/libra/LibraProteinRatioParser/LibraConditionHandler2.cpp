/*
Program       : LibraConditionHandler2
Author        : Nichole King
Date          : 10.16.05
SVN info      : $Id: LibraConditionHandler2.cpp 8748 2022-10-18 10:00:57Z real_procopio $

object holding condition attributes and methods 

Copyright (C) 2005 Nichole King

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

Nichole King
Institute for Systems Biology
401 Terry Avenue North 
Seattle, WA  98109  USA
*/

#include "LibraConditionHandler2.h"
#include <ios>
#include <string>

using namespace std;

LibraConditionHandler2::LibraConditionHandler2() {
    conditionFileName = NULL;

    // set default to -99 which means don't use a min threshhold feature
    minThresInten = (double)-99.;
    normalizationChannel = 0;
    numberOfReagentLines = 0;
    pseudoCount = 0;
}

LibraConditionHandler2::~LibraConditionHandler2() {}

void LibraConditionHandler2::setConditionFileName(const char * infile) {
    conditionFileName = infile;
}

void LibraConditionHandler2::setNormalizationChannel(int normChan) {
    normalizationChannel = normChan;
}

void LibraConditionHandler2::setMinimumThreshholdIntensity(double intensity) {
    minThresInten = intensity;
}

const char * LibraConditionHandler2::getConditionFileName() {
    return conditionFileName;
}

unsigned int LibraConditionHandler2::getPseudoCount() {
    return pseudoCount;
}

int LibraConditionHandler2::getNormalizationChannel() {
    return normalizationChannel;
}

int LibraConditionHandler2::getNumberOfReagentLines() {
    return numberOfReagentLines;
}

void LibraConditionHandler2::setQuantitationFileName(string infile) {
    quantitationFileName = infile;
}

/**
* get the name of outfile for protein quantitation.
* @return quantitationFileName
*/
string LibraConditionHandler2::getQuantitationFileName() {
    return quantitationFileName;
}


/**
* get the value to be used as a minimum integrated intensity for a reagent
* line during the outlier removal stage in LibraProteinRatioParser package
* @return minThreshInten
*/
double LibraConditionHandler2::getMinimumThreshholdIntensity() {
    return minThresInten;
}

/**
* read the condition file to parse for quantitation file name and
* minimum thrshhold intensity
* @return check is -2 for failed to open file
*/
int LibraConditionHandler2::readFile() {
    const char * infile = getConditionFileName();

    setConditionFileName(infile);

    char *nextline = new char[10000];
    char* data = NULL;
    Array<Tag*>* tags = NULL;
    Tag* tag = NULL;

    // in  Element tags for elements with children
    int inFragmentMassesElement; // 0 = not read yet; 1= in element; 2 = done reading

    std::ifstream in(conditionFileName);

    int checkread;

    try {
        if (!in) throw 100;

        inFragmentMassesElement = 0;
        numberOfReagentLines = 0;

        while( in.getline(nextline, 10000) ) {
            // check read state:
            checkread = in.rdstate();

            if (checkread & ios::eofbit) int dummyvar = 1;
            else if (checkread & ios::failbit ) throw 101;
            else if (checkread & ios::badbit ) throw 101;

            data = strstr(nextline, "<");

            while(data != NULL) {
	      tag = new Tag(data);

	      // set flag to in fragment mass element
	      if ( tag->isStart() && ! strcmp(tag->getName(), "fragmentMasses"))
		inFragmentMassesElement = 1;

	      // set flag to out of fragment element
	      if ( tag->isEnd() && ! strcmp(tag->getName(), "fragmentMasses")) {
		inFragmentMassesElement = 2;
		if (numberOfReagentLines == 0 ) throw 111;
	      }

	      if ( (inFragmentMassesElement == 1) && (tag->isStart())
		   && (! strcmp(tag->getName(), "reagent"))) {
		numberOfReagentLines++;
	      }

	      if ( tag->isEnd() && ! strcmp(tag->getName(), "pseudocount")) {
		int temp = atoi(tag->getAttributeValue("value"));
		if(temp <= 0)
		  temp = 0;
		pseudoCount = temp;
	      }

	      if ( tag->isEnd() && ! strcmp(tag->getName(), "normalization"))
		setNormalizationChannel(atoi(tag->getAttributeValue("type")));

	      if ( tag->isEnd() && ! strcmp(tag->getName(), "minimumThreshhold"))
		setMinimumThreshholdIntensity(atof(tag->getAttributeValue("value")));

	      if ( tag->isEnd() && ! strcmp(tag->getName(), "quantitationFile"))
		setQuantitationFileName(tag->getAttributeValue("name"));

	      if (tag != NULL) delete tag;
	      tag = NULL;

	      data = strstr(data+1, "<");
            }
        }
        in.close();

    } catch ( int code) {
        if( code == 100) {
	  cerr << "Error opening file " << conditionFileName << std::endl;
	  return 1; // We failed, switch to manual mode!

        } else if ( code == 111) {
	  cerr << "Error: Could not find reagent masses. " << conditionFileName << " format incorrect? \n";
        }

        in.close();
        return -2; // failure
    }
    delete []nextline;
    return 0; 
}
