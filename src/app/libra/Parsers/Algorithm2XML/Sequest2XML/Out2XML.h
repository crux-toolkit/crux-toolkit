#ifndef OUT2XML_H_
#define OUT2XML_H_


/*

Program       : Out2XML
Author        : David Shteynberg <dshteynb@systemsbiology.org>                                                       
Date          : 11.27.02 


Copyright (C) 2006 David Shteynberg

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

David Shteynberg
Institute for Systems Biology
401 Terry Avenue North 
Seattle, WA  98109  USA
akeller@systemsbiology.org

*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef _MSC_VER  // MSVC
#include <io.h>
#include <direct.h>
#include <ctype.h>
#else
#include <unistd.h>
#include <dirent.h>
#endif
#include <string.h>

#include "Parsers/Algorithm2XML/SearchResult/SequestResult.h"
#include "Common/sysdepend.h"
#include "Common/Array.h"
#include "SequestParams.h"
#include "Parsers/mzParser/mzParser.h"
#include "Common/ModificationInfo/ModificationInfo.h"
#include "Common/ModifiedResidueMass/ModifiedResidueMass.h"
#include "Parsers/Algorithm2XML/pICalculator.h"
#include "Common/Enzyme/ProteolyticEnzyme/ProteolyticEnzymeFactory/ProteolyticEnzymeFactory.h"
#include "SequestOut.h"
#include "SequestHit.h"
#include "Parsers/Parser/Parser.h"
#include "Common/constants.h"
#include "Common/util.h"
#include "gzstream.h"

// General struct to hold .out contents which
// are needed to create an INTERACT summary line.
struct SequestOutStruct
{
   char cAA1;
   char cAA2;
   char szFileName[SIZE_FILE];
   char szBaseFileName[SIZE_FILE];
   char szProt[SIZE_PEP];
   char szPlainPep[SIZE_PEP];
   char szSubPep[SIZE_PEP];
   char szDSite[SIZE_PEP];
   char szMod[SIZE_FILE];
   char szDup[SIZE_PEP];
   char szDatabase[SIZE_FILE];
   double dAMass;
   double dMass;
   double dXC;
   double dDeltCn;
   double dSp;
   double dMass1;
   double dMass2;
   double dMass3;
   double dMass4;
   double dMass5;
   double dMass6;
   int  iRankSp;
   int  iMassType;
   int  iIon;
   int  iTot;
   int  bSpecialDeltCn;
   int  bNucDb;
} ;


// General struct to hold .out contents which
// are needed to create an INTERACT header.
struct HeaderStruct
{
   char szDate[SIZE_PEP];
   char szTime[SIZE_PEP];
   char szTimeSuffix[SIZE_PEP];
   char szMassType[SIZE_PEP];
} ;


class Out2XML {
  friend class CombineOut;
 public:
  Out2XML();
  Out2XML(char* path, int topHits, char** argv, int argc);
  void processData();
  void writeOutData();
  void writeXMLHeader();
  void writeXMLFooter();
  void readOutFile(char* fileName, struct SequestOutStruct* data, struct HeaderStruct * hdr);
  void readOutFile(char* fileName, SequestOut* data, struct HeaderStruct * hdr, SequestParams* seqParams);
  void readOutFile(FILE* ppIn, SequestOut* data, struct HeaderStruct * hdr, SequestParams* seqParams);
  ~Out2XML();
  
 private:
  void init(); // ctor helper
  ogzstream* pepXmlFile_; // write out to gzip file if filename has .gz at end
  SequestParams* sequestParams_;
  char* baseName_;
  char* baseDir_;
  char* mzXmlPath_;
  char* pepXmlPath_;
  int outFileCount_;
  int numHitsReport_;
  int iAnalysisLastScan;
  mzParser::ramp_fileoffset_t *pScanIndex;
  mzParser::RAMPFILE *pFI;
  //  Array<SequestOutStruct*> outFiles_;
  Array<SequestOut*> outFiles_;
  ProteolyticEnzyme* enzyme_;
  Boolean write_all_;
  pICalculator* pi_calc_;
  Boolean maldi_;
  int mass_type_;
} ;


#endif
