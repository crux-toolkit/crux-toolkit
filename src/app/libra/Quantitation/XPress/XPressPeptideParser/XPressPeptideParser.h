#ifndef XPRESS_PEP_PARSER_H
#define XPRESS_PEP_PARSER_H

/*

Program       : XPressPeptideParser
Author        : J.Eng and Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN Info      : $Id: XPressPeptideParser.h 8025 2020-02-14 01:05:59Z mhoopmann $

added mzData support - Brian Pratt Insilicos LLC 2005


Copyright (C) 2003 Andrew Keller

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

Andrew Keller
Institute for Systems Biology
401 Terry Avenue North 
Seattle, WA  98109  USA
akeller@systemsbiology.org

*/

#include <stdio.h>
#include <math.h>
#include <time.h>

#include "Parsers/Parser/Parser.h"
#include "Parsers/Parser/TagFilter.h"
#include "Quantitation/Option.h"
//#include "Validation/MixtureModel/MixtureModel.h"
#include "Parsers/mzParser/mzParser.h"
#include "Common/constants.h"
#include "Common/ResidueMass/ResidueMass.h"
#include "Common/ModificationInfo/ModificationInfo.h"

//#define PROG_NAME "xpress"

#define PROGRAM_VERSION "2.2"   // 2.2 supports recursive expanding scan bounds as needed
#define PROGRAM_AUTHOR "Jimmy Eng"

#define HTTP_TARGET     "Win1"
#define DEFAULT_TOL      1.0
#define SIZE_BUF         8192
#define FILTER_SIZE      4
#define PROTON_MASS      1.00727646688

#define TRUE             1
#define FALSE            0


struct XpressStruct
{
   int    iLightFirstScan;
   int    iLightLastScan;
   int    iHeavyFirstScan;
   int    iHeavyLastScan;
   int    iMetabolicLabeling;
   int    bXpressLight1;
   int    iLightIntensityScan;
   int    iHeavyIntensityScan;
   double dLightPeptideMass;
   double dHeavyPeptideMass;
   double dLightArea;
   double dHeavyArea;
   double dLightIntensity;
   double dHeavyIntensity;
   double dIntensityRatio;
   double dLightIntensityRT;
   double dHeavyIntensityRT;
   double dLightFirstScanRT;
   double dLightLastScanRT;
   double dMassTol;
   char   szQuan[128];

   char   szXMLFile[SIZE_FILE];
   int    iChargeState;
   int    iMinNumChromatogramPoints;
   int    iMinNumIsotopePeaks;
   char   szOutFile[SIZE_FILE];
};


class  XPressPeptideParser:public Parser
{

 public:

   XPressPeptideParser();
   XPressPeptideParser(const char *xmlfile, const InputStruct & options, const char *testMode);
   virtual ~XPressPeptideParser();
   void   setFilter(Tag * tag);
   void validate_mzXMLfile(); // verify scan data source

 protected:
   void readAllModMasses(const char *xmlfile);
   void   parse(const char *xmlfile);
   Tag   *getRatio();

   Tag   *getSummaryTag(const InputStruct &opts);
   void   XPRESS_ANALYSIS(Boolean bLightPeptide, int iScanLocalBuffer=0);
   void   XPRESS_ANALYSIS_LABELFREE(int iScanLocalBuffer=0);

   //   void   FILTER_MS(double *dOrigMS, double *dFilteredMS, int iStart, int iEnd);

   void   FILTER_MS(const double *dOrigMS, Array<double> & dFilteredMS, int iStart, int iEnd);
   //void FILTER_MS(Array<double> *dOrigMS,
   //             Array<double> *dFilteredMS);

   void   FIND_ENDPOINTS(Array<double> & pdFiltered, //const double *pdFiltered,
                         int &iPepStartScan,
                         int &iPepEndScan,
                         int iAnalysisFirstScan,
                         int iAnalysisLastScan);

   void   FIND_ENDPOINTS_FIX(Array<double> & pdFiltered, //const double *pdFiltered,
                             int &iPepStartScan,
                             int &iPepEndScan,
                             int iAnalysisFirstScan,
                             int iAnalysisLastScan);

   void   flipRatio(char *ratio, char *flipped);
   double ELEMENT_COUNT(char *szPeptide,
                        char cElement);

   ModelOptions modelOpts_;
   ScoreOptions scoreOpts_;

   char  *testMode_;    // regression test stuff - bpratt Insilicos LLC, Nov 2005

   mzParser::RAMPFILE *fp_;
   mzParser::ramp_fileoffset_t *index_;
   int *piSequentialScan;
   int *piReverseSequentialScan;
   struct ScanCacheStruct* pCache_;
   char   mzXMLfile_[1000];
   int m_XMLfile_state;

   InputStruct pInput_;
   XpressStruct pXpress_;

   double //*dTmpFilter_,
     *dLightMS_,
     *dHeavyMS_;
   
   Array<double> * dLightFilteredMS_;
   Array<double> * dHeavyFilteredMS_;  

   Array<double> *dTmpFilter_;

#define MAX_ISOTOPES 2
   double *dIsotopeLight[MAX_ISOTOPES+1],
          *dIsotopeHeavy[MAX_ISOTOPES+1];


#ifdef USE_STD_MODS
   ModificationInfo *modinfo_;
   double  light_label_masses_[26];
   double  heavy_label_masses_[26];
   double  light_nterm_mass_;
   double  heavy_nterm_mass_;
   double  light_cterm_mass_;
   double  heavy_cterm_mass_;
   Boolean monoisotopic_;
#endif

};


#endif
