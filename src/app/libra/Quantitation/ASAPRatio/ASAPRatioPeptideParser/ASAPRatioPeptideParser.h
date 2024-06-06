/*
Program       : ASAPRatioPeptideParser
Author        : X.Li and Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN Info      : $Id: ASAPRatioPeptideParser.h 8025 2020-02-14 01:05:59Z mhoopmann $


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

#ifndef ASAP_PEP_PARSER_H
#define ASAP_PEP_PARSER_H


#include <stdio.h>
#include <math.h>
#include <time.h>

# include <stdlib.h>
# include <string.h>
# include <ctype.h>
//# include <sys/param.h>
# include <sys/stat.h> 
#include "Common/constants.h"

#include "Parsers/Parser/Parser.h"
#include "Parsers/Parser/TagFilter.h"
#include "Quantitation/Option.h"
#include "Parsers/mzParser/mzParser.h"
#include "Common/ResidueMass/ResidueMass.h"
#include "Common/ModificationInfo/ModificationInfo.h"
#include "Common/spectStrct.h"

#define PROGRAM_VERSION "3.0"
#define PROGRAM_AUTHOR "Xiaojun Li"

#define SIZE_BUF         8192
#define FILTER_SIZE      6

#define TRUE             1
#define FALSE            0

# define _MZXML_READER_MXSCANNUM_ 100000
# define _ASAPRATIO_HM_ 1.0079407539 // average mass of H
# define _ASAPRATIO_M_ 1.0033548378 // mass diff. bewteen isotopes
# define _ASAPRATIO_MZBOUND_ 0.5 // bound of an isotopic peak
# define _ASAPRATIO_ISONUM_ 3 // number of isotopes covered
# define _ASAPRATIO_EXPLCRANGE_ 50 // selected range of a LC spectrum

# define NATIVEAANUM 24 // number of native amino acids + N-end + C-end
# define MODAANUM 20 // total number of modified amino acids

# define _ASAPRATIO_CONFL_ 0.6826 // confidence level at 1 sigma

//
// Structures
//

typedef struct {
  char site;
  double lightmass;
  double heavymass;
} staticQuant;


// structure of XML index
typedef struct {
  mzParser::RAMPFILE *mzfile;
  mzParser::ramp_fileoffset_t *scanIndx;
  int totScan; // 1 .. totScan
} xmlIndxStrct;


// structure of LC spectra
class lcSpectStrct;


// structure of amino acid residue
typedef struct {
  int sz; // 1 or 2 letter represent
  char rp[3]; // represent: 1 letter for native, 1-2 for modified 
  double ms; // mono-isotopic mass
} residueStrct;


// native amino acid residues
const static residueStrct nativeAA[] = {
  {1, "n", 1.0078250321},  // n terminus
  {1, "c", 17.0027396542}, // c terminus
  {1, "J", 1.0078250321},
  {1, "U", 17.0027396542},
  {1, "G", 57.0214637236},
  {1, "A", 71.0371137878},
  {1, "V", 99.0684139162},
  {1, "L", 113.0840639804},
  {1, "I", 113.0840639804},
  {1, "S", 87.0320284099},
  {1, "C", 103.0091844778},
  {1, "T", 101.0476784741},
  {1, "M", 131.0404846062},
  {1, "P", 97.0527638520},
  {1, "F", 147.0684139162},
  {1, "Y", 163.0633285383},
  {1, "W", 186.0793129535},
  {1, "H", 137.0589118624},
  {1, "K", 128.0949630177},
  {1, "R", 156.1011110281},
  {1, "D", 115.0269430320},
  {1, "E", 129.0425930962},
  {1, "N", 114.0429274472},
  {1, "Q", 128.0585775114}
};


// structure of pairing amino acid
typedef struct {
  char prtnA[3];
  char prtnB[3];
} pairStrct;


int modAA_size_compare(void const *a, void  const *b);  
int pairStrctCmp(void const *a, void const *b); 
int comp_chars(void const *a, void const *b);


class ASAPRatioPeptideParser : public Parser {

 public:
  ASAPRatioPeptideParser(const char* xmlfile, const InputStruct &options, const char *testMode, double mzBound = _ASAPRATIO_MZBOUND_);
  ~ASAPRatioPeptideParser();
  void setFilter(Tag* tag);

 protected:

  void parse(const char* xmlfile);
  void getRatio();
  void setLabelStrings(Array<char*>* lightpartners, char* lightstring, Array<char*>* heavypartners, char* heavystring);
  Tag* getSummaryTag(const InputStruct& opts);

  Array<Tag*>* generateXML(int index, Boolean cover, Boolean data);

#ifdef USE_STD_MODS
  void evalPepDataStrct(char *pepSeq, long scan, int chrg, char *xmlFile,
			int eltn, int areaFlag, double pepmass, double computed_peptide_mass);
#endif
#ifndef USE_STD_MODS
  void evalPepDataStrct(char *pepSeq, long scan, int chrg, char *xmlFile,
			int eltn, int areaFlag, residueStrct *modAAs, int modAANum,
			pairStrct *prtnAAs, int prtnAANum);
#endif
  void getPairSequences(char *sequence, char *prtnSequence, 
			int *cidIndx, double *msLight, double *msHeavy, 
			residueStrct *modAAs, int modAANum,
			pairStrct *prtnAAs, int prtnAANum);

  double getSeqMonoMZ(char *sequence, int qState,
		      residueStrct *modifiedAA, int modAANum);

  void getPepDataStrct(pepDataStrct *data, char *xmlFile,
		       char *pngFileBase, int cgiIndx);

  int getPeakSpectra(lcSpectStrct *lcSpect, pepDataStrct *data, 
		     double timeWd, char *xmlFile);
		     
  xmlIndxStrct *createXmlIndx(const char *xmlFile); // make a new one
  xmlIndxStrct *getXmlIndx(const char *xmlFile); // return a possibly cached one
  void freeXmlIndx();

  spectStrct *getCombLCSpect(int firstScan, int lastScan, 
			     double *mzArray, int arraySize, double dMz, 
			     xmlIndxStrct *xmlIndx);

  int getScanNum(double time, xmlIndxStrct *xmlIndx);

  void getLCPeakStrct(lcPeakStrct *peak, const spectStrct &rawSpectrum, 
		      const spectStrct &fitSpectrum, int scan, int scanRange);

  void getLCSpectPeakAndValleys(int *peakScan, int *valleyScan, 
				const spectStrct &spectrum, 
				int scan, double background); 

  int spectPeak(const spectStrct &spectrum, int position, int direction);

  int findNextValley(const spectStrct &spectrum, int position, int direction);

  double getLCSpectBackground(const spectStrct &rawSpectrum, 
			      const spectStrct &fitSpectrum, 
			      int *ranges, int *valleyScan);

  int getLCSpectAreaAndTime(double *area, double *peakTime, 
			    const spectStrct &rawSpectrum, 
			    const spectStrct &fitSpectrum, 
			    int pkScan, int *valleyScan, double background);

  void plotPepPeak(const char *pngFileBase, const pepDataStrct &data, 
		   const lcSpectStrct &lcSpect, int qState);

  void getRidOfSpace(char *string);
  void smoothSpectFlx(spectStrct * spectrum, double range, int repeats);
  void smoothSpectFlx(spectStrct * spectrum, double range, int repeats, int smoothWindow);
  double smoothDataPtFlx(const spectStrct &spectrum, int dtIndx, 
			 double range, double threshold);
  double smoothDataPtFlx(const spectStrct &spectrum, int dtIndx, 
			 double range, double threshold, int smoothWindow);
  int getPeakSpectra(lcSpectStrct *lcSpect, pepDataStrct *data, 
		     double timeWd, char *xmlFile,int smoothItrs,
		     double smoothRTwindow);

 double PadeApprx(double x, double *xa, double *ya, int size);
 void DixonTest(double *data, int *outliers, int size);
 void findMeanAndStdDevWeight(double *mean, double *error,
			     double *data, double *h2l_mean, double *h2l_error,
			     double *h2l_data, double *weight, int size);
 void getDataRatio(double *ratio, double *error, double *h2l_ratio, double *h2l_error, 
		   double confL, 
		   double *data, double *dataErrs, 
		   double *dataWghs, int *dataIndx, int dataSize,
		   int testType);

 void myLUBksb(double **mtrx, int order, int *indx, double *vec);
 void myLUDcmp(double **mtrx, int order, int *indx, double *d);
 void string2AA(char *string, residueStrct *aa);

 pairStrct *collectPrtnAAStrct(int *prtnAANum, Array<char*>* lightStrings, Array<char*>* heavyStrings);
 residueStrct *collectModAAStrct(int *modAANum, Array<char*>* inputStrs);
 char **getStrSects(int *sectNum, char *string, char sep);
 void freeMtrx(void **mtrx, int size);

 int getStaticQuantStatus(char site, double mass, double error);
 int diff(double first, double second, double error);
 void enterStaticQuant(char site, double mass, double error);
 void setStaticPartners(Array<char*>* lightStrings, Array<char*>* heavyStrings);
 void setStaticQuantInfo(const char* xmlfile);
 char** getStaticModification(char site, char* symbol, int status);
 double getUserSpecifiedEquivalent(const char & aa, double avgmass);
 double getMonoisotopicEquivalent(double averagemass, double error, Tag* tag);

 ModelOptions modelOpts_;
 ScoreOptions scoreOpts_;

 FILE* fp_;
 long* index_;
 char mzXMLfile_[2000];

 char *testMode_; // regression test stuff - bpratt Insilicos LLC, Nov 2005

 long pepIdx_;

 int elution_;
 int areaFlag_;

 residueStrct* modAAs_;
 xmlIndxStrct* xmlIndx_;
 
 int modAANum_;
 int prtnAANum_;
 pairStrct* prtnAAs_;

 pepDataStrct data_;

 char paired_labels_[100]; // these are the user specified light/heavy labeled amino acids (or N/C termini)

 InputStruct pInput_;

 spectStrct* spectrum_;

 Array<staticQuant*>* static_pairs_;
 int static_status_;
 Boolean static_quant_; // whether or not quant is static (i.e. all run is either light or heavy static mod)
 double error_;
 char static_symbol_;
 Boolean static_nterm_;
 Boolean static_cterm_;
 Array<double>* average_mods_;
 Array<double>* monoisostopic_mods_;
 int last_valley_[2];
 double diff_; // error for correction of average to monoisotopic masses
 Boolean monoisotopic_;
#ifdef USE_STD_MODS
  ModificationInfo* modinfo_;
  double light_label_masses_[26];
  double heavy_label_masses_[26];
  double light_nterm_mass_;
  double heavy_nterm_mass_;
  double light_cterm_mass_;
  double heavy_cterm_mass_;
  //  double light_mass_;
  lcSpectStrct* lcSpect_;
  Boolean compute_peptide_mass_;
#endif
  double mzBound_;
  bool verbose_;
  int m_XMLfile_state;
  std::string m_lastMZXML;
};



#endif
