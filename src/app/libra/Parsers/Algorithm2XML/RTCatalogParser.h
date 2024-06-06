#ifndef _RTCATALOGPARSER_H
#define _RTCATALOGPARSER_H

/*

Program       : RTCatalog                                                       
Author        : David Shteynberg <dshteynb  AT systemsbiology.org>                                                       
Date          : 01.27.11

Primary data object holding all mixture distributions for each precursor ion charge

Copyright (C) 2011 David Shteynberg

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

*/


#include "Common/sysdepend.h"
#include <stdio.h>
#include <sys/types.h>
#include <math.h>
#include <time.h>
#include <string>
#include "Parsers/Parser/Parser.h"
#include "Common/Array.h"
#include "Parsers/Parser/Tag.h"
#include "Parsers/Algorithm2XML/RTCatalog.h"
#include "Parsers/Algorithm2XML/RTCalculator.h"

using namespace std;

typedef TPP_STDSTRING_HASHMAP(RTCalculator*) rtcalc_hash;

class RTCatalogParser : public Parser {

 public:

  RTCatalogParser(double minProb, double minPTMProb, int minRT, int maxRT, bool gradientCorr, bool tables, bool xics, bool matchedions,
		  string* ignoreFile = NULL, string * acnFile = NULL,  string * pepsbyrunFile = NULL,
		  string * gradPepsFile = NULL,  string * iRTsFile = NULL, string* worklistFile = NULL, 
		  int minGradPeps = 0);
  void addFile(const char* filename);
  ~RTCatalogParser();


  void parseSPTXT(const char* c);
  void parseTables(const char* c);
  void parse(const char* c);

  void parsePeptideQ1s(const char* c);

  void parsePeptidesByRun(const char* c);
  

  void sortRunNames();

  void writeRTCatalog(const char* file) ;

  void trackPeptidesAcrossRuns(string* chromFile, string* trackFile) ;

  void trackPeptidesReport(string* file);

  void chromPeptidesReport(string* file);
  str_hash* getPrevRTMix() { return  prev_rtmix_byrun_; }
  str_hash* getNextRTMix() { return  next_rtmix_byrun_; }
  void correctRTsByAdjacentRun();

  void correctRTsByRun();

  void calcRTsByRun();

  void setRTMixChromatograms(str_hash*);

  void gradientCorrectReport(ostream& out);


  void useXICs(bool);

  void setMods(bool, bool);
  void setTolerances(float ppm, float dal);
  void setDisco(bool disco) { disco_ = disco; rt_cat_->setDisco(disco_); }
  void setUpRTMixRuns();
  void computeRTMixRuns();
  void processRTMixChromatograms(string* , string* );
  
 protected:

  //void displayOptions(char* eng);

  RTCatalog* rt_cat_;
  
  Array<string*>* input_files_;

  str_hash* prev_rtmix_byrun_;

  str_hash* next_rtmix_byrun_;

  str_hash* rtmix_chroms_;

  //int_hash* ms_runs_;
  //  Array<RTCalculator*>* byrun_rt_calcs_;

  rtcalc_hash* byrun_rt_calcs_;

  Array<string*>* ignore_runs_;

  Array<Array<string*>*>* worklist_runs_;

  dblarr_hash_hash* rtmix_peps_q1q3_;

  dblarr_hash_hash* byrun_pep_q1s_;

  bool_hash_hash* byrun_peps_;

  Array<string*>* track_peps_;

  bool_hash* grad_peps_hash_;

  dbl_hash* irt_peps_hash_;

  bool_hash* ok_runs_; //Array<string*>* ok_runs_;
  
  dblarr_hash* track_pepq1q3_hash_;

  
  bool gradCorr_;
  bool tablesIn_;

  bool XICs_;

  bool matchedions_;

  float XICtol_;


  bool pv_mods_;
  bool cys_cam_;

  bool disco_;

  int minGradPeps_;

  double minProb_;

  double minPTMProb_;
  
  int run_idx_;
  
};











#endif
