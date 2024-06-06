/*
Program       : ASAPRatioProteinCGIDisplayParser
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: ASAPRatioProteinCGIDisplayParser.h 7996 2019-12-25 00:16:42Z real_procopio $


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

#ifndef ASAP_PRO_CGI_PARSER_H
#define ASAP_PRO_CGI_PARSER_H


#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#include "Parsers/Parser/Parser.h"
#include "Parsers/Parser/TagFilter.h"
#include "Common/constants.h"
#include "Quantitation/ASAPRatio/ASAP_structs.h"
#include "Parsers/Algorithm2XML/SearchResult/ProDataComponent.h"
#include "Parsers/Algorithm2XML/SearchResult/SearchResult.h"
#include "Parsers/Algorithm2XML/SearchResult/SequestResult.h"
#include "Parsers/Algorithm2XML/SearchResult/MascotResult.h"
#include "Parsers/Algorithm2XML/SearchResult/CometResult.h"
#include "Parsers/Algorithm2XML/SearchResult/TandemResult.h"
#include "Parsers/Algorithm2XML/SearchResult/PhenyxResult.h"
#include "Quantitation/ASAPRatio/ASAPRatio_Fns/ASAPRatio_numFns.h"
#include "Validation/MixtureDistribution/MixtureDistrFactory.h"

int OrderByXmlAndDataInds(void const *a, void const *b);
int OrderBySeqPkDataInds(void const *a, void const *b);

class ASAPRatioProteinCGIDisplayParser : public Parser {

 public:

  ASAPRatioProteinCGIDisplayParser(proDataStrct* protein, Array<char*>* xmlfiles, Boolean heavy2light, char* colored_aas);
  ~ASAPRatioProteinCGIDisplayParser();
  void setFilter(Tag* tag);

  void initMSMSRunIdx();
  void setPepDataStruct(Array<Tag*>* tags, int elution, long scan, int precursor_charge);
  proDataStrct* getProDataStruct();
  char* display(int seq, int pk, int data);
  int getResultIndex(int seq, int pk, int data, Boolean advance);
  void write(ostream& os);

 protected:

  void parse(const char* xmlfile);
  void cleanup(Array<char*>* data);
  int getTimestampIndex(const char* timestamp);
  SearchResult* getSearchResult(Array<Tag*>* tags, char* engine);
  pepDataStrct getPepDataStrct(int seq, int pk, int data);
  void makeRadioTag(char* tag, int value, char* html_name);
  void setDataRadioTag(char* tag, int seq, int pk, int data, int value);
  void setPeakRadioTag(char* tag, int seq, int pk, int value);
  void setSeqRadioTag(char* tag, int seq, int value);
  void setASAPRatioTag(char* text, double mean, double error);
  void update();
  double PadeApprx(double x, double *xa, double *ya, int size);
  void DixonTest(double *data, int *outliers, int size);
  void findMeanAndStdDevWeight(double *mean, double *error,
			       double *data, double *inv_mean, double *inv_error,
			       double *inv_data, double *weight, int size);
  void getDataRatio(double *ratio, double *error,  double *inv_ratio, double *inv_error, double confL, 
		    double *data, double *dataErrs, double *inv_data, double *inv_dataErrs,
		    double *dataWghs, int *dataIndx, int dataSize,
		    int testType);
  void updatePeakStrctRatio(int seq, int pk);
  void updateSeqStrctRatio(int seq);


  int current_index_;

  int result_index_;
  int pepdata_index_;

  proDataStrct* protein_;

  Array<char*>* xmlfiles_;

  Array<char*>* databases_;
  Array<char*>* basenames_;
  Array<char*>* pepproph_timestamps_;  
  Array<char*>* iproph_timestamps_;
  Array<char*>* asap_timestamps_;
  Array<Boolean>* asap_quantHighBGs_;
  Array<Boolean>* asap_zeroBGs_;
  Array<bool>* asap_wavelets_;
  Array<double>* asap_mzBounds_;

  Array<int>* elutions_;

  Array<char*>* aa_modifications_;
  Array<char*>* term_modifications_;
  Array<char*>* misc_run_conditions_; // only used for COMET

  Array<ProDataComponent*>* components_;

  pepDataStrct data_;
  Boolean heavy2light_;
  char* colored_aas_;
};


#endif
