/*
Program       : QuanticProteinRatioParser
Authors       : David Shteynberg
Date          : 12.02.20
SVN info      : $Id: QuanticProteinRatioParser.h 

Computes Quantic ratios and errors for proteins, then overwrites
that information onto ProteinProphet XML

Copyright (C) 2020 David Shteynberg

Based on ASAPRatioProteinRatioParser
Copyright (C) 2003 Andrew Keller, Jimmy Eng

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

#ifndef QUANTIC_PROT_RATIO_PARSER_H
#define QUANTIC_PROT_RATIO_PARSER_H

#include <sstream>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "Parsers/Parser/Parser.h"
#include "Parsers/Parser/TagFilter.h"
#include "QuanticGroupPeptideParser.h"
#include "QuanticProteinRatio.h"
#include "QuanticMSMSSeqData.h"


#define SIZE_PEPTIDE   128
#define SIZE_BUF      8192

using namespace std;
class QuanticProteinRatioParser : public Parser {

 public:

  QuanticProteinRatioParser(const char* protxmlfile, const char *testMode);
  ~QuanticProteinRatioParser();
  //void cachePepXML(); 
 protected:

  void parse(const char* protxmlfile);
  void setFilter(Tag* tag);
  char* getPeptideString(Array<const char*>* peps, const char* link);
  void getRatio(Array<UniqPeptide*>* peps, double minpepprob, double minwt);
  void setInputFiles(char** inputfiles);
  UniqPeptide* enterUnique(Array<UniqPeptide*>* uniques, const char* next, double wt, double prob, int charge);
  Array<Tag*>* getProteinRatioTags(proQuantStrct* pro_ratio);
  Array<Tag*>* getPeptideRatioTags(proQuantStrct* pro_ratio, string modpep, int chg);

  int countProteinRatioN(proQuantStrct* pro_ratio); 

  QuanticGroupPeptideParser* parser_;
  Array<const char*> input_pepxmlfiles_;

  int iNumRawData_;
  proQuantStrct* pro_ratio_;

  char *testMode_;

  bool iprophet_;
};


#endif
