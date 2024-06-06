/*
Program       : ASAPRatioProteinRatioParser
Authors       : Andrew Keller <akeller@systemsbiology.org>
                *Jimmy Eng (jeng@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: ASAPRatioProteinRatioParser.h 7996 2019-12-25 00:16:42Z real_procopio $

Computes ASAPRatio ratios and errors for proteins, then overwrites
that information onto ProteinProphet XML

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

#ifndef ASAP_PROT_RATIO_PARSER_H
#define ASAP_PROT_RATIO_PARSER_H


#include <stdio.h>
#include <math.h>
#include <time.h>

#include "Parsers/Parser/Parser.h"
#include "Parsers/Parser/TagFilter.h"
#include "Quantitation/ASAPRatio/Ratio.h"
#include "ASAPRatioGroupPeptideParser.h"
#include "Quantitation/ASAPRatio/ASAP_structs.h"
#include "Quantitation/ASAPRatio/UniquePeptide/UniquePeptide.h"
#include "Quantitation/ASAPRatio/ASAPRatioPvalueParser/ASAPRatioPvalueParser.h"

#define SIZE_PEPTIDE   128
#define SIZE_BUF      8192


class ASAPRatioProteinRatioParser : public Parser {

 public:

  ASAPRatioProteinRatioParser(const char* protxmlfile, const char *testMode);
  ~ASAPRatioProteinRatioParser();
  //void cachePepXML(); 
 protected:

  void parse(const char* protxmlfile);
  void setFilter(Tag* tag);
  char* getPeptideString(Array<const char*>* peps, const char* link);
  void getRatio(Array<UniquePeptide*>* peps, double minpepprob, double minwt, Boolean heavy2light);
  void setInputFiles(char** inputfiles);
  void enterUnique(Array<UniquePeptide*>* uniques, const char* next, double wt, double prob);
  Array<Tag*>* getProteinRatioTags(proDataStrct* pro_ratio);

  ASAPRatioGroupPeptideParser* parser_;
  Array<const char*> input_pepxmlfiles_;

  int iNumRawData_;
  proDataStrct* pro_ratio_;
  Boolean heavy2light_;
  char *testMode_;
};


#endif
