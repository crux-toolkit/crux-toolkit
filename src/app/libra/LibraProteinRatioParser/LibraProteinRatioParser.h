/*
Program       : LibraProteinRatioParser
Author        : Andrew Keller <akeller@systemsbiology.org>
                *Jimmy Eng (jeng@systemsbiology.org>
Date          : 11.27.02 
SVN info      : $Id: LibraProteinRatioParser.h 8873 2023-02-27 08:14:30Z real_procopio $

Computes LIBRA ratios and errors for proteins, then overwrites
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


#ifndef LIB_PROT_RATIO_PARSER_H
#define LIB_PROT_RATIO_PARSER_H

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string>
#include<iostream>
#include<fstream>

#include "Parsers/Parser/Parser.h"
#include "Parsers/Parser/TagFilter.h"
//#include "Quantitation/ASAPRatio/Ratio.h"
#include "LibraGroupPeptideParser.h"
#include "LibraConditionHandler2.h"
#include "StringConvertor.h"

#define SIZE_PEPTIDE   128
#define SIZE_BUF      8192


class LibraProteinRatioParser : public Parser {

 public:
  LibraProteinRatioParser(const char * xmlfile, const char * conditionFile, char *testArg);
 ~LibraProteinRatioParser();

  char* conditionFileName;
  double minimumThreshholdIntensity;
  std::string quantitationFileName;

 protected:
  void initializeQuantitationFile( int nChannels);
  void parse(const char* xmlfile);
  void setFilter(Tag* tag);
  char * getPeptideString(Array<const char *>* peps, const char * link);
  void getRatio(const char* prot, Array<const char*>* peps, double minpepprob);
  void setInputFiles(const char ** inputfiles);
  void enterUnique(Array<const char*>* uniques, const char* next);

  LibraGroupPeptideParser* libraGroupPeptideParser_;
  Array<const char*>* input_xmlfiles_;

  // user selected channel to use as reference for normalization
  // (0 means no normalization).
  int norm_channel_;
  int num_channels_;

  Array<Tag*>* summary_tags_;
  Array<Tag*>* result_tags_;

  LibraConditionHandler2* pLibraConditionHandler2;

  char *testMode_;     // regression test stuff - bpratt Insilicos LLC, Nov 2005
};

#endif
