/*
Program       : ASAPRatioPvalueParser
Authors       : Andrew Keller <akeller@systemsbiology.org>
                Xiao-jun Li (xli@systemsbiology.org>
Date          : 11.27.02
SVN Info      : $Id: ASAPRatioPvalueParser.h 7996 2019-12-25 00:16:42Z real_procopio $

Overwrites ASAPRatio protein pvalue information into ProteinProphet XML

Copyright (C) 2003 Andrew Keller, Xiao-jun Li

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

#ifndef ASAP_PVALUE_PARSER_H
#define ASAP_PVALUE_PARSER_H

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include "Parsers/Parser/Parser.h"
#include "Parsers/Parser/TagFilter.h"
#include "Quantitation/ASAPRatio/Ratio.h"
#include "Quantitation/ASAPRatio/ASAP_structs.h"
#include "Quantitation/ASAPRatio/PvalueModel/PvalueModel.h"
#include "Quantitation/ASAPRatio/Normalization/Normalization.h"

class ASAPRatioPvalueParser : public Parser {

 public:

  ASAPRatioPvalueParser();
  ASAPRatioPvalueParser(const char* xmlfile, const char *testMode);
  ASAPRatioPvalueParser(const char* xmlfile, const char* pngfile, const char *testMode);
  ~ASAPRatioPvalueParser();


 protected:

  void parse(const char* xmlfile);
  void setFilter(Tag* tag);
  void writeProteinRatio(ostream& os, proDataStrct* pro_ratio);
  pValueStrct getNormParams(double *ratios, int dataNum, const char *pngFile);
  char* getModelFile(const char* xmlfile);
  void stripWebRoot(char* filename);

  double normalize(double* ratio);
  void setModel(double* ratios); // or setModel(double** ratios);

  Boolean use_seq_ratios_; // whether to use Sequence object ratios rather than protein ratios for normalization model
  char* pngfile_;
  char *testMode_; // regression test stuff - bpratt Insilicos LLC, Nov 2005
};


#endif
