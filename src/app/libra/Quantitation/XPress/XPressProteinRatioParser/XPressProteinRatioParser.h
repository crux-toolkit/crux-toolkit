/*

Program       : XPressProteinRatioParser
Author        : Andrew Keller <akeller@systemsbiology.org>
                *Jimmy Eng (jeng@systemsbiology.org>
Date          : 11.27.02
SVN Info      : $Id: XPressProteinRatioParser.h 7705 2017-12-15 21:43:12Z real_procopio $

Computes XPRESS ratios and errors for proteins, then overwrites
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


#ifndef X_PROT_RATIO_PARSER_H
#define X_PROT_RATIO_PARSER_H

#include <stdio.h>
#include <math.h>
#include <time.h>

#include "Parsers/Parser/Parser.h"
#include "Parsers/Parser/TagFilter.h"
#include "../XPressGroupPeptideParser.h"

#define SIZE_PEPTIDE   128
#define SIZE_BUF      8192



class XPressProteinRatioParser : public Parser {

 public:

  XPressProteinRatioParser(const char * protxmlfile, const Boolean use_intensities, const char *testMode);
  XPressProteinRatioParser(Array<const char*> &input_pepxmlfiles, const peplist &peptides,double minpepprob ); // used by XPressCGIProteinDisplay

  ~XPressProteinRatioParser();
  const RatioStruct &getRatio() const { // used by XPressCGIProteinDisplay
     return pRatio_;
  }

 protected:

  void parse(const char * protxmlfile);
  void setFilter(Tag* tag);
  void cachePepXML(); // read the search hits out of the pepXML files

  char * getPeptideString(peplist* peps, const char * link);
  Boolean getRatio(const peplist* peps, double minpepprob);
  void setInputFiles(const char ** inputfiles);
  void enterUnique(peplist* uniques, const char* next);

  XPressGroupPeptideParser* parser_;
  Array<const char *> input_pepxmlfiles_;
  XPressRatioSearchHitCache *searchHits_; // caches the pepXML reads

  char *testMode_; // regression test stuff - bpratt Insilicos LLC, Nov 2005 

  int iNumRawData_;
  RatioStruct pRatio_;
  Boolean heavy2light_;
  Boolean calc_intensity_ratios_;
};


#endif
