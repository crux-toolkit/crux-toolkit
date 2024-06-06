#ifndef ANAL_SUMM_PARSER_H
#define ANAL_SUMM_PARSER_H

/*

Program       : AnalysisSummary                                                  
Author        : Andrew Keller <akeller@systemsbiology.org>                                                     
Date          : 11.27.02 

Primary data object holding all mixture distributions for each precursor ion charge

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


class AnalysisSummaryParser : public Parser {

 public:

  AnalysisSummaryParser(char* xmlfile, char* analysis, char* timestamp, Boolean timestamp_only);
  ~AnalysisSummaryParser();

  void setFilter(Tag* tag);
  void print();

 protected:

  void parse(const char* xmlfile);
  const char* getAnalysisName(const char* analysis);

  Boolean found_;

  char* analysis_;
  char* timestamp_;
  Tag* summary_;
  Boolean timestamp_only_;
};











#endif
