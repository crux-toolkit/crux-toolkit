#ifndef LIBRA_PEP_PARSER_H
#define LIBRA_PEP_PARSER_H

/*
Program       : LibraPeptideParser
Author        : Patrick Pedrioli and Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: LibraPeptideParser.h 8748 2022-10-18 10:00:57Z real_procopio $

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
#include "Parsers/mzParser/cached_ramp.h" // deals with inefficient use of RAMP file open etc
#include "Common/constants.h"

#include "LibraConditionHandler.hpp"
#include "LibraWrapper.hpp"
#include "LibraSummary.hpp"
#include "LibraResult.hpp"

#define PROGRAM_VERSION "1.0"
#define PROGRAM_AUTHOR "P.Pedrioli"


class LibraPeptideParser : public Parser {

 public:

  LibraPeptideParser(const char* xmlfile, const char* conditionFileName, const char *testMode);
  ~LibraPeptideParser();
  void setFilter(Tag* tag);

 protected:

  void parse(const char* xmlfile);
  char mzXMLfile_[10000];
  LibraConditionHandler* condition_;
  char *testMode_;
  LibraSummary* libra_summary_;
  LibraResult* libra_result_;
  LibraWrapper* libra_quantifier_;

};

#endif
