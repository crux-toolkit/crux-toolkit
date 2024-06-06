/*
Program       : ASAPRatioPeptideParser
Author        : X.Li and Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN Info      : $Id: ASAPRatioPeptideUpdateParser.h 7996 2019-12-25 00:16:42Z real_procopio $


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

#ifndef ASAP_PEP_UPDATE_PARSER_H
#define ASAP_PEP_UPDATE_PARSER_H



#include <stdio.h>
#include <math.h>
#include <time.h>

#include "Parsers/Parser/Parser.h"
#include "Parsers/Parser/TagFilter.h"
#include "Quantitation/Option.h"
//#include "Validation/MixtureModel/MixtureModel.h"


class ASAPRatioPeptideUpdateParser : public Parser {

 public:

  ASAPRatioPeptideUpdateParser(const char* xmlfile, const char* basename, int index, Array<Tag*>* replacements);
  ~ASAPRatioPeptideUpdateParser();
  void setFilter(Tag* tag);
  Boolean update();

 protected:

  void parse(const char* xmlfile);

  char* options_;
  Boolean overwrite_;
  Boolean found_;

  Array<Tag*>* replacements_;
  int index_;
  char* basename_;
};



#endif
