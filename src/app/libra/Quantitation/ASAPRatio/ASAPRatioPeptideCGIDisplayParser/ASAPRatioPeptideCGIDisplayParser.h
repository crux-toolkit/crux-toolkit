/*
Program       : ASAPRatioPeptideParser
Author        : X.Li and Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN Info      : $Id: ASAPRatioPeptideCGIDisplayParser.h 7996 2019-12-25 00:16:42Z real_procopio $


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

#ifndef ASAP_PEP_DISPLAY_PARSER_H
#define ASAP_PEP_DISPLAY_PARSER_H


#include <stdio.h>
#include <math.h>
#include <time.h>

#include "Parsers/Parser/Parser.h"
#include "Parsers/Parser/TagFilter.h"
#include "Quantitation/Option.h"


class ASAPRatioPeptideCGIDisplayParser : public Parser {

 public:

  ASAPRatioPeptideCGIDisplayParser(const char* xmlfile, const char* basename, const char* timestamp, int index, int elution, Boolean zeroBG = False, double mzBound = -1);
  ~ASAPRatioPeptideCGIDisplayParser();
  pepDataStrct getPepDataStruct();
  void setPepDataStruct(Array<Tag*>* tags, int elution, long scan, int precursor_charge);
  void setFilter(Tag* tag);
  Boolean found();

 protected:

  void parse(const char* xmlfile);

  Boolean found_;

  Boolean zeroBG_;

  int index_;
  int asap_index_;
  int elution_;
  //char* asaptimestamp_;
  pepDataStrct data_;
  char* basename_;
  char* timestamp_;
  double mzBound_;
};


#endif
