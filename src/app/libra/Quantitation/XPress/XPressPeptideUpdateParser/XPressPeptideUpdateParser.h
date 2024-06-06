#ifndef XPRESS_PEP_UPDATE_PARSER_H
#define XPRESS_PEP_UPDATE_PARSER_H

/*
Program       : XPressPeptideUpdateParser                                                   
Author        : J.Eng and Andrew Keller <akeller@systemsbiology.org>, Robert Hubley, and 
                open source code                                                       
Date          : 11.27.02 
Version       : $Id: XPressPeptideUpdateParser.h 8870 2023-02-25 06:40:57Z real_procopio $

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
#include "Validation/MixtureModel/MixtureModel.h"


class XPressPeptideUpdateParser : public Parser {

 public:

  XPressPeptideUpdateParser(const char* xmlfile, int index, Tag* replacement);
  ~XPressPeptideUpdateParser();
  void setFilter(Tag* tag);
  Boolean update();

 protected:

  void parse(const char* xmlfile);

  char* options_;
  Boolean overwrite_;
  Boolean found_;

  ModelOptions modelOpts_;
  ScoreOptions scoreOpts_;

  Tag* replacement_;
  int index_;

};

#endif
