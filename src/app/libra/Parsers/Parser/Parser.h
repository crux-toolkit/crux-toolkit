/*
Program       : XML Parser
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN Info      : $Id: Parser.h 8587 2022-04-01 00:01:29Z dshteyn $


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

#ifndef PARSER_H
#define PARSER_H


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "Common/sysdepend.h"
#include "Util/RACI/RACI.h"

#ifdef _MSC_VER  // MSVC
#include <direct.h>
#else
#include <unistd.h>
#endif
#include <time.h>

#include "Tag.h"
#include <cmath>


static const int line_width_ = LINE_WIDTH;

class Parser {

 public:

  Parser();
  Parser(const char* xmlfile);
  virtual ~Parser();
  static char* getDateTime(Boolean local = True);
  const char* getName();
  Boolean success() const { // success_ is set true in overwrite()
    return success_;
  }

 protected:

  virtual Boolean overwrite(const char* xmlfile, const char* outfile, const char* last_tag);
  Boolean tag_is_at_tail(const char* xmlfile, const char* last_tag);
  void init(const char* xmlfile);
  virtual void parse(const char* xmlfile)=0; 
  virtual void setFilter(Tag* tag);
  void parse_and_print(const char* xmlfile); // just print tags, subject to filters
  Boolean filter_;
  Boolean filter_memory_;
  char* getOrigFile(char* file);
  char* getUniqueNameSpace(Tag* tag, char* orig);
  Boolean setTagValue(Tag* tag, char* tagname, char* attr_name, int* index);

  char* xmlfile_;
  char* time_; // analysis datetime
  char* name_;
  Boolean success_; // will be set to true when overwrite succeeds
};


#endif
