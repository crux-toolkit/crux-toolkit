/*

Program       : XPressCGIParser                                                 
Author        : Andrew Keller <akeller@systemsbiology.org> 
                Jimmy Eng (jeng@systemsbiology.org>                                                      
Date          : 11.27.02 

Overwrites specified modified XPRESS protein ratio onto
ProteinProphet XML

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

#ifndef X_CGI_PARSER_H
#define X_CGI_PARSER_H

#include <stdio.h>

#include "Parsers/Parser/Parser.h"
#include "Parsers/Parser/TagFilter.h"
#include "Parsers/Parser/Tag.h"


#define SIZE_BUF      8192

class XPressCGIParser : public Parser {

 public:

  XPressCGIParser(const char* xmlfile, const char* protein, double ratio, double error, double h2l_ratio, double h2l_error, int numpeps);

 protected:

  void parse(const char* xmlfile);
  void setFilter(Tag* tag);
  char* protein_;
  char ratio_[25];
  char error_[25];
  char h2l_ratio_[25];
  char h2l_error_[25];
  char numpeps_[25];

};









#endif
