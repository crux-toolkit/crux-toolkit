#ifndef DATABASE_PARSER_H
#define DATABASE_PARSER_H

/*

Program       : DatabaseParser                                                   
Author        : Andrew Keller <akeller@systemsbiology.org>                                                       
Date          : 11.27.02 
SVN Info      : $Id: DatabaseParser.h 8800 2023-01-06 09:26:53Z real_procopio $

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
#include <string>

#include "Parsers/Parser/Parser.h"
#include "Parsers/Parser/TagFilter.h"


class DatabaseParser : public Parser {

 public:

  DatabaseParser(const char* xmlfile);
  int getNumDatabases();
  std::string getDatabases();

 protected:

  void parse(const char* xmlfile);
  Boolean enter(char* db);

  Array<char*>* databases_;

};











#endif
