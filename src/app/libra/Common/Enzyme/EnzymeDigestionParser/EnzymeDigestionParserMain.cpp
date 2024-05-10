/*

Program       : EnzymeDigestionParser
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02 
SVN Info      : $Id: EnzymeDigestionParserMain.cpp 8440 2021-04-19 22:37:03Z real_procopio $

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

#include "EnzymeDigestionParser.h"

int main(int argc, char** argv) {
  hooks_tpp handler(argc,argv); // set up install paths etc

  if(argc < 2) {
    cerr << "Error: usage: EnzymeDigestionParser <xmlfile> (<enzyme>)" << endl;
    exit(1);
  }
  EnzymeDigestionParser *e;
  if(argc == 2)
    e = new EnzymeDigestionParser(argv[1]);
  else 
    e = new EnzymeDigestionParser(argv[1], argv[2]);
  delete e;
  return 0;
}
