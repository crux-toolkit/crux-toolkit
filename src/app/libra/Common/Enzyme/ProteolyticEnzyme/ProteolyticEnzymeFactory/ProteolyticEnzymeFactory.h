#ifndef PROT_ENZ_FAC_H
#define PROT_ENZ_FAC_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <ctype.h>
#include <string>
#include "Parsers/Parser/Tag.h"
using namespace std;

#include "Common/Enzyme/ProteolyticEnzyme/ProteolyticEnzyme.h"

/*

Program       : EnzymeSpecificity for PeptideProphet
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN Info      : $Id: ProteolyticEnzymeFactory.h 8309 2020-12-03 03:17:24Z dshteyn $

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


class ProteolyticEnzymeFactory {
 public:
  ProteolyticEnzymeFactory();
  ProteolyticEnzymeFactory(int min_spacing);
  ProteolyticEnzyme* getProteolyticEnzyme(const char* name);

  void showEnzymesCatalog();

  
  ProteolyticEnzyme* getProteolyticEnzyme(const string* name);
  
  void setUserEnzyme(bool enz) {
    user_enzyme_ = enz;
  }

 private:
  ProteolyticEnzyme* getEnzymeFromCatalog(const char* name, const char* fidelity) const;
  ProteolyticEnzyme* getEnzymeFromCatalog(int index, const char* fidelity) const;
  ProteolyticEnzyme* getProteolyticEnzyme(const char *name, const char* fidelity, const char* independent, const char *description, const char *specificities) const;
  string getEnzymeStandardName(const char *names) const;
  bool isEnzymeMatched(const char *findName, const char *names) const;

 private:
  int min_spacing_;

  bool user_enzyme_;
  std::string user_enz_name_;
  
  static const char *CATALOG_ENZYMES[][4];
};


#endif
