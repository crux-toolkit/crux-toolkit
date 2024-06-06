#ifndef XPRESS_CGI_PROT_DISPLAY
#define XPRESS_CGI_PROT_DISPLAY

/*

Program       : XPressProteinDisplay                                                  
Author        : J.Eng and Andrew Keller <akeller@systemsbiology.org>, Robert Hubley, and 
                open source code                                                       
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

#include "Common/constants.h"
#include "Quantitation/XPress/XPressProteinRatioParser/XPressProteinRatioParser.h"

class XPressCGIProteinDisplay {

 public:


  XPressCGIProteinDisplay(const char* inputfiles, char* peptides, const char* protein, 
			  const char* cgihome, const char* protxmlfile, 
			  double minpepprob, Boolean heavy2light, const char* xslt, const char* mark_aas, Boolean glyc);




 protected:


  Array<const char*>* parse(const char* input, char separator);
  void writeXSLFile(const char* xslfile, const char* xmlfile, Boolean heavy2light);
  

  Array<const char*>* inputfiles_;
  Array<const char*>* inputlinks_;
  peplist* peptides_;
  RatioStruct pRatio_;
  XPressProteinRatioParser* parser_;
  char* cgihome_;
  char* xslt_;
  double minpepprob_;
  char* mark_aas_;
  Boolean glyc_;

};





#endif
