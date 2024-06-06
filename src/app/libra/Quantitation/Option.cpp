/*

Program       : Option                                                       
Author        : Andrew Keller <akeller@systemsbiology.org>                                                       
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

#include "Option.h"

/*
void setSearchEngine(ModelOptions* opts, char* engine) {
  //  if(opts->pI_)
  //    cout << "engine pI" << endl;

  if(! strcmp(engine, "SEQUEST"))
    opts->engine_ = SEQUEST;
  else if(! strcmp(engine, "MASCOT")) {
    opts->engine_ = MASCOT;
  }
  else if(! strcmp(engine, "COMET"))
    opts->engine_ = COMET;
  // still here
  else {
    cout << "error: do not recognize search engine: " << engine << endl;
    exit(1);
  }
  //  cout << "returning..." << endl;
}

char* getSearchEngine(ModelOptions opts) {
  if(opts.engine_ == SEQUEST)
    return "SEQUEST";
  if(opts.engine_ == MASCOT)
    return "MASCOT";
  if(opts.engine_ == COMET)
    return "COMET";
  if(opts.engine_ == PROBID)
    return "PROBID";
  return "";
}

void setSampleEnzyme(ModelOptions* opts, char* enzyme) {
  
  if(! strcmp(enzyme, "tryptic"))
    opts->enzyme_ = tryptic;
  else if(! strcmp(enzyme, "chymotryptic"))
    opts->enzyme_ = chymotryptic;
  else if(! strcmp(enzyme, "gluc_bicarb"))
    opts->enzyme_ = gluc_bicarb;
  else if(! strcmp(enzyme, "elastase"))
    opts->enzyme_ = elastase;
  else if(! strcmp(enzyme, "nonspecific"))
    opts->enzyme_ = nonspecific;
  else if(! strcmp(enzyme, "tca"))
    opts->enzyme_ = tca;
  else if(! strcmp(enzyme, "CNBr"))
    opts->enzyme_ = CNBr;
  else if(! strcmp(enzyme, "AspN"))
    opts->enzyme_ = AspN;
  else if(! strcmp(enzyme, "tryptic/CNBr"))
    opts->enzyme_ = tryptic_CNBr;
  else
    opts->enzyme_ = nonspecific; // default
}

char* getSampleEnzyme(ModelOptions opts) {
  if(opts.enzyme_ == tryptic)
    return "tryptic";
  if(opts.enzyme_ == chymotryptic)
    return "chymotryptic";
  if(opts.enzyme_ == gluc_bicarb)
    return "gluc_bicarb";
  if(opts.enzyme_ == elastase)
    return "elastase";
  if(opts.enzyme_ == nonspecific)
    return "nonspecific";
  if(opts.enzyme_ == tca)
    return "tca";
  if(opts.enzyme_ == CNBr)
    return "CNBr";
  if(opts.enzyme_ == AspN)
    return "AspN";
  if(opts.enzyme_ == tryptic_CNBr)
    return "tryptic/CNBr";
  return "";
}
*/
