#include "SpectraSTResult.h"

/*

Program       : SpectraSTResult for discr_calc of PeptideProphet 
Author        : Henry Lam <hlam@systemsbiology.org>                                                       
Date          : 04.10.06 


Copyright (C) 2006 Henry Lam

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

Henry Lam
Institute for Systems Biology
401 Terry Avenue North 
Seattle, WA  98109  USA
hlam@systemsbiology.org

*/


SpectraSTResult::SpectraSTResult(Array<Tag*>* tags) : SearchResult(tags) {
  if (!processed_) {
    return;
  }
  
  dot_ = -1.0;
  delta_ = -1.0;
  dotbias_ = -1.0;
  mzdiff_ = -1.0;
  charge_ = -1;
  hitnum_ = -1;
  fval_ = -1.0;
  libprob_ = 1.0;
  
  Tag* tag;
  for (int k = 0; k < tags->length(); k++) {
    tag = (*tags)[k];
    //      	tag->write(cout);
    if (!strcmp(tag->getName(), "search_score") && tag->isStart()) {
      if (!strcmp(tag->getAttributeValue("name"), "dot")) {
	dot_ = atof(tag->getAttributeValue("value"));
      } else if (!strcmp(tag->getAttributeValue("name"), "delta")) {
	delta_ = atof(tag->getAttributeValue("value"));
      } else if (!strcmp(tag->getAttributeValue("name"), "dot_bias")) {
	dotbias_ = atof(tag->getAttributeValue("value"));
      } else if (!strcmp(tag->getAttributeValue("name"), "precursor_mz_diff")) {
	mzdiff_ = atof(tag->getAttributeValue("value"));
      } else if (!strcmp(tag->getAttributeValue("name"), "charge")) {
	charge_ = atoi(tag->getAttributeValue("value"));
      } else if (!strcmp(tag->getAttributeValue("name"), "hits_num")) {
	hitnum_ = atoi(tag->getAttributeValue("value"));			
      } else if (!strcmp(tag->getAttributeValue("name"), "fval")) {
	fval_ = atof(tag->getAttributeValue("value"));
      } else if (!strcmp(tag->getAttributeValue("name"), "lib_probability")) {
        libprob_ = atof(tag->getAttributeValue("value"));
      } else if (!strcmp(tag->getAttributeValue("name"), "lib_remark")) {
	libremark_ = tag->getAttributeValue("value");
      }
    } // if score
  } // next tag
  
  if (dot_ < -0.001 || delta_ < -0.001 || dotbias_ < -0.001 || charge_ < 0 || hitnum_ < 0 || fval_ < -0.001) {
    processed_ = False;
  }
  
}


const char* SpectraSTResult::getName() {
  return "SpectraST";
}
