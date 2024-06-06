#include "PhenyxResult.h"

/*

Program       : PhenyxResult for discr_calc of PeptideProphet 
Author        : David Shteynberg <dshteynb%at%systemsbiology.org>                                                       
Date          : 05.16.07


Copyright (C) 2007 David Shteynberg and Andrew Keller

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

David Shteynberg
Institute for Systems Biology
401 Terry Avenue North 
Seattle, WA  98109  USA
akeller@systemsbiology.org

*/

PhenyxResult::PhenyxResult(char* szBuf, Boolean preexisting_probs){
  init();
  char validator[] = "http://www.ncbi.nlm.nih.gov"; //EXPECT="; //"DATABASE=";
  format_ = 0; // set

  process(szBuf, preexisting_probs, validator);
}

PhenyxResult::PhenyxResult(Array<Tag*>* tags) : SearchResult(tags) {
  char score_tag[] = "search_score";

  if(processed_) {
    //cout << "processed" << endl;
    const int num_nec_fields = 2;
    Boolean found[num_nec_fields];
    int k;
    for(k = 0; k < num_nec_fields; k++)
      found[k] = False;

    Tag* tag;
    for(k = 0; k < tags->length(); k++) {
      tag = (*tags)[k];
    //tag->write(cout);
      if(! strcmp(tag->getName(), score_tag) && tag->isStart()) {
	if(found[0] == False && ! strcmp(tag->getAttributeValue("name"), "zscore")) {
	  zscore_ = (double) atof(tag->getAttributeValue("value"));
	  found[0] = True;
	}
	else if(found[1] == False && ! strcmp(tag->getAttributeValue("name"), "origScore")) {
	  origScore_ = (double) atof(tag->getAttributeValue("value"));
	  found[1] = True;
	}

  }
    
    } // next tag
    for(k = 0; k < num_nec_fields; k++)
      if(! found[k])
        processed_ = False;


    if(zscore_ == 0)
      processed_ = False;

  } // if score proc




}

/*
char* PhenyxResult::extractDatabase(char* html) {
  char start[] = "&amp;Db=";
  char stop[] = "&amp";
  return extractDatabaseWithTags(html, start, stop);
}
*/

void PhenyxResult::process(char* szBuf, Boolean preexisting_probs, char* valid) {
}

const char* PhenyxResult::getName() {
  return "Phenyx";
}
