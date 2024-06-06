#include "ProbIDResult.h"

/*

Program       : MascotResult for discr_calc of PeptideProphet 
Author        : Andrew Keller <akeller@systemsbiology.org>                                                       
Date          : 11.27.02 


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


ProbIDResult::ProbIDResult(Array<Tag*>* tags) : SearchResult(tags) {
  char score_tag[] = "search_score";
  probability_ = -1.0; // 

  if(processed_) {
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
	if(found[0] == False && ! strcmp(tag->getAttributeValue("name"), "bays_score")) {
	  if (!strcmp(tag->getAttributeValue("value"), "inf")) {
	    bays_ = 9999;
	  }
	  else {
	    bays_ = atof(tag->getAttributeValue("value"));
	  }
	    found[0] = True;
	}
	else if(found[1] == False && ! strcmp(tag->getAttributeValue("name"), "z_score")) {
	  zscore_ = atof(tag->getAttributeValue("value"));
	  found[1] = True;
	}
	else if(! strcmp(tag->getAttributeValue("name"), "probability")) {
	  probability_ = atof(tag->getAttributeValue("value"));
	}

      } // if score
    } // next tag
    for(k = 0; k < num_nec_fields; k++)
      if(! found[k]) {
	cout << "couldnt find " << (k+1) << endl;
	processed_ = False;
      }
    /*
    if(processed_)
       cout << "processed" << endl;
      else
       cout << "NOT Processed" << endl;
    */

  } // if score proc
/*
    else if(! strcmp(tag->getName(), pepproph_tag) && tag->isStart()) {
      probability_ = atof(tag->getAttributeValue("probability"));
    }
*/


}


const char* ProbIDResult::getName() {
  return "ProbID";
}
