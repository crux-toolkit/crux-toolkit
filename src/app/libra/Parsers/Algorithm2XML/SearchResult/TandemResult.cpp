#include "TandemResult.h"

/*

Program       : TandemResult for discr_calc of PeptideProphet 
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

TandemResult::TandemResult(char* szBuf, Boolean preexisting_probs){
  init();
  char validator[] = "http://www.ncbi.nlm.nih.gov"; //EXPECT="; //"DATABASE=";
  format_ = 0; // set

  process(szBuf, preexisting_probs, validator);
}

TandemResult::TandemResult(Array<Tag*>* tags) : SearchResult(tags) {
  char score_tag[] = "search_score";
  /*
  char result_tag[] = "search_result";
  char hit_tag[] = "search_hit";
  char pepproph_tag[] = "peptideprophet_result";
  double neutral_prec_mass;
  */
  //for(int k = 0; k < tags->length(); k++) 
  // (*tags)[k]->write(cout);

  //cout << "here ready...;" << endl;


  if(processed_) {
    //cout << "processed" << endl;
    const int num_nec_fields = 3;
    Boolean found[num_nec_fields];
    int k;
    for(k = 0; k < num_nec_fields; k++)
      found[k] = False;

    Tag* tag;
    for(k = 0; k < tags->length(); k++) {
      tag = (*tags)[k];
    //tag->write(cout);
      if(found[0] == False && ! strcmp(tag->getName(), score_tag) && tag->isStart()) {
	if(! strcmp(tag->getAttributeValue("name"), "expect")) {
	  expect_ = (double) atof(tag->getAttributeValue("value"));
	  found[0] = True;
	}
	else if(found[1] == False && ! strcmp(tag->getAttributeValue("name"), "hyperscore")) {
	  hyper_ = (double) atof(tag->getAttributeValue("value"));
	  found[1] = True;
	}
	else if(found[2] == False && ! strcmp(tag->getAttributeValue("name"), "nextscore")) {
	  next_ = (double) atof(tag->getAttributeValue("value"));
	  found[2] = True;
	}
  }

/*
    else if(! strcmp(tag->getName(), pepproph_tag) && tag->isStart()) {
      probability_ = atof(tag->getAttributeValue("probability"));
    }
*/


    
    } // next tag
    for(k = 0; k < num_nec_fields; k++)
      if(! found[k])
        processed_ = False;


    if(expect_ == 0)
      processed_ = False;

  } // if score proc




}

/*
char* TandemResult::extractDatabase(char* html) {
  char start[] = "&amp;Db=";
  char stop[] = "&amp";
  return extractDatabaseWithTags(html, start, stop);
}
*/

void TandemResult::process(char* szBuf, Boolean preexisting_probs, char* valid) {
}

const char* TandemResult::getName() {
  return "X! Tandem";
}
