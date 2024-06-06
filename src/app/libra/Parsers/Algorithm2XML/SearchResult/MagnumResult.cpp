#include "MagnumResult.h"

/*

Program       : MagnumResult for discr_calc of PeptideProphet 
Author        : David Shteynberg <dshteynb%at%systemsbiology.org>                                                       
Date          : 03.28.16


Copyright (C) 2016 David Shteynberg and Andrew Keller

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
MagnumResult::~MagnumResult() {
  linkedPeps_->clear();
  delete linkedPeps_;
  
  for (int i = 0 ; i < linkedProts_->size(); i ++ ) {
    (*linkedProts_)[i]->clear();
    delete (*linkedProts_)[i];
  }
  linkedProts_->clear();
  delete linkedProts_;
  

}

MagnumResult::MagnumResult(char* szBuf, Boolean preexisting_probs){
  init();
  char validator[] = "http://www.ncbi.nlm.nih.gov"; //EXPECT="; //"DATABASE=";
  format_ = 0; // set

  process(szBuf, preexisting_probs, validator);
}

MagnumResult::MagnumResult(Array<Tag*>* tags) : SearchResult(tags) {
  char score_tag[] = "search_score";
  char xlink_score_tag[] = "xlink_score";
  bool unset = true;
  lower_score_ = 0;
  higher_score_ = 0;
  sum_rank_ = 0;
  linkedPeps_ = new vector<string>();
  linkedProts_ = new vector<vector<string>*>();
  if(processed_) {
    //cout << "processed" << endl;
    const int num_nec_fields = 3;
    Boolean found[num_nec_fields];
    int k;
    bool linked = false;
    for(k = 0; k < num_nec_fields; k++)
      found[k] = False;
    
    Tag* tag;
    for(k = 0; k < tags->length(); k++) {
      tag = (*tags)[k];
    //tag->write(cout);
      if(! strcmp(tag->getName(), score_tag) && tag->isStart()) {
	if(!found[0] &&
	   (! strcmp(tag->getAttributeValue("name"), "magnum_score") ||
	    ! strcmp(tag->getAttributeValue("name"), "Mscore")) ) {
	  magnumScore_ = (double)(atof(tag->getAttributeValue("value")));
	  found[0] = True;
	}
	else if(!found[1] &&
		(! strcmp(tag->getAttributeValue("name"), "delta_score") ||
		 ! strcmp(tag->getAttributeValue("name"), "Dscore")) )  {
	  deltaScore_ = (double)(atof(tag->getAttributeValue("value")));
	  found[1] = True;
	}
	else if(!found[2] &&
		(! strcmp(tag->getAttributeValue("name"), "e_value") ||
		 ! strcmp(tag->getAttributeValue("name"), "Evalue")) ) {
	  expectScore_ = (double)(atof(tag->getAttributeValue("value")));
	  found[2] = True;
	}
	
      }

      if(! strcmp(tag->getName(), xlink_score_tag) && tag->isStart() && ! strcmp(tag->getAttributeValue("name"), "rank")) {
	if (!sum_rank_) sum_rank_++;
	sum_rank_ += (double)(atof(tag->getAttributeValue("value")));
      }      

      if(! strcmp(tag->getName(), xlink_score_tag) && tag->isStart() && ! strcmp(tag->getAttributeValue("name"), "score")) {
	if (unset || lower_score_ > (double)(atof(tag->getAttributeValue("value")))) {
	  lower_score_ = (double)(atof(tag->getAttributeValue("value")));
	  unset = false;
	}
	if (higher_score_ < (double)(atof(tag->getAttributeValue("value")))) {
	  higher_score_ = (double)(atof(tag->getAttributeValue("value")));
	  unset = false;
	}
      }      

      if(! strcmp(tag->getName(), "linked_peptide") && tag->isStart()) {
	linked = true;
	linkedPeps_->push_back(string(tag->getAttributeValue("peptide")));
	linkedProts_->push_back(new vector<string>());
	(*linkedProts_)[linkedProts_->size()-1]->push_back(string(tag->getAttributeValue("protein")));
      }
      if(! strcmp(tag->getName(), "linked_peptide") && tag->isEnd()) {
	linked = false;
      }

      if(linked && ! strcmp(tag->getName(), "alternative_protein") && tag->isStart()) {
	(*linkedProts_)[linkedProts_->size()-1]->push_back(string(tag->getAttributeValue("protein")));
      }

    } // next tag
    for(k = 0; k < num_nec_fields; k++)
      if(! found[k])
        processed_ = False;


    //if(pvalue_ < 0)
    //  processed_ = False;

  } // if score proc




}

/*
char* MagnumResult::extractDatabase(char* html) {
  char start[] = "&amp;Db=";
  char stop[] = "&amp";
  return extractDatabaseWithTags(html, start, stop);
}
*/

void MagnumResult::process(char* szBuf, Boolean preexisting_probs, char* valid) {
}

const char* MagnumResult::getName() {
  return "Magnum";
}
