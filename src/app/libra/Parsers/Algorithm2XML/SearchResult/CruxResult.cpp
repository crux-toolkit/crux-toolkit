#include "CruxResult.h"

/*

Program       : CruxResult for discr_calc of PeptideProphet 
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


/*
format -1 uncertain
format 0 no mw column
format 1 has mw column
*/

CruxResult::CruxResult() { }

CruxResult::CruxResult(Array<Tag*>* tags) : SearchResult(tags) {
  char score_tag[] = "search_score";
  /*
  char result_tag[] = "search_result";
  char hit_tag[] = "search_hit";
  char pepproph_tag[] = "peptideprophet_result";
  double neutral_prec_mass;
  */
  if(processed_) {
    const int num_nec_fields = 2;
    Boolean found[num_nec_fields];
    int k;
    for(k = 0; k < num_nec_fields; k++)
      found[k] = False;

    Tag* tag;
    for(k = 0; k < tags->length(); k++) {
      tag = (*tags)[k];

      if(! strcmp(tag->getName(), score_tag) && tag->isStart()) {
	if(! found[0] && ! strcmp(tag->getAttributeValue("name"), "xcorr_score")) {
	  xcorr_ = atof(tag->getAttributeValue("value"));
	  found[0] = True;
	}
	else if(! found[1] && ! strcmp(tag->getAttributeValue("name"), "delta_cn")) {
	  delta_ = atof(tag->getAttributeValue("value"));
	  found[1] = True;
	}
// 	else if(! found[2] && ! strcmp(tag->getAttributeValue("name"), "deltacnstar")) {
// 	  deltastar_ = atof(tag->getAttributeValue("value"));
// 	  found[2] = True;
// 	}
// 	else if(! found[3] && ! strcmp(tag->getAttributeValue("name"), "sprank")) {
// 	  rank_ = atoi(tag->getAttributeValue("value"));
// 	  found[3] = True;
// 	}
// 	else if(! found[4] && ! strcmp(tag->getAttributeValue("name"), "spscore")) {
// 	  sp_score_ = atof(tag->getAttributeValue("value"));
// 	  found[4] = True;
// 	}

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

    if(xcorr_ == 0.0)
      processed_ = False;



  } // if main score processed

}

CruxResult::CruxResult(char* szBuf, Boolean preexisting_probs, int format) {
  init(); // initial settings

  int formats[] = {1,0};  // 1 for bHasMWColumn True, 0 for False
  char validator[] = "http://www.ncbi.nlm.nih.gov"; //EXPECT="; //"DATABASE=";
  //char db_tag[] = "&amp;Db=";

  process(szBuf, preexisting_probs, format, formats, sizeof(formats)/sizeof(int), validator);
}



void CruxResult::process(char* szBuf, Boolean preexisting_probs, int format, int* formats, int numformats, char* valid) {
  if(strstr(szBuf, valid) == NULL) {
    return;
  }
  cerr << "error parsing " << szBuf << " Crux Result" << endl;
  exit(1);

}

const char* CruxResult::getName() { return "Crux"; }
