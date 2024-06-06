#include "ProspectorResult.h"
/*

Program       : ProspectorResult for discr_calc of PeptideProphet
Author        : Peter Baker <pbaker%at%cgl.ucsf.edu>
Date          : 11.07.17
SVN Info      : $Id$

Copyright (C) 2017 Peter Baker

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

Peter Baker
Mass Spectrometry Facility
University of California San Francisco
600 16th Street,
Genentech Hall, suite N472A
San Francisco, CA 94143-2240
USA
pbaker@cgl.ucsf.edu

*/

ProspectorResult::ProspectorResult ( char* szBuf, Boolean preexisting_probs )
{
  init ();
  char validator[] = "http://www.ncbi.nlm.nih.gov"; //EXPECT="; //"DATABASE=";
  format_ = 0; // set

  process ( szBuf, preexisting_probs, validator );
}
ProspectorResult::ProspectorResult ( Array<Tag*>* tags) :
  SearchResult ( tags )
{
  char score_tag[] = "search_score";

  if ( processed_ ) {
    const int num_nec_fields = 5;
    Boolean found[num_nec_fields];
    int k;
    for ( k = 0; k < num_nec_fields; k++)
      found[k] = False;

    Tag* tag;
    for(k = 0; k < tags->length(); k++) {
      tag = (*tags)[k];
      if(! strcmp(tag->getName(), score_tag) && tag->isStart()) {
	if(found[0] == False && ! strcmp(tag->getAttributeValue("name"), "expect")) {
	  expect_ = (double) atof(tag->getAttributeValue("value"));
	  found[0] = True;
	}
	else if(found[1] == False && ! strcmp(tag->getAttributeValue("name"), "pvalue")) {
	  pvalue_ = (double) atof(tag->getAttributeValue("value"));
	  found[1] = True;
	}
	else if(found[2] == False && ! strcmp(tag->getAttributeValue("name"), "disc_score")) {
	  disc_score_ = (double) atof(tag->getAttributeValue("value"));
	  found[2] = True;
	}
	else if(found[3] == False && ! strcmp(tag->getAttributeValue("name"), "ion_score")) {
	  ion_score_ = (double) atof(tag->getAttributeValue("value"));
	  found[3] = True;
	}
	else if(found[4] == False && ! strcmp(tag->getAttributeValue("name"), "ion_score_diff")) {
	  ion_score_diff_ = (double) atof(tag->getAttributeValue("value"));
	  found[4] = True;
	}
      }
    } // next tag
    for(k = 0; k < num_nec_fields; k++)
      if(! found[k])
	processed_ = False;

    if(expect_ == 0)
      processed_ = False;

  } // if score proc
}
void ProspectorResult::process ( char* szBuf, Boolean preexisting_probs, char* valid )
{
}
const char* ProspectorResult::getName ()
{
  return "PROTEINPROSPECTOR";
}
