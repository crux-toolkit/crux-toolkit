/*
Program       : ASAPRatioPeptideParser
Author        : X.Li and Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: ASAPRatioPeptideCGIDisplayParser.cpp 8022 2020-02-12 21:47:21Z mhoopmann $


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

#include "ASAPRatioPeptideCGIDisplayParser.h"

ASAPRatioPeptideCGIDisplayParser::ASAPRatioPeptideCGIDisplayParser(const char* xmlfile, const char* basename, const char* timestamp, int index, int elution, Boolean zeroBG, double mzBound) : Parser(NULL) {
  // default settings

  index_ = index;
  //cout << "index: " << index_ << " and asapindex: " << asap_index_ << endl;

  zeroBG_ = zeroBG;
  mzBound_ = mzBound;
  found_ = False;
  elution_ = elution;
  timestamp_ = new char[strlen(timestamp)+1];
  strcpy(timestamp_, timestamp);
  basename_ = new char[strlen(basename)+1];
  strcpy(basename_, basename);

  init(xmlfile);
}

ASAPRatioPeptideCGIDisplayParser::~ASAPRatioPeptideCGIDisplayParser() {
  if(basename_ != NULL)
    delete basename_;
}

void ASAPRatioPeptideCGIDisplayParser::parse(const char* xmlfile) {
  char* engine = NULL;
  char* enzyme = NULL;
  char* massspec = NULL;
  Array<Tag*>* tags = NULL;
  Tag* tag = NULL;

  char *nextline = new char[line_width_];
  char* data = NULL;

  char match[1000];
  sprintf(match, "base_name=\"%s\"", basename_);
  Boolean analyze = False;

  char ind_match[1000];
  sprintf(ind_match, "index=\"%d\"", index_);
  
  long scan = -1;
  int precursor_charge = -1;

  int next_elution = -2;

  tags = new Array<Tag*>;
  
  RACI fin(xmlfile); // can read gzipped xml
  if(! fin) {
    cout << "ASAPRatioPeptideCGIDisplayParser: error opening " << xmlfile << endl;
    exit(1);
  }
  Boolean ready = False;
  Boolean match_found = False;
  Boolean asap_ind_match_found = False;

  //cout << "match: " << ind_match << endl;
  //cout << "elution: " << elution_ << endl;
  //cout << "asapind match: " << asap_ind_match << endl;

  while(fin.getline(nextline, line_width_)) {
    if((elution_ == -2 && strstr(nextline, "asapratio_") != NULL) ||
       strstr(nextline, match) != NULL || strstr(nextline, "analysis_") != NULL ||
       (match_found && strstr(nextline, "spectrum_query") != NULL && strstr(nextline, ind_match) != NULL) ||
       //(match_found && (index_ == 1 || strstr(nextline, asap_ind_match) != NULL)) ||
       asap_ind_match_found) {

      if(strstr(nextline, match) != NULL)
	match_found = True;
      else if(match_found && (strstr(nextline, "spectrum_query") != NULL && strstr(nextline, ind_match) != NULL))
	//(index_ == 1 || strstr(nextline, asap_ind_match) != NULL))
	asap_ind_match_found = True;

      analyze = True;
    }
    else
      analyze = False;

    if(analyze) {
      data = strstr(nextline, "<");
      while(data != NULL) {
	tag = new Tag(data);
	//tag->write(cout);
	if(tag != NULL) {
	  if(elution_ == -2 && tag->isStart() && ! strcmp(tag->getName(), "asapratio_summary")) {
	    next_elution = atoi(tag->getAttributeValue("elution"));
	  }
	  else if(elution_ == -2 && tag->isStart() && ! strcmp(tag->getName(), "analysis_timestamp") &&
		  ! strcmp(tag->getAttributeValue("time"), timestamp_)) {
	    elution_ = next_elution; // found it
	  }
	  else if(0 && ! strcmp(tag->getAttributeValue("time"), timestamp_)) {
	    elution_ = next_elution; // found it
	    //cout << "setting elution to " << elution_ << endl;
	  }

	  if(tag->isStart() && ! strcmp(tag->getName(), "spectrum_query")) {
	    scan = (long)(atoi(tag->getAttributeValue("start_scan")));
	    precursor_charge = atoi(tag->getAttributeValue("assumed_charge"));
	  }
	  else if(tag->isStart() && ! strcmp(tag->getName(), "search_hit") && ! strcmp(tag->getAttributeValue("hit_rank"), "1")) {
	    ready = True;
	  }
	  else if(tag->isEnd() && ! strcmp(tag->getName(), "search_hit")) {
	    ready = False;
	    if(found_) {
	      setPepDataStruct(tags, elution_, scan, precursor_charge);
	      // clean up
	      for(int k = 0; k < tags->length(); k++)
		if((*tags)[k] != NULL)
		  delete (*tags)[k];

	      fin.close();
	      delete tags;
	      delete[] nextline;
	      return;
	    }
	  }

	  if(ready) {
	    if(tag->isStart() && ! strcmp(tag->getName(), "asapratio_result"))
	      found_ = True;
	    
	    if(found_) {
	      tags->insertAtEnd(tag);
	    }
	    else
	      if(tag != NULL)
		delete tag;
	  }
	  else
	    if(tag != NULL)
	      delete tag;

	}

	data = strstr(data+1, "<");
      }
    }
  }
  fin.close();
  delete tags;
  delete[] nextline;
}


Boolean ASAPRatioPeptideCGIDisplayParser::found() {
  return found_;
}

void ASAPRatioPeptideCGIDisplayParser::setFilter(Tag* tag) {
  if(tag == NULL)
    return;

  if(filter_memory_) {
    filter_memory_ = False;
    filter_ = False;
  }

  if(! strcmp(tag->getName(), "spectrum_query")) {
    if(tag->isStart()) {
      //tag->print();
      filter_ = True;
    }
    else {
      filter_memory_ = True;
    }
  }
}

pepDataStrct ASAPRatioPeptideCGIDisplayParser::getPepDataStruct() {
  return data_;
}

void ASAPRatioPeptideCGIDisplayParser::setPepDataStruct(Array<Tag*>* tags, int elution, long scan, int precursor_charge) {
  Tag* next;
  int charge=0;
  data_.scan = scan;
  data_.chrg = precursor_charge;
  data_.eltn = elution;

  for(int k = 0; k < tags->length(); k++) {
    next = (*tags)[k];
    if(next->isStart()) {
      if(! strcmp(next->getName(), "asapratio_result")) {
	data_.pepRatio[0] = atof(next->getAttributeValue("mean"));
	data_.pepRatio[1] = atof(next->getAttributeValue("error"));
	data_.pepH2LRatio[0] = atof(next->getAttributeValue("heavy2light_mean"));
	data_.pepH2LRatio[1] = atof(next->getAttributeValue("heavy2light_error"));
	// make change here
	if(data_.pepRatio[0] == -1) {
	  data_.pepRatio[0] = -2;
	  data_.pepH2LRatio[0] = -2;
	}
	else if(data_.pepRatio[0] >= 999.0) {
	  data_.pepRatio[0] = -1;
	  data_.pepH2LRatio[0] = -1;
	}
      }
      else if(! strcmp(next->getName(), "asapratio_peptide_data")) {
	data_.indx = atoi(next->getAttributeValue("status"));
	data_.cidIndx = atoi(next->getAttributeValue("cidIndex"));
	data_.msLight = atof(next->getAttributeValue("light_mass"));
	data_.msHeavy = atof(next->getAttributeValue("heavy_mass"));
	data_.areaFlag = atoi(next->getAttributeValue("area_flag"));
      }
      else if(! strcmp(next->getName(), "asapratio_contribution")) {
	charge = atoi(next->getAttributeValue("charge"));
	data_.pkRatio[charge-1] = atof(next->getAttributeValue("ratio"));
	data_.pkError[charge-1] = atof(next->getAttributeValue("error"));
	data_.pkCount[charge-1] = atoi(next->getAttributeValue("use"));
      }
      else if(! strcmp(next->getName(), "asapratio_lc_lightpeak")) {
	int label = 0;
	data_.peaks[charge-1][label].indx = atoi(next->getAttributeValue("status"));
	data_.peaks[charge-1][label].valley[0] = atoi(next->getAttributeValue("left_valley"));
	data_.peaks[charge-1][label].valley[1] = atoi(next->getAttributeValue("right_valley"));
	if (zeroBG_) {
	  data_.peaks[charge-1][label].bckgrnd = 0.;
	}
	else {
	  data_.peaks[charge-1][label].bckgrnd = atof(next->getAttributeValue("background"));
	}
	data_.peaks[charge-1][label].area[0] = atof(next->getAttributeValue("area"));
	data_.peaks[charge-1][label].area[1] = atof(next->getAttributeValue("area_error"));
	data_.peaks[charge-1][label].time[0] = atof(next->getAttributeValue("time"));
	data_.peaks[charge-1][label].time[1] = atof(next->getAttributeValue("time_width"));
	data_.peaks[charge-1][label].peak = atoi(next->getAttributeValue("is_heavy"));
      }
      else if(! strcmp(next->getName(), "asapratio_lc_heavypeak")) {
	int label = 1;
	data_.peaks[charge-1][label].indx = atoi(next->getAttributeValue("status"));
	data_.peaks[charge-1][label].valley[0] = atoi(next->getAttributeValue("left_valley"));
	data_.peaks[charge-1][label].valley[1] = atoi(next->getAttributeValue("right_valley"));
	if (zeroBG_) {
	  data_.peaks[charge-1][label].bckgrnd = 0.;
	}
	else {
	  data_.peaks[charge-1][label].bckgrnd = atof(next->getAttributeValue("background"));
	}
	data_.peaks[charge-1][label].area[0] = atof(next->getAttributeValue("area"));
	data_.peaks[charge-1][label].area[1] = atof(next->getAttributeValue("area_error"));
	data_.peaks[charge-1][label].time[0] = atof(next->getAttributeValue("time"));
	data_.peaks[charge-1][label].time[1] = atof(next->getAttributeValue("time_width"));
	data_.peaks[charge-1][label].peak = atoi(next->getAttributeValue("is_heavy"));
      }
    }
  }
}
