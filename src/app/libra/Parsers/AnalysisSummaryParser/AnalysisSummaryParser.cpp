/*

Program       : AnalysisSummary                                                  
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

#include "AnalysisSummaryParser.h"
#include "Util/RACI/RACI.h" // for reading gzipped files with efficient seeks

AnalysisSummaryParser::AnalysisSummaryParser(char* xmlfile, char* analysis, char* timestamp, Boolean timestamp_only) : Parser(NULL) {

  analysis_ = new char[strlen(analysis)+1];
  strcpy(analysis_, analysis);
  timestamp_ = new char[strlen(timestamp)+1];
  strcpy(timestamp_, timestamp);
  summary_ = NULL;
  found_ = False;
  timestamp_only_ = timestamp_only;

  init(xmlfile);
}

AnalysisSummaryParser::~AnalysisSummaryParser() {
  if(analysis_ != NULL)
    delete analysis_;
  if(timestamp_ != NULL)
    delete timestamp_;
  if(summary_ != NULL)
    delete summary_;
}


void AnalysisSummaryParser::parse(const char* xmlfile) {

  Tag* tag = NULL;

  //  int line_width = 10000;
  char *nextline=new char[line_width_];
  char* data = NULL;


  char summary_suff[] = "_summary";
  char timestamp_suff[] = "_timestamp";

  char* summary_name = new char[strlen(analysis_) + strlen(summary_suff) + 1];
  strcpy(summary_name, analysis_);
  strcat(summary_name, summary_suff);

  char* timestamp_name = new char[strlen(analysis_) + strlen(timestamp_suff) + 1];
  strcpy(timestamp_name, analysis_);
  strcat(timestamp_name, timestamp_suff);

  //cout << "summary name: " << summary_name << " timestamp name: " << timestamp_name << endl;

  const char* analysis_name = getAnalysisName(analysis_);
  if(analysis_name == NULL) {
    cout << " error: no analysis name" << endl;
    exit(1);
  }

  RACI fin(xmlfile); // can read gzipped xml
  if(! fin) {
    cerr << "AnalysisSummaryParser: error opening " << xmlfile << endl;
    exit(1);
  }
  Boolean found = False;
  while(fin.getline(nextline, line_width_)) {
    //cout << "next: " << nextline << endl;

    //   if(strstr(nextline, "msms_run_summary") != NULL || strstr(nextline, "refresh_timestamp") != NULL) {
    
      data = strstr(nextline, "<");
      while(data != NULL) {
	tag = new Tag(data);
	//setFilter(tag);
	//tag->write(cout);

	if(tag != NULL) {
	  if(! timestamp_only_ && tag->isStart() && ! strcmp(tag->getName(), "analysis_summary") &&
	     ! strcmp(tag->getAttributeValue("analysis"), analysis_name) &&
	     ! strcmp(tag->getAttributeValue("time"), timestamp_)) {
	    found = True;
	  }

	  else if(found && ! timestamp_only_ && tag->isStart() && ! strcmp(tag->getName(), summary_name)) {
	    if(summary_ != NULL) {
	      delete summary_;
	      summary_ = NULL;
	    }
	    summary_ = tag;
	    found_ = True;
	    return;
	    //strcpy(current_database, tag->getAttributeValue("database"));
	  }
	  else if(timestamp_only_ && tag->isStart() && ! strcmp(tag->getName(), "analysis_timestamp") &&
	     ! strcmp(tag->getAttributeValue("analysis"), analysis_name) &&
	     ! strcmp(tag->getAttributeValue("time"), timestamp_)) {
	    found = True;
	  }
	  else if(found && timestamp_only_ && tag->isStart() && ! strcmp(tag->getName(), timestamp_name)) {
	    if(summary_ != NULL) {
	      delete summary_;
	      summary_ = NULL;
	    }
	    found_ = True;
	    summary_ = tag;
	    return;
	  }
	  else if(! timestamp_only_ && tag->isStart() && ! strcmp(tag->getName(), "msms_run_summary")) { // process
	    fin.close();
	    delete tag;
	    return; // no chance now
	  }
	  else
	    delete tag;
	} // if not null

	data = strstr(data+1, "<");
      } // next tag
    
      //   } // if have reson to parse tags
  } // next line
  fin.close();

  delete[] nextline;

}


void AnalysisSummaryParser::setFilter(Tag* tag) {
  if(tag == NULL)
    return;

  if(filter_memory_) {
    filter_memory_ = False;
    filter_ = False;
  }

  if(! strcmp(tag->getName(), "search_result")) {
    if(tag->isStart()) {
      //tag->print();
      filter_ = True;
    }else{
      filter_memory_ = True;
    }
  }

}

const char* AnalysisSummaryParser::getAnalysisName(const char* analysis) {
  if(! strcmp(analysis, "xpressratio"))
    return "xpress";
  // default
  return analysis;

}

void AnalysisSummaryParser::print() {
  if(summary_ != NULL)
    summary_->print();
}
