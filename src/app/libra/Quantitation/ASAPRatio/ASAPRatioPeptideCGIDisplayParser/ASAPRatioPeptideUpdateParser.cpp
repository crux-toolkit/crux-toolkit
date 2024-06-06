/*
Program       : ASAPRatioPeptideParser
Author        : X.Li and Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: ASAPRatioPeptideUpdateParser.cpp 8022 2020-02-12 21:47:21Z mhoopmann $



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

#include "ASAPRatioPeptideUpdateParser.h"


ASAPRatioPeptideUpdateParser::ASAPRatioPeptideUpdateParser(const char* xmlfile, const char* basename, int index, Array<Tag*>* replacements) : Parser(NULL) {
  // default settings
  index_ = index;
  replacements_ = replacements;
  overwrite_ = False;
  found_ = False;
  basename_ = new char[strlen(basename)+1];
  strcpy(basename_, basename);

  if(replacements_ == NULL) {
    cout << "error: null replacments for index " << index << " in " << xmlfile << endl;
    exit(1);
  }

  init(xmlfile);
}

ASAPRatioPeptideUpdateParser::~ASAPRatioPeptideUpdateParser() {
  if(replacements_ != NULL) {
    for(int k = 0; k < replacements_->length(); k++)
      if((*replacements_)[k] != NULL)
	delete (*replacements_)[k];
    delete replacements_;
  }
  if(basename_ != NULL)
    delete basename_;
}

void ASAPRatioPeptideUpdateParser::parse(const char* xmlfile) {
  char* engine = NULL;
  char* enzyme = NULL;
  char* massspec = NULL;
  Array<Tag*>* tags = NULL;
  Tag* tag = NULL;

  char *nextline = new char[line_width_];
  char* data = NULL;

  // construct a tempfile name, possibly in the tmp dir if so configured
  std::string outfile = make_tmpfile_name(xmlfile);

  //cerr << "writing data to " << outfile << endl;
  ofstream fout(outfile.c_str());
  if(! fout) {
    cerr << "cannot write output to file " << outfile << endl;
    exit(1);
  }

  Boolean done = False;
  Boolean replace = False;

  char match[1000];
  sprintf(match, "base_name=\"%s\"", basename_); // not really needed anymore
  Boolean analyze = False;

  char index_match[1000];
  sprintf(index_match, "index=\"%d\"", index_);

  RACI fin(xmlfile); // can read gzipped xml
  if(! fin) {
    cout << "ASAPRatioPeptideUpdateParser: error opening " << xmlfile << endl;
    exit(1);
  }

  Boolean ready = False;

  while(fin.getline(nextline, line_width_)) {
    if(! done && (analyze || (strstr(nextline, "spectrum_query") != NULL && strstr(nextline, index_match) != NULL))) {
      //if(! done && (analyze  || strstr(nextline, match) != NULL)) {
      analyze = True;
    
      data = strstr(nextline, "<");
      while(data != NULL) {
	tag = new Tag(data);

	if(tag != NULL) {
	  if(tag->isStart() && ! strcmp(tag->getName(), "spectrum_query") && atoi(tag->getAttributeValue("index")) != index_)
	    analyze = False;
	  
	  else if(tag->isStart() && ! strcmp(tag->getName(), "search_hit") && ! strcmp(tag->getAttributeValue("hit_rank"), "1")) {
	    ready = True;
	  }
	  else if(tag->isEnd() && ! strcmp(tag->getName(), "asapratio_result")) {
	    ready = False;
	    if(found_) 
	      done = True; // done
	  }

	  if(ready && tag->isStart() && ! strcmp(tag->getName(), "asapratio_result")) { // && 
	    //atoi(tag->getAttributeValue("index")) == index_) {
	    found_ = True;
	    // time to write out new stuff
	    for(int k = 0; k < replacements_->length(); k++)
	      if((*replacements_)[k] != NULL)
		if(! (*replacements_)[k]->isEnd() || strcmp((*replacements_)[k]->getName(), "asapratio_result"))
		  (*replacements_)[k]->write(fout);
	  }

	  if(done || ! found_) {
	    //if(! tag->isEnd() || strcmp(tag->getName(), "asapratio_result"))
	    tag->write(fout);
	  } 
	  delete tag;

	}
	data = strstr(data+1, "<");
      }
    }
    else
      fout << nextline << endl;

  }
  fin.close();
  fout.close();

  if(! overwrite(xmlfile, outfile.c_str(), "</msms_pipeline_analysis>")) {
    cerr << "error: no ASAPRatio data written to file " << xmlfile << endl;
    overwrite_ = False;
  }
  else
    overwrite_ = True;
  delete[] nextline;
}


Boolean ASAPRatioPeptideUpdateParser::update() {
  return found_ && overwrite_;
}

void ASAPRatioPeptideUpdateParser::setFilter(Tag* tag) {
  if(tag == NULL)
    return;

  if(filter_memory_) {
    filter_memory_ = False;
    filter_ = False;
  }

  if(! strcmp(tag->getName(), "spectrum_query") && tag->isStart() && atoi(tag->getAttributeValue("index")) == index_){ 
    if(tag->isStart()) {
      //tag->print();
      filter_ = True;
    }
    else {
      filter_memory_ = True;
    }
  }
}
