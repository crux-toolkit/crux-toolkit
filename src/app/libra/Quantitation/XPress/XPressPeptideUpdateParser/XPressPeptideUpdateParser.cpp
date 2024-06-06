/*
Program       : XPressPeptideUpdateParser                                                   
Author        : J.Eng and Andrew Keller <akeller@systemsbiology.org>, Robert Hubley, and 
                open source code                                                       
Date          : 11.27.02 
Version       : $Id: XPressPeptideUpdateParser.cpp 8870 2023-02-25 06:40:57Z real_procopio $

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

#include "XPressPeptideUpdateParser.h"


XPressPeptideUpdateParser::XPressPeptideUpdateParser(const char* xmlfile, int index, Tag* replacement) : Parser(NULL) {
  // default settings

  index_ = index;
  replacement_ = replacement;
  overwrite_ = False;
  found_ = False;
  if(replacement_ == NULL) {
    cout << "error: null replacment for index " << index << " in " << xmlfile << std::endl;
    exit(1);
  }

  init(xmlfile);
}

XPressPeptideUpdateParser::~XPressPeptideUpdateParser() {
  if(replacement_ != NULL)
    delete replacement_;

}

void XPressPeptideUpdateParser::parse(const char* xmlfile) {
  char* engine = NULL;
  char* enzyme = NULL;
  char* massspec = NULL;
  Array<Tag*>* tags = NULL;
  Tag* tag = NULL;
  char* data = NULL;

  // construct a tmpfile name based on xmlfile
  std::string outfile = make_tmpfile_name(xmlfile);
  ofstream fout(outfile.c_str());
  if(! fout) {
    cerr << "cannot write output to file " << outfile << endl;
    exit(1);
  }

  Boolean done = False;
  Boolean replace = False;

  char match[1000];
  sprintf(match, "index=\"%d\"", index_);
  Boolean analyze = False;

  
  RACI fin(xmlfile); // can read gzipped xml
  if(! fin) {
    cout << "XPressPeptideUpdateParser: error opening " << xmlfile << endl;
    exit(1);
  }
  char *nextline = new char[line_width_];
  while(fin.getline(nextline, line_width_)) {

    if(done) {
      fout << nextline << endl;
    }
    else if(! analyze && strstr(nextline, match) == NULL) {
      fout << nextline << endl;
    }

    else { // not yet done
      analyze = True;

      data = strstr(nextline, "<");
      while(data != NULL) {
	tag = new Tag(data);

	if(tag != NULL) {
	  setFilter(tag);

	  if(filter_) {
	    if(replace && tag->isStart() && strcmp(tag->getName(), replacement_->getName())==0 ) {
	      replacement_->write(fout);
	      done = True;
	      found_ = True;
	      analyze = False;
	    }
	    else {
	      if(tag->isStart() && ! strcmp(tag->getName(), "search_hit")) {
		if(! strcmp(tag->getAttributeValue("hit_rank"), "1"))
		  replace = True;
		else
		  replace = False;
	      }
	      tag->write(fout);
	    }
	  } // if filter
	  else {
	    tag->write(fout);
	  }
	  delete tag;
	} // if not null

	data = strstr(data+1, "<");
      } // next tag
    } // not done

  } // next line
  fin.close();
  fout.close();

  delete [] nextline;

  if(! overwrite(xmlfile, outfile.c_str(), "</msms_pipeline_analysis>")) {
    cerr << "error: no xpress data written to file " << xmlfile << endl;
    overwrite_ = False;
  }
  else
    overwrite_ = True;
}


Boolean XPressPeptideUpdateParser::update() {
  return found_ && overwrite_;
}

void XPressPeptideUpdateParser::setFilter(Tag* tag) {
  if(tag == NULL)
    return;

  if(filter_memory_) {
    filter_memory_ = False;
    filter_ = False;
  }

  if(! strcmp(tag->getName(), "spectrum_query") && tag->isStart() && atoi(tag->getAttributeValue("index")) == index_) {
    if(tag->isStart())
      filter_ = True;
    else
      filter_memory_ = True;

  }
}

