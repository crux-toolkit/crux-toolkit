/*
Program       : EnzymeDigestionParser
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02 
SVN Info      : $Id: EnzymeDigestionParser.cpp 8440 2021-04-19 22:37:03Z real_procopio $

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

#include "EnzymeDigestionParser.h"

EnzymeDigestionParser::EnzymeDigestionParser(char* xmlfile, char* enzyme) : Parser(NULL) {
  enzyme_ = new char[strlen(enzyme)+1];
  strcpy(enzyme_, enzyme);
  enzyme_spec_ = new ProteolyticEnzymeFactory();

  enzyme_digestion_ = enzyme_spec_->getProteolyticEnzyme(enzyme_);
  if(enzyme_digestion_ == NULL) {
    cerr << "Error: enzyme " << enzyme_ << " not recognized" << endl;
    exit(1);
  }
  cerr << "Enzyme specified: " << enzyme_ << endl;
  enzyme_tags_ = enzyme_digestion_->getPepXMLTags();

  init(xmlfile);
}

EnzymeDigestionParser::EnzymeDigestionParser(char* xmlfile) : Parser(NULL) {
  enzyme_ = NULL;
  enzyme_spec_ = new ProteolyticEnzymeFactory();

  enzyme_digestion_ = NULL;  
  enzyme_tags_ = NULL;

  init(xmlfile);
}

EnzymeDigestionParser::~EnzymeDigestionParser() {
  if(enzyme_ != NULL)
    delete enzyme_;
  if(enzyme_spec_ != NULL)
    delete enzyme_spec_;
  if(enzyme_digestion_ != NULL)
    delete enzyme_digestion_;
}


void EnzymeDigestionParser::parse(const char* xmlfile) {
  Tag* tag = NULL;
  char* data = NULL;

  RACI fin(xmlfile); // can read gzipped xml
  if(! fin) {
    cerr << "EnzymeDigestionParser: error opening " << xmlfile << endl;
    exit(1);
  }
  // construct a tmpfile name based on xmlfile
  std::string outfile = make_tmpfile_name(xmlfile);
  ofstream fout(outfile.c_str());
  if(! fout) {
    cerr << "cannot write output to file " << outfile << endl;
    exit(1);
  }
  char enzyme_match[] = "<sample_enzyme";
  char msms_match[]   = "<msms_run_summary";
  char hit_match[]    = "<search_hit";
  Boolean analyze = True;
  char text[100];
  Boolean enzyme_on = False;
  char *nextline = new char[line_width_];
  while(fin.getline(nextline, line_width_)) {

    if(strstr(nextline, enzyme_match) != NULL ||
       strstr(nextline, msms_match)   != NULL ||
       strstr(nextline, hit_match)    != NULL ||
       enzyme_on) {

      data = strstr(nextline, "<");
      while(data != NULL) {
	tag = new Tag(data);

	if(tag != NULL) {
	  if(0 && tag->isStart() && ! strcmp(tag->getName(), "msms_run_summary")) {
	    if(enzyme_ == NULL) {
	      if(enzyme_digestion_ != NULL)
		delete enzyme_digestion_;
	      /*
		enzyme_digestion_ = enzyme_spec_->getEnzymeDigestion(tag->getAttributeValue("sample_enzyme"));
		if(enzyme_digestion_ == NULL) {
		cerr << "Error: enzyme " << tag->getAttributeValue("sample_enzyme") << " not recognized" << endl;
		exit(1);
		}
	      */
	      analyze = 1;
	    } // if no enzyme specified on command line
	    else {
	      analyze = strcmp(enzyme_, tag->getAttributeValue("sample_enzyme"));
	      if(analyze)
		tag->setAttributeValue("sample_enzyme", enzyme_);

	      //strcpy(current_database, tag->getAttributeValue("database"));
	    }
	  }
	  else if(tag->isStart() && ! strcmp(tag->getName(), "sample_enzyme")) {
	    if(enzyme_ == NULL) {
	      if(enzyme_digestion_ != NULL)
		delete enzyme_digestion_;
	      enzyme_digestion_ = new ProteolyticEnzyme(tag);
	      if(enzyme_digestion_ == NULL) {
		cerr << "Error: enzyme ";
                tag->write(cerr);
                cerr << " not recognized" << endl;
		exit(1);
	      }
	      analyze = 1;
	    } // if no enzyme specified on command line
	    else { // only analyze if specified enz differs from one already present
	      analyze = strcmp(enzyme_, tag->getAttributeValue("name"));
	    }
	    enzyme_on = True;
	  } // if sample enzyme start
	  else if(enzyme_on) {
	    if(tag->isEnd() && ! strcmp(tag->getName(), "sample_enzyme")) {
	      if(enzyme_ == NULL && enzyme_digestion_ != NULL) {
		enzyme_digestion_->fixSpecificity();
		enzyme_tags_ = enzyme_digestion_->getPepXMLTags();
	      } // ino specified enzyme
	      if(enzyme_tags_ != NULL) {
		for(int k = 0; k < enzyme_tags_->length(); k++) {
		  if((*enzyme_tags_)[k] != NULL) {
		    (*enzyme_tags_)[k]->write(fout);
		    if(enzyme_ == NULL)
		      delete (*enzyme_tags_)[k];
		  }
		}
		if(enzyme_ == NULL)
		  delete enzyme_tags_;
	      }
	      enzyme_on = False;		
	    } // last one
	    else if(enzyme_ == NULL && enzyme_digestion_ != NULL) {
	      enzyme_digestion_->enterSpecificity(tag);
	      cerr << "here!!" << endl;
	    }

	  } // if enzyme on
	  if(tag->isStart() && ! strcmp(tag->getName(), "database_refresh_timestamp")) {
	    cout << "Error: EnzymeDigestionParser cannot be applied to data with database refreshment" << endl;
	    exit(1);
	  }
	  if(analyze && tag->isStart() && ! strcmp(tag->getName(), "search_hit") 
	     && ! strcmp(tag->getAttributeValue("hit_rank"), "1")) { // process
	    if(tag->getAttributeValue("peptide_prev_aa") != NULL && tag->getAttributeValue("peptide_next_aa") != NULL) {
	      sprintf(text, "%d", enzyme_digestion_->getNumTolTerm((tag->getAttributeValue("peptide_prev_aa"))[0], 
#ifdef USE_STD_MODS
								   tag->getAttributeValue("peptide"), 
#endif
#ifndef USE_STD_MODS
								   tag->getAttributeValue("stripped_peptide"),
#endif 
								   (tag->getAttributeValue("peptide_next_aa"))[0]));
	      tag->setAttributeValue("num_tol_term", text);
	    }
	    else
	      tag->setAttributeValue("num_tol_term", "2"); // default

	    sprintf(text, "%d", enzyme_digestion_->getNumMissedCleavages(tag->getAttributeValue("peptide"))); 
	    tag->setAttributeValue("num_missed_cleavages", text); // default
	  }
	  if(! enzyme_on && (tag->isStart() || strcmp(tag->getName(), "sample_enzyme")))
	    tag->write(fout);

	  delete tag;
	} // if not null

	data = strstr(data+1, "<");
      } // next tag
    } // if have reson to parse tags
    else
      fout << nextline << endl;
  } // next line
  fin.close();
  fout.close();

  delete [] nextline;

  if(! overwrite(xmlfile, outfile.c_str(), "</msms_pipeline_analysis>")) {
    cerr << "error: no enzyme data written to file " << xmlfile << endl;
  }

}

void EnzymeDigestionParser::setFilter(Tag* tag) {
  if(tag == NULL)
    return;

  if(filter_memory_) {
    filter_memory_ = False;
    filter_ = False;
  }

  if(! strcmp(tag->getName(), "spectrum_query")){ 
    if(tag->isStart())
      filter_ = True;
    else
      filter_memory_ = True;
  }

}
