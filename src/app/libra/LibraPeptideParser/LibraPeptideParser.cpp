/*
Program       : LibraPeptideParser
Author        : Patrick Pedrioli and Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: LibraPeptideParser.cpp 8868 2023-02-24 18:46:05Z real_procopio $

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

#include "Common/sysdepend.h"
#include "LibraPeptideParser.h"
#include "Common/util.h"
#include "Parsers/Parser/TagListComparator.h" // regression test stuff - bpratt Insilicos LLC, Nov 2005

using namespace mzParser;

// will probably have to define LibraInputStruct instead with relevant Libra parameters
LibraPeptideParser::LibraPeptideParser(const char* xmlfile, const char* conditionFileName, const char *testMode) : Parser("libra") {
  condition_ = new LibraConditionHandler();
  condition_->setFileName( conditionFileName );
  if (condition_->readFile() != 0) {
    exit(1);
  }

  libra_summary_ = NULL;
  libra_result_ = NULL;
  libra_quantifier_ = NULL;
  testMode_ = testMode?strdup(testMode):NULL; // regression test activity?

  mzXMLfile_[0] = 0;

  init(xmlfile);
}


LibraPeptideParser::~LibraPeptideParser() {
  if(libra_quantifier_ != NULL)
    delete libra_quantifier_;

  if(libra_summary_ != NULL)
    delete libra_summary_;

  if(libra_result_ != NULL)
    delete libra_result_;

  if(condition_ != NULL)
    delete condition_;

  free(testMode_);
}


void LibraPeptideParser::parse(const char* xmlfile) {
  Array<Tag*>* tags = NULL;
  Tag* tag = NULL;

  char *nextline = new char[line_width_];
  char* data = NULL;

  //
  // regression test stuff - bpratt Insilicos LLC, Nov 2005
  //
  Array<Tag*> test_tags;
  eTagListFilePurpose testType;
  char *testFileName=NULL;
  checkRegressionTestArgs(testMode_,testType);
  if (testType!=NO_TEST) {
    testFileName = constructTagListFilename(xmlfile, // input file
					    testMode_, // program args
					    "libraPeptideParser", // program name
					    testType); // user info output
  }
#define RECORD(tag) {(tag)->write(fout);if (testType!=NO_TEST) {test_tags.insertAtEnd(new Tag(*(tag)));}}

  Tag* timestamp = new Tag("analysis_timestamp", True, True);
  timestamp->setAttributeValue("analysis", getName());
  timestamp->setAttributeValue("time", time_);
  timestamp->setAttributeValue("id", "1");

  Tag* summary_start = new Tag("analysis_summary", True, False);
  summary_start->setAttributeValue("analysis", getName());
  summary_start->setAttributeValue("time", time_);

  //  Tag* summary_term = new Tag("libra_summary", False, True);
  Tag* summary_stop = new Tag("analysis_summary", False, True);
  Array<Tag*>* summary_tags = condition_->getPepXMLTags();

  Tag* result_start = new Tag("analysis_result", True, False);
  result_start->setAttributeValue("analysis", getName());

  Tag* result_stop = new Tag("analysis_result", False, True);

  TagFilter* libra_filter = new TagFilter("analysis_result");
  libra_filter->enterRequiredAttributeVal("analysis", getName());
  TagFilter* libra_time_filter = new TagFilter("analysis_timestamp");
  libra_time_filter->enterRequiredAttributeVal("analysis", getName());
  TagFilter* libra_summ_filter = new TagFilter("analysis_summary");
  libra_summ_filter->enterRequiredAttributeVal("analysis", getName());

  // construct a tempfile name, possibly in the tmp dir if so configured
  std::string outfile = make_tmpfile_name(xmlfile);
  ofstream fout(outfile.c_str());

  if(! fout) {
    cerr << "cannot write output to file " << outfile << endl;
    exit(1);
  }

  Boolean first = False;
  Array<Tag*>* libra_tags = NULL;
  Boolean collected = False;

  RAMPFILE *pFI = NULL;
  
  RACI fin(xmlfile); // can read gzipped XML

  if(! fin) {
    cerr << "LibraPeptideParser:error opening " << xmlfile << endl;
    exit(1);
  }
  else {
    cout << "Calculating and updating peptide ratios in " << xmlfile << "..." << endl;
  }


  while(fin.getline(nextline, line_width_)) {
    //cout << "next: " << nextline << endl;

    data = strstr(nextline, "<");

    while(data != NULL) {
      tag = new Tag(data);

      //tag->write(cout);
      collected = False;

      setFilter(tag);

      if((! libra_filter->filter(tag) && ! libra_summ_filter->filter(tag) 
	  && ! libra_time_filter->filter(tag))) {

	if(tag->isStart() && ! strcmp(tag->getName(), "msms_pipeline_analysis")) {
	  RECORD(tag);
	  RECORD(summary_start);

	  if(summary_tags != NULL) {
	    for(int k = 0; k < summary_tags->length(); k++) {
	      if((*summary_tags)[k] != NULL) {
		RECORD((*summary_tags)[k]);
		delete (*summary_tags)[k];
	      }
	    }
	    delete summary_tags;
	  }

	  RECORD(summary_stop);
	  delete summary_start;
	  delete summary_stop;

	} else if(tag->isStart() && ! strcmp(tag->getName(), "msms_run_summary")) {
	  if(rampConstructInputFileName(mzXMLfile_, 
					sizeof(mzXMLfile_), tag->getAttributeValue("base_name"))==NULL){
	    strcpy(mzXMLfile_,tag->getAttributeValue("base_name"));
	    strcat(mzXMLfile_,".mzXML");
	  }
	  unCygwinify(mzXMLfile_); // no effect in cygwin builds

	  if (pFI) {
	    rampCloseFile(pFI);
	  }
	  if ( (pFI = rampOpenFile(mzXMLfile_)) == NULL) {
	    printf("could not open input file %s\n", mzXMLfile_);
	    exit(1);
	  }
	  
	  libra_quantifier_ = new LibraWrapper(condition_, mzXMLfile_, pFI);

	  if (libra_quantifier_==NULL) {
	    printf("Error - cannot quantify with %s\n", mzXMLfile_);
	    exit(1);
	  }

	  first = True;
	  RECORD(tag);

	} else if(tag->isEnd() && ! strcmp(tag->getName(), "search_summary")) {
	  // msms_summ
	  RECORD(tag);
	  RECORD(timestamp);

	} else if(filter_) {

	  if(tag->isStart() && ! strcmp("spectrum_query", tag->getName())) {

	    libra_result_ = libra_quantifier_->getLibraResult(atoi(tag->getAttributeValue("start_scan")));
	    libra_tags = libra_result_->getPepXMLTags();

	    // progress bar  (assumes index attribute increases monotonically)
	    const char* idx = tag->getAttributeValue("index");
	    if (atoi(idx) % 1000 == 0)
	      cout << " " << idx << endl;
	    else if (atoi(idx) % 100 == 0)
	      cout << ":";
	    else if (atoi(idx) % 10 == 0)
	      cout << ".";

	    flush(cout);

	  } else if(0 && tag->isStart() && ! strcmp("search_hit", tag->getName()) 
                    && ! strcmp("1", tag->getAttributeValue("hit_rank"))) {
     	  }

	  if(tags == NULL)
	    tags = new Array<Tag*>;

	  tags->insertAtEnd(tag);
	  collected = True;

	} else {

	  if(tag->isEnd() && ! strcmp(tag->getName(), "msms_run_summary")) {
	    rampCloseFile(pFI);
	    pFI = NULL;
	    if(libra_quantifier_ != NULL) {
	      delete libra_quantifier_;
	      libra_quantifier_ = NULL;
	    }
	  }

	  if(tag != NULL) {
	    RECORD(tag);
	  }

	}

	if(filter_memory_) { // process
	  if(tags != NULL) {
	    for(int k = 0; k < tags->length(); k++) {
	      if((*tags)[k] != NULL) {
		if(! libra_filter->filter((*tags)[k])) {
		  // here check for correct time to write xpress tag
		  if((*tags)[k]->isEnd() && 
		     ! strcmp((*tags)[k]->getName(), "search_hit") 
		     && libra_tags != NULL) {

		    RECORD(result_start);

		    for(int j = 0; j < libra_tags->length(); j++) {
		      if((*libra_tags)[j] != NULL) {
			RECORD((*libra_tags)[j]);
			delete (*libra_tags)[j];
		      }
		    } // end for loop over libra_tags

		    delete libra_tags;
		    libra_tags = NULL;
		    RECORD(result_stop);
		  }
		  RECORD((*tags)[k]);
		}

		delete (*tags)[k];
	      }

	    } // next tag

	    delete tags;
	    tags = NULL;
	  }

	  if(libra_result_ != NULL) {
	    delete libra_result_;
	    libra_result_ = NULL;
	  }

	} // end filter_memory 

      } else {
	// if not filtered

      } // filtered

      if(! collected && tag != NULL)
	delete tag;

      data = strstr(data+1, "<");

    } // next tag
  } // next line

  fin.close();
  fout.close();

  if (testType!=NO_TEST) {
    //
    // regression test stuff - bpratt Insilicos LLC, Nov 2005
    //

    TagListComparator("LibraPeptideParser",testType,test_tags,testFileName);
    delete[] testFileName;
    for(int k = test_tags.length(); k--;) {
      delete test_tags[k];
    }
  }

  if(! overwrite(xmlfile, outfile.c_str(), "</msms_pipeline_analysis>"))
    cerr << "error: no libra data written to file " << xmlfile << endl;
  else
    cout << "...done! " << endl;

  delete timestamp;
  delete result_start;
  delete result_stop;
  delete libra_filter;
  delete libra_time_filter;
  delete libra_summ_filter;

  delete[] nextline;
}


void LibraPeptideParser::setFilter(Tag* tag) {
  if(tag == NULL)
    return;

  if(filter_memory_) {
    filter_memory_ = False;
    filter_ = False;
  }

  if(! strcmp(tag->getName(), "spectrum_query")){ 
    if(tag->isStart()) {
      //tag->print();
      filter_ = True;
    } else {
      filter_memory_ = True;
    }
  }
}
