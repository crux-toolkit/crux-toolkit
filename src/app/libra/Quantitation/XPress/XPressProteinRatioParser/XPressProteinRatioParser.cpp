/*

Program       : XPressProteinRatioParser
Author        : Andrew Keller <akeller@systemsbiology.org>
                *Jimmy Eng (jeng@systemsbiology.org>
Date          : 11.27.02
SVN Info      : $Id: XPressProteinRatioParser.cpp 8022 2020-02-12 21:47:21Z mhoopmann $

Computes XPRESS ratios and errors for proteins, then overwrites
that information onto ProteinProphet XML

Copyright (C) 2003 Andrew Keller, Jimmy Eng

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

#include "XPressProteinRatioParser.h"
#include "Parsers/Parser/TagListComparator.h" // regression test stuff - bpratt Insilicos LLC, Nov 2005

XPressProteinRatioParser::XPressProteinRatioParser(const char* protxmlfile, const Boolean use_intensities, const char *testMode) : Parser("xpress") { 
  iNumRawData_=0;
  heavy2light_ = False;
  calc_intensity_ratios_ = use_intensities;

  parser_ = NULL;
  searchHits_ = NULL;
  testMode_ = testMode?strdup(testMode):NULL;

  init(protxmlfile);
}


XPressProteinRatioParser::XPressProteinRatioParser(Array<const char*> &input_pepxmlfiles, 
                                                   const peplist &peptides, 
                                                   double min_pep_prob ) : Parser("xpress") { 
  iNumRawData_=0;
  heavy2light_ = False;
  calc_intensity_ratios_ = False;  // not implemented for XPressCGIProtein

  parser_ = NULL;
  searchHits_ = NULL;
  testMode_ = NULL;
  for (int i=0;i<input_pepxmlfiles.size();i++) {
     input_pepxmlfiles_.insertAtEnd(strCopy(input_pepxmlfiles[i]));
  }
  cachePepXML();
  getRatio(&peptides, min_pep_prob);
}


XPressProteinRatioParser::~XPressProteinRatioParser() {
  for(int k = 0; k < input_pepxmlfiles_.length(); k++)
    if(input_pepxmlfiles_[k] != NULL)
      delete[] input_pepxmlfiles_[k];

  delete searchHits_;
  free(testMode_);
}

void XPressProteinRatioParser::parse(const char* protxmlfile) {
  //open file and pass along
  char* data = NULL;
  Tag* tag;
  Boolean heavy2light = False;

  RACI fin(protxmlfile); // can read gzipped xml
  if(! fin) {
     cout << "XPressProteinRatioParser: error opening " << protxmlfile << std::endl;
    exit(1);
  }

  //
  // regression test stuff - bpratt Insilicos LLC, Nov 2005
  //
  Array<Tag*> test_tags;
  eTagListFilePurpose testType;
  char *testFileName=NULL;
  checkRegressionTestArgs(testMode_,testType);
  if (testType !=NO_TEST) {
     std::string options;
     testFileName = constructTagListFilename(protxmlfile, // input file
        testMode_, // program args
        "XPressProteinRatioParser", // program name
        testType); // user info output
  }

#define RECORD(tag) {(tag)->write(fout);if (testType !=NO_TEST) {test_tags.insertAtEnd(new Tag(*(tag)));}}


  double MIN_PEP_WT = 0.5;
  double MIN_PEP_PROB = 0.8;
  Array<Tag*>* tags = NULL;
  peplist* peps = NULL;
  double next_prot_prob = 2.0;
  double MIN_PROB = 0.2;
  Boolean done = False;
  Boolean old_summary = False;

  TagFilter* summary_filter = new TagFilter("analysis_summary", 1);
  summary_filter->enterRequiredAttributeVal("analysis", getName());

  TagFilter* ratio_filter = new TagFilter("analysis_result");
  ratio_filter->enterRequiredAttributeVal("analysis", getName());


  Tag* result_start = new Tag("analysis_result", True, False);
  result_start->setAttributeValue("analysis", getName());
  Tag* result_stop = new Tag("analysis_result", False, True);
  Tag* summary_start = new Tag("analysis_summary", True, False);
  summary_start->setAttributeValue("analysis", getName());
  summary_start->setAttributeValue("time", time_);
  summary_start->setAttributeValue("id", "1");
  Tag* summary_stop = new Tag("analysis_summary", False, True);

  char** altpeps = NULL;

  // construct a tmpfile name based on xmlfile
  std::string outfile = make_tmpfile_name(protxmlfile);
  //cerr << "writing data to " << outfile << endl;
  ofstream fout(outfile.c_str());
  if(! fout) {
    cout << "cannot write output to file " << outfile << endl;
    exit(1);
  }



  char *nextline = new char[line_width_];
  while(fin.getline(nextline, line_width_)) {
    data = strchr(nextline, '<');
    while(data != NULL) {
      tag = new Tag(data);

      //tag->write(cout);

      setFilter(tag);
      Boolean stored = False;
      // parse input file names
      if(tag->isStart() && ! strcmp(tag->getName(), "protein_summary_header")) {
	const char* files = tag->getAttributeValue("source_files_alt");
	// parse through
	int i = 0;
	int last_i = 0;
	char* nextfile;
	while(files[i]) {
	  i++;
	  if(files[i] == '+' || !files[i]) {
	    nextfile = new char[i - last_i + 1];
	    for(int z = last_i; z < i; z++)
	      nextfile[z-last_i] = files[z];
	    nextfile[i - last_i] = 0;
	    unCygwinify(nextfile); // no effect in Cygwin builds
	    input_pepxmlfiles_.insertAtEnd(nextfile);
	    last_i = i+1;
	  }
	} // while

	// now call to set


      } // if protein summary header

      if(! done && tag->isStart() && ! strcmp(tag->getName(), "protein_group")) {
	const char* next = tag->getAttributeValue("probability");
	done = (next == NULL || atof(next) < MIN_PROB);
      }
      if(! done && tag->isStart() && ! strcmp(tag->getName(), "protein")) {
	const char* next = tag->getAttributeValue("probability");
	next_prot_prob = next == NULL ? 0.0 : atof(next);
      }

      // filter out entries below min prob, and exclude all previous ASAPRatio calculations
      if(! done && filter_ && next_prot_prob >= MIN_PROB) {
	// new protein

	if(tags == NULL)
	  tags = new Array<Tag*>;

	if(! ratio_filter->filter(tag)) {
	  tags->insertAtEnd(tag);
	  stored = True;
	}

	if(peps == NULL)
	  peps = new peplist;

	if(tag->isStart() && ! strcmp(tag->getName(), "peptide")) {

	  // check that weight and prob are above minimum.....
	  if(atof(tag->getAttributeValue("weight")) > MIN_PEP_WT && atof(tag->getAttributeValue("nsp_adjusted_probability")) >= MIN_PEP_PROB)
	    enterUnique(peps, tag->getAttributeValue("peptide_sequence"));
	}

	else if(filter_memory_ && tags != NULL && peps != NULL) {
	  // add the asap pro
	  Tag* protag = NULL;
	  if(peps->size() > 0  && 
	     getRatio(peps, MIN_PEP_PROB)) {

	    //cout << "ready to compute ratio..." << endl; 
	    //getRatio(peps, MIN_PEP_PROB);

	    char* peptidestring = getPeptideString(peps, "+");

	    protag = new Tag("XPressRatio", True, True);
	    if(protag == NULL) {
	      cerr << "error with protag" << endl;
	      exit(1);
	    }

	    char next[50];
	    sprintf(next, "%0.2f", pRatio_.dRatio);
	    protag->setAttributeValue("ratio_mean", next);
	    sprintf(next, "%0.2f", pRatio_.dStdDev);
	    protag->setAttributeValue("ratio_standard_dev", next);
	    sprintf(next, "%d", pRatio_.iNumPeptides);
	    protag->setAttributeValue("ratio_number_peptides", next);

	    // now the heavy2light
	    if(pRatio_.dRatio == 0.0) {
	      protag->setAttributeValue("heavy2light_ratio_mean", "999.");
	      protag->setAttributeValue("heavy2light_ratio_standard_dev", "0.00");
	    }
	    else if(pRatio_.dRatio >= 999.0) {
	      protag->setAttributeValue("heavy2light_ratio_mean", "0.00");
	      protag->setAttributeValue("heavy2light_ratio_standard_dev", "0.00");
	    }
	    else {
	      sprintf(next, "%0.2f", pRatio_.dh2lRatio);
	      protag->setAttributeValue("heavy2light_ratio_mean", next);
	      sprintf(next, "%0.2f", pRatio_.dh2lStdDev);
	      protag->setAttributeValue("heavy2light_ratio_standard_dev", next);
	    }
	    protag->setAttributeValue("peptide_string", peptidestring);

	    if(peptidestring != NULL)
	      delete [] peptidestring;

	  } // only if have enough peps
	  int k;
	  for(k = 0; k < tags->length(); k++) {
	    RECORD((*tags)[k]); //print();
	    if(k == 0 && protag != NULL) {
	      //protag->write(cout);
	      RECORD(result_start);
	      RECORD(protag); //print();
	      RECORD(result_stop);
	      delete protag;
	    }
	  }

	  if(tags != NULL) {
	    for(k = 0; k < tags->length(); k++)
	      if((*tags)[k] != NULL)
		delete (*tags)[k];

	    delete tags;
	    tags = NULL;
	  }

	  if(peps != NULL) {
	    delete peps;
	    peps = NULL;
	  }

	}
      }
      else {
	// add after first tag the additional asap_pro tag...
	if(! summary_filter->filter(tag))
	  RECORD(tag); //print();
	if(tag->isEnd() && ! strcmp(tag->getName(), "protein_summary_header")) {
	  Tag* summary = new Tag("XPress_analysis_summary", True, True);
	  // get time info
	  char next[20];
	  sprintf(next, "%0.2f", MIN_PEP_PROB);
	  summary->setAttributeValue("min_peptide_probability", next);
	  sprintf(next, "%0.2f", MIN_PEP_WT);
	  summary->setAttributeValue("min_peptide_weight", next);
	  sprintf(next, "%0.2f", MIN_PROB);
	  summary->setAttributeValue("min_protein_probability", next);

	  if(heavy2light_)
	    summary->setAttributeValue("reference_isotope", "light");
	  else
	    summary->setAttributeValue("reference_isotope", "heavy");

	  if(calc_intensity_ratios_)
	    summary->setAttributeValue("ratio_basis", "intensity");
	  else
	    summary->setAttributeValue("ratio_basis", "area");

	  RECORD(summary_start);
	  RECORD(summary); //print();
	  RECORD(summary_stop);
	  delete summary;
	}
      }

      if(! stored) {
	delete tag;
	tag = NULL;
      }

      data = strchr(data+1, '<');
    }
  }
  fin.close();
  fout.close();

  if(! overwrite(protxmlfile, outfile.c_str(), "</protein_summary>")) {
    cout << "error: no xpress data written to file " << protxmlfile << endl;
  }

  if (testType!=NO_TEST) {
     //
     // regression test stuff - bpratt Insilicos LLC, Nov 2005
     //
     TagListComparator("XPressProteinRatioParser",testType,test_tags,testFileName);
     delete[] testFileName;
     for(int k = test_tags.length(); k--;) {
        delete test_tags[k];
     }
  }

  delete [] nextline;

  delete summary_start;
  delete summary_stop;
  delete result_start;
  delete result_stop;
  delete summary_filter;
  delete ratio_filter;
}

void XPressProteinRatioParser::setFilter(Tag* tag) {
  if(tag == NULL)
    return;

  if(filter_memory_) {
    filter_memory_ = False;
    filter_ = False;
  }

  if(! strcmp(tag->getName(), "protein")){ 
    if(tag->isStart()) {
      filter_ = True;
    }else{
      filter_memory_ = True;
    }
  }

}


void XPressProteinRatioParser::enterUnique(peplist* uniques, const char* next) {
   if (!uniques->find(next)) {
      uniques->add(next);
   }
}


// must substitute # -> ~
char* XPressProteinRatioParser::getPeptideString(peplist* peps, const char* link) {

  char* output = NULL;
  if(peps == NULL || peps->size() == 0)
    return output;

  std::string str;
  for(int k = 0; k < peps->size(); k++) {
    if(k) {
      str += link;
    }
    str += peps->getNthCharptr(k); // write out in order read in
  }
  output = new char[str.length()+1];
  strcpy(output,str.c_str());

  for (char *c=output;*c;c++) {
     if ('#'==*c) {
        *c = '~';
     }
  }
  return output;
}

Boolean XPressProteinRatioParser::getRatio(const peplist* peptides, double dProbability) { //char** peps, double minpepprob) {

  double mean = 0.0;
  double meansq = 0.0;
  int num = 0;

  if (!searchHits_) {
     // create a pepXML cache, ignoring all but the interesting tags and attributes
     cachePepXML();
  }

  parser_ = new XPressGroupPeptideParser(*searchHits_, peptides, dProbability);

  if(parser_ != NULL) {
    pRatio_ = parser_->getRatio();
    delete parser_;
    if (pRatio_.iNumPeptides == 0) {
      return False;
    }
  }
  else {
    cout << "Error: null parser for inputfiles: ";
    int k;
    for(k = 0; k < input_pepxmlfiles_.length(); k++)
      cout << input_pepxmlfiles_[k] << " ";
    cout << " and peptides: ";
    for(k = 0; k < peptides->size(); k++) {
       if(k) {
          cout << " ";
       }
       cout << peptides->getNthCharptr(k); // write out in order read in
    }
    cout << endl;
    exit(1);
  }

  return True;

}

//
// NOTE:
// if you find you are interested in tags named other than
//    "search_hit","xpressratio_result","peptideprophet_result"
// or attributes other than
//    "hit_rank","peptide","heavy_area","light_area","probability",
//  you will need to alter these lists:
//
static const char *XPressProteinRatioParser_interesting_tags[] = {
   "search_hit",
   "xpressratio_result",
   "peptideprophet_result",
   NULL
};
static const char *XPressProteinRatioParser_interesting_attrs[] = 
   {"hit_rank","peptide","heavy_area","light_area","heavy_intensity","light_intensity","probability",NULL};



void XPressProteinRatioParser::cachePepXML() {

   searchHits_ = new XPressRatioSearchHitCache;

   Tag* tag = NULL;

   char* data = NULL;

   Boolean analyze = False;
   Boolean ratio_found = False;
   Boolean search_score_found = False;

   XPressRatioSearchHit hit;
   TagInclusionTest tagtest(XPressProteinRatioParser_interesting_tags,
                            XPressProteinRatioParser_interesting_attrs);

   char *nextline = new char[line_width_];
   std::string pepname;

   for(int k = 0; k < input_pepxmlfiles_.length(); k++) {
      RACI fin(input_pepxmlfiles_[k]); // can read gzipped xml
      if(! fin) {
         cout << "error opening " << input_pepxmlfiles_[k] << endl;
         exit(1);
      }
      while(fin.getline(nextline, line_width_)) {

            data = strstr(nextline, "<");
            while(data != NULL) {
               tag = new Tag(data,&tagtest);

               if(tag && tag->getName()) {

                  bool isSearchHitTag = ! strcmp(tag->getName(), "search_hit");
                  if (tag->isStart() && isSearchHitTag) {
                     ratio_found = False;
                     pepname = tag->getAttributeValue("peptide");
                  }

                  if(tag->isStart() && isSearchHitTag && ! strcmp(tag->getAttributeValue("hit_rank"), "1")) {
                     analyze = True;
                  }
                  else if(analyze) {                    
                     if(tag->isStart() && ! strcmp(tag->getName(), "xpressratio_result")) {

		       if (calc_intensity_ratios_) {
			 if (tag->getAttributeValue("heavy_intensity") == NULL ||
			     tag->getAttributeValue("light_intensity") == NULL) {
			   cerr << "ERROR: Could not find XPRESS intensities in pepXML source file(s). Exiting without calculating protein ratios." << endl;
			   exit(1);
			 }
			 hit.heavy_ = atof(tag->getAttributeValue("heavy_intensity"));
			 hit.light_ = atof(tag->getAttributeValue("light_intensity"));
		       }
		       else {
			 hit.heavy_ = atof(tag->getAttributeValue("heavy_area"));
			 hit.light_ = atof(tag->getAttributeValue("light_area"));
		       }

                        ratio_found = True;
                     }
                     else if(tag->isStart() && ! strcmp(tag->getName(), "peptideprophet_result")) {
                        hit.probability_ = atof(tag->getAttributeValue("probability"));
                     }
                     else if(tag->isEnd() && isSearchHitTag && ratio_found) { // process
                        searchHits_->addHit(pepname.c_str(),hit);                          
                     } // if process

                     if (tag->isEnd() && isSearchHitTag) { 
                        hit.probability_ = -1.0;
                        analyze = False;
                        ratio_found = False;
                     }
                  }

               } //  if not null
               delete tag;
               data = strstr(data+1, "<");
            } // next tag
         
      } // next line
      fin.close();
   } // next inputfile
   
   //fout.close();
   
   delete [] nextline;  
}

