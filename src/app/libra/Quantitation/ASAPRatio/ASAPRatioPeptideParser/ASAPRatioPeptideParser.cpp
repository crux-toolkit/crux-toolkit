/*
Program       : ASAPRatioPeptideParser
Author        : X.Li and Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: ASAPRatioPeptideParser.cpp 9064 2023-11-29 02:08:59Z dshteyn $


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

#include "ASAPRatioPeptideParser.h"

#include "Common/TPPVersion.h" // contains version number, name, revision

#include "Parsers/Parser/TagListComparator.h" // regression test stuff - bpratt Insilicos LLC, Nov 2005
#include "ASAPRatioPeptideParserTagListComparator.h" // regression test stuff - bpratt Insilicos LLC, Nov 2005

#include <string>
#include <math.h>
#define CACHED_RAMP_HOME
#include "Parsers/mzParser/cached_ramp.h"

using namespace mzParser;

//#define DEBUG_NOFREE

int comp_chars(const void*a, const void*b) {
  char* num1_ = (char*)a;
  char* num2_ = (char*)b;
  if(*num1_ < *num2_) return -1;
  if(*num1_ == *num2_) return 0;
  /* if(*num1_ > *num2_) */ return 1;
}

ASAPRatioPeptideParser::ASAPRatioPeptideParser(const char* xmlfile, const InputStruct &options, const char *testMode, double mzBound) : Parser("asapratio") {
  // default settings
  mzBound_ = mzBound;
  modAAs_ = NULL;
  prtnAAs_ = NULL;

  spectrum_ = new spectStrct();
  spectrum_->size = 10;
  spectrum_->xval = (double *) calloc(_MZXML_READER_MXSCANNUM_, sizeof(double));
  spectrum_->yval = (double *) calloc(_MZXML_READER_MXSCANNUM_, sizeof(double));
  pInput_ = options;
  verbose_ = false;

  xmlIndx_ = NULL;
  m_XMLfile_state = -1; // no mzxml file yet

  strcpy(paired_labels_, pInput_.szXpressResidues);

  testMode_ = testMode?strdup(testMode):NULL; // regression test stuff - bpratt Insilicos LLC, Nov 2005

  // now alphabetize....
  qsort(paired_labels_, strlen(paired_labels_), sizeof(char), (int(*)(const void*, const void*)) comp_chars);

  elution_ = pInput_.bUseSameScanRange ? 0 : -1; //1; // 0 for elute at same time, -1 if heavy elutes first, 1 if light elutes first
  areaFlag_ = pInput_.bXpressLight1; // how to present the ratio (not yet in use)

  mzXMLfile_[0] = 0;

  // HERE ARE THE CORRECTIONS TO DERIVE MONOISOTOPIC MASSES FROM AVERAGE MASSES
  // cl icat light, cl icat heavy, D+14, D+17, c+14, c+17, E+14, E+17
  double averagevalues[] = {330.26, 339.27, 129.08859, 132.08859, 31.0027, 34.0027, 143.11549, 146.11549};
  double monoisotopicvalues[] = {330.136, 339.166, 129.0429, 132.0429, 31.0187, 34.0187, 143.0586, 146.0586};
  int num_vals = sizeof(averagevalues)/sizeof(double);
  diff_ = 0.0005; // mass error allowed for recognizing average mass to be corrected

  average_mods_ = new Array<double>;
  monoisostopic_mods_ = new Array<double>;
  for(int k = 0; k < num_vals; k++) {
    average_mods_->insertAtEnd(averagevalues[k]); // cl icat light
    monoisostopic_mods_->insertAtEnd(monoisotopicvalues[k]);
  }

  error_ = 0.3;
  static_quant_ =  pInput_.staticQuant;

  static_status_ = 0;
  static_pairs_ = NULL;
  static_symbol_ = '\'';
  static_nterm_ = False;
  static_cterm_ = False;
 
  if(static_quant_) {
    static_pairs_ = new Array<staticQuant*>;

    if(strlen(pInput_.labelMasses) > 1) { // must enter user specified label masses
      char nextmass[100];
      int index = 1;
      int start = 0;

      while(index < (int) strlen(pInput_.labelMasses)) {
	while(index < (int) strlen(pInput_.labelMasses) &&
	      (pInput_.labelMasses[index] < 'a' || pInput_.labelMasses[index] > 'z') && 
	      (pInput_.labelMasses[index] < 'A' || pInput_.labelMasses[index] > 'Z')) 
	  index++;
	// process
	strncpy(nextmass, pInput_.labelMasses + start + 1, index - start - 1);
	nextmass[index - start - 1] = 0;
	if (verbose_) {
	  cout << "entering " << pInput_.labelMasses[start] << " with mass " << nextmass << endl;
	}
	enterStaticQuant(pInput_.labelMasses[start], atof(nextmass), error_);
	start = index;
	index++;
      }
    }
  }

#ifdef USE_STD_MODS
  modinfo_ = NULL;
  compute_peptide_mass_ = True;
#endif

  init(xmlfile);
}

ASAPRatioPeptideParser::~ASAPRatioPeptideParser() {
  if(static_quant_ && static_pairs_ != NULL) 
    delete static_pairs_;

  if(average_mods_ != NULL)
    delete average_mods_;

  if(monoisostopic_mods_ != NULL)
    delete monoisostopic_mods_;

  delete spectrum_;

  freeXmlIndx();
  free(testMode_);
}


void ASAPRatioPeptideParser::parse(const char* xmlfile) {
  pepIdx_ = 0;
  if(static_quant_)
    setStaticQuantInfo(xmlfile);

  char* engine = NULL;
  char* enzyme = NULL;
  char* massspec = NULL;
  Array<Tag*>* tags = NULL;
  Tag* tag = NULL;

  char *nextline = new char[line_width_];
  char* data = NULL;

  Boolean first_hit = False;
  
  double prob = -100;
  double iprob = -100;

  //
  // regression test stuff - bpratt Insilicos LLC, Nov 2005
  //
  Array<Tag*> test_tags;
  eTagListFilePurpose testType;
  char *testFileName=NULL;
  checkRegressionTestArgs(testMode_,testType);
  if (testType!=NO_TEST) {
     std::string options;
     testFileName = constructTagListFilename(xmlfile, // input file
        testMode_, // program args
        "ASAPRatioPeptideParser", // program name
        testType); // user info output
  }

#define RECORD(tag) {(tag)->write(fout);if (testType!=NO_TEST) {test_tags.insertAtEnd(new Tag(*(tag)));}}

  Tag* summary_start = new Tag("analysis_summary", True, False);
  summary_start->setAttributeValue("analysis", getName());
  summary_start->setAttributeValue("time", time_);
  Tag* summary_stop = new Tag("analysis_summary", False, True);

  Tag* timestamp_start = new Tag("analysis_timestamp", True, False);
  timestamp_start->setAttributeValue("analysis", getName());
  timestamp_start->setAttributeValue("time", time_);
  timestamp_start->setAttributeValue("id", "1");
  Tag* timestamp_stop = new Tag("analysis_timestamp", False, True);

  Tag* timestamp = new Tag("asapratio_timestamp", True, True);

  Tag* summary = NULL;

  summary = getSummaryTag(pInput_); //new Tag("asapratio_summary", True, False);

  Array<Tag*>* paramTags = new Array<Tag*>;
  if (pInput_.bQuantHighBackGrnd) {
    Tag* param = new Tag("parameter", True, True);
    param->setAttributeValue("name", "quantHighBG");
    param->setAttributeValue("value", "True");
    paramTags->insertAtEnd(param);
  }
  if (pInput_.bUseWaveletSmoothing) {
    Tag* param = new Tag("parameter", True, True);
    param->setAttributeValue("name", "wavelet");
    param->setAttributeValue("value", "True");
    paramTags->insertAtEnd(param);
  }
  else {
    Tag* param = new Tag("parameter", True, True);
    param->setAttributeValue("name", "wavelet");
    param->setAttributeValue("value", "False");
    paramTags->insertAtEnd(param);
  }
  if (pInput_.bZeroAllBackGrnd) {
    Tag* param = new Tag("parameter", True, True);
    param->setAttributeValue("name", "zeroBG");
    param->setAttributeValue("value", "True");
    paramTags->insertAtEnd(param);
  }
  if (pInput_.dMassTol > 0 && pInput_.dMassTol < 1) {
    Tag* param = new Tag("parameter", True, True);
    param->setAttributeValue("name", "mzBound");
    char text[100];
    sprintf(text, "%0.4f", pInput_.dMassTol);
    param->setAttributeValue("value", text);
    paramTags->insertAtEnd(param);
  }

  Tag* result_start = new Tag("analysis_result", True, False);
  result_start->setAttributeValue("analysis", getName());
  Tag* result_stop = new Tag("analysis_result", False, True);

  // other global stuff here....
  TagFilter* asap_filter = new TagFilter("analysis_result");
  asap_filter->enterRequiredAttributeVal("analysis", getName());
  TagFilter* asap_summ_filter = new TagFilter("analysis_summary");
  asap_summ_filter->enterRequiredAttributeVal("analysis", getName());
  TagFilter* asap_time_filter = new TagFilter("analysis_timestamp");
  asap_time_filter->enterRequiredAttributeVal("analysis", getName());

  // construct a tmpfile name based on xmlfile
  std::string outfile = make_tmpfile_name(xmlfile);
  ofstream fout(outfile.c_str());
  if(! fout) {
    cerr << "cannot write output to file " << outfile << endl;
    exit(1);
  }

  Boolean first = False;
  int ratio_tags_written = 0;

#ifndef USE_STD_MODS
  Array<char*>* modifications = NULL;
  Array<char*>* lightpartners = NULL;
  Array<char*>* heavypartners = NULL;
#endif

  Array<Tag*>* asap_tags = NULL;
  int asap_index = 0;

#ifndef USE_STD_MODS
  char lightstring[500];
  char heavystring[500];
#endif

  Boolean collected = False;
  monoisotopic_ = True; 

#ifdef USE_STD_MODS
  Boolean mod_on = False;
  Array<Tag*>* modification_tags = NULL;
#endif

  RACI fin(xmlfile); // can read gzipped xml
  if(! fin) {
    cerr << "ASAPRatioPeptideParser: error opening " << xmlfile << endl;
    exit(1);
  }
  while(fin.getline(nextline, line_width_)) {
    //cout << "next: " << nextline << endl;

    data = strstr(nextline, "<");
    while(data != NULL) {
      tag = new Tag(data);
      setFilter(tag);
      //tag->write(cout);
      collected = False;
      if(! asap_filter->filter(tag) && ! asap_summ_filter->filter(tag) && ! asap_time_filter->filter(tag)) {

	if(tag->isStart() && ! strcmp(tag->getName(), "asapratio_timestamp")) {
	  if (verbose_) {
	    tag->write(cout);
	  }
	  if(asap_summ_filter->filter(tag)) {
	    if (verbose_) {
	      cout << "filtered" << endl;
	    }
	  } else {
	    if (verbose_) {
	      cout << "not filtered" << endl;
	    }
	  }
	}

	if(tag->isStart() && ! strcmp(tag->getName(), "msms_pipeline_analysis")) {
	  RECORD(tag);
	  RECORD(summary_start);
	  RECORD(summary);
	  for (int i=0; i < paramTags->size(); i++) {
	    RECORD((*paramTags)[i]);
	    delete (*paramTags)[i];
	  }
	  delete paramTags;
	  RECORD(summary_stop);
#ifndef DEBUG_NOFREE
	  delete summary_start;
	  delete summary;
	  delete summary_stop;
#endif
	}
	else if(tag->isStart() && ! strcmp(tag->getName(), "msms_run_summary")) {
	  rampConstructInputFileName(mzXMLfile_, sizeof(mzXMLfile_), tag->getAttributeValue("base_name"));
	  freeXmlIndx(); // starting on a new mzXML file
#ifdef USE_STD_MODS
	  if(! static_quant_) {
	    for(int k = 0; k < 26; k++) {
	      light_label_masses_[k] = 0.0;
	      heavy_label_masses_[k] = 0.0;
	    }

	    int k = 0;
	    int start = 0;
	    char nextmass[100];
	    while(k < (int) strlen(pInput_.labelMasses)) {
	      char next_res = pInput_.labelMasses[start];
	      k++;
	      while(k < (int) strlen(pInput_.labelMasses) &&
		    (pInput_.labelMasses[k] < 'a' || pInput_.labelMasses[k] > 'z') && 
		    (pInput_.labelMasses[k] < 'A' || pInput_.labelMasses[k] > 'Z')) {
		k++;
	      }

	      strncpy(nextmass, pInput_.labelMasses + start + 1, k - start - 1);
	      nextmass[k - start - 1] = 0;
	      start = k;

	      double usermass = atof(nextmass);
	      double monomass = ResidueMass::getMass(next_res, monoisotopic_);
	      if (next_res == 'n') {
		if (usermass > monomass) {
		  light_nterm_mass_ = monomass;
		  heavy_nterm_mass_ = usermass;
		}
		else {
		  heavy_nterm_mass_ = monomass;
		  light_nterm_mass_ = usermass;
		}
	      }
	      else if (next_res == 'c') {
		if (usermass > monomass) {
		  light_cterm_mass_ = monomass;
		  heavy_cterm_mass_ = usermass;
		}
		else {
		  heavy_cterm_mass_ = monomass;
		  light_cterm_mass_ = usermass;
		}
	      }
	      else {
		if (usermass > monomass) {
		  light_label_masses_[next_res - 'A'] = monomass;
		  heavy_label_masses_[next_res - 'A'] = usermass;
		}
		else {
		  heavy_label_masses_[next_res - 'A'] = monomass;
		  light_label_masses_[next_res - 'A'] = usermass;
		}
	      }
	    }
	  }
#endif

#ifndef USE_STD_MODS
	  modifications = new Array<char*>;
	  lightpartners = new Array<char*>;
	  heavypartners = new Array<char*>;
	  // fill them up with pointers for each user specified paired label
	  for(int k = 0; k < strlen(paired_labels_); k++) {
	    lightpartners->insertAtEnd(NULL);
	    heavypartners->insertAtEnd(NULL);
	  }
#endif
	  pInput_.iAnalysisFirstScan = 1;

	  first = True;
	  RECORD(tag);
	}
	else if(tag->isEnd() && ! strcmp(tag->getName(), "msms_run_summary")) {
	  RECORD(tag);
#ifndef USE_STD_MODS
	  if(modifications != NULL) {
#ifndef DEBUG_NOFREE
	    for(int k = 0; k < modifications->length(); k++)
	      if((*modifications)[k] != NULL)
		delete (*modifications)[k];
	    delete modifications;
#endif
	    modifications = NULL;
	  }

	  if(lightpartners != NULL) {
#ifndef DEBUG_NOFREE
	    for(int k = 0; k < lightpartners->length(); k++)
	      if((*lightpartners)[k] != NULL)
		delete (*lightpartners)[k];
	    delete lightpartners;
#endif
	    lightpartners = NULL;
	  }
	  if(heavypartners != NULL) {
#ifndef DEBUG_NOFREE
	    for(int k = 0; k < heavypartners->length(); k++)
	      if((*heavypartners)[k] != NULL)
		delete (*heavypartners)[k];
	    delete heavypartners;
#endif
	    heavypartners = NULL;
	  }
#endif
#ifndef DEBUG_NOFREE
	  if(static_quant_) {
	    static_status_ = 0; // reset
	    static_nterm_ = False;
	    static_cterm_ = False;
	  }
	  if(modAAs_ != NULL)
	    delete modAAs_;
	  if(prtnAAs_ != NULL)
	    delete prtnAAs_;
#endif
	}

#ifdef USE_STD_MODS
	// have a modification worth recording here
	else if(tag->isStart() && ! strcmp(tag->getName(), "terminal_modification") && 
		strchr(pInput_.szXpressResidues, tag->getAttributeValue("terminus")[0]) != NULL) { 
	  char next_res = tag->getAttributeValue("terminus")[0];
	  double nextmass = atof(tag->getAttributeValue("mass"));
	  // adjust average masses to monoisotpic

	  if(monoisotopic_)
	    nextmass = getMonoisotopicEquivalent(nextmass, diff_, NULL);

	  if(next_res == 'n') {
	    if(! strcmp(tag->getAttributeValue("variable"), "Y")) { // must be heavy

	      if(static_quant_ && strstr(paired_labels_, tag->getAttributeValue("terminus")) != NULL) {
		cout << "error: cannot have variable modified " << tag->getAttributeValue("terminus") << " in static mode" << endl;
		exit(1);
	      }
	      //if(heavy_nterm_mass_ > 0.0) {
	      // error
	      //exit(1);
	      //}

	      heavy_nterm_mass_ = nextmass;
	    } // variable
	    else {
	      // static check here.....?
	      static_nterm_ = True;
	      if(strchr(paired_labels_, 'n') != NULL && static_quant_) {
		int nextstatus = getStaticQuantStatus('n', atof(tag->getAttributeValue("mass")), error_);
		if(static_status_ && nextstatus != static_status_) {
		  cout << "error: have mixture of light and heavy static modifications" << endl;
		  cout << "n: " << tag->getAttributeValue("mass") << " " << nextstatus << " " << static_status_ << endl;
		  exit(1);
		}
		static_status_ = nextstatus;
		if(! static_status_)
		  static_status_ = -1; // the default is light
		if(static_status_ == -1)  // light mod
		  light_nterm_mass_ = nextmass;
		else
		  heavy_nterm_mass_ = nextmass;
	      } // if n terminus is a quant label and static quant
	      else {
		// place in light for now
		light_nterm_mass_ = nextmass;
	      }
	    }

	    if (light_nterm_mass_ > heavy_nterm_mass_) {
	      double tmp = light_nterm_mass_;
	      light_nterm_mass_ = heavy_nterm_mass_;
	      heavy_nterm_mass_ = tmp;
	    }
	  } // n terminal

	  else if(next_res == 'c') {
	    if(! strcmp(tag->getAttributeValue("variable"), "Y")) { // must be heavy

	      if(static_quant_ && strstr(paired_labels_, tag->getAttributeValue("terminus")) != NULL) {
		cout << "error: cannot have variable modified " << tag->getAttributeValue("terminus") << " in static mode" << endl;
		exit(1);
	      }
	      // if(heavy_cterm_mass_ > 0.0) {
	      // error
	      //exit(1);
	      //}

	      heavy_cterm_mass_ = nextmass;
	    } // variable
	    else {
	      // static check here?
	      static_cterm_ = True;
	      if(strchr(paired_labels_, 'c') != NULL && static_quant_) {
		int nextstatus = getStaticQuantStatus('c', atof(tag->getAttributeValue("mass")), error_);
		if(static_status_ && nextstatus != static_status_) {
		  cout << "error: have mixture of light and heavy static modifications" << endl;
		  cout << "c: " << tag->getAttributeValue("mass") << " " << nextstatus << " " << static_status_ << endl;
		  exit(1);
		}
		static_status_ = nextstatus;
		if(! static_status_)
		  static_status_ = -1; // the default is light
		if(static_status_ == -1)  // light mod
		  light_cterm_mass_ = nextmass;
		else
		  heavy_cterm_mass_ = nextmass;
	      } // if n terminus is a quant label and static quant
	      else {
		// place in light for now
		light_cterm_mass_ = atof(tag->getAttributeValue("mass"));
	      }
	    }

	    if (light_cterm_mass_ > heavy_cterm_mass_) {
	      double tmp = light_cterm_mass_;
	      light_cterm_mass_ = heavy_cterm_mass_;
	      heavy_cterm_mass_ = tmp;
	    }

	  } // c
	  RECORD(tag);
	} // term modification

	else if(tag->isStart() && ! strcmp(tag->getName(), "aminoacid_modification") && 
		strchr(paired_labels_, tag->getAttributeValue("aminoacid")[0]) != NULL) { 
    
	  if (verbose_) {
	    tag->write(cout);
	  }
	  char next_res = tag->getAttributeValue("aminoacid")[0];
	  double nextmass = atof(tag->getAttributeValue("mass"));
	  double massd = atof(tag->getAttributeValue("massdiff"));
	  // adjust average masses to monoisotpic
	  if(monoisotopic_)
	    nextmass = getMonoisotopicEquivalent(nextmass, diff_, NULL);
	  if(! strcmp(tag->getAttributeValue("variable"), "Y")) { // must be heavy if only one
	    if(static_quant_ && strstr(paired_labels_, tag->getAttributeValue("aminoacid")) != NULL) {
	      cout << "ERROR: cannot have variable modified " << tag->getAttributeValue("aminoacid") << " in static mode" << endl;
	      exit(1);
	    }

	    if (massd < 0) {
	      if (heavy_label_masses_[next_res-'A'] > nextmass) {
		light_label_masses_[next_res-'A']  = nextmass;
	      }
	    }
	    else {
	      if (heavy_label_masses_[next_res-'A'] > 0) {
		if (fabs(heavy_label_masses_[next_res-'A']-nextmass) > 0.01) {
		  cout << "WARNING: Found more than one variable mod on \'" << next_res << "\'. Please make sure to specify a heavy mass for this residue." << endl;
		}
	      }
	      else {
		heavy_label_masses_[next_res-'A']  = nextmass;
	      }
	    }
	  } // variable
	  else {
	    // static check here?
	    if(static_quant_) {
	      int nextstatus = getStaticQuantStatus(next_res, atof(tag->getAttributeValue("mass")), error_);
	      if(static_status_ && nextstatus != static_status_) {
		cout << "ERROR: have mixture of light and heavy static modifications" << endl;
		cout << next_res << " " << tag->getAttributeValue("mass") << " " << nextstatus << " " << static_status_ << endl;
		exit(1);
	      }
	      static_status_ = nextstatus;
	      if(! static_status_)
		static_status_ = -1; // the default is light
	      if(static_status_ == -1) { // light mod
		light_label_masses_[next_res-'A'] = nextmass;
	      }
	      else {
		heavy_label_masses_[next_res-'A'] = nextmass;
	      }
	    } // static quant
	    else {
	      // place in light for now
	      light_label_masses_[next_res-'A'] = nextmass;
	    }
	  } // static mod
	  RECORD(tag);
	} // aa modification
#endif
#ifndef USE_STD_MODS

	// old format (to be replaced)
	else if(tag->isStart() && ! strcmp(tag->getName(), "search_modification")) {
	  RECORD(tag);
	  char* nextmod = NULL;
	  if(! strcmp(tag->getAttributeValue("variable"), "Y")) {
	    nextmod = new char[strlen(tag->getAttributeValue("aminoacid")) + 
			       strlen(tag->getAttributeValue("symbol")) +
			       strlen(tag->getAttributeValue("mass")) + 3];
	    strcpy(nextmod, tag->getAttributeValue("aminoacid"));
	    strcat(nextmod, tag->getAttributeValue("symbol"));
	    strcat(nextmod, "[");

	    double nextmass = atof(tag->getAttributeValue("mass"));

	    // adjust average masses to monoisotpic
	    if(monoisotopic_)
	      nextmass = getMonoisotopicEquivalent(nextmass, diff_, tag);
	    //	    getMonoisotopicEquivalent(nextmass, diff_, tag);
	    strcat(nextmod, tag->getAttributeValue("mass"));
	    strcat(nextmod, "]");
	    // check whether on list of paired aa's
	    for(int k = 0; k < strlen(paired_labels_); k++)
	      if(paired_labels_[k] == (tag->getAttributeValue("aminoacid"))[0]) { // set heavy
		(*heavypartners)[k] = new char[strlen(tag->getAttributeValue("aminoacid")) + 
					       strlen(tag->getAttributeValue("symbol")) + 1];
		strcpy((*heavypartners)[k], tag->getAttributeValue("aminoacid"));
		strcat((*heavypartners)[k], tag->getAttributeValue("symbol"));
		break; // done
	      }
	  }
	  else { // static
	    nextmod = new char[strlen(tag->getAttributeValue("aminoacid")) + 
			       strlen(tag->getAttributeValue("mass")) + 3];
	    strcpy(nextmod, tag->getAttributeValue("aminoacid"));
	    strcat(nextmod, "[");
	    // make correction here if nec
	    double nextmass = atof(tag->getAttributeValue("mass"));

	    // adjust average masses to monoisotpic
	    if(monoisotopic_)
	      nextmass = getMonoisotopicEquivalent(nextmass, diff_, tag);
	    //	    getMonoisotopicEquivalent(nextmass, diff_, tag);
	    strcat(nextmod, tag->getAttributeValue("mass"));
	    strcat(nextmod, "]");
	    // check whether on list of paired aa's
	    for(int k = 0; k < strlen(paired_labels_); k++)
	      if(paired_labels_[k] == (tag->getAttributeValue("aminoacid"))[0]) { // set heavy
		(*lightpartners)[k] = new char[strlen(tag->getAttributeValue("aminoacid")) + 1];
		strcpy((*lightpartners)[k], tag->getAttributeValue("aminoacid"));
		break; // done
	      }
	  }	  
	  modifications->insertAtEnd(nextmod);
	}

	else if(tag->isStart() && ! strcmp(tag->getName(), "aminoacid_modification")) {
	  RECORD(tag);
	  //tag->write(cout);
	  char* nextmod = NULL;

	  if(! strcmp(tag->getAttributeValue("variable"), "Y")) {
	    // check to make sure not on list (if static mode)
	    if(static_quant_ && strstr(paired_labels_, tag->getAttributeValue("aminoacid")) != NULL) {
	      cout << "error: cannot have variable modified " << tag->getAttributeValue("aminoacid") << " in static mode" << endl;
	      exit(1);
	    }

	    nextmod = new char[strlen(tag->getAttributeValue("aminoacid")) + 
			       strlen(tag->getAttributeValue("symbol")) +
			       strlen(tag->getAttributeValue("mass")) + 3];
	    strcpy(nextmod, tag->getAttributeValue("aminoacid"));
	    strcat(nextmod, tag->getAttributeValue("symbol"));
	    strcat(nextmod, "[");
	    double nextmass = atof(tag->getAttributeValue("mass"));
	    // adjust average masses to monoisotpic
	    if(monoisotopic_)
	      nextmass = getMonoisotopicEquivalent(nextmass, diff_, tag);
	    strcat(nextmod, tag->getAttributeValue("mass"));
	    strcat(nextmod, "]");
	    // check whether on list of paired aa's
	    for(int k = 0; k < strlen(paired_labels_); k++)
	      if(paired_labels_[k] == (tag->getAttributeValue("aminoacid"))[0]) { // set heavy
		(*heavypartners)[k] = new char[strlen(tag->getAttributeValue("aminoacid")) + 
					       strlen(tag->getAttributeValue("symbol")) + 1];
		strcpy((*heavypartners)[k], tag->getAttributeValue("aminoacid"));
		strcat((*heavypartners)[k], tag->getAttributeValue("symbol"));
		break; // done
	      }
	    modifications->insertAtEnd(nextmod);
	  }
	  else { // static
	    char** nextmodifs = NULL;

	    // check whether on list of paired aa's
	    for(int k = 0; k < strlen(paired_labels_); k++)
	      if(paired_labels_[k] == (tag->getAttributeValue("aminoacid"))[0]) { // set heavy

		if(static_quant_) {
		  int nextstatus = getStaticQuantStatus(paired_labels_[k], atof(tag->getAttributeValue("mass")), error_);
		  if(static_status_ && nextstatus != static_status_) {
		    cout << "error: have mixture of light and heavy static modifications" << endl;
		    cout << paired_labels_[k] << " " << tag->getAttributeValue("mass") << " " << nextstatus << " " << static_status_ << endl;
		    exit(1);
		  }
		  static_status_ = nextstatus;
		  if(! static_status_)
		    static_status_ = -1; // the default is light

		  nextmodifs = getStaticModification((tag->getAttributeValue("aminoacid"))[0], tag->getAttributeValue("symbol"), static_status_);
		  int nextsymbol = tag->getAttributeValue("symbol") != NULL ? strlen(tag->getAttributeValue("symbol")) : 0;

		  if(static_status_ == -1) { // light mod
		    (*lightpartners)[k] = new char[2 + nextsymbol];
		    (*lightpartners)[k][0] = paired_labels_[k];
		    (*lightpartners)[k][1] = 0;
		    if(nextsymbol)
		      strcat((*lightpartners)[k], tag->getAttributeValue("symbol"));
		    (*heavypartners)[k] = new char[3];
		    (*heavypartners)[k][0] = paired_labels_[k];
		    (*heavypartners)[k][1] = static_symbol_;
		    (*heavypartners)[k][2] = 0;
		  }
		  else { // heavy
		    (*lightpartners)[k] = new char[3];
		    (*lightpartners)[k][0] = paired_labels_[k];
		    (*lightpartners)[k][1] = static_symbol_;
		    (*lightpartners)[k][2] = 0;
		    (*heavypartners)[k] = new char[2 + nextsymbol];
		    (*heavypartners)[k][0] = paired_labels_[k];
		    (*heavypartners)[k][1] = 0;
		    if(nextsymbol)
		      strcat((*heavypartners)[k], tag->getAttributeValue("symbol"));
		  }
		  if(nextmodifs != NULL) {
		    modifications->insertAtEnd(nextmodifs[0]);
		    modifications->insertAtEnd(nextmodifs[1]);
		  }

		} // if static quant
		else {
		  (*lightpartners)[k] = new char[strlen(tag->getAttributeValue("aminoacid")) + 1];
		  strcpy((*lightpartners)[k], tag->getAttributeValue("aminoacid"));

		  nextmod = new char[strlen(tag->getAttributeValue("aminoacid")) + 
				     strlen(tag->getAttributeValue("mass")) + 3];
		  strcpy(nextmod, tag->getAttributeValue("aminoacid"));
		  strcat(nextmod, "[");
		  // make correction here if nec
		  double nextmass = atof(tag->getAttributeValue("mass"));
		  // adjust average masses to monoisotpic
		  if(monoisotopic_)
		    nextmass = getMonoisotopicEquivalent(nextmass, diff_, tag);
		  //		  getMonoisotopicEquivalent(nextmass, diff_, tag);
		  strcat(nextmod, tag->getAttributeValue("mass"));
		  strcat(nextmod, "]");

		  modifications->insertAtEnd(nextmod);
		}
		break; // done
	      }
	  }
	}

	else if(tag->isStart() && ! strcmp(tag->getName(), "terminal_modification")) {
	  RECORD(tag);

	  char* nextmod = NULL;
	  if(! strcmp(tag->getAttributeValue("variable"), "Y")) {

	    if(static_quant_ && strstr(paired_labels_, tag->getAttributeValue("terminus")) != NULL) {
	      cout << "error: cannot have variable modified " << tag->getAttributeValue("terminus") << " in static mode" << endl;
	      exit(1);
	    }

	    nextmod = new char[strlen(tag->getAttributeValue("terminus")) + 
			       strlen(tag->getAttributeValue("symbol")) +
			       strlen(tag->getAttributeValue("mass")) + 3];

	    strcpy(nextmod, tag->getAttributeValue("terminus"));
	    strcat(nextmod, tag->getAttributeValue("symbol"));
	    strcat(nextmod, "[");
	    double nextmass = atof(tag->getAttributeValue("mass"));
	    // adjust average masses to monoisotpic
	    if(monoisotopic_)
	      nextmass = getMonoisotopicEquivalent(nextmass, diff_, tag);
	    strcat(nextmod, tag->getAttributeValue("mass"));
	    strcat(nextmod, "]");
	    // check whether on list of paired aa's
	    for(int k = 0; k < strlen(paired_labels_); k++)
	      if(paired_labels_[k] == (tag->getAttributeValue("terminus"))[0]) { // set heavy
		(*heavypartners)[k] = new char[strlen(tag->getAttributeValue("aminoacid")) + 
					       strlen(tag->getAttributeValue("symbol")) + 1];
		strcpy((*heavypartners)[k], tag->getAttributeValue("aminoacid"));
		strcat((*heavypartners)[k], tag->getAttributeValue("symbol"));
		break; // done
	      }
	    modifications->insertAtEnd(nextmod);
	  }
	  else { // static
	    if(! strcmp(tag->getAttributeValue("terminus"), "n"))
	      static_nterm_ = True;
	    else if(! strcmp(tag->getAttributeValue("terminus"), "c"))
	      static_cterm_ = True;

	    char** nextmodifs = NULL;

	    // check whether on list of paired aa's
	    for(int k = 0; k < strlen(paired_labels_); k++)
	      if(paired_labels_[k] == (tag->getAttributeValue("terminus"))[0]) { // set heavy

		if(static_quant_) {
		  int nextstatus = getStaticQuantStatus(paired_labels_[k], atof(tag->getAttributeValue("mass")), error_);
		  if(static_status_ && nextstatus != static_status_) {
		    cout << "error: have mixture of light and heavy static modifications" << endl;
		    cout << paired_labels_[k] << " " << tag->getAttributeValue("mass") << " " << nextstatus << " " << static_status_ << endl;
		    exit(1);
		  }
		  static_status_ = nextstatus;
		  if(! static_status_)
		    static_status_ = -1; // the default is light

		  nextmodifs = getStaticModification((tag->getAttributeValue("terminus"))[0], tag->getAttributeValue("symbol"), static_status_);

		  int nextsymbol = tag->getAttributeValue("symbol") != NULL ? strlen(tag->getAttributeValue("symbol")) : 0;
		  char adjsite = paired_labels_[k];

		  if(static_status_ == -1) { // light mod
		    (*lightpartners)[k] = new char[2 + nextsymbol];
		    (*lightpartners)[k][0] = paired_labels_[k];
		    (*lightpartners)[k][1] = 0;
		    if(nextsymbol)
		      strcat((*lightpartners)[k], tag->getAttributeValue("symbol"));
		    (*heavypartners)[k] = new char[3];
		    (*heavypartners)[k][0] = paired_labels_[k];
		    (*heavypartners)[k][1] = static_symbol_;
		    (*heavypartners)[k][2] = 0;
		  }
		  else { // heavy
		    (*lightpartners)[k] = new char[3];
		    (*lightpartners)[k][0] = paired_labels_[k];
		    (*lightpartners)[k][1] = static_symbol_;
		    (*lightpartners)[k][2] = 0;
		    (*heavypartners)[k] = new char[2 + nextsymbol];
		    (*heavypartners)[k][0] = paired_labels_[k];
		    (*heavypartners)[k][1] = 0;
		    if(nextsymbol)
		      strcat((*heavypartners)[k], tag->getAttributeValue("symbol"));
		  }

		  if(nextmodifs != NULL) {
		    modifications->insertAtEnd(nextmodifs[0]);
		    modifications->insertAtEnd(nextmodifs[1]);
		  }
		} // if static quant
		else {
		  (*lightpartners)[k] = new char[strlen(tag->getAttributeValue("terminus")) + 1];
		  strcpy((*lightpartners)[k], tag->getAttributeValue("terminus"));

		  nextmod = new char[strlen(tag->getAttributeValue("terminus")) + 
				     strlen(tag->getAttributeValue("mass")) + 3];
		  strcpy(nextmod, tag->getAttributeValue("terminus"));
		  strcat(nextmod, "[");
		  // make correction here if nec
		  double nextmass = atof(tag->getAttributeValue("mass"));
		  // adjust average masses to monoisotpic
		  if(monoisotopic_)
		    nextmass = getMonoisotopicEquivalent(nextmass, diff_, tag);
		  strcat(nextmod, tag->getAttributeValue("mass"));
		  strcat(nextmod, "]");

		  modifications->insertAtEnd(nextmod);
		}
		break; // done
	      }
	  }
	}

#endif
	//else if(tag->isStart() && ! strcmp(tag->getName(), "search_summary")) {
	//  monoisotopic_ = strcmp(tag->getAttributeValue("precursor_mass_type"), "monoisotopic") == 0 ? True : False ;
	//}
	else if(tag->isEnd() && ! strcmp(tag->getName(), "search_summary")) {
	
#ifdef USE_STD_MODS
	  // check label tags here
	  int k;
	  for(k = 0; paired_labels_[k]; k++) {
	    if(paired_labels_[k] == 'n') {
	      if(heavy_nterm_mass_ > 0.0 && light_nterm_mass_ == 0.0)
		light_nterm_mass_ = ResidueMass::getMass('n', monoisotopic_);
	      else if(heavy_nterm_mass_ == 0.0 && light_nterm_mass_ > 0.0) { //must switch them around
		heavy_nterm_mass_ = light_nterm_mass_;
		light_nterm_mass_ = ResidueMass::getMass('n', monoisotopic_);
	      }
	      else if(heavy_nterm_mass_ == 0.0 && light_nterm_mass_ == 0.0) {
		// error
		cout << "nterm: " << light_nterm_mass_ << " vs " << heavy_nterm_mass_ << endl;
		exit(1);
	      }
	    } // n case
	    else if(paired_labels_[k] == 'c') {
	      if(heavy_cterm_mass_ > 0.0 && light_cterm_mass_ == 0.0)
		light_cterm_mass_ = ResidueMass::getMass('c', monoisotopic_);
	      else if(heavy_cterm_mass_ == 0.0 && light_cterm_mass_ > 0.0) { //must switch them around
		heavy_cterm_mass_ = light_cterm_mass_;
		light_cterm_mass_ = ResidueMass::getMass('c', monoisotopic_);
	      }
	      else if(heavy_cterm_mass_ == 0.0 && light_cterm_mass_ == 0.0) {
		// error
		cout << "cterm: " << light_cterm_mass_ << " vs " << heavy_cterm_mass_ << endl;
		exit(1);
	      }
	    } // n case
	    else {
	      char next_res = paired_labels_[k];
	      if(heavy_label_masses_[next_res-'A'] > 0.0 && light_label_masses_[next_res-'A'] == 0.0)
		light_label_masses_[next_res-'A'] = ResidueMass::getMass(next_res, monoisotopic_);
	      else if(heavy_label_masses_[next_res-'A'] == 0.0 && light_label_masses_[next_res-'A'] > 0) { 
		heavy_label_masses_[next_res-'A'] = light_label_masses_[next_res-'A'];
		light_label_masses_[next_res-'A'] = ResidueMass::getMass(next_res, monoisotopic_);
	      }
	      else if(heavy_label_masses_[next_res-'A'] == 0.0 && light_label_masses_[next_res-'A'] == 0) { 
		// error
		cout << "label " << next_res << ": " << light_label_masses_[next_res-'A'] << " vs " << heavy_label_masses_[next_res-'A'] << endl;
	      }
	    } // aa case
	  } // next label
	  if(static_quant_) {
	    cout << "static ";
	    if(static_status_ == -1)
	      cout << "light";
	    else if(static_status_ == 1)
	      cout << "heavy";
	    else {
	      static_status_ = -1;
	      cout << "unknown (assuming light)";
	    }
	    cout << " label for file " << mzXMLfile_ << endl;
	  }
#endif
#ifndef USE_STD_MODS
	  // time to set prtnAAs_ (using paired_labels_) and modAAs_
	  if(static_quant_) {
	    cout << "static ";
	    if(static_status_ == -1)
	      cout << "light";
	    else if(static_status_ == 1)
	      cout << "heavy";
	    else
	      cout << "unknown";
	    cout << " label for file " << mzXMLfile_ << endl;
	    setStaticPartners(lightpartners, heavypartners);
	  }

	  for(int k = 0; k < modifications->length(); k++) 
	    cout << (k+1) << ": " << (*modifications)[k] << endl;

	  for(int k = 0; k < lightpartners->length(); k++)
	    if((*lightpartners)[k] != NULL && (*heavypartners)[k] != NULL)
	      cout << (k+1) << ": " << (*lightpartners)[k] << " vs " << (*heavypartners)[k] << endl;

	  modAAs_ = collectModAAStrct(&modAANum_, modifications);  //new residueStrct(); // pass modifications
	  prtnAAs_ = collectPrtnAAStrct(&prtnAANum_, lightpartners, heavypartners); //new pairStrct(); // pass lightpartners and heavypartners

	  // print out results here....
	  cout << "paired info" << endl;
	  if(prtnAAs_ != NULL)
	    for(int k = 0; k < prtnAANum_; k++)
	      cout << prtnAAs_[k].prtnA << "<->" << prtnAAs_[k].prtnB << endl;
	  cout << "residue info" << endl;
	  if(modAAs_ != NULL)
	    for(int k = 0; k < modAANum_; k++)
	      cout << modAAs_[k].sz << " " << modAAs_[k].rp << " " << modAAs_[k].ms << endl;

	  setLabelStrings(lightpartners, lightstring, heavypartners, heavystring);
#endif
	  RECORD(tag);
	  // now write asaptimestamp

#ifdef USE_STD_MODS
	  char text[100];
	  std::string label_masses; // static buffer has overflow risk 
	  for(k = 0; paired_labels_[k]; k++) {
	    label_masses += paired_labels_[k];
	    label_masses += '[';
	    if(paired_labels_[k] == 'n') {
	      sprintf(text, "%0.3f,%0.3f]", light_nterm_mass_, heavy_nterm_mass_);
	    }
	    else if(paired_labels_[k] == 'c') {
	      sprintf(text, "%0.3f,%0.3f]", light_cterm_mass_, heavy_cterm_mass_);
	    }
	    else {
	      sprintf(text, "%0.3f,%0.3f]", light_label_masses_[paired_labels_[k]-'A'], heavy_label_masses_[paired_labels_[k]-'A']);
	    }
	    label_masses += text;
	  } // next label
	  timestamp->setAttributeValue("quant_label_masses", label_masses.c_str());  
#endif

	  if(static_quant_) {
	    if(static_status_ == -1)
	      timestamp->setAttributeValue("static_quant_label", "light");
	    else if(static_status_ == 1)
	      timestamp->setAttributeValue("static_quant_label", "heavy");
	  }
	  RECORD(timestamp_start);
	  RECORD(timestamp);
	  RECORD(timestamp_stop);
	}	
	else if(filter_) {
	  if(tag->isStart() && ! strcmp("spectrum_query", tag->getName())) {
	    strncpy(pInput_.szSpectrumName,tag->getAttributeValue("spectrum")?tag->getAttributeValue("spectrum"):"",
		    sizeof(pInput_.szSpectrumName)); // useful for deriving mzxml filename if need be
	    pInput_.iFirstScan = atoi(tag->getAttributeValue("start_scan"));
	    pInput_.iLastScan = atoi(tag->getAttributeValue("end_scan"));
	    pInput_.iChargeState = atoi(tag->getAttributeValue("assumed_charge"));
	    pInput_.dPeptideMass = (double)(atof(tag->getAttributeValue("precursor_neutral_mass")) + 1.008);
	  }
	  else if(tag->isStart() && ! strcmp("search_hit", tag->getName()) && ! strcmp("1", tag->getAttributeValue("hit_rank"))) {
	    strcpy(pInput_.szPeptide, tag->getAttributeValue("peptide"));
	    pInput_.dCalcPeptideMass = atof(tag->getAttributeValue("calc_neutral_pep_mass")) + 1.008;
	    first_hit = True;
	    prob = -100;
	    iprob = -100;
	  }
	  else if(! strcmp("search_hit", tag->getName()) ) {
	    first_hit = False;
	  }
	  else if ( tag->isStart() && ! strcmp("peptideprophet_result", tag->getName()) ) {
	    prob = atof(tag->getAttributeValue("probability"));
	  }
	  else if ( tag->isStart() && ! strcmp("interprophet_result", tag->getName()) ) {
	    iprob = atof(tag->getAttributeValue("probability"));
	  }

#ifdef USE_STD_MODS
	  else if(first_hit && tag->isStart() && ! strcmp("modification_info", tag->getName())) {
	    if(modification_tags == NULL)
	      modification_tags = new Array<Tag*>;
	    modification_tags->insertAtEnd(tag);
	    mod_on = !tag->isEnd();
	    if (!mod_on) { // tag already closed, process it now
	      modinfo_ = new ModificationInfo(modification_tags);
	    }
	  }
	  else if(first_hit && mod_on && tag->isEnd() && ! strcmp("modification_info", tag->getName())) {
	    modification_tags->insertAtEnd(tag);
	    modinfo_ = new ModificationInfo(modification_tags);
	    mod_on = False;
	  }
	  else if(first_hit && mod_on) {
	    modification_tags->insertAtEnd(tag);
	  }
#endif
	  if(tags == NULL)
	    tags = new Array<Tag*>;
	  tags->insertAtEnd(tag);
	  collected = True;
	}
	else {
#ifdef USE_STD_MODS
	  if(tag->isEnd() && ! strcmp(tag->getName(), "msms_run_summary") && ! static_quant_) {
	    // reset label masses
	    for(int k = 0; pInput_.szXpressResidues[k]; k++) {
	      if(pInput_.szXpressResidues[k] == 'n') {
		light_nterm_mass_ = 0.0;
		heavy_nterm_mass_ = 0.0;
	      }
	      else if(pInput_.szXpressResidues[k] == 'c') {
		light_cterm_mass_ = 0.0;
		heavy_cterm_mass_ = 0.0;
	      }
	      else {
		char next_res = pInput_.szXpressResidues[k];
		light_label_masses_[next_res-'A'] = 0.0;
		heavy_label_masses_[next_res-'A'] = 0.0;
	      }
	    } // next label
	  }
#endif
	  if(tag != NULL) {
	    RECORD(tag);
	  }
	}

	if(filter_memory_) { // process
	  if(tags != NULL) {
	    if((pInput_.iChargeState >= 0) && (pInput_.iChargeState <= _ASAPRATIO_MXQ_)){
	      // set the input
	      if ((prob < -10 || prob > pInput_.dMinPprob) && (iprob < -10 || iprob > pInput_.dMinIprob)) {
		getRatio();
	      }
	      pepIdx_++;
	      for(int k = 0; k < tags->length(); k++)
		if((*tags)[k] != NULL) {
		  if(! asap_filter->filter((*tags)[k])) {

		    // here check for correct time to write asapratio tag
		    if((prob < -10 || prob > pInput_.dMinPprob) && (iprob < -10 || iprob > pInput_.dMinIprob) && first_hit && (*tags)[k]->isEnd() && ! strcmp((*tags)[k]->getName(), "search_hit") ) {
		      asap_tags = generateXML(asap_index, True, True);
		      if(asap_tags != NULL) {
			if(asap_tags->length() > 0) {
			  ratio_tags_written++;
			  if (ratio_tags_written % 1000 == 0)
			    cout << " " << ratio_tags_written/1000 << "k" << endl;
			  else if (ratio_tags_written % 100 == 0)
			    cout << ":";
			  else if (ratio_tags_written % 10 == 0) {
			    cout << ".";
			    fflush(stdout);
			  }

			  RECORD(result_start);
			  for(int j = 0;  j < asap_tags->length(); j++) {
			    if((*asap_tags)[j] != NULL) {
			      RECORD((*asap_tags)[j]);
#ifndef DEBUG_NOFREE
			      delete (*asap_tags)[j];
#endif
			    } // if
			  } // next
			  RECORD(result_stop);
			} // have asaptags
#ifndef DEBUG_NOFREE
			delete asap_tags;
#endif
		      }
		    }
		    else if ((*tags)[k]->isStart() && ! strcmp((*tags)[k]->getName(), "search_hit")  && ! strcmp("1", (*tags)[k]->getAttributeValue("hit_rank"))) {
		      first_hit = True;
		    }
		    else if ((*tags)[k]->isStart() && ! strcmp((*tags)[k]->getName(), "search_hit")) {
		      first_hit = False;
		    }
		    RECORD((*tags)[k]);
		  }
#ifndef DEBUG_NOFREE
		  delete (*tags)[k];
#endif
		}
#ifndef DEBUG_NOFREE
	      delete tags;
#endif
	      tags = NULL;
	    }
#ifdef USE_STD_MODS
#ifndef DEBUG_NOFREE
	    if(modification_tags != NULL) {
	      delete modification_tags;
	      modification_tags = NULL;
	    }
#endif
#endif
	    pInput_.iChargeState = -1; // reset
#ifdef USE_STD_MODS
#ifndef DEBUG_NOFREE
	    if(modinfo_ != NULL)
	      delete modinfo_;
	    modinfo_ = NULL;
#endif
#endif

	  }
	}
      } // if not filtered
#ifndef DEBUG_NOFREE
      if(! collected && tag != NULL)
	delete tag;
#endif
      data = strstr(data+1, "<");
    } // next tag
  } // next line
  fin.close();
  fout.close();
  freeXmlIndx(); //DDS
  
  cout << " Total: " << ratio_tags_written << endl;

  if(! overwrite(xmlfile, outfile.c_str(), "</msms_pipeline_analysis>")) {
    cerr << "error: no ASAPRatio data written to file " << xmlfile << endl;
  }
  
  if (testType!=NO_TEST) {
    //
    // regression test stuff - bpratt Insilicos LLC, Nov 2005
    //
    ASAPRatioPeptideParserTagListComparator("ASAPRatioPeptideParser",testType,test_tags,testFileName);
    delete[] testFileName;
    for(int k = test_tags.length(); k--;) {
      delete test_tags[k];
    }
  }

  delete timestamp_start;
  delete timestamp_stop;
  delete timestamp;
  delete result_start;
  delete result_stop;
  delete asap_filter;
  delete asap_summ_filter;
  delete asap_time_filter;
  delete[] nextline;
}


double ASAPRatioPeptideParser::getUserSpecifiedEquivalent(const char & aa, double avgmass) {
  double rtn = -1;
  double diff = -1;
  if(strlen(pInput_.labelMasses) > 1) {
    int k = 0;
    int start = 0;
    char nextmass[100];
    while(k < (int) strlen(pInput_.labelMasses)) {
      char next_res = pInput_.labelMasses[start];
      k++;

      if (next_res == aa) {
	while(k < (int) strlen(pInput_.labelMasses) &&
	      (pInput_.labelMasses[k] < 'a' || pInput_.labelMasses[k] > 'z') && 
	      (pInput_.labelMasses[k] < 'A' || pInput_.labelMasses[k] > 'Z')) {
	  k++;
	}

	strncpy(nextmass, pInput_.labelMasses + start + 1, k - start - 1);
	nextmass[k - start - 1] = 0;
	rtn = atof(nextmass);
	if (diff < 0 || fabs(rtn - avgmass) < diff) {
	  diff = fabs(rtn - avgmass);
	}
      }
      else {
	while(k < (int) strlen(pInput_.labelMasses) &&
	      (pInput_.labelMasses[k] < 'a' || pInput_.labelMasses[k] > 'z') && 
	      (pInput_.labelMasses[k] < 'A' || pInput_.labelMasses[k] > 'Z')) {
	  k++;
	}
	start = k;
      }
    }
  }
  if (rtn < 0) {
    rtn = avgmass;
  }
  return rtn;
}


double ASAPRatioPeptideParser::getMonoisotopicEquivalent(double averagemass, double error, Tag* tag) {
  char text[100];
  //Use user mass first
  if(strlen(pInput_.labelMasses) > 1) { // must enter user specified label masses   
    char nextmass[100];
    int index = 1;
    int start = 0;

    while(index < (int) strlen(pInput_.labelMasses)) {
      while(index < (int) strlen(pInput_.labelMasses) &&
	    (pInput_.labelMasses[index] < 'a' || pInput_.labelMasses[index] > 'z') && 
	    (pInput_.labelMasses[index] < 'A' || pInput_.labelMasses[index] > 'Z')) 
	index++;
      // process
      strncpy(nextmass, pInput_.labelMasses + start + 1, index - start - 1);
      nextmass[index - start - 1] = 0;
      double nmass = atof(nextmass);
      if (averagemass >= nmass - error && averagemass <= nmass + error) {
	//cerr << "INFO: Adjusting Average Mass " << averagemass << " to Mono Mass " <<  nmass << endl;
	return nmass;
      }
      start = index;
      index++;
    }
  }

  for(int z = 0; z < average_mods_->length(); z++) {
    if(averagemass >= (*average_mods_)[z] - error && averagemass <= (*average_mods_)[z] + error) { // make change
      cout << "adjusting " << averagemass << " to " <<  (*monoisostopic_mods_)[z] << endl;
      if(tag != NULL) {
	sprintf(text, "%0.4f", (*monoisostopic_mods_)[z]);
	tag->setAttributeValue("mass", text);
      }
      return (*monoisostopic_mods_)[z];
    }
  }

  return averagemass;
}

// returns light and heavy
char** ASAPRatioPeptideParser::getStaticModification(char site, char* symbol, int status) {
  if(! status)
    return NULL;
//  char** output = new char* [2];
  char text[100];
  int symbol_len = symbol != NULL ? (int)strlen(symbol) : 0;

  for(int k = 0; k < static_pairs_->length(); k++) {
    if((*static_pairs_)[k]->site == site) {
      char** output = new char* [2];
      if(status == -1) { // light

	// put in corrections for mono masses here...
	sprintf(text, "%0.5f", getMonoisotopicEquivalent((*static_pairs_)[k]->lightmass, diff_, NULL));
	output[0] = new char[strlen(text) + 1 + symbol_len + 3];
	output[0][0] = site;
	output[0][1] = 0;
	if(symbol_len)
	  strcat(output[0], symbol);
	strcat(output[0], "[");
	strcat(output[0], text);
	strcat(output[0], "]");

	sprintf(text, "%0.5f", getMonoisotopicEquivalent((*static_pairs_)[k]->heavymass, diff_, NULL));
	output[1] = new char[strlen(text) + 1 + 1 + 3];
	output[1][0] = site;
	output[1][1] = static_symbol_;
	output[1][2] = 0;
	strcat(output[1], "[");
	strcat(output[1], text);
	strcat(output[1], "]");
      }
      else { // heavy
	sprintf(text, "%0.5f", getMonoisotopicEquivalent((*static_pairs_)[k]->lightmass, diff_, NULL));
	output[0] = new char[strlen(text) + 1 + 1 + 3];
	output[0][0] = site;
	output[0][1] = static_symbol_;
	output[0][2] = 0;
	strcat(output[0], "[");
	strcat(output[0], text);
	strcat(output[0], "]");

	sprintf(text, "%0.5f", getMonoisotopicEquivalent((*static_pairs_)[k]->heavymass, diff_, NULL));
	output[1] = new char[strlen(text) + 1 + symbol_len  + 3];
	output[1][0] = site;
	output[1][1] = 0;
	if(symbol_len)
	  strcat(output[1], symbol);
	strcat(output[1], "[");
	strcat(output[1], text);
	strcat(output[1], "]");

      }
      return output;
    }
  }
  return NULL;
}


void ASAPRatioPeptideParser::setStaticQuantInfo(const char* xmlfile) {

  //  int line_width = 10000;
  char *nextline = new char[line_width_];
  char* data = NULL;
  Tag* tag = NULL;
  int k;

#ifdef USE_STD_MODS
  for(k = 0; k < 26; k++) {
    light_label_masses_[k] = 0.0;
    heavy_label_masses_[k] = 0.0;
  }
  light_nterm_mass_ = 0.0;
  heavy_nterm_mass_ = 0.0;
  light_cterm_mass_ = 0.0;
  heavy_cterm_mass_ = 0.0;
#endif
  RACI fin(xmlfile); // can read gzipped xml
  if(! fin) {
    cerr << "error opening " << xmlfile << endl;
    exit(1);
  }
  while(fin.getline(nextline, line_width_)) {
    if(strstr(nextline, "modification") != NULL) {
      data = strstr(nextline, "<");
      while(data != NULL) {
	tag = new Tag(data);
	if(tag != NULL) {
	  if(tag->isStart() && ! strcmp(tag->getName(), "aminoacid_modification") &&
	     ! strcmp(tag->getAttributeValue("variable"), "N")) {
	    enterStaticQuant((tag->getAttributeValue("aminoacid"))[0], atof(tag->getAttributeValue("mass")), error_);
	  }
	  else if(tag->isStart() && ! strcmp(tag->getName(), "terminal_modification") &&
		  ! strcmp(tag->getAttributeValue("variable"), "N")) {
	    enterStaticQuant((tag->getAttributeValue("terminus"))[0], atof(tag->getAttributeValue("mass")), error_);
	  }
	}
#ifndef DEBUG
	delete tag;
#endif
	data = strstr(data+1, "<");

      }

    } // if have modification info
  } // next line

  fin.close();

  for(k = 0; k < static_pairs_->length(); k++) {
    cout << (*static_pairs_)[k]->site << " " << (*static_pairs_)[k]->lightmass << " " << (*static_pairs_)[k]->heavymass << endl;
#ifdef USE_STD_MODS
    if((*static_pairs_)[k]->site == 'n') {
      light_nterm_mass_ = (*static_pairs_)[k]->lightmass;
      heavy_nterm_mass_ = (*static_pairs_)[k]->heavymass;
    }
    else if((*static_pairs_)[k]->site == 'c') {
      light_cterm_mass_ = (*static_pairs_)[k]->lightmass;
      heavy_cterm_mass_ = (*static_pairs_)[k]->heavymass;
    }
    else { // aa
      light_label_masses_[(*static_pairs_)[k]->site - 'A'] = (*static_pairs_)[k]->lightmass;
      heavy_label_masses_[(*static_pairs_)[k]->site - 'A'] = (*static_pairs_)[k]->heavymass;
    }
#endif

  }
  //exit(1);
  delete [] nextline;
}


Tag* ASAPRatioPeptideParser::getSummaryTag(const InputStruct& opts) {
  Tag* output = new Tag("asapratio_summary", True, True);
  char text[1000];
  snprintf(text,sizeof(text),"%s (%s)",PROGRAM_VERSION,szTPPVersionInfo);
  output->setAttributeValue("version", text);
  output->setAttributeValue("author", PROGRAM_AUTHOR);
  
  if(opts.bUseSameScanRange)
    output->setAttributeValue("elution", "0");
  else
    output->setAttributeValue("elution", "1");


  // use paired_lables_ after alphabetizing....
  //  output->setAttributeValue("labeled_residues", opts.szXpressResidues);  

  output->setAttributeValue("labeled_residues", paired_labels_);  


  sprintf(text, "%d", opts.bXpressLight1);
  output->setAttributeValue("area_flag", text);
  if(static_quant_)
    output->setAttributeValue("static_quant", "Y");
  else
    output->setAttributeValue("static_quant", "N");
  if(strlen(pInput_.labelMasses) > 1) 
    output->setAttributeValue("specified_residue_masses", pInput_.labelMasses);
 

  return output;
}

void ASAPRatioPeptideParser::setFilter(Tag* tag) {
  if(tag == NULL)
    return;

  if(filter_memory_) {
    filter_memory_ = False;
    filter_ = False;
  }

  if(! strcmp(tag->getName(), "spectrum_query")){
    if(tag->isStart()) {
      filter_ = True;
    }else{
      filter_memory_ = True;
    }
  }
}

// only have to deal with null ones
void ASAPRatioPeptideParser::setStaticPartners(Array<char*>* lightpartners, Array<char*>* heavypartners) {
  if(! static_quant_ || ! static_status_)
    return;
  for(int k = 0; k < static_pairs_->length(); k++) {
    for(int j = 0; paired_labels_[j]; j++) {
      if(((*lightpartners)[j] == NULL  || (*heavypartners)[j] == NULL) && (*static_pairs_)[k]->site == paired_labels_[j]) {
	if(static_status_ == -1) { // light mod
	  (*lightpartners)[j] = new char[2];
	  (*lightpartners)[j][0] = paired_labels_[j];
	  (*lightpartners)[j][1] = 0;
	  (*heavypartners)[j] = new char[3];
	  (*heavypartners)[j][0] = paired_labels_[j];
	  (*heavypartners)[j][1] = static_symbol_;
	  (*heavypartners)[j][2] = 0;
	}
	else { // heavy
	  (*lightpartners)[j] = new char[3];
	  (*lightpartners)[j][0] = paired_labels_[j];
	  (*lightpartners)[j][1] = static_symbol_;
	  (*lightpartners)[j][2] = 0;
	  (*heavypartners)[j] = new char[2];
	  (*heavypartners)[j][0] = paired_labels_[j];
	  (*heavypartners)[j][1] = 0;
	}

      } // have a match

    } // next paired quant label

  } // next static pair identified
}

// -1 for light, 0 for not applicable, 1 for heavy
int ASAPRatioPeptideParser::getStaticQuantStatus(char site, double mass, double error) {
  if(strchr(paired_labels_, site) == NULL)
    return 0;
  for(int k = 0; k < static_pairs_->length(); k++) {
    if((*static_pairs_)[k]->site == site) {
      if(! diff((*static_pairs_)[k]->lightmass, mass, error))
	return -1;
      else if(! diff((*static_pairs_)[k]->heavymass, mass, error))
	return 1;
      else return 0;
    }
  }
  return 0;
}



// -1 for less, 0 for equal, 1 for greater
int ASAPRatioPeptideParser::diff(double first, double second, double error) {
  if(first == second)
    return 0;
  if(first > second) {
    if(first > second + error)
      return 1;
    return 0;
  }
  if(second > first + error)
    return -1;
  return 0;
}


void ASAPRatioPeptideParser::enterStaticQuant(char site, double mass, double error) {
  if(strchr(paired_labels_, site) == NULL)
    return; // nothing to do

  for(int k = 0; k < static_pairs_->length(); k++){
    if((*static_pairs_)[k]->site == site){
      if(! diff((*static_pairs_)[k]->lightmass, (*static_pairs_)[k]->heavymass, error)) { // check whether new mass is lighter or heavier
	int difference = diff(mass, (*static_pairs_)[k]->lightmass, error);
	if(! difference) // same as last one
	  return;
	if(difference == 1) {
	  (*static_pairs_)[k]->heavymass = mass;
	}else{
	  (*static_pairs_)[k]->lightmass = mass;
        }
	return;
       
      } // only one value seen thus far
      else { // check whether current is third
	if(diff((*static_pairs_)[k]->lightmass, mass, error) && diff((*static_pairs_)[k]->heavymass, mass, error)) {
	  cout << "error: have 3 values for " << site << ": " << (*static_pairs_)[k]->lightmass << ", " << (*static_pairs_)[k]->heavymass << ", " << mass << endl;
	  exit(1);
	}
	return;
      }
    }
  }

  // still here, add new one
  staticQuant* nextQ = new staticQuant();
  nextQ->site = site;
  nextQ->lightmass = mass;
  nextQ->heavymass = mass;
  static_pairs_->insertAtEnd(nextQ);

}

void ASAPRatioPeptideParser::setLabelStrings(Array<char*>* lights, char* lightstring, Array<char*>* heavys, char* heavystring) {
  lightstring[0] = 0;
  heavystring[0] = 0;
  Boolean first = True;
  for(int k = 0; k < lights->length(); k++)
    if((*lights)[k] != NULL && (*heavys)[k] != NULL) {
      if(first) {
	first = False;
      }
      else {
	strcat(lightstring, ",");
	strcat(heavystring, ",");
      }
      strcat(lightstring, (*lights)[k]);
      strcat(heavystring, (*heavys)[k]);
    } // if have both labels
}

// This function frees a matrix.
void ASAPRatioPeptideParser::freeMtrx(void **mtrx, int size)
{
#ifndef DEBUG_NOFREE
  int i;
  
  for (i = 0; i < size; ++i)
    free(mtrx[i]);
  free(mtrx);
#endif
  return;
}

// This function converts "char *string" to "residueStrct AA". aa->sz = -1 if invalid input.
void ASAPRatioPeptideParser::string2AA(char *string, residueStrct *aa)
{
  char **sects;
  int sectNum;
  char *tmpStr;
  int lngth;
  char rp[3];
  double mass;
  char test;
  int indx;

  int i;

  // check input
  getRidOfSpace(string);
  // andy here
  //cnvtUpper(string);
  lngth = (int)strlen(string);
  if(lngth < 1) {
    printf("Invalid input \"%s\" for modified amino acid.\n", string);
    aa->sz = -1;
    return;
  }

  // size and representation
  sects = getStrSects(&sectNum, string, '[');
  if(sectNum != 2 
     || (lngth = (int)strlen(sects[0])) > 2
     || (lngth == 2 
	 && (isalpha(sects[0][1]) != 0
	     || sects[0][1] == ']'
	     || sects[0][1] == '['))) {
    printf("Invalid input \"%s\" for modified amino acid.\n", string);
    aa->sz = -1;
    return;
  }
  aa->sz = lngth;
  strcpy(aa->rp, sects[0]);

  // mass

  // text for mass 
  if((tmpStr = strchr(sects[1], ']')) == NULL
     || (lngth = (int)strlen(sects[1])-(int)strlen(tmpStr)) < 1) {
    printf("Invalid input \"%s\" for modified amino acid.\n", string);
    aa->sz = -1;
    return;
  }
  else {
    sects[1][lngth] = '\0';
    getRidOfSpace(sects[1]);
    lngth = (int)strlen(sects[1]);
    tmpStr = (char *) calloc(lngth+1, sizeof(char));
    strcpy(tmpStr, sects[1]);
    freeMtrx((void **)sects, sectNum);
  }

  // format
  if(strchr(tmpStr, '+') != NULL) { // format [R+m]
    sects = getStrSects(&sectNum, tmpStr, '+');
    if(sectNum != 2 
       || (lngth = (int)strlen(sects[0])) > 2
       || sscanf(sects[1], "%lf%c", &mass, &test) != 1) {
      printf("Invalid input \"%s\" for modified amino acid.\n", string);
      aa->sz = -1;
#ifndef DEBUG_NOFREE
      free(tmpStr);
      freeMtrx((void **)sects, sectNum);
#endif
      return;
    }
    else {
      strcpy(rp, sects[0]);
      aa->ms = mass;
#ifndef DEBUG_NOFREE
      free(tmpStr);
#endif
      freeMtrx((void **)sects, sectNum);
    }
  } // if(strchr(tmpStr, '+') != NULL) { // format [R+m]
  else if(strchr(tmpStr, '-') != NULL) { // format [R-m]
    sects = getStrSects(&sectNum, tmpStr, '-');
    if(sectNum != 2 
       || (lngth = (int)strlen(sects[0])) > 2
       || sscanf(sects[1], "%lf%c", &mass, &test) != 1) {
      printf("Invalid input \"%s\" for modified amino acid.\n", string);
      aa->sz = -1;
#ifndef DEBUG_NOFREE
      free(tmpStr);
      freeMtrx((void **)sects, sectNum);
#endif
      return;
    }
    else {
      strcpy(rp, sects[0]);
      aa->ms = -mass;
#ifndef DEBUG_NOFREE
      free(tmpStr);
#endif
      freeMtrx((void **)sects, sectNum);
    }
  } // else if(strchr(tmpStr, '-') != NULL) { // format [R-m]
  else if (isalpha(tmpStr[0]) != 0) { // format [R] 
    if(lngth > 2) {
      printf("Invalid input \"%s\" for modified amino acid.\n", string);
      aa->sz = -1;

#ifndef DEBUG_NOFREE
      free(tmpStr);
#endif
      return;      
    }
    else {
      aa->ms = 0.;
      strcpy(rp, tmpStr);
#ifndef DEBUG_NOFREE
      free(tmpStr);
#endif
    }
  } // else if (isalpha(tmpStr[0]) != 0) { // format [R] 
  else { // format [m] 
    if(sscanf(tmpStr, "%lf%c", &mass, &test) != 1 
       || mass < 0.){
      printf("Invalid input \"%s\" for modified amino acid.\n", string);
      aa->sz = -1;
#ifndef DEBUG_NOFREE
      free(tmpStr);
#endif
      return;      
    }
    else {
      aa->ms = mass;
#ifndef DEBUG_NOFREE
      free(tmpStr);
#endif
      return;      
    }
  } // else { // format [m] 

  // mass of native amino acid
  indx = -1;
  for (i = 0; i < NATIVEAANUM; ++i) {
    if(strcmp(nativeAA[i].rp, rp) == 0) {
      aa->ms += nativeAA[i].ms;
      indx = i;
      break;
    }
  }
  if(indx < 0 || aa->ms < 0.) {
    printf("Invalid input \"%s\" for modified amino acid.\n", string);
    aa->sz = -1;      
  }

  return;
}


// This function gets consecutive sections of a string, separated by "sep".
char **ASAPRatioPeptideParser::getStrSects(int *sectNum, char *string, char sep)
{
  int lngth = (int)strlen(string);
  char **sects;
  int startIndx, endIndx;
  int i;

  // sectNum
  *sectNum = 1;
  for (i = 0; i < lngth; ++i) {
    if (string[i] == sep) {
      ++(*sectNum);
    }
  }

  // sects
  sects = (char **) calloc(*sectNum, sizeof(char *));
  *sectNum = 0;
  startIndx = 0;
  for (i = 0; i < lngth; ++i) {
    if (string[i] == sep) {
      endIndx = i;
      sects[*sectNum] = (char *) calloc(endIndx-startIndx+1, sizeof(char));
      strncpy(sects[*sectNum], string+startIndx, endIndx-startIndx);
      getRidOfSpace(sects[*sectNum]);
      ++(*sectNum);
      startIndx = endIndx + 1;
    }
  }
  endIndx = i;
  sects[*sectNum] = (char *) calloc(endIndx-startIndx+1, sizeof(char));
  strncpy(sects[*sectNum], string+startIndx, endIndx-startIndx);
  getRidOfSpace(sects[*sectNum]);
  ++(*sectNum);
  
  return sects;
}


// This function coverts a string array into a list of modified amino acids.
// input string:  "R*[mass]" or "R*[N+dmass]".
int modAA_size_compare(void const *a, void  const *b) {  
  return -((residueStrct *)a)->sz + ((residueStrct *)b)->sz;
}

residueStrct *ASAPRatioPeptideParser::collectModAAStrct(int *modAANum, Array<char*>* inputStrs)
{
  residueStrct *modAAs;
  int tmpNum;
  int i;

  modAAs = (residueStrct *) calloc(inputStrs->length(), sizeof(residueStrct));

  tmpNum = 0;
  for (i = 0; i < inputStrs->length(); ++i) {
    string2AA((*inputStrs)[i], modAAs+tmpNum);
    if(modAAs[tmpNum].sz > 0)
      ++tmpNum;
    else {
      printf("Warning!\n");
      fflush(stdout);
    }
  }
  qsort(modAAs, tmpNum, sizeof(residueStrct), modAA_size_compare);

  if(tmpNum > 0){
    *modAANum = tmpNum;
    return modAAs;
  }
  else {
#ifndef DEBUG_NOFREE
    free(modAAs);
#endif
    return NULL;
  }
}


// This function coverts two strings into a list of amino acid partners. It returns NULL on failure.
// lightString: L1, L2, L3, ...
// heavyString: H1, H2, H3, ...
// Pairs (L1, H1), (L2, H2), (L3, H3), ..., must appear in order.
int pairStrctCmp(void const *a, void const *b) {  
  int cmp = (int)strlen(((pairStrct *)b)->prtnA) - (int)strlen(((pairStrct *)a)->prtnA);
  if(cmp != 0)
    return cmp;
  else
    cmp = (int)strlen(((pairStrct *)b)->prtnB) - (int)strlen(((pairStrct *)a)->prtnB);
  return (cmp?cmp:(char *)a-(char *)b); // stabilize sort
}

pairStrct *ASAPRatioPeptideParser::collectPrtnAAStrct(int *prtnAANum, Array<char*>* lightStrings, Array<char*>* heavyStrings)
{
  pairStrct *prtnAAs;
  //char **tmpInputs;
  char tmpString[100];
  int tmpNum = 0;
  int i,k;

  *prtnAANum = 0;
  // light
  for(k = 0; k < lightStrings->length(); k++)
    if((*lightStrings)[k] != NULL) {
      //cout << "here incr" << endl;
      // andy here
      //cnvtUpper((*lightStrings)[k]);
      (*prtnAANum)++;
    }
  cout << "total: " << *prtnAANum << " partners found" << endl;
  //cnvtUpper(lightString);
  //tmpInputs = getStrSects(prtnAANum, lightString, ',');
  // replaced by lightStrings....

  if(*prtnAANum < 1) {
    printf("Invalid input for LIGHT isotopes: ");
    for(int k = 0; k < lightStrings->length(); k++)
      if((*lightStrings)[k] != NULL)
	printf("\"%s\"\n", (*lightStrings)[k]);
    fflush(stdout);
    //freeMtrx((void **)tmpInputs, *prtnAANum);
    return NULL;
  }
  prtnAAs = (pairStrct *) calloc(*prtnAANum, sizeof(pairStrct));
  for (i = 0; i < *prtnAANum; ++i) 
    if((*lightStrings)[i] != NULL)
      strcpy(prtnAAs[i].prtnA, (*lightStrings)[i]);
  //freeMtrx((void **)tmpInputs, *prtnAANum);

  // heavy
  for(k = 0; k < heavyStrings->length(); k++)
    if((*heavyStrings)[k] != NULL) {
      // andy here
      //cnvtUpper((*heavyStrings)[k]);
      tmpNum++;
    }


  //cnvtUpper(heavyString);
  //tmpInputs = getStrSects(&tmpNum, heavyString, ',');
  if(*prtnAANum != tmpNum) {
    printf("Number of LIGHT isotopes (%d) doesn't match that of HEAVY (%d). \n", 
	   *prtnAANum, tmpNum);
    fflush(stdout);
#ifndef DEBUG_NOFREE
    free(prtnAAs);
#endif
    //freeMtrx((void **)tmpInputs, tmpNum);
    return NULL;
  }
  for (i = 0; i < *prtnAANum; ++i) 
    if((*heavyStrings)[i] != NULL)
      strcpy(prtnAAs[i].prtnB, (*heavyStrings)[i]);
  //freeMtrx((void **)tmpInputs, *prtnAANum);

  // sort prtnAAs
  for (i = 0; i < *prtnAANum; ++i) {
    if(strlen(prtnAAs[i].prtnA) < strlen(prtnAAs[i].prtnB)) {
      strcpy(tmpString, prtnAAs[i].prtnA);
      strcpy(prtnAAs[i].prtnA, prtnAAs[i].prtnB);
      strcpy(prtnAAs[i].prtnB, tmpString);
    }
  }
  qsort(prtnAAs, *prtnAANum, sizeof(pairStrct), pairStrctCmp);

  return prtnAAs;
}

#ifdef USE_STD_MODS
void ASAPRatioPeptideParser::evalPepDataStrct(char *pepSeq, long scan, int chrg, char *xmlFile,
					      int eltn, int areaFlag, double pepMass, double computed_peptide_mass) {


#endif
#ifndef USE_STD_MODS
void ASAPRatioPeptideParser::evalPepDataStrct(char *pepSeq, long scan, int chrg, char *xmlFile,
					      int eltn, int areaFlag, residueStrct *modAAs, int modAANum,
					      pairStrct *prtnAAs, int prtnAANum) {
#endif
#ifndef USE_STD_MODS
  char prtnSeq[1000];
#endif
  double msLight = -1, msHeavy=-1;
  int cidIndx;

  //cout << "Pep: " << pepSeq << " "; 
  
  // initialize
  //data = (pepDataStrct *)calloc(1, sizeof(pepDataStrct));
  memset(&data_, 0, sizeof(data_));;

#ifdef USE_STD_MODS

  // MUST MAKE SURE THIS CAN WORK IN STATIC MODE

  //  Boolean compute_peptide_mass = (pepMass > 0.0);

  // check modifications only to see if heavy or light or none
  Boolean light = False;
  Boolean heavy = False;
  Boolean unmod = False; // could be counted as light

  // calculate if correct, and mass diff between heavy and light
  double massdiff = 0.0;
  double peptide_mass = 0.0; // calculate this as well
  peptide_mass += _ASAPRATIO_HM_; //Store the M+H
  // MUST CHECK FOR N AND C TERMINAL MODS HERE FIRST............
  if(modinfo_ != NULL) {
    if(strchr(paired_labels_, 'n') != NULL) {

      if(modinfo_->getNtermModMass() > 0.0) {

	double nextmass = modinfo_->getNtermModMass();
	if(monoisotopic_)
	  nextmass = getMonoisotopicEquivalent(nextmass, diff_, NULL);

	if(nextmass - light_nterm_mass_ <= MOD_ERROR &&
	   light_nterm_mass_ - nextmass <= MOD_ERROR) 
	  light = True;
	else if(nextmass - heavy_nterm_mass_ <= MOD_ERROR &&
	   heavy_nterm_mass_ - nextmass <= MOD_ERROR)
	  heavy = True;
	if(compute_peptide_mass_)
	  peptide_mass += nextmass;
      } // have modified terminus
      else {
	if(fabs(light_nterm_mass_ - ResidueMass::getMass('n', monoisotopic_)) < 0.01)
	  light = True;
	else unmod = True; // illegal
	if(compute_peptide_mass_)
	  peptide_mass += ResidueMass::getMass('n', monoisotopic_);
      }
      massdiff += heavy_nterm_mass_ - light_nterm_mass_;
    } // n
    else if(compute_peptide_mass_) { // whether modified anyway
      if(modinfo_->getNtermModMass() > 0.0) {
	double nextmass = modinfo_->getNtermModMass();
	if(monoisotopic_)
	  nextmass = getUserSpecifiedEquivalent('n', nextmass);
	peptide_mass += nextmass;

      }
      else
	peptide_mass += ResidueMass::getMass('n', monoisotopic_);
    }
    if(strchr(paired_labels_, 'c') != NULL) {

      if(modinfo_->getCtermModMass() > 0.0) {
	double nextmass = modinfo_->getCtermModMass();
	if(monoisotopic_)
	  nextmass = getMonoisotopicEquivalent(nextmass, diff_, NULL);
	if(nextmass - light_cterm_mass_ <= MOD_ERROR &&
	   light_cterm_mass_ - nextmass <= MOD_ERROR) 
	  light = True;
	else if(nextmass - heavy_cterm_mass_ <= MOD_ERROR &&
	   heavy_cterm_mass_ - nextmass <= MOD_ERROR)
	  heavy = True;
	if(compute_peptide_mass_)
	  peptide_mass += nextmass;
      } // have modified terminus
      else {
	if(fabs(light_cterm_mass_ - ResidueMass::getMass('c', monoisotopic_)) < 0.01)
	  light = True;
	else unmod = True; // illegal
	if(compute_peptide_mass_)
	  peptide_mass += ResidueMass::getMass('c', monoisotopic_);
      }
      massdiff += heavy_cterm_mass_ - light_cterm_mass_;
    } // c
    else if(compute_peptide_mass_) {
      if(modinfo_->getCtermModMass() > 0.0) {
	double nextmass = modinfo_->getCtermModMass();
	if(monoisotopic_)
	  nextmass = getUserSpecifiedEquivalent('c', nextmass);
	peptide_mass += nextmass;
      }
      else
	peptide_mass += ResidueMass::getMass('c', monoisotopic_);
    }
  } // some mods to look at
  else if(compute_peptide_mass_) { // no mods
      peptide_mass += ResidueMass::getMass('n', monoisotopic_) + ResidueMass::getMass('c', monoisotopic_);
  }

  if (compute_peptide_mass_) {
    for(int k = 0; pepSeq[k]; k++) {
      if(strchr(paired_labels_, pepSeq[k]) != NULL) { // have a labeled aa
	//      cout << k << ": " << pepSeq[k] << endl; exit(1);
	double nextmass = modinfo_ == NULL ? 0.0 : modinfo_->getModifiedResidueMass(k);

	if(monoisotopic_)
	  nextmass = getMonoisotopicEquivalent(nextmass, diff_, NULL);

	//     cout << "next mass: " << nextmass << endl;
	if(nextmass > 0.0) { // really is modified, find out if light or heavy
	    
	  //		cout << pepSeq[k] << " masses: " << light_label_masses_[pepSeq[k]-'A'] << " vs " << heavy_label_masses_[pepSeq[k]-'A'] << endl;
	    
	  if(nextmass - light_label_masses_[pepSeq[k]-'A'] <= MOD_ERROR &&
	     light_label_masses_[pepSeq[k]-'A'] - nextmass <= MOD_ERROR) 
	    light = True;
	  else if(nextmass - heavy_label_masses_[pepSeq[k]-'A'] <= MOD_ERROR &&
		  heavy_label_masses_[pepSeq[k]-'A'] - nextmass <= MOD_ERROR)
	    heavy = True;
	  if(compute_peptide_mass_)
	    peptide_mass += nextmass;
	    
	} // have modified aa
	else {
	  if(fabs(light_label_masses_[pepSeq[k] - 'A'] - ResidueMass::getMass(pepSeq[k], monoisotopic_)) < 0.01)
	    light = True;
	  else unmod = True; // illegal
	  if(compute_peptide_mass_)
	    peptide_mass += ResidueMass::getMass(pepSeq[k], monoisotopic_);
	}
	massdiff += heavy_label_masses_[pepSeq[k]-'A'] - light_label_masses_[pepSeq[k]-'A'];

      } // modified res
      else if(compute_peptide_mass_) {
	if(modinfo_ != NULL && modinfo_->getModifiedResidueMass(k) > 0.0) {
	  double nextmass =  modinfo_->getModifiedResidueMass(k);

	  if(monoisotopic_)
	    nextmass = getUserSpecifiedEquivalent(pepSeq[k],nextmass);

	  peptide_mass += nextmass;
	}
	else
	  peptide_mass += ResidueMass::getMass(pepSeq[k], monoisotopic_);
      }
    }
  }
  else {
    int lablLen = (int)strlen(paired_labels_);
    unmod = True;
    for(int k = 0; k < lablLen; k++) {
      char * tmp = strchr(pepSeq, paired_labels_[k]);
      if ( tmp != NULL) {
	unmod = False;
	int idx = (int)(tmp-pepSeq);
	double nextmass = modinfo_ == NULL ? 0.0 : modinfo_->getModifiedResidueMass(idx);
	if(monoisotopic_)
	  nextmass = getMonoisotopicEquivalent(nextmass, diff_, NULL);
	//     cout << "next mass: " << nextmass << endl;
	if(nextmass > 0.0) { // really is modified, find out if light or heavy
	  
	  //		cout << pepSeq[k] << " masses: " << light_label_masses_[pepSeq[k]-'A'] << " vs " << heavy_label_masses_[pepSeq[k]-'A'] << endl;

	  if(nextmass - light_label_masses_[pepSeq[idx]-'A'] <= MOD_ERROR &&
	     light_label_masses_[pepSeq[idx]-'A'] - nextmass <= MOD_ERROR) 
	    light = True;
	  else if(nextmass - heavy_label_masses_[pepSeq[(int)(tmp-pepSeq)]-'A'] <= MOD_ERROR &&
		  heavy_label_masses_[pepSeq[idx]-'A'] - nextmass <= MOD_ERROR)
	    heavy = True;

	} // have modified aa
	else {
	  if(fabs(light_label_masses_[pepSeq[idx] - 'A'] - ResidueMass::getMass(pepSeq[idx], monoisotopic_)) < 0.01)
	    light = True;
	  else unmod = True; // illegal
	}
	massdiff += heavy_label_masses_[pepSeq[idx]-'A'] - light_label_masses_[pepSeq[idx]-'A'];

      } // modified res
    }
    peptide_mass = computed_peptide_mass;
  }


  //  pInput_.dPeptideMass = peptide_mass;
  //  cout << "here" << endl;
  //DDS: Check if the light mod is on the terminals
  if (!heavy && !light) {
    if (strchr(pInput_.szXpressResidues, 'c') != NULL) {
      massdiff += heavy_cterm_mass_ - light_cterm_mass_;
      light = True;
    }
    if (strchr(pInput_.szXpressResidues, 'n') != NULL) {
      massdiff += heavy_nterm_mass_ - light_nterm_mass_;
      light = True;
    } 
  }
  
  if(! unmod && ((light && !heavy) || (! light && heavy)) && massdiff > 0.0) { // valid quant data to process
    if(compute_peptide_mass_) {
      msLight = light ? peptide_mass : peptide_mass - massdiff;
      msHeavy = heavy ? peptide_mass : peptide_mass + massdiff;
    }
#ifndef USE_STD_MODS
    else {
      msLight = light ? pepMass : pInput_.dPeptideMass - massdiff;
      msHeavy = heavy ? pepMass : pInput_.dPeptideMass + massdiff;
    }
#else
    else {
      msLight = light ? pInput_.dCalcPeptideMass : pInput_.dCalcPeptideMass - massdiff;
      msHeavy = heavy ? pInput_.dCalcPeptideMass : pInput_.dCalcPeptideMass + massdiff;
    }
#endif
    if(light)
      cidIndx = 0;
    else
      cidIndx = 1;
    //     cout << "light: " << msLight << " vs heavy: " << msHeavy << endl;
    // how about cidIndx????
  }
  else if (light || heavy) {
    //msLight = 0.0;
    //msHeavy = 0.0;
    //DDS:ASSUME Light
     if(compute_peptide_mass_) {
      msLight = light ? peptide_mass : peptide_mass - massdiff;
      msHeavy = heavy ? peptide_mass : peptide_mass + massdiff;
    }
#ifndef USE_STD_MODS
    else {
      msLight = light ? pInput_.dPeptideMass : pInput_.dPeptideMass - massdiff;
      msHeavy = heavy ? pInput_.dPeptideMass : pInput_.dPeptideMass + massdiff;
    }
#else
    else {
      msLight = light ? pInput_.dCalcPeptideMass : pInput_.dCalcPeptideMass - massdiff;;
      msHeavy = heavy ? pInput_.dCalcPeptideMass : pInput_.dCalcPeptideMass + massdiff;
    }
#endif
    //cout << "no massdiff" << endl;
    cidIndx = 0; // assume light
  }
  // light_mass_ = msLight; // zero if not valid quant peptide

#endif
#ifndef USE_STD_MODS
  // light and heavy sequences
  getPairSequences(pepSeq, prtnSeq, &cidIndx, &msLight, &msHeavy, 
		   modAAs, modAANum, prtnAAs, prtnAANum);
  
#endif
  // get peptide data
  if(msLight <= 0. || msHeavy <= 0. 
     || msHeavy <= msLight) {
    //free(data);
    //cout << "warning1: light " << msLight << " and heavy: " << msHeavy << " for " << pepSeq << endl;
    data_.indx = -1; // set as warning
    //cout << "returning here" << endl;
    return;
  }

  // initialize data
  data_.scan = scan;
  data_.chrg = chrg;
  data_.cidIndx = cidIndx;
  data_.msLight = msLight;
  data_.msHeavy = msHeavy;
  data_.eltn = eltn;
  data_.areaFlag = areaFlag;

  //cout << "reqady for pepdatastrct" << endl;
  // calculate pepDataStrct
  getPepDataStrct(&data_, xmlFile, NULL, 0);

  //  cout << "got struct with data_.index: " << data_.indx << endl;

  return;
}


// This function smoothes rough spectrum within a x range.
void ASAPRatioPeptideParser::smoothSpectFlx(spectStrct * spectrum, double range, int repeats) {
  smoothSpectFlx(spectrum, range, repeats, 10);
}


void ASAPRatioPeptideParser::smoothSpectFlx(spectStrct * spectrum, double range, int repeats, int smoothWindow)
{

  if (pInput_.bUseWaveletSmoothing) {
    // use wavelet smoothing from WaveletQuant implementation

    // begin code from WaveletQuant
    //------------------------------wavelet-----------------------------
    double *temp;
    int size = spectrum->size;
    double *C0, *h, *g, *C1,*C2,*C4;                //sqrt(2)/2,-sqrt(2)/2,0,0
    double **d;
    int i,j,m,loop=0;
    double x,y,k=0;
    if (size<16) {
      m=size;
      while(m>1) {
	m=m/2;
	loop++;
      }
    }
    else loop=4;

    d=(double **)malloc(loop*sizeof(double *));
    C0=(double *)calloc(size, sizeof(double));

    for(i=0;i<size;i++) {
      C0[i]=spectrum->yval[i];
    }
    h=(double *)calloc(size, sizeof(double));
    g=(double *)calloc(size, sizeof(double));
    for(i=0;i<size;i++) {
      h[i]=0;
    }
    for(i=0;i<size;i++) {
      g[i]=0;
    }	
    //db4
    h[0]=0.2304;h[1]=0.7148;h[2]=0.6309;h[3]=-0.0280;h[4]=-0.1870;h[5]=0.0308;h[6]=0.0329;h[7]=-0.0106;
    int t=8;
    int flag=-1;
    for(i=0;i<t;i++) {
      g[i]=flag*h[t-1-i];
      flag*=(-1);
    }
    //loop (unreadable)
    C2=(double *)calloc(size,sizeof(double));
    C1=(double *)calloc(size,sizeof(double));
    for(m=0;m<loop;m++) {
      for(i=0;i<size;i++) {
	C1[i]=0;
	C2[i]=0;
      }
      for(i=0; i<size; i++) {
	for(j=0; j<size; j++) {
	  if(j-2*i<0)
	    x=0;
	  else
	    x=C0[j]*h[j-2*i];

	  C1[i]+=x;

	  if(j-2*i<0) 
	    y=0;
	  else
	    y=C0[j]*g[j-2*i];

	  C2[i]+=y;
	}
      }

      for(i=0; i<size; i++) {
	C0[i]=C1[i];
      }
      d[m]=(double *)malloc(size*sizeof(double));
      for(i=0; i<size; i++) {
	d[m][i]=C2[i];
      }
    }
    // (unreadable)
    double MAD,T,P_temp,avaY;
    double sum=0.0;
    long N,P;
    N=spectrum->size;
    temp=(double *)malloc(N*sizeof(double));
    for(i=0;i<N;i++) {
      temp[i]=d[0][i];
    }
    for(i=0;i<N-1;i++) {
      P=i;
      for(j=i+1;j<N;j++) {
	if(temp[j]<temp[P])	P=j;
      }
      P_temp=temp[P];
      temp[P]=temp[i];
      temp[i]=P_temp;
    }
    avaY=temp[int(N/2)];
    for(i=0;i<N;i++) {
      temp[i]=0;
    }
    for(i=0;i<N;i++) {
      if((d[0][i]-avaY)>=0)
	temp[i]=d[0][i]-avaY;
      else
	temp[i]=avaY-d[0][i];
    }
    for(i=0;i<N-1;i++) {
      P=i;
      for(j=i+1;j<N;j++){
	if(temp[j]<temp[P])	P=j;
      }
      P_temp=temp[P];
      temp[P]=temp[i];
      temp[i]=P_temp;
    }
    MAD=temp[int(N/2)];
    T=MAD/0.6745*sqrt(2*log((double)N));
    //free(temp);
    // (unreadable)
    for(m=0;m<loop;m++)
      for(i=0; i<size; i++) {
	if(fabs(d[m][i])<=T)
	  d[m][i]=0;
      }
    // (unreadable)
    C4=(double *)calloc(size,sizeof(double));

    for(m=loop-1;m>=0;m--) {
      for(i=0;i<size;i++)
	C4[i]=0;
      for(i=0; i<size; i++)
	for(j=0; j<size; j++) {
	  if(i-2*j<0) {
	    x=0;
	    y=0;
	  }
	  else {
	    x=C0[j]*h[i-2*j];
	    y=d[m][j]*g[i-2*j];
	  }
	  C4[i]+=x+y;
	}
      for(i=0; i<size; i++)
	C0[i]=C4[i];
    }
    for(i=0;i<size;i++) {
      spectrum->yval[i]=C0[i];
    }
    free(C0);
    free(h);
    free(g);
    free(C1);
    free(C2);
    free(C4);
    for(i=0;i<loop;i++)
      free(d[i]);
    free(d);
    //------------------------------------------------------------------------------------------
    // end code from WaveletQuant

#ifndef DEBUG_NOFREE
    free(temp);
#endif

  }
  else {
    // use S-G smoothing

    int size = spectrum->size;
    double *temp;
    double threshold;
    int i, j;
    //cout << "here and ready" << endl;

    temp = (double *) calloc(size, sizeof(double));

    threshold = spectrum->yval[0];
    for (i = 0; i < size; ++i) {
      if(spectrum->yval[i] < threshold) threshold = spectrum->yval[i];
    }

    for (j = 0; j < repeats; ++j) {
      for (i = 0; i < size; ++i) {
	//cout << "repeat " << j << ", size: " << i << " of " << size << endl;
	temp[i] = smoothDataPtFlx(*spectrum, i, range, threshold, smoothWindow);
      }
      spectrum->set_yvals(temp);
    }
    //cout << "there and tready" << endl;

#ifndef DEBUG_NOFREE
    free(temp);
#endif

    // end S-G smoothing
  }
  return;
}


// This function smoothes rough spectrum at a specific point.
double ASAPRatioPeptideParser::smoothDataPtFlx(const spectStrct &spectrum, int dtIndx, 
					       double range, double threshold) { 
  return smoothDataPtFlx(spectrum, dtIndx, 
			 range, threshold, 10);
}

double ASAPRatioPeptideParser::smoothDataPtFlx(const spectStrct &spectrum, int dtIndx, 
					       double range, double threshold, int smoothWindow)
{
  int order = 4;
  int lower, upper;
  double value;
  
  // get boundary
  lower = dtIndx;

  while(lower > 0 
	&& spectrum.xval[dtIndx]-spectrum.xval[lower] < range) 
    --lower;
  upper = dtIndx;
  while(upper < spectrum.size-1
	&& spectrum.xval[upper]-spectrum.xval[dtIndx] < range) 
    ++upper;

  //cout << "upper: " << upper << " and lower: " << lower << endl;

  // get filter value
  while(upper-lower < smoothWindow) {
    if(lower > 0)
      --lower;
    if(upper < spectrum.size-1) 
      ++upper;
    if(lower <= 0 
       && upper >= spectrum.size-1) 
      break;
  }
  order = order < upper-lower-1 ? order : upper-lower-1;
  if(order < 1)
    return spectrum.yval[dtIndx];
  else
    value = spectrum.dataFilter(dtIndx, lower, upper, order); 

  if(value > threshold) 
    return value;
  else
    return threshold;

}  

// This function performs LU decomposition.
void ASAPRatioPeptideParser::myLUDcmp(double **mtrx, int order, int *indx, double *d)
{
  double big, dum, sum, temp;
  double *vv;
  int i, j, k, imax=-1;

  vv = (double *) calloc(order, sizeof(double));

  *d = 1.0;
  for (i = 0; i < order; i++) {
    big = 0.0;
    for (j = 0; j < order; j++)
      if ((temp=fabs(mtrx[i][j])) > big) big = temp;
    vv[i] = 1.0/big;
  } 

  for (j = 0; j < order; j++) {
    for (i = 0; i < j; i++) {
      sum = mtrx[i][j];
      for (k = 0; k < i; k++) sum -= mtrx[i][k]*mtrx[k][j];
      mtrx[i][j] = sum;
    }
    big = 0.0;
    for (i = j; i < order; i++) {
      sum = mtrx[i][j];
      for (k = 0; k < j; k++)
	sum -= mtrx[i][k]*mtrx[k][j];
      mtrx[i][j] = sum;
      if ((dum=vv[i]*fabs(sum)) >= big) {
	big = dum;
	imax = i;
      }
    }
    if (j != imax) {
      for (k = 0; k < order; k++) {
	dum = mtrx[imax][k];
	mtrx[imax][k] = mtrx[j][k];
	mtrx[j][k] = dum;
      }
      *d *= -1;
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (mtrx[j][j] == 0.0) mtrx[j][j] = 1.e-20;
    if (j != order) {
      dum = 1.0/(mtrx[j][j]);
      for (i = j+1; i < order; i++) mtrx[i][j] *= dum;
    }
  }
#ifndef DEBUG_NOFREE
  free(vv);
#endif 
  return;
}

// This function performs LU backsubstition.
void ASAPRatioPeptideParser::myLUBksb(double **mtrx, int order, int *indx, double *vec)
{
  double sum;
  int ii, ip;
  int i, j;

  ii = -1;
  for (i = 0; i < order; i++) {
    ip = indx[i];
    sum = vec[ip];
    vec[ip] = vec[i];
    if (ii != -1)
      for (j = ii; j < i; j++) sum -= mtrx[i][j]*vec[j];
    else if (sum) ii = i;
    vec[i] = sum;
  }
  for (i = order-1; i >= 0; i--) {
    sum = vec[i];
    for (j = i+1; j < order; j++) sum -= mtrx[i][j]*vec[j];
    vec[i] = sum/mtrx[i][i];
  }

  return;
}



// For a given set of data, dataErrs, and dataWghs, this function identifies 
// any outliers and gets the mean and confidence interval of ratio.
void ASAPRatioPeptideParser::getDataRatio(double *ratio, double *error, 
					  double *h2l_ratio, double *h2l_error, 
					  double confL, 
					  double *data, double *dataErrs, 
					  double *dataWghs, int *dataIndx, int dataSize,
					  int testType)
{

  int counts[4] = {0, 0, 0, 0};
  double *dtRatios, *dtErrors, *dtWeights;
  double *dt_h2l_Ratios,  *dt_h2l_Errors;
  int *dtIndx;
  int pass;
  int count, vldNum;
  double sum;
  double tmpError;
  double acc = 0.01;
  int i, j;


  // check whether there are valid data
  for(i = 0; i < dataSize; ++i) {
    if(data[i] == 0.) 
      ++counts[0];
    else if(data[i] == -1.) 
      ++counts[1];
    else if(data[i] == -2.) 
      ++counts[2];
    else
      ++counts[3];
  }

  // easy output for invalid data set
  if(counts[3] < 1) {
    if(counts[0] > counts[1]) {
      *ratio = 0.;
      *error = 0.;
      *h2l_ratio = 0.;
      *h2l_error = 0.;
      if(testType == 1) {
	for (i = 0; i < dataSize; ++i) 
	  if(data[i] != 0.)
	    dataIndx[i] = 1;
	  else
	    dataIndx[i] = 0;
      }
    }
    else if(counts[0] < counts[1]) {
      *ratio = -1.;
      *error = 0.;
      *h2l_ratio = -1.;
      *h2l_error = 0.;
      if(testType == 1) {
	for (i = 0; i < dataSize; ++i) 
	  if(data[i] != -1.)
	    dataIndx[i] = 1;
	  else
	    dataIndx[i] = 0;
      }
    }
    else {
      *ratio = -2.;
      *error = 0.;
      *h2l_ratio = -2.;
      *h2l_error = 0.;
      if(testType == 1) {
	for (i = 0; i < dataSize; ++i) 
	  if(data[i] != -2.)
	    dataIndx[i] = 1;
	  else
	    dataIndx[i] = 0;
      }
    }
    return;
  } // if(count[3] < 1) {

  // allocate memory
  dtRatios = (double *) calloc(dataSize, sizeof(double));
  dtErrors = (double *) calloc(dataSize, sizeof(double));

  dt_h2l_Ratios = (double *) calloc(dataSize, sizeof(double));
  dt_h2l_Errors = (double *) calloc(dataSize, sizeof(double));

  dtWeights = (double *) calloc(dataSize, sizeof(double));
  dtIndx = (int *) calloc(dataSize, sizeof(int));


  // identify outliers

  // collect valid data, transform into log(ratio)
  for (i = 0; i < dataSize; ++i) 
    if(data[i] > 0.)
      dataIndx[i] = 0;
    else
      dataIndx[i] = 1;

  count = 0;
  for (i = 0; i < dataSize; ++i) {
    if(dataIndx[i] == 0) {
      pass = 1;
      for (j = 0; pass == 1 && j < count; ++j) {
	if(fabs(log(data[i])-dtRatios[j]) < acc*dtRatios[j]){
	  pass = 0;
	}
      }
      if(pass == 1) {
	dtRatios[count] = log(data[i]);
	dt_h2l_Ratios[count] = log(1/data[i]); 
	++count;
      }
    }
  }

  // identify any outliers
  if(testType != 0)
    DixonTest(dtRatios, dtIndx, count);
  else
    for (i = 0; i < count; ++i)
      dtIndx[i] = 0;

  for (i = 0; i < dataSize; ++i) {
    if(dataIndx[i] == 0) {
      pass = 1;
      for (j = 0; pass == 1 && j < count; ++j) {
	if(dtIndx[j] == 1 
	   && fabs(log(data[i])-dtRatios[j]) < acc*dtRatios[j]){
	  pass = 0;
	}
      }
      if(pass == 0) 
	dataIndx[i] = 1;
    } // if(dataIndx[i] == 0) {
  } //for (i = 0; i < dataSize; ++i) {


  // get ratio and error

  // collect valid date
  count = 0; 
  for (i = 0; i < dataSize; ++i) {
    if(dataIndx[i] == 0) { // not an outlier
      dtRatios[count] = data[i];
      dt_h2l_Ratios[count] = 1 / data[i];
      dtErrors[count] = dataErrs[i];
      dt_h2l_Errors[count] = dataErrs[i] / (data[i] * data[i]);
      dtWeights[count] = dataWghs[i];
      ++count;
    }
  }

  // calculate ratio and error
  if(count < 1) { // no valid data
    *ratio = -2.;
    *error = 0.;
    *h2l_ratio = -2.;
    *h2l_error = 0.;
  }
  else if(count == 1) { // only one valid data
    *ratio = dtRatios[0];
    *error = dtErrors[0];
    *h2l_ratio =  dt_h2l_Ratios[0];
    *h2l_error =  dt_h2l_Errors[0];
  }
  else {
    // transform into log(ratio)
    for (i = 0; i < count; ++i) {
      dtErrors[i] /= dtRatios[i];
      dtRatios[i] = log(dtRatios[i]);
      dt_h2l_Errors[i] /= dt_h2l_Ratios[i];
      dt_h2l_Ratios[i] = log(dt_h2l_Ratios[i]);
    }
    // calculate the light:heavy ratio by weight 
    findMeanAndStdDevWeight(ratio, error, dtRatios, h2l_ratio, h2l_error, dt_h2l_Ratios, dtWeights, count);

    sum = 0.;
    vldNum = 0;
    for (i = 0; i < count; ++i) {
      if(dtErrors[i] > 0.) {
	sum += 1./dtErrors[i]/dtErrors[i];
	++vldNum;
      }
    }
    if(vldNum > 0) {
      tmpError = 1./sqrt(sum);
    }
    else
      tmpError = 0.;

    *error = sqrt((*error)*(*error)+tmpError*tmpError);

    sum = 0.;
    vldNum = 0;
    for (i = 0; i < count; ++i) {
      if(dt_h2l_Errors[i] > 0.) {
	sum += 1./dt_h2l_Errors[i]/dt_h2l_Errors[i];
	++vldNum;
      }
    }
    if(vldNum > 0) {
      tmpError = 1./sqrt(sum);
    } // if(vldNum > 0) {
    else
      tmpError = 0.;

    *h2l_error = sqrt((*h2l_error)*(*h2l_error)+tmpError*tmpError);

    // transform back to ratio
    *ratio = exp(*ratio);
    *error *= (*ratio);
    *h2l_ratio = exp(*h2l_ratio);
    *h2l_error *= (*h2l_ratio);
  }
  
  // free memory
#ifndef DEBUG_NOFREE
  free(dtRatios);
  free(dtErrors);
  free(dt_h2l_Ratios);
  free(dt_h2l_Errors);
  free(dtWeights);
  free(dtIndx);
#endif
  return;
}


// For a set of data and weight, this function finds the mean and standard deviation.
void ASAPRatioPeptideParser::findMeanAndStdDevWeight(double *mean, double *error, double *data, double *h2l_mean, double *h2l_error,
						     double *h2l_data, double *weight, int size)
{
  double sum0, sum1, sum2, h2l_sum1, h2l_sum2, sumW;
  double nEff;
  double mnValue, mxValue, mn_inv_Value, mx_inv_Value;
  int count, i;
  
  if (size < 2) {
    *mean = *data;
    *error = 0.;
    *h2l_mean = *h2l_data;
    *h2l_error = 0.;
    return;
  }

  // ensure weight is valid
  count = 0;
  sum0 = 0.;
  mnValue = data[0];
  mxValue = data[0];
  mn_inv_Value = h2l_data[0];
  mx_inv_Value = h2l_data[0];
  for (i = 0; i < size; ++i) {
    if(weight[i] >= 0.) {
      ++count;
      sum0 += weight[i];
    }
    mnValue = mnValue < data[i] ? mnValue : data[i];
    mxValue = mxValue > data[i] ? mxValue : data[i];

    mn_inv_Value = mn_inv_Value < h2l_data[i] ? mn_inv_Value : h2l_data[i];
    mx_inv_Value = mx_inv_Value > h2l_data[i] ? mx_inv_Value : h2l_data[i];
  }

  if(mnValue >= mxValue) {
    *mean = mnValue;
    *error = 0.;
    *h2l_mean = mn_inv_Value;
    *h2l_error = 0.;
    return;
  }

  if(count < size || sum0 == 0.) {   // no all have valid weight
    if(count < 1 || sum0 == 0.) {   // if no data has weight
      for (i = 0; i < size; ++i) {
	weight[i] = 1.;
      }
    }
    else {
      sum0 /= count;
      for (i = 0; i < size; ++i) {
	if(weight[i] < 0.) {
	  weight[i] = sum0;
	}
      }
    }
  }

  // get mean and std. dev.
  sum0 = 0.;
  sum1 = 0.;
  sum2 = 0.;
  h2l_sum1 = 0.;
  h2l_sum2 = 0.;
  sumW = 0.;
  for (i = 0; i < size; ++i) {
    sum0 += weight[i];
    sum1 += data[i]*weight[i];
    sum2 += data[i]*data[i]*weight[i];
    h2l_sum1 += h2l_data[i]*weight[i];
    h2l_sum2 += h2l_data[i]*h2l_data[i]*weight[i];
    sumW += weight[i]*weight[i];
  }

  // get mean
  *mean = sum1/sum0;
  *h2l_mean = h2l_sum1/sum0;

  // get std. dev.
  if(sum2*sum0-sum1*sum1 > 0.) {
    nEff = sum0*sum0/sumW;
    if(nEff > 2.) {
      *error = sqrt((sum2*sum0-sum1*sum1)*nEff/(nEff-1.))/sum0;
      *h2l_error = sqrt((h2l_sum2*sum0-h2l_sum1*h2l_sum1)*nEff/(nEff-1.))/sum0;
    }
    else {
      *error = sqrt(2.*(sum2*sum0-sum1*sum1))/sum0;
      *h2l_error = sqrt(2.*(h2l_sum2*sum0-h2l_sum1*h2l_sum1))/sum0;
    }
  }
  else {
    *error = 0.;
    *h2l_error = 0.;
  }    

  return;
}


// This function uses Dixon's test with alpha = 0.05 to identify any outliers.
void ASAPRatioPeptideParser::DixonTest(double *data, int *outliers, int size)
{
  // cutoff values in Dixon's test: n = 3, ..., 30, INF.
  double ya[29] = {0.941, 0.765, 0.642, 0.560, 0.507, 0.554,
		   0.512, 0.477, 0.576, 0.546, 0.521, 0.546,
		   0.525, 0.507, 0.490, 0.475, 0.462, 0.450,
		   0.440, 0.430, 0.421, 0.413, 0.406, 0.399,
		   0.393, 0.387, 0.381, 0.376, 0.};

  // values of 1/n: 1/3, ..., 1/30, 1/INF.
  double xa[29] = {0.333333, 0.250000, 0.200000, 0.166667, 0.142857, 
		   0.125000, 0.111111, 0.100000, 0.090909, 0.083333, 
		   0.076923, 0.071429, 0.066667, 0.062500, 0.058824, 
		   0.055556, 0.052632, 0.050000, 0.047619, 0.045455, 
		   0.043478, 0.041667, 0.040000, 0.038462, 0.037037, 
		   0.035714, 0.034483, 0.033333, 0.};
  
  int cnstSize = 29;
  int *dataIndx;
  int startIndx, endIndx;
  int count;
  double ratio1, ratio2;
  double cutoff, x;
  int i, j;

  // assume none is an outlier
  for (i = 0; i < size; ++i)
    outliers[i] = 0;
  if (size < 3) // not enough data for checking
    return;
  
  // get dataIndx for ordered data
  dataIndx = (int *) calloc(size, sizeof(int));
  for (i = 0; i < size; ++i)
    dataIndx[i] = i;
  for(i = 0; i < size; ++i) {
    for(j = 0; j < size-i-1; ++j) {
      if(data[dataIndx[j]] > data[dataIndx[j+1]]) {
	count = dataIndx[j];
	dataIndx[j] = dataIndx[j+1];
	dataIndx[j+1] = count;
      }
    }
  }
  
  // check for outliers
  count = 0;
  startIndx = 0;
  endIndx = size;
  while(size > 2  // enough data for checking
	&& count != size // look for more when an outlier is identified
	&& data[dataIndx[startIndx]] != data[dataIndx[endIndx-1]]) { 

    // restore size
    count = size;

    // get cutoff
    if (size < 3)
      cutoff = 1.;
    else if(size <= cnstSize+1)
      cutoff = ya[size-3];
    else {
      x = 1./((double) size);
      cutoff = PadeApprx(x, xa, ya, cnstSize);
    }

    // get ratio
    if(size < 8) {
      ratio1 = (data[dataIndx[startIndx+1]]-data[dataIndx[startIndx]])
	/(data[dataIndx[endIndx-1]]-data[dataIndx[startIndx]]);
      ratio2 = (data[dataIndx[endIndx-1]]-data[dataIndx[endIndx-2]])
	/(data[dataIndx[endIndx-1]]-data[dataIndx[startIndx]]);
    }
    else if(size < 11) {
      if(data[dataIndx[startIndx]] != data[dataIndx[endIndx-2]]) { 
	ratio1 = (data[dataIndx[startIndx+1]]-data[dataIndx[startIndx]])
	  /(data[dataIndx[endIndx-2]]-data[dataIndx[startIndx]]);
      }
      else
	ratio1 = 0.;
      if(data[dataIndx[startIndx+1]] != data[dataIndx[endIndx-1]]) { 
	ratio2 = (data[dataIndx[endIndx-1]]-data[dataIndx[endIndx-2]])
	  /(data[dataIndx[endIndx-1]]-data[dataIndx[startIndx+1]]);
      }
      else
	ratio2 = 0.;  
    }
    else if(size < 14) {
      if(data[dataIndx[startIndx]] != data[dataIndx[endIndx-2]]) { 
	ratio1 = (data[dataIndx[startIndx+2]]-data[dataIndx[startIndx]])
	  /(data[dataIndx[endIndx-2]]-data[dataIndx[startIndx]]);
      }
      else
	ratio1 = 0.;
      if(data[dataIndx[startIndx+1]] != data[dataIndx[endIndx-1]]) { 
	ratio2 = (data[dataIndx[endIndx-1]]-data[dataIndx[endIndx-3]])
	  /(data[dataIndx[endIndx-1]]-data[dataIndx[startIndx+1]]);
      }
      else
	ratio2 = 0.;  
    }
    else {
      if(data[dataIndx[startIndx]] != data[dataIndx[endIndx-3]]) { 
	ratio1 = (data[dataIndx[startIndx+2]]-data[dataIndx[startIndx]])
	  /(data[dataIndx[endIndx-3]]-data[dataIndx[startIndx]]);
      }
      else
	ratio1 = 0.;
      if(data[dataIndx[startIndx+2]] != data[dataIndx[endIndx-1]]) { 
	ratio2 = (data[dataIndx[endIndx-1]]-data[dataIndx[endIndx-3]])
	  /(data[dataIndx[endIndx-1]]-data[dataIndx[startIndx+2]]);
      }
      else
	ratio2 = 0.;  
    }

    // check ratio
    if(ratio1 > ratio2) {
      if(ratio1 > cutoff) { // an outlier
	outliers[dataIndx[startIndx]] = 1;
	--size;
	++startIndx;
      }
    }
    else {
      if(ratio2 > cutoff) { // an outlier
	outliers[dataIndx[endIndx-1]] = 1;
	--size;
	--endIndx;
      }
    }
  } // while(size > 2  // enough data for checking

#ifndef DEBUG_NOFREE
  free(dataIndx);
#endif
  return;
}


// This function returns the value of Pade Approximation.
double ASAPRatioPeptideParser::PadeApprx(double x, double *xa, double *ya, int size)
{
  double y, dy;
  double tiny = 1.e-25;
  int m,i,ns=1;
  double w,t,hh,h,dd,*c,*d;
  double *xb, *yb;
  int n = size;

  // convert into 1 ... n
  xb = (double *) calloc(n+1, sizeof(double));
  yb = (double *) calloc(n+1, sizeof(double));
  for (i = 0; i < size; ++i) {
    xb[i+1] = xa[i];
    yb[i+1] = ya[i];
  }

  // use ratint
  c = (double *) calloc(n+1, sizeof(double));
  d = (double *) calloc(n+1, sizeof(double));

  hh=fabs(x-xb[1]);
  for (i=1;i<=n;i++) {
    h=fabs(x-xb[i]);
    if (h == 0.0) {
      y=yb[i];
      dy=0.0;
#ifndef DEBUG_NOFREE
      free(c);
      free(d);
      free(xb);
      free(yb);
#endif
      return y;
    } else if (h < hh) {
      ns=i;
      hh=h;
    }
    c[i]=yb[i];
    d[i]=yb[i]+tiny;
  }
  y=yb[ns--];
  for (m=1;m<n;m++) {
    for (i=1;i<=n-m;i++) {
      w=c[i+1]-d[i];
      h=xb[i+m]-x;
      t=(xb[i]-x)*d[i]/h;
      dd=t-c[i+1];
      if (dd == 0.0) {
	printf("Error in routine PadeApprx\n");
#ifndef DEBUG_NOFREE
	free(c);
	free(d);
	free(xb);
	free(yb);
#endif
	return y;
      }
      dd=w/dd;
      d[i]=c[i+1]*dd;
      c[i]=t*dd;
    }
    y += (dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
#ifndef DEBUG_NOFREE
  free(c);
  free(d);
  free(xb);
  free(yb);
#endif
  return y;
}


// This function returns a light sequence and its heavy partner.
void ASAPRatioPeptideParser::getPairSequences(char *sequence, char *prtnSequence, 
					      int *cidIndx, double *msLight, double *msHeavy, 
					      residueStrct *modAAs, int modAANum,
					      pairStrct *prtnAAs, int prtnAANum)
{
  char tmpSequence[1000];
  int seqIndx, prtnIndx;
  int lngth = (int)strlen(sequence);
  //char end[3];
  int flag;
  double msTh, msPt;
  int i;

  // get partner sequence

  // check for N-end
  seqIndx = 0;
  prtnIndx = 0;

  // get any n term label at outset
  for (i = 0; i < prtnAANum; ++i) {
    if(strlen(prtnAAs[i].prtnA) > 0 && prtnAAs[i].prtnA[0] == 'n') { // n term
      if(strlen(prtnAAs[i].prtnA) > 1 && sequence[0] == prtnAAs[i].prtnA[1]) { // got a match
	if(strlen(prtnAAs[i].prtnB) > 1)  // use its symbol at front end
	  prtnSequence[prtnIndx++] = prtnAAs[i].prtnB[1];
	seqIndx++;
      }
      else if(strlen(prtnAAs[i].prtnB) > 1 && sequence[0] == prtnAAs[i].prtnB[1]) { // got a match
	if(strlen(prtnAAs[i].prtnA) > 1)  // use its symbol at front end
	  prtnSequence[prtnIndx++] = prtnAAs[i].prtnA[1];
	seqIndx++;
      }
      else { // choose the first heavy one for the ending
	if(strlen(prtnAAs[i].prtnA) > 1)
	  prtnSequence[prtnIndx++] = prtnAAs[i].prtnA[1];
	else if(strlen(prtnAAs[i].prtnB) > 1)
	  prtnSequence[prtnIndx++] = prtnAAs[i].prtnB[1];
      }
      i = prtnAANum;
    } // n term
   }

  // check sequence one by one
  while(seqIndx < lngth) {
    flag = -1; // assume having no partner
    // check on prtnA
    for (i = 0; i < prtnAANum && flag == -1; ++i) {
      if(strncmp(sequence+seqIndx, prtnAAs[i].prtnA, 
		 strlen(prtnAAs[i].prtnA)) == 0) {
	flag = i;
      }
    }
    if (flag != -1) { // has a partner
      strcpy(prtnSequence+prtnIndx, prtnAAs[flag].prtnB);
      seqIndx += (int)strlen(prtnAAs[flag].prtnA);
      prtnIndx += (int)strlen(prtnAAs[flag].prtnB);
      continue;
    }
    // check on prtnB
    for (i = 0; i < prtnAANum && flag == -1; ++i) {
      if(strncmp(sequence+seqIndx, prtnAAs[i].prtnB, 
		 strlen(prtnAAs[i].prtnB)) == 0) {
	flag = i;
      }
    }
    if (flag == -1) { // has no partner
      prtnSequence[prtnIndx++] = sequence[seqIndx++];
    }
    else { // has a partner
      strcpy(prtnSequence+prtnIndx, prtnAAs[flag].prtnA);
      seqIndx += (int)strlen(prtnAAs[flag].prtnB);
      prtnIndx += (int)strlen(prtnAAs[flag].prtnA);
    }
  } // while(seqIndx < lngth) {


  // get any c term label now
  for (i = 0; i < prtnAANum; ++i) {
    //cout << prtnAAs[i].prtnA << " vs " << prtnAAs[i].prtnB << endl;
    if(strlen(prtnAAs[i].prtnA) > 0 && prtnAAs[i].prtnA[0] == 'c') { // c term
      if(strlen(prtnAAs[i].prtnA) > 1 && sequence[lngth-1] == prtnAAs[i].prtnA[1]) { // got a match
	if(strlen(prtnAAs[i].prtnB) > 1)  // use its symbol at front end
	  prtnSequence[prtnIndx++] = prtnAAs[i].prtnB[1];
      }
      else if(strlen(prtnAAs[i].prtnB) > 1 && sequence[lngth-1] == prtnAAs[i].prtnB[1]) { // got a match
	if(strlen(prtnAAs[i].prtnA) > 1)  // use its symbol at front end
	  prtnSequence[prtnIndx++] = prtnAAs[i].prtnA[1];
      }
      else { // choose the first heavy one for the ending
	if(strlen(prtnAAs[i].prtnA) > 1)
	  prtnSequence[prtnIndx++] = prtnAAs[i].prtnA[1];
	else if(strlen(prtnAAs[i].prtnB) > 1)
	  prtnSequence[prtnIndx++] = prtnAAs[i].prtnB[1];
      }
      i = prtnAANum;
    } // c term
  }


  prtnSequence[prtnIndx] = '\0';
  //cout << "partner of " << sequence << " is " << prtnSequence << endl;


  // decide sequence

  // theoretical mass
  msTh = getSeqMonoMZ(sequence, 0, modAAs, modAANum);
  msPt = getSeqMonoMZ(prtnSequence, 0, modAAs, modAANum);
  if(msTh <= 0. || msPt <= 0. || msTh == msPt){
    *msLight = -1.;
    *msHeavy = -1.;
    //cout <<  "Error, mass for " << sequence << " is " << msTh << " and for " << prtnSequence << " is " << msPt << endl; 
    //cout << modAAs_[0].sz << " " << modAAs_[0].rp << " " << modAAs_[0].ms << endl;

    if(strcmp(sequence, prtnSequence) > 0) {
      *cidIndx = 1;
      strcpy(tmpSequence, sequence);
      strcpy(sequence, prtnSequence);
      strcpy(prtnSequence, tmpSequence);
    }
    else
      *cidIndx = 0;    
  }
  else if(msTh < msPt){
    *cidIndx = 0;
    *msLight = msTh;
    *msHeavy = msPt;
  }
  else if(msTh > msPt){
    *cidIndx = 1;
    *msLight = msPt;
    *msHeavy = msTh;
    strcpy(tmpSequence, sequence);
    strcpy(sequence, prtnSequence);
    strcpy(prtnSequence, tmpSequence);
  }

  return;
}


// For a given sequence "char * sequence" of charge "int qState", this function gives its (mono-isotopic mass)/charge.
double ASAPRatioPeptideParser::getSeqMonoMZ(char *sequence, int qState,
					    residueStrct *modifiedAA, int modAANum)
{
  double mass;
  int length = (int)strlen(sequence);
  int count = 0; 
  int flag = 0;
  char end[3];
  char *tmpSeq;
  int i;

  Boolean verbose = False;    

  // only positive or neutral charge states allowed
  if (qState < 0) {
    printf("Cannot handle negative charge state.\n");
    return -1.;
  }

  // check input
  getRidOfSpace(sequence);
  cnvtUpper(sequence);

  // mass of "qState*H"
  mass = qState*nativeAA[0].ms;


  // mass of N-end
  strcpy(end, "n");

  if(static_quant_) {
    for (i = 0; i < modAANum && flag == 0; ++i) {
      //cout << modifiedAA[i].rp << " 9878 " << modifiedAA[i].sz << endl;
      if(modifiedAA[i].sz > 1 && end[0] == modifiedAA[i].rp[0] && sequence[0] == modifiedAA[i].rp[1]) {
	flag = 1;
	mass += modifiedAA[i].ms;
	count++;
      }
    }
    for (i = 0; i < modAANum && flag == 0; ++i) {
      //cout << modifiedAA[i].rp << " 9878 " << modifiedAA[i].sz << endl;
      if(modifiedAA[i].sz == 1 && end[0] == modifiedAA[i].rp[0]) {
	flag = 1;
	mass += modifiedAA[i].ms;
      }
    }
    if(flag == 0) { // && static_status_ == -1) { // light label
      for (i = 0; i < NATIVEAANUM && flag == 0; ++i) {
	if(strncmp(end, nativeAA[i].rp, nativeAA[i].sz) == 0) {
	  flag = 1;
	  mass += nativeAA[i].ms;
	  if(verbose)
	    cout << "adding..." << nativeAA[i].ms << " ";
	}
      }
    }
    // if no match
    if(flag == 0) {
      printf("N-end unspecified.\n");
      return -1.;
    }
  }
  else { // non-static
    for (i = 0; i < modAANum && flag == 0; ++i) {
      //cout << modifiedAA[i].rp << " 9878 " << modifiedAA[i].sz << endl;
      if(modifiedAA[i].sz > 1 && sequence[0] == modifiedAA[i].rp[1]) {
	flag = 1;
	mass += modifiedAA[i].ms;
	count++;
      }
    }
    // check for native Amino Acids
    for (i = 0; i < NATIVEAANUM && flag == 0; ++i) {
      if(strncmp(end, nativeAA[i].rp, nativeAA[i].sz) == 0) {
	flag = 1;
	mass += nativeAA[i].ms;
      }
    }
    // if no match
    if(flag == 0) {
      printf("N-end unspecified.\n");
      return -1.;
    }
  } // non static


  // mass of c-end
  flag = 0;
  strcpy(end, "c");

  if(static_quant_) {
    for (i = 0; i < modAANum && flag == 0; ++i) {
      //cout << modifiedAA[i].rp << " 9878 " << modifiedAA[i].sz << " and seq term: " << sequence + length -1  << endl;
      if(modifiedAA[i].sz > 1 && end[0] == modifiedAA[i].rp[0] && sequence[length-1] == modifiedAA[i].rp[1]) {
	//cout << "YES!" << endl;
	flag = 1;
	mass += modifiedAA[i].ms;
	if(verbose)
	  cout << modifiedAA[i].ms << " ";
	length--;
      }
    }
    for (i = 0; i < modAANum && flag == 0; ++i) {
      if(modifiedAA[i].sz == 1 && end[0] == modifiedAA[i].rp[0]) {
	flag = 1;
	mass += modifiedAA[i].ms;
	  if(verbose)
	cout << modifiedAA[i].ms << " ";
      }
    }
    if(flag == 0) { // && static_status_ == -1) { // light label
      for (i = 0; i < NATIVEAANUM && flag == 0; ++i) {
	if(strncmp(end, nativeAA[i].rp, nativeAA[i].sz) == 0) {
	  flag = 1;
	  mass += nativeAA[i].ms;
	  if(verbose)
	    cout << nativeAA[i].ms << " ";
	}
      }
    }
    // if no match
    if(flag == 0) {
      printf("C-end unspecified.\n");
      return -1.;
    }

  }
  else { // non-static
    for (i = 0; i < modAANum && flag == 0; ++i) {
      //cout << modifiedAA[i].rp << " 9878 " << modifiedAA[i].sz << endl;
      if(modifiedAA[i].sz > 1 && sequence[length-1] == modifiedAA[i].rp[1]) {
	flag = 1;
	mass += modifiedAA[i].ms;
	length--;
      }
    }
    // check for native Amino Acids
    for (i = 0; i < NATIVEAANUM && flag == 0; ++i) {
      if(strncmp(end, nativeAA[i].rp, nativeAA[i].sz) == 0) {
	flag = 1;
	mass += nativeAA[i].ms;
      }
    }
    // if no match
    if(flag == 0) {
      printf("C-end unspecified.\n");
      return -1.;
    }
  } // non static

  //cout << "seq: " << sequence << endl;

  // add sequence mass
  while(count < length) {
    flag = 0;

    // check for modified Amino Acids first (with extra symbol)
    for (i = 0; i < modAANum && flag == 0; ++i) {
      //cout << modifiedAA[i].rp << " " << modifiedAA[i].sz << endl;
      if(modifiedAA[i].sz > 1 && strncmp(sequence+count, modifiedAA[i].rp, modifiedAA[i].sz) == 0) {
	flag = 1;
	mass += modifiedAA[i].ms;
	if(verbose)
	  cout << modifiedAA[i].ms << " ";
	count += modifiedAA[i].sz;
      }
    }
    // check for modified Amino Acids first (without extra symbol)
    for (i = 0; i < modAANum && flag == 0; ++i) {
      //cout << modifiedAA[i].rp << " " << modifiedAA[i].sz << endl;
      if(modifiedAA[i].sz == 1 && strncmp(sequence+count, modifiedAA[i].rp, modifiedAA[i].sz) == 0) {
	flag = 1;
	mass += modifiedAA[i].ms;
	if(verbose)
	  cout << modifiedAA[i].ms << " ";
	count += modifiedAA[i].sz;
      }
    }

    // check for native Amino Acids
    for (i = 0; i < NATIVEAANUM && flag == 0; ++i) {
      if(strncmp(sequence+count, nativeAA[i].rp, nativeAA[i].sz) == 0) {
	flag = 1;
	mass += nativeAA[i].ms;
	if(verbose)
	  cout << nativeAA[i].ms << " ";
	count += nativeAA[i].sz;
      }
    }
    // if no match
    if(flag == 0) {
      tmpSeq = (char *) calloc(count+1, sizeof(count));
      strncpy(tmpSeq, sequence, count);
      tmpSeq[count] = '\0';
      printf("Unknown amino acid at \"%s'%c'%s\". \n", 
	     tmpSeq, sequence[count], sequence+count+1);
#ifndef DEBUG_NOFREE
      free(tmpSeq);
#endif
      return -1.;
    }
  } // while(count < length) {

  if(verbose)
    printf("\ncomputed mass of %s: %0.2f (ch: %d)\n", sequence, mass, qState);
  //cout << "computed mass of " << sequence << ": " << mass << " (ch: " << qState << ")" << endl;

  if(qState > 0)
    return mass/qState;
  else
    return mass;
}


// This function evaluates the light:heavy ratio of a peptide.
// For cgi display: cgiIndx = 1 and pngFileBase != NULL.
void ASAPRatioPeptideParser::getPepDataStrct(pepDataStrct *data, char *xmlFile,
					     char *pngFileBase, int cgiIndx)
{
  lcSpectStrct lcSpect;
  int lcSize, cidScan;
  int cidIndx = data->cidIndx;
  int count, qState;
  double value[_ASAPRATIO_MXQ_];
  double error[_ASAPRATIO_MXQ_];
  double weight[_ASAPRATIO_MXQ_];
  int outliers[_ASAPRATIO_MXQ_];
  double timeWd = 3.; // time range for LC spetra in minutes
  int scanRange = 50;
  int pkIndx[2];
  int valleyScan;
  double area;
  int startZ = 0;
  int endZ = _ASAPRATIO_MXQ_;
  if (pInput_.bQuantCIDChrgOnly) {
    startZ = pInput_.iChargeState-1;
    endZ = pInput_.iChargeState;
  }
  int i, j;

  if(data->indx < 0) {
    cout << "returning with dataindx 0" << endl;
    return;
  }


  // get LC spectrum

  //cout << "goping to getpeak spectra with time " << timeWd << " file " << xmlFile << endl;

  if(getPeakSpectra(&lcSpect, data, timeWd, xmlFile) == 0){
    data->indx = -1;
    //cout << "returning with index -1" << endl;
    return;
  }

  lcSize = lcSpect.rawSpectrum[data->chrg-1][cidIndx].size;
  cidScan = 0;
  while(cidScan < lcSize-1 
	&& lcSpect.expScanNums[cidScan] < data->scan) 
    ++cidScan;

  // convert valleys into scan numbers in LC/MS spectrum
  if(data->indx > 0) {
    for (i = startZ; i < endZ; ++i){
      for (j = 0; j < 2; ++j) {
	if(data->peaks[i][j].indx >= 0) {
	  if(lcSpect.expScanNums[0] > data->peaks[i][j].valley[1]
	     || lcSpect.expScanNums[lcSize-1] < data->peaks[i][j].valley[0]) {
	    data->peaks[i][j].peak = lcSize/2;
	    data->peaks[i][j].valley[0] = lcSize/2;
	    data->peaks[i][j].valley[1] = lcSize/2;
	  }
	  else {
	    valleyScan = 0;
	    while(valleyScan < lcSize-1
		  && lcSpect.expScanNums[valleyScan] 
		  < data->peaks[i][j].valley[0]) 
	      ++valleyScan;
	    data->peaks[i][j].valley[0] = valleyScan;
	    while(valleyScan < lcSize-1
		  && lcSpect.expScanNums[valleyScan] 
		  < data->peaks[i][j].peak) 
	      ++valleyScan;
	    data->peaks[i][j].peak = valleyScan;
	    while(valleyScan < lcSize-1
		  && lcSpect.expScanNums[valleyScan] 
		  < data->peaks[i][j].valley[1]) 
	      ++valleyScan;
	    data->peaks[i][j].valley[1] = valleyScan;
	    if(data->peaks[i][j].valley[0] >= data->peaks[i][j].valley[1]) {
	      data->peaks[i][j].peak = lcSize/2;
	      data->peaks[i][j].valley[0] = lcSize/2;
	      data->peaks[i][j].valley[1] = lcSize/2;
	    }
	  } // else {
	} //if(data->peaks[i][j].indx >= 0) {
      } //for (j = 0; j < 2; ++j) {
    } //for (i = 0; i < _ASAPRATIO_MXQ_; ++i){
  } // if(data->indx > 0) {

  // get peak scan
  last_valley_[0] = -1;
  last_valley_[1] = -1;
  if(data->indx == 0) {
    // identified peptide 
    pkIndx[cidIndx] = -1;
    qState = data->chrg;
    if(data->peaks[qState-1][cidIndx].indx >= 0) {
      getLCPeakStrct(&(data->peaks[qState-1][cidIndx]), 
		     lcSpect.rawSpectrum[qState-1][cidIndx],
		     lcSpect.fitSpectrum[qState-1][cidIndx],
		     cidScan, scanRange);
      pkIndx[cidIndx] = qState - 1;
    } // if(data->peaks[qState-1][cidIndx].indx >= 0) {
    else {
      for (i = 1; i <  _ASAPRATIO_MXQ_; ++i){
	if(data->chrg > _ASAPRATIO_MXQ_/2) {
	  qState = data->chrg - i;
	  if(qState <= 0) qState += _ASAPRATIO_MXQ_;
	}	    
	else {
	  qState = data->chrg + i;
	  if(qState > _ASAPRATIO_MXQ_) qState -= _ASAPRATIO_MXQ_;
	}
	if (pInput_.bQuantCIDChrgOnly) {
	  if (pInput_.iChargeState != qState)
	    break;
	}
	if(data->peaks[qState-1][cidIndx].indx >= 0) {
	  getLCPeakStrct(&(data->peaks[qState-1][cidIndx]), 
			 lcSpect.rawSpectrum[qState-1][cidIndx],
			 lcSpect.fitSpectrum[qState-1][cidIndx],
			 cidScan, scanRange);
	  if(data->peaks[qState-1][cidIndx].indx > 1) {
	    pkIndx[cidIndx] = qState - 1;
	    break;
	  }
	} // if(data->peaks[qState-1][cidIndx].indx >= 0) {
      } // for (i = 1; i < _ASAPRATIO_MXQ_; ++i){
      if(pkIndx[cidIndx] < 0) {
	area = 0.;
	for (i = startZ; i < endZ; ++i){
	  if(data->peaks[i][cidIndx].indx > 0
	     && data->peaks[i][cidIndx].area[0] > area) {
	    area = data->peaks[i][cidIndx].area[0];
	    pkIndx[cidIndx] = i;
	  }
	}
      }
      if(pkIndx[cidIndx] < 0) 
	pkIndx[cidIndx] = data->chrg - 1;

      for (i = startZ; i < endZ; ++i){
	if(data->peaks[i][cidIndx].indx > 0
	   && pkIndx[cidIndx] != i) {
	  data->peaks[i][cidIndx].indx = 0;
	}
      }
    } // else {

    // peptide partner
    if(cidIndx == 0) { // light
      cidScan = data->peaks[pkIndx[cidIndx]][cidIndx].peak 
	+ data->eltn*(data->peaks[pkIndx[cidIndx]][cidIndx].valley[1]
		      -data->peaks[pkIndx[cidIndx]][cidIndx].valley[0])/4;
    }
    else { // heavy
      cidScan = data->peaks[pkIndx[cidIndx]][cidIndx].peak 
	- data->eltn*(data->peaks[pkIndx[cidIndx]][cidIndx].valley[1]
		      -data->peaks[pkIndx[cidIndx]][cidIndx].valley[0])/4;
    }
    cidScan = cidScan > 0 ? cidScan : 0;
    cidScan = cidScan < lcSize-1 ? cidScan : lcSize-1;

    cidIndx = (cidIndx+1)%2;
    pkIndx[cidIndx] = -1;
    qState = pkIndx[(cidIndx+1)%2] + 1; 
    if(data->peaks[qState-1][cidIndx].indx >= 0) {
      getLCPeakStrct(&(data->peaks[qState-1][cidIndx]), 
		     lcSpect.rawSpectrum[qState-1][cidIndx],
		     lcSpect.fitSpectrum[qState-1][cidIndx],
		     cidScan, scanRange);
      if(data->peaks[qState-1][cidIndx].indx > 1) {
	pkIndx[cidIndx] = qState - 1;
      }
    } // if(data->peaks[qState-1][cidIndx].indx >= 0) {
    if(pkIndx[cidIndx] < 0) {
      for (i = 1; i < _ASAPRATIO_MXQ_; ++i){
	if(data->chrg > _ASAPRATIO_MXQ_/2) {
	  qState = data->chrg - i;
	  
	  if(qState <= 0) qState += _ASAPRATIO_MXQ_;
	}	    
	else {
	  qState = data->chrg + i;

	  if(qState > _ASAPRATIO_MXQ_) qState -= _ASAPRATIO_MXQ_;
	}
	if (pInput_.bQuantCIDChrgOnly) {
	  if (pInput_.iChargeState != qState)
	    break;
	}
	if(data->peaks[qState-1][cidIndx].indx >= 0) {
	  getLCPeakStrct(&(data->peaks[qState-1][cidIndx]), 
			 lcSpect.rawSpectrum[qState-1][cidIndx],
			 lcSpect.fitSpectrum[qState-1][cidIndx],
			 cidScan, scanRange);
	  if(data->peaks[qState-1][cidIndx].indx > 1) {
	    pkIndx[cidIndx] = qState - 1;
	    break;
	  }
	} // if(data->peaks[qState-1][cidIndx].indx >= 0) {
      } // for (i = 1; i < _ASAPRATIO_MXQ_; ++i){
    } // if(pkIndx[cidIndx] < 0) {
    if(pkIndx[cidIndx] < 0) {
      area = 0.;
      for (i = startZ; i < endZ; ++i){
	if(data->peaks[i][cidIndx].indx > 0
	   && data->peaks[i][cidIndx].area[0] > area) {
	  area = data->peaks[i][cidIndx].area[0];
	  pkIndx[cidIndx] = i;
	}
      }
    }
    if(pkIndx[cidIndx] < 0) 
      pkIndx[cidIndx] = pkIndx[(cidIndx+1)%2];

    for (i = startZ; i < endZ; ++i){
      if(data->peaks[i][cidIndx].indx > 0
	 && pkIndx[cidIndx] != i) {
	data->peaks[i][cidIndx].indx = 0;
      }
    }
  } //   if(data->indx == 0) {


  // get individual ratios
  area = 0.;

    
  for (i = startZ; i < endZ; ++i){
    if(data->indx == 0) 
      data->pkCount[i] = 1;

    //last_valley_[0] = -1;
    //last_valley_[1] = -1;
    for (j = 0; j < 2; ++j) {
      if(data->peaks[i][j].indx == 0) {
	getLCPeakStrct(&(data->peaks[i][j]), lcSpect.rawSpectrum[i][j], 
		       lcSpect.fitSpectrum[i][j], 
		       data->peaks[pkIndx[j]][j].peak, scanRange);
      }
      else if(data->peaks[i][j].indx > 0
	      && data->pkCount[i] == 1) { 
	data->peaks[i][j].indx 
	  = getLCSpectAreaAndTime(data->peaks[i][j].area, data->peaks[i][j].time,
				  lcSpect.rawSpectrum[i][j], lcSpect.fitSpectrum[i][j], 
				  data->peaks[i][j].peak, data->peaks[i][j].valley,
				  data->peaks[i][j].bckgrnd);
      }
      if(data->indx == 0
	 && data->peaks[i][j].indx >= 0
	 && fabs((double)data->peaks[i][j].peak-data->peaks[pkIndx[j]][j].peak)
	 > (double)(data->peaks[i][j].valley[1]-data->peaks[i][j].valley[0])/4.0
	 + (double)(data->peaks[pkIndx[j]][j].valley[1]-data->peaks[pkIndx[j]][j].valley[0])/4.0) 
	data->pkCount[i] = 0;
    }
    if(data->peaks[i][0].indx < 0 // at least 1 invalid peak
       || data->peaks[i][1].indx < 0
       || (data->peaks[i][0].indx == 1 // no valid area 
	   && data->peaks[i][1].indx == 1)) { 
      data->pkRatio[i] = -2;
      data->pkError[i] = 0.;      
      data->pkCount[i] = 0;      
    }
    else if(data->peaks[i][0].indx == 1) { // light invalid
      data->pkRatio[i] = 0.;
      data->pkError[i] = 0.;      
    }
    else if(data->peaks[i][1].indx == 1) { // heavy invalid
      data->pkRatio[i] = -1.;
      data->pkError[i] = 0.;
    }
    else { // both valid
      data->pkRatio[i] = data->peaks[i][0].area[0]/data->peaks[i][1].area[0];
      data->pkError[i] = data->pkRatio[i]
	*sqrt(data->peaks[i][0].area[1]*data->peaks[i][0].area[1]
	      /data->peaks[i][0].area[0]/data->peaks[i][0].area[0]
	      + data->peaks[i][1].area[1]*data->peaks[i][1].area[1]
	      /data->peaks[i][1].area[0]/data->peaks[i][1].area[0]);
    }
    if(data->pkCount[i] == 1
       && area < data->peaks[i][0].area[0]+data->peaks[i][1].area[0])
      area = data->peaks[i][0].area[0] + data->peaks[i][1].area[0];      
  } //for (i = 0; i < _ASAPRATIO_MXQ_; ++i){

  // strong data only
  for (i = startZ; i < endZ; ++i){
    if ((pInput_.bQuantCIDChrgOnly && data->chrg!=i+1) ||
	(data->indx == 0
	 && data->pkCount[i] == 1
	 && data->peaks[i][0].area[0]+data->peaks[i][1].area[0] < 0.1*area)) {
      data->pkCount[i] = 0;
    }
  }

  // get peptide output

  // collect valid date
  count = 0;
  for (i = startZ; i < endZ; ++i){
    if (data->pkCount[i] == 1) {  
      value[count] = data->pkRatio[i];
      error[count] = data->pkError[i];
      weight[count] = data->peaks[i][0].area[0] 
	+ data->peaks[i][1].area[0];
      outliers[count] = 0; // assume not an outlier
      ++count;
    }
  }

  // calculate ratio and error
  if(count > 0) {
    // get ratio and error
    if(data->indx == 0) 
      getDataRatio(&(data->pepRatio[0]), &(data->pepRatio[1]), 
		   &(data->pepH2LRatio[0]), &(data->pepH2LRatio[1]), 
		   _ASAPRATIO_CONFL_,
		   value, error, weight, outliers, count, 1);
    else
      getDataRatio(&(data->pepRatio[0]), &(data->pepRatio[1]), 
		   &(data->pepH2LRatio[0]), &(data->pepH2LRatio[1]), 
		   _ASAPRATIO_CONFL_,
		   value, error, weight, outliers, count, 0);

    // set pepCount
    count = 0;
    for (i = startZ; i < endZ; ++i){
      if (data->pkCount[i] == 1) {  
	data->pkCount[i] = 1 - outliers[count];            
	++count;
      }
    }

    // get time and error
    for (i = 0; i < 2; ++i) {
      area = 0.;
      qState = data->chrg-1;      
      for (j = startZ; j < endZ; ++j){
	if(data->peaks[j][i].indx == 2
	   && data->pkCount[j] == 1
	   && data->peaks[j][i].area[0] > area) {
	  qState = j;
	  area = data->peaks[j][i].area[0];
	}
      }
      data->pepTime[i][0] = data->peaks[qState][i].time[0];
      data->pepTime[i][1] = data->peaks[qState][i].time[1];
    }

    // get area
    data->pepArea = 0.;
    for (i = startZ; i < endZ; ++i){
      if (data->pkCount[i] == 1
	  && data->pepArea 
	  < data->peaks[i][0].area[0] + data->peaks[i][1].area[0])
	data->pepArea 
	  = data->peaks[i][0].area[0] + data->peaks[i][1].area[0];
    }
  } //  if(count > 0) 
  else {
    data->pepRatio[0] = -2.;
    data->pepRatio[1] = 0.;
    data->pepTime[0][0] = -1.;
    data->pepTime[0][1] = -1.;
    data->pepTime[1][0] = -1.;
    data->pepTime[1][1] = -1.;
    data->pepArea = -1.;
  }


  // generate .pngFiles
  if(cgiIndx == 1) {
    for (i = startZ; i < endZ; ++i){
      if(data->peaks[i][0].indx > 0
	 || data->peaks[i][1].indx > 0) {
	plotPepPeak(pngFileBase, *data, lcSpect, i+1);
      }
    }
  }


  // convert valleys into scan numbers in .dat file
  for (i = startZ; i < endZ; ++i){
    for (j = 0; j < 2; ++j) {
      if(data->peaks[i][j].indx >= 0) {
	data->peaks[i][j].peak
	  = lcSpect.expScanNums[data->peaks[i][j].peak];
	data->peaks[i][j].valley[0] 
	  = lcSpect.expScanNums[data->peaks[i][j].valley[0]];
	data->peaks[i][j].valley[1] 
	  = lcSpect.expScanNums[data->peaks[i][j].valley[1]];
      }
    }
  }

  if(data->indx == 0)
    data->indx = 1;

  return;
}


// This function gets rid of any space at the beginning or end of a string.
void ASAPRatioPeptideParser::getRidOfSpace(char *string) 
{
  int lngth;
  int i;

  // start
  for (i = 0; string[i]; ++i) {
    if (isspace(string[i]) == 0) // not a space
      break;
  }
  // BSP: The behavior of strcpy is undefined if the source and destination strings overlap.
  // strcpy(string, string+i);
  memmove(string,string+i,strlen(string+i)+1);

  // end
  lngth = (int)strlen(string);
  for (i = lngth-1; i >= 0; --i) {
    if (isspace(string[i]) == 0) // not a space
      break;
  }
  string[i+1] = '\0';

  return;
}

// This function collects LC spectra for a pepDataStrct.  On success, it returns 1; on failure, 0.
int ASAPRatioPeptideParser::getPeakSpectra(lcSpectStrct *lcSpect, pepDataStrct *data, 
					   double timeWd, char *xmlFile) {
  return getPeakSpectra(lcSpect, data, timeWd, xmlFile, 10, 0.5);
}

int ASAPRatioPeptideParser::getPeakSpectra(lcSpectStrct *lcSpect, pepDataStrct *data, 
					   double timeWd, char *xmlFile,int smoothItrs,double smoothRTwindow)
{
  //xmlIndxStrct *xmlIndx;
  const struct ScanHeaderStruct* header;
  spectStrct *tmpSpect = NULL;
  double scanTime;
  long startScan, endScan;
  double mzArray[_ASAPRATIO_ISONUM_];
  double dMz;
  double mass[2];
  int lcScan, qState, size;
  double fitTime = smoothRTwindow;

  int repeat = smoothItrs;

  long scan;

  int i, j, k;

  // initial steps
  if((xmlIndx_ = getXmlIndx(xmlFile)) == NULL){
    printf("Error In Reading %s.\n", xmlFile);
    fflush(stdout);
    return 0;
  }

  if(data->scan > xmlIndx_->totScan) {
    data->scan = xmlIndx_->totScan / 2;
  }

  // get scanTime
  header = cached_ramp_readHeader(xmlIndx_->mzfile, xmlIndx_->scanIndx[data->scan]);
  scanTime = header->retentionTime/60.;

  // find scan range
  startScan = (long)getScanNum(scanTime-timeWd, xmlIndx_);
  endScan = (long)getScanNum(scanTime+timeWd, xmlIndx_);

  // find mass range
  mass[0] = data->msLight;
  mass[1] = data->msHeavy;

  // get spectra
  int startZ = 0;
  int endZ = _ASAPRATIO_MXQ_;
  if (pInput_.bQuantCIDChrgOnly) {
    startZ = pInput_.iChargeState-1;
    endZ = pInput_.iChargeState;
  }
    
  for (i = startZ; i < endZ; ++i){
    qState = i + 1;
    for (j = 0; j < 2; ++j) {

      if(data->indx != 0 && data->peaks[i][j].indx == -1)
	continue;

      // get mzArray
      for (k = 0; k < _ASAPRATIO_ISONUM_; ++k) { 
	mzArray[k] = (mass[j]+k*_ASAPRATIO_M_)/qState + _ASAPRATIO_HM_*(qState-1)/(qState);
      }
      dMz = mzBound_ < _ASAPRATIO_M_/qState/2. ? 
		       mzBound_ : _ASAPRATIO_M_/qState/2.;

      if((tmpSpect = getCombLCSpect(startScan, endScan, mzArray, 
				    _ASAPRATIO_ISONUM_, dMz, xmlIndx_)) == NULL) {
	//cout << "null lc spect" << endl; 
	data->peaks[i][j].indx = -1;
	//exit(1);
	continue;
      }
      else if(data->indx == 0) {
	//cout << "NON_NULL speect" << endl;
	data->peaks[i][j].indx = 0;
      }
      //  cout << "got fit" << endl;

      // get fitting spectrum
      lcSpect->rawSpectrum[i][j] = *tmpSpect;
      lcSpect->fitSpectrum[i][j] = *tmpSpect;
      smoothSpectFlx(&(lcSpect->fitSpectrum[i][j]), fitTime, 10);

    }//for (j = 0; j < 2; ++j) {
  } //  for (i = 0; i < _ASAPRATIO_MXQ_; ++i){

  // get expScanNums in mzXML file and cidScan in spectrum
  size = endScan - startScan + 1;
  lcSpect->expScanNums = (long *) calloc(size, sizeof(long));

  lcScan = 0;
  for (scan = startScan; scan <= endScan; ++scan) {
    header = cached_ramp_readHeader(xmlIndx_->mzfile, xmlIndx_->scanIndx[scan]);
    if (header->msLevel == 1) {
      lcSpect->expScanNums[lcScan] = scan;
      ++lcScan;
    }
  }

#ifndef DEBUG_NOFREE
  //freeXmlIndx(xmlIndx_);
  //free(xmlIndx_); // added by AK....taken away by AK 7.15.04
#endif

  return 1;
}


// This function opens an xmlFile, gets scanIndx, and total scan number.
xmlIndxStrct *ASAPRatioPeptideParser::createXmlIndx(const char *xmlFile)
{
  if (xmlIndx_ == NULL) {
    ramp_fileoffset_t indexOffset;
    xmlIndx_ = (xmlIndxStrct *) calloc(1, sizeof(xmlIndxStrct));

    // open file
    if ((xmlIndx_->mzfile = cached_ramp_rampOpenFile(xmlFile)) == NULL) {
      freeXmlIndx();
      return NULL;
    }

    // Read the offset of the index
    indexOffset = cached_ramp_getIndexOffset(xmlIndx_->mzfile);

    // Read the scan index into a vector, get LastScan
    xmlIndx_->scanIndx = cached_ramp_readIndex(xmlIndx_->mzfile, indexOffset, &(xmlIndx_->totScan));
  }
  return xmlIndx_;
}


xmlIndxStrct *ASAPRatioPeptideParser::getXmlIndx(const char *xmlfile)
{
  xmlIndxStrct *result;
  if (-1 == m_XMLfile_state) {
    // first pass
    result = createXmlIndx(xmlfile);
    if (result) { // basename-derived file opened OK
      m_XMLfile_state = 1; // normal operation
    } else {
      m_XMLfile_state = 0; // custom operation
    }
  } 
  if (!m_XMLfile_state) { // deriving filenames from spectrum names      
    bool bChange = false;
    // construct a plausible mzxml filename from spectrum name
    int len = (int)(strlen(xmlfile)+strlen(pInput_.szSpectrumName));
    char *tryName = (char *)malloc(len);
    strcpy(tryName,xmlfile);
    char *slash=findRightmostPathSeperator(tryName);
    if (!slash) {
      slash = tryName-1;
    }
    strcpy(slash+1,pInput_.szSpectrumName);
    char *dot = strchr(slash+1,'.');
    if (dot) {
      *dot = 0;
    }
    rampConstructInputFileName(tryName,len,tryName); // .mzXML or .mzData
    if (strcmp(tryName,m_lastMZXML.c_str())) {
      // switching files
      m_lastMZXML = tryName;
      bChange = true;
      if (verbose_) {
	cerr << "mzXML file \"" << xmlfile << "\" not found, using filename \"" <<
	  tryName << "\" constructed from spectrum name \"" <<
	  pInput_.szSpectrumName << "\"" << endl;           	        
      }
    }
    free(tryName);

    if (bChange) {
      freeXmlIndx();
      ramp_fileoffset_t indexOffset;

      result = (xmlIndxStrct *) calloc(1, sizeof(xmlIndxStrct));
      result->mzfile = cached_ramp_rampOpenFile(m_lastMZXML.c_str());

      if (!result->mzfile) {
	cerr << "Could not open input file " << xmlfile << endl;
	fflush(stdout);
	free(result);
	return NULL;
      }

      // Read the offset of the index
      indexOffset = cached_ramp_getIndexOffset(result->mzfile);

      // Read the scan index into a vector, get LastScan
      result->scanIndx = cached_ramp_readIndex(result->mzfile, indexOffset, &(result->totScan));
    } else {
      result = xmlIndx_; // use same as last time
    }
  } else { // standard operation
    result = createXmlIndx(xmlfile);
  }
  return result;
}


// This function frees a xmlIndxStrct.
void ASAPRatioPeptideParser::freeXmlIndx()
{
#ifndef DEBUG_NOFREE
  if (xmlIndx_) {
    free(xmlIndx_->scanIndx);
    cached_ramp_rampCloseFile(xmlIndx_->mzfile);
    free(xmlIndx_);
    xmlIndx_ = NULL;
  }
#endif

  return;
}


// This function gets LC spectrum by summing up ion intensities in MS spectra within a series of mass windows.
spectStrct *ASAPRatioPeptideParser::getCombLCSpect(int firstScan, int lastScan, 
						   double *mzArray, int arraySize, double dMz, 
						   xmlIndxStrct *xmlIndx)
{
  //  spectStrct *spectrum;
  const struct ScanHeaderStruct* scanHeader;
  int iSt, iNd;
  int size;
  int dataNum;
  const RAMPREAL *pPeaks;
  int peakCount;

  int i, j;

  // check on scan
  iSt = firstScan > 1 ? firstScan : 1;
  iNd = lastScan < xmlIndx->totScan ? lastScan : xmlIndx->totScan;

  if(iSt > iNd)
    return NULL;

  spectrum_->size = _MZXML_READER_MXSCANNUM_;

  size = 0;
  dataNum = 0;
  for (i = iSt; i <= iNd; ++i) {
    scanHeader = cached_ramp_readHeader(xmlIndx->mzfile, xmlIndx->scanIndx[i]);

    if (scanHeader->msLevel == 1 && size < spectrum_->size) { // LC/MS scan
      spectrum_->xval[size] = scanHeader->retentionTime/60.;  // in minute
      spectrum_->yval[size] = 0.;

      peakCount = 0;
 
      pPeaks = cached_ramp_readPeaks_const(xmlIndx->mzfile, xmlIndx->scanIndx[i]);

      if(pPeaks == NULL) {
	return NULL;
      }

      int l = 0;
      int r = scanHeader->peaksCount;
      RAMPREAL fMass;
      while (l < r) {
	int m = (l + r) / 2;
	fMass = pPeaks[m*2];
	if (fMass == mzArray[0]-dMz) {
	  l = r = m;
	  break;
	}
	else if (fMass > mzArray[0]-dMz) {
	  if (r == m) {
	    r = m - 1;
	  }
	  else {
	    r = m ;
	  }
	}
	else {
	  if (l == m) {
	    l = m + 1;
	  }
	  else {
	    l = m ;
	  }
	}
      }

      int peakIx = r * 2; //We should have the smallest m/z in the spectrum matching to the target range

      if(peakIx <= scanHeader->peaksCount*2 && pPeaks[peakIx] >= mzArray[0]-dMz) {
	// Sum the signal over all peaks falling within the target range around each isotope
	j = 0;
	while (j < arraySize && pPeaks[peakIx] != -1 && pPeaks[peakIx] <= mzArray[arraySize-1]+dMz ) {
	  if (pPeaks[peakIx] <= mzArray[j]+dMz
	      && pPeaks[peakIx] >= mzArray[j]-dMz) {
	    spectrum_->yval[size] += pPeaks[peakIx+1];
	    peakIx += 2;
	  }
	  else if (pPeaks[peakIx] > mzArray[j]+dMz) {
	    j++;
	  }
	  else if (pPeaks[peakIx] < mzArray[j]-dMz) {
	    peakIx += 2;
	  }
	  else {
	    peakIx += 2; // Should never get here
	  }
	}
      }

      if(size < spectrum_->size && spectrum_->yval[size] > 0.)
	++dataNum;
      ++size;
    } //     if (scanHeader->msLevel == 1)  // LC/MS scan
  } // for (i = iSt; i <= iNd; ++i) 

  if(dataNum < 1) {
    spectrum_->size = 0;
    return NULL;
  }
  else {
    spectrum_->size = size;
    return spectrum_;
  }
}


// For a given time, this function returns the closest scan number with fabs(scanTime-time) minimium.
int ASAPRatioPeptideParser::getScanNum(double time, xmlIndxStrct *xmlIndx)
{
  const struct ScanHeaderStruct* header;
  int scanNum;
  int bnd[2];
  double scanTime[2];
  double tmpTime;

  // 1st spect
  bnd[0] = 1;
  header = cached_ramp_readHeader(xmlIndx->mzfile, xmlIndx->scanIndx[bnd[0]]);
  scanTime[0] = header->retentionTime/60.;

  if(scanTime[0] >= time)
    return bnd[0];

  // last spect
  bnd[1] = xmlIndx->totScan;
  header = cached_ramp_readHeader(xmlIndx->mzfile, xmlIndx->scanIndx[bnd[1]]);
  scanTime[1] = header->retentionTime/60.;

  if(scanTime[1] <= time)
    return bnd[1];

  // search for scanNum
  while(bnd[1]-bnd[0] > 1) {
    scanNum = (bnd[0]+bnd[1])/2;
    header = cached_ramp_readHeader(xmlIndx->mzfile, xmlIndx->scanIndx[scanNum]);
    tmpTime = header->retentionTime/60.;
    if(tmpTime == time)
      return scanNum;
    else if(tmpTime < time) {
      bnd[0] = scanNum;
      scanTime[0] = tmpTime;
    }
    else {
      bnd[1] = scanNum;
      scanTime[1] = tmpTime;
    }

  }

  if(fabs(scanTime[0]-time) <= fabs(scanTime[1]-time))
    return bnd[0];
  else
    return bnd[1];
}


// This function gets a lcPeakStrct for a spectrum
void ASAPRatioPeptideParser::getLCPeakStrct(lcPeakStrct *peak, const spectStrct &rawSpectrum, 
					    const spectStrct &fitSpectrum, int scan, int scanRange)
{
  int ranges[2];

  // get peakScan and valley
  getLCSpectPeakAndValleys(&(peak->peak), peak->valley, fitSpectrum, scan, 0.); 

  // get background
  ranges[0] = peak->valley[0] - scanRange;
  ranges[0] = ranges[0] > 0 ? ranges[0] : 0;
  ranges[1] = peak->valley[1] + scanRange;
  ranges[1] = ranges[1] < rawSpectrum.size-1 ? ranges[1] : rawSpectrum.size-1;
  peak->bckgrnd = getLCSpectBackground(rawSpectrum, fitSpectrum, ranges, peak->valley);

  // reset valleys
  getLCSpectPeakAndValleys(&(peak->peak), peak->valley, fitSpectrum, 
			   peak->peak, peak->bckgrnd); 

  // get area and time
  peak->indx = getLCSpectAreaAndTime(peak->area, peak->time, 
				     rawSpectrum, fitSpectrum, 
				     peak->peak, peak->valley, 
				     peak->bckgrnd);

  if (last_valley_[0] == -1 && pInput_.bUseFixedScanRange) {
    last_valley_[0] = peak->valley[0];
  }
  if (last_valley_[1] == -1 && pInput_.bUseFixedScanRange) {
    last_valley_[1] = peak->valley[1];
  }
  return;
}


// This function finds the peak and valleys of a LC spectrum.
void ASAPRatioPeptideParser::getLCSpectPeakAndValleys(int *peakScan, int *valleyScan, 
						      const spectStrct &spectrum, 
						      int scan, double background) 
{
  int i;

  // find peak scan
  *peakScan = spectPeak(spectrum, scan, -1);

  // find valley positions
  if (last_valley_[0] == -1 || !pInput_.bUseFixedScanRange) {
    valleyScan[0] = findNextValley(spectrum, *peakScan, -1);
  }
  else {
    valleyScan[0] = last_valley_[0];
  }

  if (valleyScan[0] == -1) { // if no valley found, set to boundary
    valleyScan[0] = 0;
  }

  if (last_valley_[1] == -1 || !pInput_.bUseFixedScanRange) {
    valleyScan[1] = findNextValley(spectrum, *peakScan, 1);
  }
  else {
    valleyScan[1] = last_valley_[1];
  }

  if (valleyScan[1] == -1) { // if no valley found, set to boundary
    valleyScan[1] = spectrum.size-1;
  }

  // find positions at background level
  if (last_valley_[1] == -1 && last_valley_[0] == -1) {
    for (i = *peakScan; i >= valleyScan[0]; --i){
      if(spectrum.yval[i] < background) {
	break;
      }
    }
    valleyScan[0] = i + 1;
    for (i = *peakScan; i <= valleyScan[1]; ++i){
      if(spectrum.yval[i] < background)
	break;
    }
    valleyScan[1] = i - 1;
  }

  // check valley
  if(valleyScan[0] > valleyScan[1]) {
    valleyScan[0] = *peakScan;
    valleyScan[1] = *peakScan;
  }

  return;
}


// This function finds the peak at "position" in a spectrum.  
//  If "position" is a valley, then find the left peak (when 
//  "direction = -1") or the right one (when "direction = 1").
int ASAPRatioPeptideParser::spectPeak(const spectStrct &spectrum, int position, int direction)
{
  int right, left;
  int i;
  // determine slope on the right hand side 
  i = 1;
  while(position+i < spectrum.size &&
	spectrum.yval[position] == spectrum.yval[position+i]){
    ++i;
  }
  if (position+i >= spectrum.size) {
    right = 0;
  }
  else if (spectrum.yval[position] < spectrum.yval[position+i]){
    right = 1;
  }
  else
    right = -1;

  // determine slope on the left hand side 
  i = -1;
  while(position+i >= 0 &&
	spectrum.yval[position] == spectrum.yval[position+i]){
    --i;
  }
  if (position+i < 0) {
    left = 0;
  }
  else if (spectrum.yval[position] < spectrum.yval[position+i]){
    left = 1;
  }
  else
    left = -1;

  // determine "direction"
  if((right == -1 && left == -1) || (right == 0 && left == 0)) {
    return position;
  }
  else if (right == -1 && left == 0) {
    return 0;
  }
  else if (right == 0 && left == -1) {
    return spectrum.size-1;
  }
  else if (right <= 0 && left == 1) {
    direction = -1;
  }
  else if (right == 1 && left <= 0) {
    direction = 1;
  }

  // search for peak
  while (position+direction >= 0 && position+direction < spectrum.size
	 && spectrum.yval[position] <= spectrum.yval[position+direction]) {
    position += direction;
  }

  // return index for peak
  if (position+direction < 0) 
    return 0;
  else if (position+direction >= spectrum.size) 
    return spectrum.size-1;
  else 
    return position;
}


// This function finds the next valley from left (when "direction = -1") 
// or right (when "direction = 1"), starting at "position" in a spectrum.  
int ASAPRatioPeptideParser::findNextValley(const spectStrct &spectrum, int position, int direction)
{
  position += direction;

  // search for valley
  while (position-direction >= 0 && position-direction < spectrum.size 
	 && position+direction >= 0 && position+direction < spectrum.size 
	 && (spectrum.yval[position] > spectrum.yval[position-direction] ||
	     spectrum.yval[position] >= spectrum.yval[position+direction])) {
    position += direction;
  }

  // return index for valley
  if (position-direction < 0 || position-direction >= spectrum.size 
      || position+direction < 0 || position+direction >= spectrum.size) {
    return -1;
  } 
  else if(spectrum.yval[position] <= spectrum.yval[position-direction] 
	  && spectrum.yval[position] < spectrum.yval[position+direction]){
    return position;
  }
  else {
    return -1;
  }
}


// This function gets the background of a LC spectrum.
double ASAPRatioPeptideParser::getLCSpectBackground(const spectStrct &rawSpectrum, 
						    const spectStrct &fitSpectrum, 
						    int *ranges, int *valleyScan)
{
  double noise, background, ave;
  int count, count2;
  int i;

  // check on ranges
  ranges[0] = ranges[0] > 0 ? ranges[0] : 0;
  ranges[1] = ranges[1] < rawSpectrum.size-1 ?
			  ranges[1] : rawSpectrum.size-1;

  // get noise level
  noise = 0.;
  count = 0;
  for (i = ranges[0]; i <= ranges[1]; ++i){
    if (rawSpectrum.yval[i] > 0.) {
      noise += (rawSpectrum.yval[i]-fitSpectrum.yval[i])
	*(rawSpectrum.yval[i]-fitSpectrum.yval[i]);
      ++count;
    }
  }
  if(count > 0)
    noise = sqrt(noise/count);
  else
    noise = 1.;

  // get background
  background = 0.;
  ave = noise;
  while(fabs(background-ave) > 0.01*noise) {
    background = ave;
    ave = 0.;
    count = 0;
    count2 = 0;
    for (i = ranges[0]; i <= ranges[1]; ++i){
      if(i < valleyScan[0] || i > valleyScan[1]) {
	if (rawSpectrum.yval[i] > 0.
	    && rawSpectrum.yval[i] < background + noise) {
	  ave += rawSpectrum.yval[i];
	  ++count;
	}
	else if (rawSpectrum.yval[i] >= background + noise) 
	  ++count2;
      }
    }
    if(count > 0)
      ave /= count;
    else if(count2 > 0)
      ave = background + 0.2*noise;
    else
      ave = 0.;
  }

  if (pInput_.bZeroAllBackGrnd) {
    background = 0.;
  }
  else {
    background = ave;
  }
  return background;
}


// This function finds the area and time of a LC peak. It returns 2 if the peak passes certain tests, 1 if not.
int ASAPRatioPeptideParser::getLCSpectAreaAndTime(double *area, double *peakTime, 
						  const spectStrct &rawSpectrum, 
						  const spectStrct &fitSpectrum, 
						  int pkScan, int *valleyScan, double background)
{
  int areaIndx = 2;
  double rArea, fArea;
  double rVal, fVal;
  double mxVal, ave, areaErr;
  int pkTime1, pkTime2;
  int count;
  int i;

  // get mxVal
  mxVal = 0.;
  for (i = valleyScan[0]; i <= valleyScan[1]; ++i) {
    mxVal = mxVal > fitSpectrum.yval[i] ?
      mxVal : fitSpectrum.yval[i];
  }

  // get peak area and error
  rArea = 0.;
  fArea = 0.;
  areaErr = 0.;
  int imax = Min(valleyScan[1],rawSpectrum.size-2); // don't run off end BSP Dec 2005
  for (i = valleyScan[0]; i <= imax; ++i) {
    rVal = 0.5*(rawSpectrum.yval[i+1]+rawSpectrum.yval[i]) - background;
    rVal = rVal > 0. ? rVal : 0.;
    rArea += rVal*(rawSpectrum.xval[i+1]-rawSpectrum.xval[i]);
    fVal = 0.5*(fitSpectrum.yval[i+1]+fitSpectrum.yval[i]) - background;
    fVal = fVal > 0. ? fVal : 0.;
    fArea += fVal*(fitSpectrum.xval[i+1]-fitSpectrum.xval[i]);
    areaErr += 0.5*(rVal-fVal)*(rVal-fVal)
      *(fitSpectrum.xval[i+1]-fitSpectrum.xval[i])
      *(fitSpectrum.xval[i+1]-fitSpectrum.xval[i]);
  }

  // check error
  if (data_.areaFlag == 1) {
    area[0] = rArea;
  }
  else if (data_.areaFlag == 2) {
    area[0] = fArea;
  }
  else {
    area[0] = (fArea+rArea)/2.;
  }
  area[1] = sqrt(areaErr);
  if(rArea <= 0. || fArea <= 0. ) {
    areaIndx = 1;
  }
  else if (pInput_.bQuantHighBackGrnd) {
    areaIndx = 2;
  }
  else if(area[0] < area[1] || mxVal < 2.*background) {
    areaIndx = 1;
  }

  // find peakTime 
  peakTime[0] = fitSpectrum.xval[pkScan];

  // find error in peakTime
  pkTime1 = pkScan;
  while(pkTime1 >= valleyScan[0] // width at half peak height: left
	&& fitSpectrum.yval[pkTime1]-background 
	> 0.5*(fitSpectrum.yval[pkScan]-background)) {
    --pkTime1;
  }
  pkTime2 = pkScan;
  while(pkTime2 < valleyScan[1] // width at half peak height: right
	&& fitSpectrum.yval[pkTime2]-background 
	> 0.5*(fitSpectrum.yval[pkScan]-background)) {
    ++pkTime2;
  }
  pkTime1=Max(0,pkTime1);  // don't run off end BSP Dec 2005
  pkTime2=Min(pkTime2,fitSpectrum.size-1);  // don't run off end BSP Dec 2005
  if(pkTime2 <= pkTime1)
    peakTime[1] = 0.;
  else {
    peakTime[1] = 0.5*(fitSpectrum.xval[pkTime2]-fitSpectrum.xval[pkTime1]);
  }

  // double check on bad data
  if(area[0] > 0. && mxVal < 5.*background) {
    count = 0;
    for (i = pkTime1; i <= pkTime2; ++i){
      if(rawSpectrum.yval[i] < background) {
	++count;
      }
    }
    if(count > (pkTime2-pkTime1+1)/4) {
      count = 0;
      ave = fitSpectrum.yval[pkTime1] < fitSpectrum.yval[pkTime2]?
	fitSpectrum.yval[pkTime1] : fitSpectrum.yval[pkTime2];
      for (i = pkTime1; i <= pkTime2; ++i)
	if(rawSpectrum.yval[i] > ave)
	  ++count;
      if (pInput_.bQuantHighBackGrnd) {
	areaIndx = 2;
      }
      else if(count < (pkTime2-pkTime1+1)/2) {
	areaIndx = 1;
      }
    }
  }

  return areaIndx;
}


// This function generates pngFiles for LC/MS spectra.
void ASAPRatioPeptideParser::plotPepPeak(const char *pngFileBaseIn, const pepDataStrct &data, 
					 const lcSpectStrct &lcSpect, int qState)
{
  FILE *file, *file2;

  char *pngFile;
  char *gnuFile;
  char *datFiles[3][2];
  char *gnuCommand;
  int plotRange[2];
  long expRange[2];
  int cidIndx = data.cidIndx;
  int size = lcSpect.rawSpectrum[data.chrg-1][cidIndx].size;
  long cidTime = data.scan;
  double mxPeak, yMx, yOrder;
  int lcScan;
  
  int i, j;

  char *pngFileBase = strdup(pngFileBaseIn);
  fixPath(pngFileBase,0); // pretty up the path seperators etc
  int startZ = 0;
  int endZ = _ASAPRATIO_MXQ_;
  if (pInput_.bQuantCIDChrgOnly) {
    startZ = pInput_.iChargeState-1;
    endZ = pInput_.iChargeState;
  }
  // get timeRange
  plotRange[0] = data.peaks[data.chrg-1][cidIndx].valley[0];
  plotRange[1] = data.peaks[data.chrg-1][cidIndx].valley[1];
  for (i = startZ; i < endZ; ++i){
    for (j = 0; j < 2; ++j) {
      if(data.peaks[i][j].indx >= 0) {
	if(plotRange[0] > data.peaks[i][j].valley[0])
	  plotRange[0] = data.peaks[i][j].valley[0];
	if(plotRange[1] < data.peaks[i][j].valley[1])
	  plotRange[1] = data.peaks[i][j].valley[1];
      }
    }
  }
  lcScan = plotRange[0];
  while(lcScan > 0 
	&& lcSpect.expScanNums[lcScan] 
	> lcSpect.expScanNums[plotRange[0]]-_ASAPRATIO_EXPLCRANGE_/5)
    --lcScan;
  plotRange[0] = lcScan;
  lcScan = plotRange[1];
  while(lcScan < size-1
	&& lcSpect.expScanNums[lcScan] 
	< lcSpect.expScanNums[plotRange[1]]+_ASAPRATIO_EXPLCRANGE_/5)
    ++lcScan;
  plotRange[1] = lcScan;

  // get expRange
  expRange[0] = 5*(lcSpect.expScanNums[plotRange[0]]/5);
  expRange[1] = 5*(lcSpect.expScanNums[plotRange[1]]/5+1);

  // get datFiles
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 2; ++j) {
      datFiles[i][j] = (char *) calloc(strlen(pngFileBase)+10, sizeof(char));
      sprintf(datFiles[i][j], "%s_%d_%d_%d.dat", pngFileBase, qState, i, j);
    }
  }

  // get pngFile
  pngFile = (char *) calloc(strlen(pngFileBase)+10, sizeof(char));
  sprintf(pngFile, "%s_%d.png", pngFileBase, qState);

  // get gnuFile
  gnuFile = (char *) calloc(strlen(pngFileBase)+10, sizeof(char));
  sprintf(gnuFile, "%s_%d.gnu", pngFileBase, qState);

  if((file = fopen(gnuFile, "w")) == NULL) {
    printf("Cannot Write to File \"%s\"!\n", gnuFile);
    fflush(stdout);
    exit(0);
  }
  fprintf(file, "set termnal png\n");
  fprintf(file, "set output \'%s\'\n", pngFile);
  fprintf(file, "set format y \"%%7.2g\"\n");
  fprintf(file, "set key samplen -1\n");
  fprintf(file, "set multiplot\n");
  for (j = 0; j < 2; ++j) { 
    if(data.peaks[qState-1][j].indx >= 0) {
      fprintf(file, "set size 1.,0.5\n");
      if(j == 0) { // light
	fprintf(file, "set origin 0,0.5\n");
      }
      else { // heavy
	fprintf(file, "set origin 0,0\n");
      }
      fprintf(file, "set xrange [%ld:%ld]\n", expRange[0], expRange[1]);
      fprintf(file, "set yrange [0:*]\n");

      // datFile: rawSpectrum
      yMx = 0.;
      mxPeak = 0.;
      if((file2 = fopen(datFiles[0][j], "w")) == NULL) {
	printf("Cannot Write to File \"%s\"!\n", datFiles[0][j]);
	fflush(stdout);
	exit(0);
      }
      for (i = plotRange[0]-10 > 0 ? plotRange[0]-10 : 0;
	   i < plotRange[1]+10 && i < lcSpect.rawSpectrum[qState-1][j].size; 
	   ++i) {
	if(i > data.peaks[qState-1][j].valley[0]
	   && i < data.peaks[qState-1][j].valley[1]
	   && mxPeak < lcSpect.rawSpectrum[qState-1][j].yval[i])
	  mxPeak = lcSpect.rawSpectrum[qState-1][j].yval[i];
	if(i > plotRange[0] && i < plotRange[1] 
	   && yMx < lcSpect.rawSpectrum[qState-1][j].yval[i])
	  yMx = lcSpect.rawSpectrum[qState-1][j].yval[i];
	fprintf(file2, "%ld %f\n", lcSpect.expScanNums[i], 
		lcSpect.rawSpectrum[qState-1][j].yval[i]);
      }
      fclose(file2);
      if(yMx > 4.*mxPeak && mxPeak > 0.) {
	yOrder = pow(10., (int)log10(mxPeak));
	mxPeak = 2.*ceil(mxPeak/yOrder)*yOrder;
	fprintf(file, "set yrange [0:%f]\n", mxPeak);
      }

      // datFile: fitSpectrum
      if((file2 = fopen(datFiles[1][j], "w")) == NULL) {
	printf("Cannot Write to File \"%s\"!\n", datFiles[1][j]);	
	fflush(stdout);
	exit(0);
      }
      for (i = plotRange[0]-10 > 0 ? plotRange[0]-10 : 0;
	   i < plotRange[1]+10 && i < lcSpect.fitSpectrum[qState-1][j].size; 
	   ++i)
	fprintf(file2, "%ld %f\n", lcSpect.expScanNums[i], 
		lcSpect.fitSpectrum[qState-1][j].yval[i]);
      fclose(file2);

      // datFile: peakSpectrum
      if(data.peaks[qState-1][j].indx > 1 
	 && data.pkCount[qState-1] == 1) {
	if((file2 = fopen(datFiles[2][j], "w")) == NULL) {
	  printf("Cannot Write to File \"%s\"!\n", datFiles[2][j]);
	  fflush(stdout);
	  exit(0);
	}
	for (i = data.peaks[qState-1][j].valley[0];
	     i <= data.peaks[qState-1][j].valley[1]; ++i) {
	  fprintf(file2, "%ld %f\n", lcSpect.expScanNums[i], 
		  data.peaks[qState-1][j].bckgrnd);
	  if(data.areaFlag == 1) {
	    if (lcSpect.rawSpectrum[qState-1][j].yval[i] 
		> data.peaks[qState-1][j].bckgrnd)
	      fprintf(file2, "%ld %f\n", lcSpect.expScanNums[i], 
		      lcSpect.rawSpectrum[qState-1][j].yval[i]);
	    else
	      fprintf(file2, "%ld %f\n", lcSpect.expScanNums[i], 
		      data.peaks[qState-1][j].bckgrnd);
	  }
	  else if(data.areaFlag == 2) {
	    if (lcSpect.fitSpectrum[qState-1][j].yval[i] 
		> data.peaks[qState-1][j].bckgrnd)
	      fprintf(file2, "%ld %f\n", lcSpect.expScanNums[i], 
		      lcSpect.fitSpectrum[qState-1][j].yval[i]);
	    else
	      fprintf(file2, "%ld %f\n", lcSpect.expScanNums[i], 
		      data.peaks[qState-1][j].bckgrnd);
	  }
	  else {
	    if (0.5*(lcSpect.rawSpectrum[qState-1][j].yval[i]
		     +lcSpect.fitSpectrum[qState-1][j].yval[i])
		> data.peaks[qState-1][j].bckgrnd) 
	      fprintf(file2, "%ld %f\n", lcSpect.expScanNums[i], 
		      0.5*(lcSpect.rawSpectrum[qState-1][j].yval[i]
			   +lcSpect.fitSpectrum[qState-1][j].yval[i]));
	    else
	      fprintf(file2, "%ld %f\n", lcSpect.expScanNums[i], 
		      data.peaks[qState-1][j].bckgrnd);
	  }
	}
	fclose(file2);
      } // if(data.peaks[qState-1][j].area[0] > 0.) 

      // cid indx
      if(qState == data.chrg && j == cidIndx) {
	fprintf(file, "set arrow from  %ld, graph -0.2 to %ld, 0 lt 1 lw 2\n",
		cidTime, cidTime);
      }
      else 
	fprintf(file, "set noarrow\n");

      // peakSpectrum
      if(data.peaks[qState-1][j].indx > 1
	 && data.pkCount[qState-1] == 1) 
	fprintf(file, "pl '%s't \"\"w l lt 2 lw 2, ", datFiles[2][j]);
      else 
      	fprintf(file, "pl ");

      // rawSpectrum 
      if(j == 0) 
	fprintf(file, "'%s't \"light +%d\"w l lt 1 lw 2, ",
		datFiles[0][j], qState);
      else
	fprintf(file, "'%s't \"heavy +%d\"w l lt 1 lw 2, ",
		datFiles[0][j], qState);

      // fitSpectrum
      fprintf(file, "'%s't \"\"w l lt 3 lw 2, ", datFiles[1][j]);

      // back Ground
      fprintf(file, "%f t \"\"w l lt 4 lw 2\n", 
	      data.peaks[qState-1][j].bckgrnd);
    } // if(data.peaks[qState-1][j].indx >= 0) {
  } // for (j = 0; j < 2; ++j) 

  fprintf(file, "set nomultipl\n");
  fprintf(file, "unset output\n");
  fprintf(file, "quit\n");  
  fflush(file);
  fclose(file);

  // get image
  gnuCommand = (char *) calloc(strlen(gnuFile)+10, sizeof(char));
  sprintf(gnuCommand, "%s %s", GNUPLOT_BINARY, gnuFile);
  verified_unlink(pngFile);
  verified_system(gnuCommand);
#ifndef DEBUG_NOFREE
  free(gnuCommand);
#endif

  //Sleep Max 30 secs
  FILE* test;
  int count = 0;
  while (count++ < 30 && (test = fopen(pngFile,"r"))==NULL) {
    sleep(1);
  }

  if (count >= 30) {
    fprintf(stderr,"WARNING: Max Timeout reached waiting for gnuplot ... ");
  }
  else if ( test!=NULL ) {
    fclose(test);
  }
		
  // clear
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 2; ++j) {
      verified_unlink(datFiles[i][j]);
#ifndef DEBUG_NOFREE
      free(datFiles[i][j]);
#endif
    }
  }

  verified_unlink(gnuFile);
#ifndef DEBUG_NOFREE
  free(gnuFile);
  free(pngFile);
  free(pngFileBase);
#endif
  return;
}


Array<Tag*>* ASAPRatioPeptideParser::generateXML(int index, Boolean cover, Boolean data) {
  Array<Tag*>* output = new Array<Tag*>;

  if(data_.indx == -1) {
    //cout << "returning -1" << endl;
    return output; // done
  }

  Tag* next;
  char text[100];
  const char* lcpeakNames[] = {"asapratio_lc_lightpeak", "asapratio_lc_heavypeak"};

  if(cover) {
    if(data)
      next = new Tag("asapratio_result", True, False);
    else
      next = new Tag("asapratio_result", True, True);
    // fill it up
    if(index) {
      sprintf(text, "%d", index);
      next->setAttributeValue("index", text);
    }
    if(data_.pepRatio[0] == -2) {
      next->setAttributeValue("mean", "-1");
      next->setAttributeValue("error", "-1");
    }    
    else if(data_.pepRatio[0] == -1) {
      next->setAttributeValue("mean", "999.");
      next->setAttributeValue("error", "0.00");
    }
    else {
      sprintf(text, "%0.2f", data_.pepRatio[0]);
      next->setAttributeValue("mean", text);
      sprintf(text, "%0.2f", data_.pepRatio[1]);
      next->setAttributeValue("error", text);
    }

    // now the heavy2light.....
    if(data_.pepRatio[0] == -2) {
      next->setAttributeValue("heavy2light_mean", "-1");
      next->setAttributeValue("heavy2light_error", "-1");
    }
    else if(data_.pepRatio[0] == 0.0) {
      next->setAttributeValue("heavy2light_mean", "999.");
      next->setAttributeValue("heavy2light_error", "0.00");
    }
    else if(data_.pepRatio[0] == -1) {
      next->setAttributeValue("heavy2light_mean", "0.00");
      next->setAttributeValue("heavy2light_error", "0.00");
    }
    else {
      sprintf(text, "%0.2f", data_.pepH2LRatio[0]);
      next->setAttributeValue("heavy2light_mean", text);
      sprintf(text, "%0.2f", data_.pepH2LRatio[1]);
      next->setAttributeValue("heavy2light_error", text);
    }
#ifdef USE_STD_MODS    
    // now the lightmass
    //    sprintf(text, "%0.2f", light_mass_);
    //    next->setAttributeValue("light_mass", text);

    //    // now get the quant label positions in the peptide
    //    next->setAttributeValue("quant_label_posns", quant_label_positions_);
#endif
    if (verbose_) {
      next->write(cout);
    }
    output->insertAtEnd(next);
  }
  int startZ = 0;
  int endZ = _ASAPRATIO_MXQ_;
  if (pInput_.bQuantCIDChrgOnly) {
    startZ = pInput_.iChargeState-1;
    endZ = pInput_.iChargeState;
  }
  if(data) {
    next = new Tag("asapratio_peptide_data", True, False);
    if(index) {
      sprintf(text, "%d", index);
      next->setAttributeValue("index", text);
    }
    sprintf(text, "%d", data_.indx);
    next->setAttributeValue("status", text);
    sprintf(text, "%d", data_.cidIndx);
    next->setAttributeValue("cidIndex", text);
    sprintf(text, "%0.4f", data_.msLight);
    next->setAttributeValue("light_mass", text);
    sprintf(text, "%0.4f", data_.msHeavy);
    next->setAttributeValue("heavy_mass", text);
    sprintf(text, "%d", data_.areaFlag);
    next->setAttributeValue("area_flag", text);
    output->insertAtEnd(next);
    for(int k = startZ; k < endZ; k++) {
      //cout << "K: " << k;
      next = new Tag("asapratio_contribution", True, False);
      sprintf(text, "%0.4f", data_.pkRatio[k]);
      next->setAttributeValue("ratio", text);
      sprintf(text, "%0.4f", data_.pkError[k]);
      next->setAttributeValue("error", text);
      sprintf(text, "%d", k+1);
      next->setAttributeValue("charge", text);
      sprintf(text, "%d", data_.pkCount[k]);
      next->setAttributeValue("use", text);
      output->insertAtEnd(next);
      for(int j = 0; j < 2; j++) {
	//cout << " J: " << j << endl;
	next = new Tag(lcpeakNames[j], True, True);
	sprintf(text, "%d", data_.peaks[k][j].indx);
	next->setAttributeValue("status", text);
	sprintf(text, "%d", data_.peaks[k][j].valley[0]);
	next->setAttributeValue("left_valley", text);
	sprintf(text, "%d", data_.peaks[k][j].valley[1]);
	next->setAttributeValue("right_valley", text);
	sprintf(text, "%0.2e", data_.peaks[k][j].bckgrnd);
	next->setAttributeValue("background", text);
	sprintf(text, "%0.2e", data_.peaks[k][j].area[0]);
	next->setAttributeValue("area", text);
	sprintf(text, "%0.2e", data_.peaks[k][j].area[1]);
	next->setAttributeValue("area_error", text);
	sprintf(text, "%0.4f", data_.peaks[k][j].time[0]);
	next->setAttributeValue("time", text);
	sprintf(text, "%0.4f", data_.peaks[k][j].time[1]);
	next->setAttributeValue("time_width", text);
	sprintf(text, "%d", data_.peaks[k][j].peak);
	next->setAttributeValue("is_heavy", text);
	output->insertAtEnd(next);
      } // light/heavy
      next = new Tag("asapratio_contribution", False, True);
      output->insertAtEnd(next);
    } // next precursor ion charge

    next = new Tag("asapratio_peptide_data", False, True);
    output->insertAtEnd(next);

    if(cover) {
      next = new Tag("asapratio_result", False, True);
      output->insertAtEnd(next);
    } // if also cover
  } // if data

  return output;
}


void ASAPRatioPeptideParser::getRatio() {
#ifdef USE_STD_MODS
  evalPepDataStrct(pInput_.szPeptide, (long)(pInput_.iFirstScan), pInput_.iChargeState, mzXMLfile_,
		   elution_, areaFlag_, pInput_.dPeptideMass, pInput_.dCalcPeptideMass);
#else
  if(modAAs_ == NULL) {
    cout << "null modaas" << endl;
    exit(1);
  }

  if(prtnAAs_ == NULL) {
    cout << "null prtnaas" << endl;
    exit(1);
  }
  evalPepDataStrct(pInput_.szPeptide, (long)(pInput_.iFirstScan), pInput_.iChargeState, mzXMLfile_,
		   elution_, areaFlag_, modAAs_, modAANum_, prtnAAs_, prtnAANum_);
#endif
}
