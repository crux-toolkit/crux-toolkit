/*
Program       : ASAPRatioProteinRatioParser
Author        : Andrew Keller <akeller@systemsbiology.org>
                *Xiao-jun Li (xli@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: ASAPRatioProteinRatioParser.cpp 8022 2020-02-12 21:47:21Z mhoopmann $

Computes ASAP ratios and errors for proteins, then overwrites
that information into ProteinProphet XML

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

#include "ASAPRatioProteinRatioParser.h"
#include "Parsers/Parser/TagListComparator.h" // regression test stuff - bpratt Insilicos LLC, Nov 2005
#include "Common/TPPVersion.h" // contains version number, name, revision

static bool verbose = false;

ASAPRatioProteinRatioParser::ASAPRatioProteinRatioParser(const char* protxmlfile, const char *testMode) : Parser("asapratio") { 
  iNumRawData_=0;
  heavy2light_ = False;

  parser_ = NULL;
  pro_ratio_ = NULL;
  testMode_ = testMode?strdup(testMode):NULL; // regression test stuff - bpratt Insilicos LLC, Nov 2005

  init(protxmlfile);
}

ASAPRatioProteinRatioParser::~ASAPRatioProteinRatioParser() { 
  for(int k = 0; k < input_pepxmlfiles_.length(); k++)
    if(input_pepxmlfiles_[k] != NULL)
      delete[] input_pepxmlfiles_[k];

  if(parser_ != NULL)
    delete parser_;
  free(testMode_);
}

void ASAPRatioProteinRatioParser::parse(const char* protxmlfile) {
  //open file and pass along
  char *nextline = new char[line_width_];
  char* data = NULL;
  Tag* tag;
  Boolean heavy2light = False;

  RACI fin(protxmlfile); // can read gzipped xml
  if(! fin) {
    cout << "ASAPRatioProteinRatioParser: error opening " << protxmlfile << endl;
    exit(1);
  }

  double MIN_PEP_WT = 0.5;
  double MIN_PEP_PROB = 0.5;
  Array<Tag*>* tags = NULL;
  Array<UniquePeptide*>* peps = NULL;
  double next_prot_prob = 2.0;
  double MIN_PROB = 0.2;
  Boolean done = False;
  int ratio_tags_written = 0;

  ASAPRatioPvalueParser* pvalue_parser = new ASAPRatioPvalueParser();
  const char* pvalue_name = NULL;
  if(pvalue_parser != NULL)
    pvalue_name = pvalue_parser->getName();
  else {
    pvalue_name = "asapratio_pvalue";
  }

  TagFilter* summary_filter = new TagFilter("analysis_summary");
  summary_filter->enterRequiredAttributeVal("analysis", getName());
  TagFilter* pval_summary_filter = new TagFilter("analysis_summary");
  pval_summary_filter->enterRequiredAttributeVal("analysis", pvalue_name);

  TagFilter* result_filter = new TagFilter("analysis_result");
  result_filter->enterRequiredAttributeVal("analysis", getName());
  TagFilter* pval_result_filter = new TagFilter("analysis_result");
  pval_result_filter->enterRequiredAttributeVal("analysis", pvalue_name);

  delete pvalue_parser;

  TagFilter* ratio_filter = new TagFilter("ASAPRatio");
  TagFilter* seq_filter = new TagFilter("ASAP_Seq");
  TagFilter* peak_filter = new TagFilter("ASAP_Peak");
  TagFilter* dta_filter = new TagFilter("ASAP_Dta");

  Tag* result_start = new Tag("analysis_result", True, False);
  result_start->setAttributeValue("analysis", getName());
  Tag* result_stop = new Tag("analysis_result", False, True);

  Tag* summary_start = new Tag("analysis_summary", True, False);
  summary_start->setAttributeValue("analysis", getName());
  summary_start->setAttributeValue("time", time_);
  summary_start->setAttributeValue("id", "1");
  Tag* summary_stop = new Tag("analysis_summary", False, True);

  Tag* summary = NULL;

  char** altpeps = NULL;

  // construct a tempfile name, possibly in the tmp dir if so configured
  std::string outfile = make_tmpfile_name(protxmlfile);

  ofstream fout(outfile.c_str());
  if(! fout) {
    cerr << "cannot write output to file " << outfile << endl;
    exit(1);
  }

  //
  // regression test stuff - bpratt Insilicos LLC, Nov 2005
  //
  Array<Tag*> test_tags;
  eTagListFilePurpose testType;
  char *testFileName=NULL;
  checkRegressionTestArgs(testMode_,testType);
  if (testType!=NO_TEST) {
    testFileName = constructTagListFilename(protxmlfile, // input file
					    testMode_, // program args
					    "ASAPRatioProteinRatioParser", // program name
					    testType); // user info output
  }

#define RECORD(tag) {(tag)->write(fout);if (testType!=NO_TEST) {test_tags.insertAtEnd(new Tag(*(tag)));}}

  while(fin.getline(nextline, line_width_)) {
    data = strstr(nextline, "<");
    while(data != NULL) {
      tag = new Tag(data);

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
	    //cout << "next file: " << nextfile << endl;
	    input_pepxmlfiles_.insertAtEnd(nextfile);
	    last_i = i+1;

	  }
	} // while


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

	if(! ratio_filter->filter(tag) && ! seq_filter->filter(tag) && ! peak_filter->filter(tag) && ! dta_filter->filter(tag) &&
	   ! result_filter->filter(tag) && ! pval_result_filter->filter(tag)) {
	  tags->insertAtEnd(tag);
	  stored = True;
	}

	if(peps == NULL)
	  peps = new Array<UniquePeptide*>;
	if(tag->isStart() && ! strcmp(tag->getName(), "peptide")) {
	  // check that weight and prob are above minimum.....
	  enterUnique(peps, tag->getAttributeValue("peptide_sequence"), 
		      atof(tag->getAttributeValue("weight")), 
		      atof(tag->getAttributeValue("nsp_adjusted_probability")));
	}

	else if(filter_memory_ && tags != NULL && peps != NULL) {
	  // add the asap pro
	  Tag* protag = NULL;
	  if(peps->length() > 0) {
	    getRatio(peps, MIN_PEP_PROB, MIN_PEP_WT, heavy2light);
	  } // only if have enough peps

	  int k;
	  for(k = 0; k < tags->length(); k++) {
	    RECORD((*tags)[k]); //print();

	    if(
	       (k == 0 && (tags->length() < 2 || strcmp((*tags)[k+1]->getName(), "XPressRatio"))) ||
	       (k == 1 && ! strcmp((*tags)[k]->getName(), "XPressRatio"))) {
	      Array<Tag*>* asap_tags = getProteinRatioTags(pro_ratio_);
	      if(asap_tags != NULL) {
		ratio_tags_written++;
		const int FEEDBACK = 10;
		if (ratio_tags_written % FEEDBACK == 0) {
		  cout << "..." << ratio_tags_written;
		  fflush(stdout);
		}
		if (ratio_tags_written % (FEEDBACK * 10) == 0)
		  cout << endl;

		RECORD(result_start);
		for(int k = 0; k < asap_tags->length(); k++)
		  if((*asap_tags)[k] != NULL) {
		    RECORD((*asap_tags)[k]);
		    delete (*asap_tags)[k];
		  }
		RECORD(result_stop);
		delete asap_tags;
	      }
	    }
	  }
	  for(k = 0; k < tags->length(); k++)
	    if((*tags)[k] != NULL)
	      delete (*tags)[k];

	  if(tags != NULL) {
	    delete tags;
	    tags = NULL;
	  }

	  if(peps != NULL) {
	    // now done with peps
	    for(int k = 0; k < peps->length(); k++)
	      if((*peps)[k] != NULL)
		delete (*peps)[k];
	    delete peps;
	    peps = NULL;
	  }
	}
      } // if done
      else {
	// add after first tag the additional asap_pro tag...

	// add after header and XPress summary (if exists)
	if(summary != NULL && strcmp(tag->getName(), "XPress_analysis_summary")) {
	  RECORD(summary_start);
	  RECORD(summary); // do this first
	  RECORD(summary_stop);
	  delete summary;
	  summary = NULL;
	}
	if(! summary_filter->filter(tag) && ! pval_summary_filter->filter(tag) && 
	   ! result_filter->filter(tag) && ! pval_result_filter->filter(tag)) {
	  RECORD(tag); //print();
	}
	if(summary != NULL && ! strcmp(tag->getName(), "XPress_analysis_summary")) {
	  RECORD(summary_start);
	  RECORD(summary); // do this first
	  RECORD(summary_stop);

	  delete summary;
	  summary = NULL;
	}

	if(tag->isEnd() && ! strcmp(tag->getName(), "protein_summary_header")) {
	  summary = new Tag("ASAP_prot_analysis_summary", True, True);
	  // get time info

	  char next[20];
	  summary->setAttributeValue("version", szTPPVersionInfo);
	  sprintf(next, "%0.2f", MIN_PEP_PROB);
	  summary->setAttributeValue("min_peptide_probability", next);
	  sprintf(next, "%0.2f", MIN_PEP_WT);
	  summary->setAttributeValue("min_peptide_weight", next);
	  sprintf(next, "%0.2f", MIN_PROB);
	  summary->setAttributeValue("min_protein_probability", next);

	  // add here whether ratio is H/L
	  if(heavy2light_)
	    summary->setAttributeValue("reference_isotope", "light");
	  else
	    summary->setAttributeValue("reference_isotope", "heavy");

	  RECORD(summary_start);
	  RECORD(summary); //print();
	  RECORD(summary_stop);
	  delete summary;
	  summary = NULL;
	} // if end of header
      }

      if(! stored) {
	delete tag;
	tag = NULL;
      }

      data = strstr(data+1, "<");
    }
  }
  fin.close();
  fout.close();

  if (testType!=NO_TEST) {
    //
    // regression test stuff - bpratt Insilicos LLC, Nov 2005
    //
    TagListComparator("ASAPRatioProteinRatioParser",testType,test_tags,testFileName);
    delete[] testFileName;
    for(int k = test_tags.length(); k--;) {
      delete test_tags[k];
    }
  }
  if(! overwrite(protxmlfile, outfile.c_str(), "</protein_summary>")) {
    cout << "error: no asapratio data written to file " << protxmlfile << endl;
  }

  delete summary_filter;
  delete pval_summary_filter;
  delete result_filter;
  delete pval_result_filter;
  delete ratio_filter;
  delete seq_filter;
  delete peak_filter;
  delete dta_filter;
  delete result_start;
  delete result_stop;
  delete summary_start;
  delete summary_stop;
  delete summary;
  delete[] nextline;
}


void ASAPRatioProteinRatioParser::setFilter(Tag* tag) {
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


void ASAPRatioProteinRatioParser::enterUnique(Array<UniquePeptide*>* uniques, const char* next, double wt, double prob) {
  for(int k = 0; k < uniques->length(); k++)
    if(! strcmp((*uniques)[k]->peptide_, next)) {
      if(prob > (*uniques)[k]->probability_)
	(*uniques)[k]->probability_ = prob;
      return;
    }

  // still here
  uniques->insertAtEnd(new UniquePeptide(next, wt, prob));
}


// must substitute # -> ~
char* ASAPRatioProteinRatioParser::getPeptideString(Array<const char*>* peps, const char* link) {
  const int max_pep_length = 500;
  char encoded[max_pep_length];

  char* output = NULL;
  if(peps == NULL || peps->length() == 0)
    return output;

  int k,totlen= 0;
  for(k = 0; k < peps->length(); k++) {
    totlen += (int)strlen((*peps)[k]);
    if(k > 0)
      totlen += (int)strlen(link);
  }
  output = new char[totlen+1];

  if(strstr((*peps)[0], "#") == NULL)
    strcpy(output, (*peps)[0]);
  else {
    int j;
    for(j = 0; ((*peps)[0])[j]; j++) {
      if(((*peps)[0])[j] == '#')
	encoded[j] = '~';
      else 
	encoded[j] = ((*peps)[0])[j];
    }
    encoded[j] = 0;
    strcpy(output, encoded);
  }

  for(k = 1; k < peps->length(); k++) {
    strcat(output, link);
    if(strstr((*peps)[k], "#") == NULL)
      strcat(output, (*peps)[k]);
    else {
      int j;
      for(j = 0; ((*peps)[k])[j]; j++) {
	if(((*peps)[k])[j] == '#')
	  encoded[j] = '~';
	else 
	  encoded[j] = ((*peps)[k])[j];
      }
      encoded[j] = 0;
      strcat(output, encoded);
    }
  }
  output[totlen] = 0;

  return output;
}


void ASAPRatioProteinRatioParser::getRatio(Array<UniquePeptide*>* peptides,
					   double dProbability, double dMinWt, Boolean heavy2light) {

  double mean = 0.0;
  double meansq = 0.0;
  int num = 0;

  // LM: use a single parser instance to store all peptide quant info; boost speed by >100x
  // calculate individual protein ratios by passing peptides list to getProDataStruct
  if(parser_ == NULL) {
    parser_ = new ASAPRatioGroupPeptideParser(&input_pepxmlfiles_, dProbability, dMinWt, heavy2light);
  }
  pro_ratio_ = NULL;

  if(parser_ != NULL) {
    pro_ratio_ = parser_->getProDataStruct(peptides);
    if (verbose) {
      cout << pro_ratio_->ratio[0] << " +- " << pro_ratio_->ratio[1] << " (" << pro_ratio_->dataNum << ")" << endl;
    }
  }
  else {
    cout << "Error: null parser for inputfiles: ";
    int k;
    for(k = 0; k < input_pepxmlfiles_.length(); k++)
      cout << input_pepxmlfiles_[k] << " ";
    cout << " and peptides: ";
    for(k = 0; k < peptides->length(); k++)
      cout << ((*peptides)[k])->peptide_ << " ";
    cout << endl;
    exit(1);
  }
}


Array<Tag*>* ASAPRatioProteinRatioParser::getProteinRatioTags(proDataStrct* pro_ratio) {
  Array<Tag*>* output = new Array<Tag*>;

  Tag* protag = new Tag("ASAPRatio", True, False);
  char next[50];

  double ratiomean = (double)pro_ratio->ratio[0];
  if(ratiomean == -2.0) 
    sprintf(next, "%0.2f", -1.0);
  else if(ratiomean == -1.0) 
    strcpy(next, "999.");
  else 
    sprintf(next, "%0.2f", ratiomean);

  protag->setAttributeValue("ratio_mean", next);
  sprintf(next, "%0.2f", (double)pro_ratio->ratio[1]);
  protag->setAttributeValue("ratio_standard_dev", next);
  sprintf(next, "%d", pro_ratio->dataNum);
  protag->setAttributeValue("ratio_number_peptides", next);

  // now heavy2light
  if(ratiomean == -2.0) {
    protag->setAttributeValue("heavy2light_ratio_mean", "-1.00");
    protag->setAttributeValue("heavy2light_ratio_standard_dev", "0.00");
  }
  else if(ratiomean == -1.0) {
    protag->setAttributeValue("heavy2light_ratio_mean", "0.00");
    protag->setAttributeValue("heavy2light_ratio_standard_dev", "0.00");
  }
  else if(ratiomean == 0.0) {
    protag->setAttributeValue("heavy2light_ratio_mean", "999.");
    protag->setAttributeValue("heavy2light_ratio_standard_dev", "0.00");
  }
  else {
    sprintf(next, "%0.2f", (double)pro_ratio->inv_ratio[0]);
    protag->setAttributeValue("heavy2light_ratio_mean", next);
    if (isnan((double)pro_ratio->inv_ratio[1])) {
      protag->setAttributeValue("heavy2light_ratio_standard_dev", "NaN");
    }
    else {
      sprintf(next, "%0.2f", (double)pro_ratio->inv_ratio[1]);
      protag->setAttributeValue("heavy2light_ratio_standard_dev", next);
    }
  }

  sprintf(next, "%d", pro_ratio->indx);
  protag->setAttributeValue("status", next);

  output->insertAtEnd(protag);
  for(int k = 0; k < pro_ratio->dataNum; k++) {
    Tag* seqtag = new Tag("ASAP_Seq", True, False);
    sprintf(next, "%d", pro_ratio->sequences[k].indx);
    seqtag->setAttributeValue("status", next);
    sprintf(next, "%d", pro_ratio->dataCnts[k]);
    seqtag->setAttributeValue("include", next);
    sprintf(next, "%d", pro_ratio->sequences[k].dataNum);
    seqtag->setAttributeValue("datanum", next);
    sprintf(next, "%0.4f", (double)pro_ratio->sequences[k].ratio[0]);
    seqtag->setAttributeValue("ratio_mean", next);
    sprintf(next, "%0.4f", (double)pro_ratio->sequences[k].ratio[1]);
    seqtag->setAttributeValue("ratio_standard_dev", next);
    sprintf(next, "%0.4f", (double)pro_ratio->sequences[k].inv_ratio[0]);
    seqtag->setAttributeValue("heavy2light_ratio_mean", next);
    if (isnan((double)pro_ratio->sequences[k].inv_ratio[1])) {
      seqtag->setAttributeValue("heavy2light_ratio_standard_dev", "NaN");
    }
    else {
      sprintf(next, "%0.4f", (double)pro_ratio->sequences[k].inv_ratio[1]);
      seqtag->setAttributeValue("heavy2light_ratio_standard_dev", next);
    }
    sprintf(next, "%0.2f", pro_ratio->sequences[k].weight);
    seqtag->setAttributeValue("weight", next);
    output->insertAtEnd(seqtag);
    seqtag->setAttributeValue("light_sequence", pro_ratio->sequences[k].lightSeq);

    pro_ratio->sequences[k].sort_for_output(); // place in a logical output order
    for(int j = 0; j < pro_ratio->sequences[k].dataNum; j++) {
      Tag* peaktag = new Tag("ASAP_Peak", True, False);
      sprintf(next, "%d", pro_ratio->sequences[k].peaks[j].indx);
      peaktag->setAttributeValue("status", next);
      sprintf(next, "%d", pro_ratio->sequences[k].dataCnts[j]);
      peaktag->setAttributeValue("include", next);
      sprintf(next, "%d", pro_ratio->sequences[k].peaks[j].dataNum);
      peaktag->setAttributeValue("datanum", next);
      sprintf(next, "%0.4f", (double)pro_ratio->sequences[k].peaks[j].ratio[0]);
      peaktag->setAttributeValue("ratio_mean", next);
      sprintf(next, "%0.4f", (double)pro_ratio->sequences[k].peaks[j].ratio[1]);
      peaktag->setAttributeValue("ratio_standard_dev", next);
      sprintf(next, "%0.4f", (double)pro_ratio->sequences[k].peaks[j].inv_ratio[0]);
      peaktag->setAttributeValue("heavy2light_ratio_mean", next);
      if (isnan((double)pro_ratio->sequences[k].peaks[j].inv_ratio[1])) {
	peaktag->setAttributeValue("heavy2light_ratio_standard_dev", "NaN");
      }
      else {
	sprintf(next, "%0.4f", (double)pro_ratio->sequences[k].peaks[j].inv_ratio[1]);
	peaktag->setAttributeValue("heavy2light_ratio_standard_dev", next);
      }
      sprintf(next, "%0.1f", (double)pro_ratio->sequences[k].peaks[j].weight);
      peaktag->setAttributeValue("weight", next);

      //DDS:
      sprintf(next, "%d", pro_ratio->sequences[k].peaks[j].bofIndx);
      peaktag->setAttributeValue("peptide_binary_ind", next);     
      //sprintf(next, "%d", pro_ratio->sequences[k].peaks[j].msms_run_idx);
      //peaktag->setAttributeValue("peptide_binary_ind", next);

      output->insertAtEnd(peaktag);
      for(int i = 0; i < pro_ratio->sequences[k].peaks[j].dataNum; i++) {
	Tag* dtatag = new Tag("ASAP_Dta", True, True);
	sprintf(next, "%d", pro_ratio->sequences[k].peaks[j].dataIndx[i]);
	dtatag->setAttributeValue("peptide_index", next);
	sprintf(next, "%d", pro_ratio->sequences[k].peaks[j].dataCnts[i]);
	dtatag->setAttributeValue("include", next);
	output->insertAtEnd(dtatag);
      } // next dta

      // close it
      peaktag = new Tag("ASAP_Peak", False, True);
      output->insertAtEnd(peaktag);
    } // next peak

    // close it
    seqtag = new Tag("ASAP_Seq", False, True);
    output->insertAtEnd(seqtag);

  } // next seq

  protag = new Tag("ASAPRatio", False, True);
  output->insertAtEnd(protag);
  return output;
}
