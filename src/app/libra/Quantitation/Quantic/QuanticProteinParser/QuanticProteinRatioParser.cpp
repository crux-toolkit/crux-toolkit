/*
Program       : QuanticProteinRatioParser
Author        : David Shteynberg 
Date          : 11.30.20
SVN info      : $Id: QuanticProteinRatioParser.cpp 

Computes Quantic based quantitation for proteins, then overwrites
that information into ProteinProphet XML.  Quantic quantitative values
are based on cleavable linker neutral losses quantitation.


Copyright (C) 2020 David Shteynberg

Based on ASAPRatioProteinRatioParser
Copyright (C)  2003 Andrew Keller, Jimmy Eng


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

David Shteynberg
Institute for Systems Biology
401 Terry Avenue North
Seattle, WA  98109  USA
akeller@systemsbiology.org
*/

#include "QuanticProteinRatioParser.h"
#include "Parsers/Parser/TagListComparator.h" // regression test stuff - bpratt Insilicos LLC, Nov 2005
#include "Common/TPPVersion.h" // contains version number, name, revision

static bool verbose = false;

QuanticProteinRatioParser::QuanticProteinRatioParser(const char* protxmlfile, const char *testMode) : Parser("quantic") { 
  iNumRawData_=0;

  parser_ = NULL;
  pro_ratio_ = NULL;
  testMode_ = testMode?strdup(testMode):NULL; // regression test stuff - bpratt Insilicos LLC, Nov 2005

  iprophet_ = false;
  
  init(protxmlfile);
}

QuanticProteinRatioParser::~QuanticProteinRatioParser() { 
  for(int k = 0; k < input_pepxmlfiles_.length(); k++)
    if(input_pepxmlfiles_[k] != NULL)
      delete[] input_pepxmlfiles_[k];

  if(parser_ != NULL)
    delete parser_;
  free(testMode_);
}

void QuanticProteinRatioParser::parse(const char* protxmlfile) {
  //open file and pass along
  char *nextline = new char[line_width_];
  char* data = NULL;
  Tag* tag;

  RACI fin(protxmlfile); // can read gzipped xml
  if(! fin) {
    cout << "QuanticProteinRatioParser: error opening " << protxmlfile << endl;
    exit(1);
  }

  double MIN_PEP_WT = 0.5;
  double MIN_PEP_PROB = 0.5;
  Array<Tag*>* tags = NULL;
  Array<UniqPeptide*>* peps = NULL;
  double next_prot_prob = 2.0;
  double MIN_PROB = 0.2;
  Boolean done = False;
  int ratio_tags_written = 0;

  //  QuanticPvalueParser* pvalue_parser = new QuanticPvalueParser();
  //  const char* pvalue_name = NULL;
  //  if(pvalue_parser != NULL)
  //    pvalue_name = pvalue_parser->getName();
  //  else {
  //    pvalue_name = "quantic_pvalue";
  //  }

  TagFilter* summary_filter = new TagFilter("analysis_summary");
  summary_filter->enterRequiredAttributeVal("analysis", getName());
  //  TagFilter* pval_summary_filter = new TagFilter("analysis_summary");
  //  pval_summary_filter->enterRequiredAttributeVal("analysis", pvalue_name);

  TagFilter* result_filter = new TagFilter("analysis_result");
  result_filter->enterRequiredAttributeVal("analysis", getName());
  //TagFilter* pval_result_filter = new TagFilter("analysis_result");
  //  pval_result_filter->enterRequiredAttributeVal("analysis", pvalue_name);

  //  delete pvalue_parser;

  TagFilter* ratio_filter = new TagFilter("Quantic");
  TagFilter* seq_filter = new TagFilter("Quantic_Seq");

  TagFilter* psm_filter = new TagFilter("Quantic_Psm");

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


  UniqPeptide* upep;

  // construct a tempfile name, possibly in the tmp dir if so configured
  std::string outfile = make_tmpfile_name(protxmlfile);
  std::string mod_pep = "";
  std::string peptide_mod_pep = "";
  int charge = -1;
  int peptide_charge = -1;
  double wt;
  double nspPr;

  bool have_indist = false;


  bool in_peptide = false;
  
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
					    "QuanticProteinRatioParser", // program name
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

      if(tag->isStart() && ! strcmp(tag->getName(), "proteinprophet_details")) {
	const char* next = tag->getAttributeValue("runoptions");

	if (next != NULL)
	  next = strstr(next, "IPROPHET");

	if (next != NULL) {
	  iprophet_ = true;
	}

      }

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

      // filter out entries below min prob, and exclude all previous Quantic calculations
      if(! done && filter_ && next_prot_prob >= MIN_PROB) {
	// new protein
	if(tags == NULL)
	  tags = new Array<Tag*>;

	if(! ratio_filter->filter(tag) && ! seq_filter->filter(tag) && ! psm_filter->filter(tag) &&
	   ! result_filter->filter(tag)) {
	  tags->insertAtEnd(tag);
	  stored = True;
	}

	if(peps == NULL)
	  peps = new Array<UniqPeptide*>;
	if(tag->isStart() && ( ! strcmp(tag->getName(), "peptide") || ! strcmp(tag->getName(), "indistinguishable_peptide") )  ) {
	  upep = NULL;

	  mod_pep = string(tag->getAttributeValue("peptide_sequence"));
	  charge = atoi(tag->getAttributeValue("charge"));
	  
	  
	  if (! strcmp(tag->getName(), "peptide")) {
	    peptide_charge = charge;
	    peptide_mod_pep = mod_pep;
	    have_indist = false;
	    in_peptide = true;
	    wt =  atof(tag->getAttributeValue("weight"));
	    nspPr = atof(tag->getAttributeValue("nsp_adjusted_probability"));
	  }
	  if (! strcmp(tag->getName(), "indistinguishable_peptide")) {
	    in_peptide = false;
	    have_indist = true;
	  }
	}
	else if ( tag->isEnd() && ( ! strcmp(tag->getName(), "peptide") || ! strcmp(tag->getName(), "indistinguishable_peptide") )  ) {
	  if( (! strcmp(tag->getName(), "peptide") && ! have_indist) ||
	      (have_indist && ! strcmp(tag->getName(), "indistinguishable_peptide") )) {

	    upep = enterUnique(peps, mod_pep.c_str(), 
			       wt, nspPr, (double)charge);

	    upep->insertModPep(mod_pep);
	    
	  }

	  if( ! strcmp(tag->getName(), "peptide") && have_indist ) {

	    if (have_indist) {
	      upep = enterUnique(peps, peptide_mod_pep.c_str(), 
				 wt, nspPr, (double)peptide_charge);
	      
	      upep->insertModPep(mod_pep);
	    }
	    
	    have_indist = false;
	    charge = -1;
	    wt = -1;
	    nspPr = -1;
	  }
	}
	else if(tag->isStart() && ! strcmp(tag->getName(), "modification_info")) {
	  // check that weight and prob are above minimum.....	  
	  //
	  mod_pep = string(tag->getAttributeValue("modified_peptide"));
	  if (in_peptide)
	     peptide_mod_pep = mod_pep;

	}

	else if(filter_memory_ && tags != NULL && peps != NULL) {
	  // add the  pro
	  Tag* protag = NULL;
	  if(peps->length() > 0) {
	    getRatio(peps, MIN_PEP_PROB, MIN_PEP_WT);
	  } // only if have enough peps

	  int k;
	  have_indist = false;
	  for(k = 0; k < tags->length(); k++) {
	    if((*tags)[k]->isStart() &&
	       ( ! strcmp((*tags)[k]->getName(), "peptide") || ! strcmp((*tags)[k]->getName(), "indistinguishable_peptide") )) {
	      mod_pep = string((*tags)[k]->getAttributeValue("peptide_sequence"));
	      charge = atoi((*tags)[k]->getAttributeValue("charge"));
	      if (! strcmp((*tags)[k]->getName(), "peptide") ) {
		peptide_charge = charge;
		peptide_mod_pep = mod_pep;
		in_peptide = true;
	      }
	      if (! strcmp((*tags)[k]->getName(), "indistinguishable_peptide") ) {
		have_indist = true;
		in_peptide = false;
	      }
	      RECORD((*tags)[k]); //print();
	    }
	    else if ((*tags)[k]->isStart() && ! strcmp((*tags)[k]->getName(), "modification_info")) {
	      mod_pep = string((*tags)[k]->getAttributeValue("modified_peptide"));
	      if (in_peptide)
		peptide_mod_pep = mod_pep;
	      
	      RECORD((*tags)[k]); //print();
	    }
	    else if((*tags)[k]->isEnd() &&
		    ( ! strcmp((*tags)[k]->getName(), "peptide") || ! strcmp((*tags)[k]->getName(), "indistinguishable_peptide") ) ) {

	      if ( (! strcmp((*tags)[k]->getName(), "peptide") && ! have_indist) ||
		   (have_indist && ! strcmp((*tags)[k]->getName(), "indistinguishable_peptide") )) {
		
		Array<Tag*>* quant_tags = NULL;
		quant_tags = getPeptideRatioTags(pro_ratio_, mod_pep, charge);
		if(quant_tags != NULL && quant_tags->length()) {
		  // ratio_tags_written++;
		  //const int FEEDBACK = 10;
		  //if (ratio_tags_written % FEEDBACK == 0) {
		  //  cout << "..." << ratio_tags_written;
		  //  fflush(stdout);
		  // }
		  //if (ratio_tags_written % (FEEDBACK * 10) == 0)
		  // cout << endl;
		  
		  RECORD(result_start);
		  for(int t = 0; t < quant_tags->length(); t++)
		    if((*quant_tags)[t] != NULL) {
		      RECORD((*quant_tags)[t]);
		      delete (*quant_tags)[t];
		    }
		  RECORD(result_stop);
		  delete quant_tags;
		  
		  
		}
		mod_pep = "";
		charge = -1;
		
		
	      }
	    
              if (! strcmp((*tags)[k]->getName(), "peptide")) {
		if (have_indist) {
		  	
		  Array<Tag*>* quant_tags = NULL;
		  quant_tags = getPeptideRatioTags(pro_ratio_, peptide_mod_pep, peptide_charge);
		  if(quant_tags != NULL && quant_tags->length()) {
		    // ratio_tags_written++;
		    //const int FEEDBACK = 10;
		    //if (ratio_tags_written % FEEDBACK == 0) {
		    //  cout << "..." << ratio_tags_written;
		    //  fflush(stdout);
		    // }
		    //if (ratio_tags_written % (FEEDBACK * 10) == 0)
		    // cout << endl;
		    
		    RECORD(result_start);
		    for(int t = 0; t < quant_tags->length(); t++)
		      if((*quant_tags)[t] != NULL) {
			RECORD((*quant_tags)[t]);
			delete (*quant_tags)[t];
		      }
		    RECORD(result_stop);
		    delete quant_tags;
		    
		    
		  }
		  
		  peptide_mod_pep = "";
		  peptide_charge = -1;
		
		}
		
                have_indist = false;
              }
	      RECORD((*tags)[k]); //print();
	    }
	    else if((k == 0 && (tags->length() < 2 || strcmp((*tags)[k+1]->getName(), "XPressRatio"))) ||
		    (k == 1 && ! strcmp((*tags)[k]->getName(), "XPressRatio"))) {

	      RECORD((*tags)[k]); //print();

	      Array<Tag*>* quant_tags = getProteinRatioTags(pro_ratio_);
	      if(quant_tags != NULL) {
		ratio_tags_written++;
		const int FEEDBACK = 10;
		if (ratio_tags_written % FEEDBACK == 0) {
		  cout << "..." << ratio_tags_written;
		  fflush(stdout);
		}
		if (ratio_tags_written % (FEEDBACK * 10) == 0)
		  cout << endl;

		RECORD(result_start);
		for(int k = 0; k < quant_tags->length(); k++)
		  if((*quant_tags)[k] != NULL) {
		    RECORD((*quant_tags)[k]);
		    delete (*quant_tags)[k];
		  }
		RECORD(result_stop);
		delete quant_tags;
	      }
	    }
	    else {
	      RECORD((*tags)[k]); //print();

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
      }
      else {
	// add after first tag the additional quant_pro tag...

	// add after header and XPress summary (if exists)
	if(summary != NULL && strcmp(tag->getName(), "XPress_analysis_summary")) {
	  RECORD(summary_start);
	  RECORD(summary); // do this first
	  RECORD(summary_stop);
	  delete summary;
	  summary = NULL;
	}
	if(! summary_filter->filter(tag) && // ! pval_summary_filter->filter(tag) && 
	   ! result_filter->filter(tag) ) {    // &&! pval_result_filter->filter(tag)) {
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
	  summary = new Tag("Quantic_prot_analysis_summary", True, True);
	  // get time info

	  char next[20];
	  summary->setAttributeValue("version", szTPPVersionInfo);
	  sprintf(next, "%0.2f", MIN_PEP_PROB);
	  summary->setAttributeValue("min_peptide_probability", next);
	  sprintf(next, "%0.2f", MIN_PEP_WT);
	  summary->setAttributeValue("min_peptide_weight", next);
	  sprintf(next, "%0.2f", MIN_PROB);
	  summary->setAttributeValue("min_protein_probability", next);



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
    TagListComparator("QuanticProteinRatioParser",testType,test_tags,testFileName);
    delete[] testFileName;
    for(int k = test_tags.length(); k--;) {
      delete test_tags[k];
    }
  }
  if(! overwrite(protxmlfile, outfile.c_str(), "</protein_summary>")) {
    cout << "error: no data written to file " << protxmlfile << endl;
  }

  delete summary_filter;
  delete result_filter;
  delete ratio_filter;
  delete seq_filter;
  //delete peak_filter;
  delete psm_filter;
  delete result_start;
  delete result_stop;
  delete summary_start;
  delete summary_stop;
  delete summary;
  delete[] nextline;
}


void QuanticProteinRatioParser::setFilter(Tag* tag) {
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


UniqPeptide* QuanticProteinRatioParser::enterUnique(Array<UniqPeptide*>* uniques, const char* next,
					    double wt, double prob, int charge) {
  for(int k = 0; k < uniques->length(); k++)
    if((*uniques)[k]->charge_ == charge && ! strcmp((*uniques)[k]->peptide_.c_str(), next)) {
      if(prob > (*uniques)[k]->probability_)
	(*uniques)[k]->probability_ = prob;

      return (*uniques)[k];
    }

  // still here
  uniques->insertAtEnd(new UniqPeptide(next, wt, prob, charge));

  return (*uniques)[uniques->length()-1];
	
}


// must substitute # -> ~
char* QuanticProteinRatioParser::getPeptideString(Array<const char*>* peps, const char* link) {
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


void QuanticProteinRatioParser::getRatio(Array<UniqPeptide*>* peptides,
					   double dProbability, double dMinWt) {

  double mean = 0.0;
  double meansq = 0.0;
  int num = 0;

  // LM: use a single parser instance to store all peptide quant info; boost speed by >100x
  // calculate individual protein ratios by passing peptides list to getProDataStruct
  if(parser_ == NULL) {
    bool usech = true;
    if (peptides->size()) {
      usech = ((*peptides)[0]->charge_ != 0);
    }

    parser_ = new QuanticGroupPeptideParser(&input_pepxmlfiles_, dProbability, dMinWt, usech);  

  }

  pro_ratio_ = NULL;

  if(parser_ != NULL) {
    pro_ratio_ = parser_->getProQuantStruct(peptides);
    //    pro_ratio_ = parser_->getProQuantStruct();
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

  //delete parser_;
  //parser_ = NULL;
}


Array<Tag*>* QuanticProteinRatioParser::getPeptideRatioTags(proQuantStrct* pro_ratio, string modpep, int chg) {
  Array<Tag*>* output = new Array<Tag*>;


  string out_str = "";
  ostringstream rep;
  if (!pro_ratio) {
    delete output;
    return NULL;
  }
  
  for(int k = 0; k < pro_ratio->peps_.size(); k++) {
    if ( pro_ratio->peps_[k]->pep_ == modpep && pro_ratio->peps_[k]->chg_ == chg) {
      Tag* seqtag = new Tag("Quantic_Seq", True, False);
      
      rep.str("");  rep.clear();
      rep << std::fixed << std::setprecision(0)
      	  << pro_ratio->peps_[k]->chg_;
      
      seqtag->setAttributeValue("charge",  rep.str().c_str());
      
      //rep.str("");  rep.clear();
      //rep << std::fixed << std::setprecision(0)
      //	  << pro_ratio->peps_[k]->dataInc_ ? 1 : 0;
      
      //seqtag->setAttributeValue("include",  rep.str().c_str());
      
      rep.str("");  rep.clear();
      rep << std::fixed << std::setprecision(0)
	  << pro_ratio->peps_[k]->psms_.size();
      
      
      seqtag->setAttributeValue("datanum", rep.str().c_str());
      
      rep.str("");  rep.clear();
      
      for (int q=0; q<pro_ratio->peps_[k]->quantmeans_.size(); q++) {
	rep << std::fixed << std::setprecision(3)
	    << pro_ratio->peps_[k]->quantmeans_[q];
	if (q<pro_ratio->peps_[k]->quantmeans_.size()-1) {
	  rep << ",";
	}
      }
      
      seqtag->setAttributeValue("quant_means", rep.str().c_str());
      
      
      rep.str("");  rep.clear();
      
      for (int q=0; q<pro_ratio->peps_[k]->quantstdvs_.size(); q++) {
	rep << std::fixed << std::setprecision(3)
	    << pro_ratio->peps_[k]->quantstdvs_[q];
	if (q<pro_ratio->peps_[k]->quantstdvs_.size()-1) {
	  rep << ",";
	}
      }
      
      seqtag->setAttributeValue("quant_stdvs", rep.str().c_str());
      
      rep.str("");  rep.clear();
      rep << std::scientific << std::setprecision(3)
	  << pro_ratio->peps_[k]->weight_;
      
      seqtag->setAttributeValue("weight", rep.str().c_str());
      
      rep.str("");  rep.clear();
      rep << pro_ratio->peps_[k]->pep_;
      
      seqtag->setAttributeValue("pepseq", rep.str().c_str());
      
      output->insertAtEnd(seqtag);
    //pro_ratio->peps_[k].sort_for_output(); // place in a logical output order

      for(int j = 0; j < pro_ratio->peps_[k]->psms_.size(); j++) {
	Tag* psmtag = new Tag("Quantic_Psm", True, True);
	
	//rep.str("");  rep.clear();
	//rep << std::fixed << std::setprecision(0)
	//    << pro_ratio->peps_[k]->psms_[j]->indx_;
	
	//psmtag->setAttributeValue("status",  rep.str().c_str());
	
	//rep.str("");  rep.clear();
	//rep << std::fixed << std::setprecision(0)
	 //   << pro_ratio->peps_[k]->psms_[j]->dataInc_;
	
	//psmtag->setAttributeValue("include",  rep.str().c_str());
	
	rep.str("");  rep.clear();
	rep << std::fixed << std::setprecision(0)
	    << pro_ratio->peps_[k]->psms_[j]->specName_;
	
	psmtag->setAttributeValue("name",  rep.str().c_str());

	rep.str("");  rep.clear();
	
	for (int q=0; q<pro_ratio->peps_[k]->psms_[j]->quantmeans_.size(); q++) {
	  rep << std::fixed << std::setprecision(3)
	      << pro_ratio->peps_[k]->psms_[j]->quantmeans_[q];
	  if (q<pro_ratio->peps_[k]->psms_[j]->quantmeans_.size()-1) {
	    rep << ",";
	  }
	}
	
	psmtag->setAttributeValue("quant_means", rep.str().c_str());
	
	
	rep.str("");  rep.clear();
	
	for (int q=0; q<pro_ratio->peps_[k]->psms_[j]->quantstdvs_.size(); q++) {
	  rep << std::fixed << std::setprecision(3)
	      << pro_ratio->peps_[k]->psms_[j]->quantstdvs_[q];
	  if (q<pro_ratio->peps_[k]->psms_[j]->quantstdvs_.size()-1) {
	    rep << ",";
	  }
	}
	
	psmtag->setAttributeValue("quant_stdvs", rep.str().c_str());
	
	
	
	rep.str("");  rep.clear();
	rep << std::fixed << std::setprecision(0)
	    << pro_ratio->peps_[k]->psms_[j]->weight_;
	
	psmtag->setAttributeValue("weight",rep.str().c_str());
	
	
	// rep.str("");  rep.clear();  
	// rep << std::fixed << std::setprecision(0)
	//     << pro_ratio->peps_[k]->psms_[j]->bofIndx_;
	
	// psmtag->setAttributeValue("peptide_binary_ind", rep.str().c_str());     
	
	output->insertAtEnd(psmtag);
	
	
	// close it
	//psmtag = new Tag("Quantic_Psm", False, True);
	//output->insertAtEnd(psmtag);
      } // next peak
      
      // close it
      seqtag = new Tag("Quantic_Seq", False, True);
      output->insertAtEnd(seqtag);
    }
  }// next seq
 return output;
}

Array<Tag*>* QuanticProteinRatioParser::getProteinRatioTags(proQuantStrct* pro_ratio) {
  Array<Tag*>* output = new Array<Tag*>;

  Tag* protag = new Tag("Quantic", True, True);

  string out_str = "";
  ostringstream rep;


  for (int q=0; q<pro_ratio->quantmeans_.size(); q++) {
    rep << std::fixed << std::setprecision(3)
	<< pro_ratio->quantmeans_[q];
    if (q<pro_ratio->quantmeans_.size()-1) {
      rep << ",";
    }
  }

  protag->setAttributeValue("quant_means", rep.str().c_str());

  rep.str("");
  rep.clear();
  
  for (int q=0; q<pro_ratio->quantstdvs_.size(); q++) {
    rep << std::fixed << std::setprecision(3)
	<< pro_ratio->quantstdvs_[q];
    if (q<pro_ratio->quantstdvs_.size()-1) {
      rep << ",";
    }
  }

  protag->setAttributeValue("quant_stdvs", rep.str().c_str());

  rep.str("");
  rep.clear();

  rep << std::fixed << std::setprecision(0)
      << countProteinRatioN(pro_ratio);

  protag->setAttributeValue("quant_npeps", rep.str().c_str());

  rep.str("");
  rep.clear();

  rep << std::fixed << std::setprecision(0)
      << pro_ratio->indx_;
  
  protag->setAttributeValue("status",  rep.str().c_str());

  output->insertAtEnd(protag);

  return output;
}

int QuanticProteinRatioParser::countProteinRatioN(proQuantStrct* pro_ratio) {
  int rtn = 0;

  for (int i=0; i<pro_ratio->peps_.size(); i++) {
    if (pro_ratio->peps_[i]->weight_ > 0.0000) {
      rtn++;
    }

  }
  return rtn;
}
  
