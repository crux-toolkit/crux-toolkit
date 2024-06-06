/*
Program       : QuanticGroupPeptideParser
Author        : David Shteynberg
Date          : 11.30.20
SVN info      : $Id: QuanticGroupPeptideParser.cpp 8245 2020-08-21 04:55:02Z real_procopio $

Copyright (C) 2020 David Shteynberg

Based on ASAPRatioGroupPeptideParser
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

#include "QuanticGroupPeptideParser.h"

QuanticGroupPeptideParser::QuanticGroupPeptideParser(const char* pepxmlfile, Array<UniqPeptide*>* peptides, double minprob, double minwt) : Parser(NULL) {
  parse_all_ = False;
  //  memset(&data_,0,sizeof(data_));
  peptides_ = peptides;
  min_probability_ = minprob;
 
  data_ = new QuanticPSMData();
 
  pepxmlfiles_ = new Array<const char*>;
  char* next = new char[strlen(pepxmlfile)+1];
  strcpy(next, pepxmlfile);
  pepxmlfiles_->insertAtEnd(next);
  single_input_ = True;
  ratio_ = NULL;

  use_charge_ = true;

  delete[] next;
  init(NULL);
}


QuanticGroupPeptideParser::QuanticGroupPeptideParser(Array<const char*>* pepxmlfiles, Array<UniqPeptide*>* peptides, double minprob, double minwt) : Parser(NULL) {
  parse_all_ = False;
  //  memset(&data_,0,sizeof(data_));
  peptides_ = peptides;
  min_probability_ = minprob;

  min_weight_ = minwt;

  data_ = new QuanticPSMData();
  
  use_charge_ = true;

  pepxmlfiles_ = pepxmlfiles;
  single_input_ = False;

  ratio_ = new QuanticProteinRatio(min_probability_, min_weight_);

  init(NULL);
}


QuanticGroupPeptideParser::QuanticGroupPeptideParser(Array<const char*>* pepxmlfiles, double minprob, double minwt, bool ch) : Parser(NULL) {
  parse_all_ = True;
  //  memset(&data_,0,sizeof(data_));
  peptides_ = NULL;
  min_probability_ = minprob;
  min_weight_ = minwt;

  data_ = NULL;
  
  use_charge_ = ch;

  pepxmlfiles_ = pepxmlfiles;
  single_input_ = False;

  ratio_ = NULL;

  init(NULL);
}


QuanticGroupPeptideParser::~QuanticGroupPeptideParser() {
  for (map <string, vector<psm> >::const_iterator it = all_peptides_map_.begin(); it != all_peptides_map_.end(); it++) {
    for (vector<psm>::const_iterator p = it->second.begin(); p != it->second.end(); p++) {
      delete p->modinfo;
      delete p->data;
    }
  }

  if(ratio_ != NULL)
    delete ratio_;
}

void QuanticGroupPeptideParser::setUseCharge(bool use) {
    use_charge_ = use;
}

void QuanticGroupPeptideParser::parse(const char* xmlfile) {
  Tag* tag = NULL;

  //  int line_width = 10000;
  char *nextline = new char[line_width_];
  char* data = NULL;

  Boolean analyze = False;
  Boolean collection = False;

  double probability = -1.0;
  double weight = -1.0;

  Array<Tag*>* tags = NULL;
  //int elution;
  long scan;
  int precursor_charge;
  char peptide[200];
  int index;
  string scan_name;
  char quant_sum_match[]      = "quantic_summary";
  char msms_sum_match[]      = "msms_run_summary";
  char search_result_match[] = "spectrum_query";

  Boolean setpepstruct = False;

  //for(int k = 0; k < peptides_->length(); k++)
  //  cout << (*peptides_)[k] << " ";
  // cout << endl;

  ModificationInfo* modinfo = NULL;
  char quant_labels[500];
  quant_labels[0] = 0;
  Array<Tag*>* modification_tags = NULL;
  Boolean mod_on = False;
  double lightmass = 0.0;


  int msms_run_idx=-1;

  for(int k = 0; k < pepxmlfiles_->length(); k++) {
    RACI fin((*pepxmlfiles_)[k]); // can read gzipped xml
    if(! fin) {
      cout << "error opening " << (*pepxmlfiles_)[k] << endl;
      exit(1);
    }
 
    while(fin.getline(nextline, line_width_)) {

      if(analyze ||
	 strstr(nextline, quant_sum_match)      != NULL ||
	 strstr(nextline, msms_sum_match)      != NULL ||
	 strstr(nextline, search_result_match) != NULL ||
	 possiblePeptideListMember(nextline)) {

	data = strstr(nextline, "<");
	while(data != NULL) {
	  tag = new Tag(data);
	  if(tag != NULL) {

	    if(tag->isStart() && ! strcmp(tag->getName(), "quantic_summary")) {
	      //	      elution = atoi(tag->getAttributeValue("elution"));
	      strcpy(quant_labels, tag->getAttributeValue("options"));
	      delete tag;
	      tag = NULL;
	    }
	    else if(tag->isStart() && ! strcmp(tag->getName(), "msms_run_summary")) {
	      msms_run_idx++;
	      delete tag;
	      tag = NULL;
	    }
	    else if(tag->isStart() && ! strcmp(tag->getName(), "spectrum_query")) {
	      precursor_charge = atoi(tag->getAttributeValue("assumed_charge"));
	      scan_name = string(tag->getAttributeValue("spectrum"));
	      scan = (long)(atoi(tag->getAttributeValue("start_scan")));
	      index = atoi(tag->getAttributeValue("index"));
	      delete tag;
	      tag = NULL;
	    }

	    else if(tag->isStart() && ! strcmp(tag->getName(), "search_hit") && ! strcmp(tag->getAttributeValue("hit_rank"), "1")
	       && peptideListMember(tag->getAttributeValue("peptide"), &weight)) {
	      analyze = True;
	      tags = new Array<Tag*>;
	      strcpy(peptide, tag->getAttributeValue("peptide"));
	      delete tag;
	      tag = NULL;
	    }

	    else if(tag->isStart() && ! strcmp("modification_info", tag->getName())) {
	      if(modification_tags == NULL)
		modification_tags = new Array<Tag*>;
	      modification_tags->insertAtEnd(tag);
	      mod_on = !tag->isEnd();
	      if (!mod_on) { // tag already closed, process it now
		modinfo = new ModificationInfo(modification_tags);
	      }
	      strcpy(peptide, tag->getAttributeValue("modified_peptide"));
	    }
	    else if(mod_on && tag->isEnd() && ! strcmp("modification_info", tag->getName())) {
	      modification_tags->insertAtEnd(tag);
	      modinfo = new ModificationInfo(modification_tags);
	      mod_on = False;
	    }
	    else if(mod_on) {
	      modification_tags->insertAtEnd(tag);
	    }


	    // must add probability here......
	    else if(tag->isStart() && ! strcmp(tag->getName(), "peptideprophet_result")) {
	      probability = atof(tag->getAttributeValue("probability"));
	      delete tag;
	      tag = NULL;
	    }
	    else if(tag->isStart() && ! strcmp(tag->getName(), "interprophet_result")) {
	      probability = atof(tag->getAttributeValue("probability"));
	      delete tag;
	      tag = NULL;
	    }
	    else if(analyze) {
	      if(tag->isStart() && ! strcmp(tag->getName(), "quantic_result")) {
		if (data_ != NULL) {
		  delete data_;
		  data_ = NULL;
		}
		//data_ = new QuanticPSMData();
		
		tags->insertAtEnd(tag);
		collection = True;
	      }

	      if(tag->isEnd() && ! strcmp(tag->getName(), "quantic_result")) {
		// process
		if(!tag->isStart()) {
		   tags->insertAtEnd(tag);
		}
		
		//
		setPepDataStruct(tags, scan_name, scan, precursor_charge);
		setpepstruct = True;

		// cleanup
		for(int k = 0; k < tags->length(); k++)
		  if((*tags)[k] != NULL)
		    delete (*tags)[k];
		delete tags;
		tags = NULL;

		collection = False;
	      }
	      else if(collection)
		tags->insertAtEnd(tag);


	      else if(tag->isEnd() && ! strcmp(tag->getName(), "search_hit")) {

		if(index >= 0 && weight > -1.0 && strlen(peptide) && setpepstruct) {

		  if (parse_all_) {
		    psm pep;
		    pep.sequence = string(peptide);
		    pep.data = data_;
		    pep.index = index;
		    pep.xml_index = k;
		    pep.wt = weight;
		    pep.prob = probability;
		    pep.msms_run_idx = msms_run_idx;
		    if (use_charge_) 
		      pep.charge = precursor_charge;
		    else 
		      pep.charge = 0;
		    
		    pep.modinfo = modinfo;
		    
		    all_peptides_map_[pep.sequence].push_back(pep);
		    data_ = NULL;
		  }
		  else {
		    ratio_->enter(peptide,
				  modinfo,
				  data_,
				  index,
				  k,
				  weight,
				  probability,
				  msms_run_idx,
				  precursor_charge);
		  }
		}

		// reset
		lightmass = 0.0;
		if(modification_tags != NULL) {
		  for(int k = 0; k  < modification_tags->length(); k++)
		    if((*modification_tags)[k] != NULL)
		      delete (*modification_tags)[k];
		  delete modification_tags;
		  modification_tags = NULL;
		}
		quant_labels[0] = 0;
		//if (modinfo)
		//  delete modinfo;
		modinfo = NULL;
		analyze = False;
		peptide[0] = 0;
		index = -1;
		weight = -1.0;
		probability = -1.0;
		setpepstruct = False;

		delete tag;
		tag = NULL;

		if(tags != NULL) {
		  for(int k = 0; k < tags->length(); k++)
		    if((*tags)[k] != NULL)
		      delete (*tags)[k];
		  delete tags;
		  tags = NULL;
		}

	      }
	      else {
		delete tag; // don't keep
		tag = NULL;


		if(modification_tags != NULL) {
		  for(int k = 0; k < modification_tags->length(); k++)
		    if((*modification_tags)[k] != NULL)
		      delete (*modification_tags)[k];
		  delete modification_tags;
		  modification_tags = NULL;
		}
		//delete modinfo;
		//modinfo = NULL;

	      }

	    } // if analyze
	    else {
	      delete tag;
	      tag = NULL;
	    }

	  } //  if not null
	  data = strstr(data+1, "<");
	} // next tag
      } // if possible peptide present

    } // next line
    fin.close();
  } // next inputfile


  delete modinfo;


  delete[] nextline;
}


void QuanticGroupPeptideParser::setFilter(Tag* tag) {
  if(tag == NULL)
    return;

  if(filter_memory_) {
    filter_memory_ = False;
    filter_ = False;
  }

  if(! strcmp(tag->getName(), "search_hit")) {
    if(tag->isStart() && ! strcmp(tag->getAttributeValue("hit_rank"), "1")) {
      filter_ = True;
    }else{
       if(filter_ && tag->isEnd())
        filter_memory_ = True;
    }
  }

}


void QuanticGroupPeptideParser::clearRatio() {
  if(ratio_ != NULL)
    delete ratio_;
}


proQuantStrct* QuanticGroupPeptideParser::getProQuantStruct(Array<UniqPeptide*>* peptides) {
  Boolean verbose = False;
  peptides_ = NULL;
  peptides_ = peptides;

  clearRatio();
  ratio_ = new QuanticProteinRatio(min_probability_, min_weight_);

  for(int k = 0; k < peptides_->length(); k++) {
    std::string pep_k = (*peptides_)[k]->peptide_;

    if (all_peptides_map_.find(pep_k) == all_peptides_map_.end()) {
      if(verbose)
	cout << endl << "unmatched: " << pep_k << endl;
    }
    else {
      for (int i=0; i < all_peptides_map_.at(pep_k).size(); i++) {
	if(verbose)
	  cout << "found one psm: " << all_peptides_map_.at(pep_k)[i].sequence << endl;

	const char* peptideseq = all_peptides_map_.at(pep_k)[i].sequence.c_str();
	const ModificationInfo* modinfo = all_peptides_map_.at(pep_k)[i].modinfo;
	//const char* quant_labels = all_peptides_map_.at(pep_k)[i].quant_labels.c_str();
	//const pepQuantStrct data = all_peptides_map_.at(pep_k)[i].data;


	if ( ( (*peptides_)[k]->charge_ == 0 || (*peptides_)[k]->charge_ == all_peptides_map_.at(pep_k)[i].charge ) &&
	     (*peptides_)[k]->findModPep(modinfo->getModifiedPeptide())) {

	  //IF equal charges (or charge is 0 and ignored) and mods match
	  
	  QuanticPSMData* data = all_peptides_map_.at(pep_k)[i].data;
	  double wt = all_peptides_map_.at(pep_k)[i].wt;
	  
	  
	  ratio_->enter(modinfo->getModifiedPeptide(), //peptideseq,
			modinfo,
			data,
			all_peptides_map_.at(pep_k)[i].index,
			all_peptides_map_.at(pep_k)[i].xml_index,
			wt,
			all_peptides_map_.at(pep_k)[i].prob, 
			all_peptides_map_.at(pep_k)[i].msms_run_idx,
			all_peptides_map_.at(pep_k)[i].charge);
	}
      }
    }
  }


  if(ratio_ == NULL) {
    cout << "error: null ratio" << endl;
    return NULL;
  }

  return ratio_->getProQuantStrct();
}


proQuantStrct* QuanticGroupPeptideParser::getProQuantStruct() {
  if(ratio_ == NULL) {
    cout << "error: null ratio" << endl;
    return NULL;
  }
  return ratio_->getProQuantStrct();
}



Boolean QuanticGroupPeptideParser::peptideListMember(const char* pep, double* wt) {
  Boolean verbose = False; //strstr(pep, "YEFCTILKK") != NULL;
  if(verbose) {
    cout << "comparing " << pep << " with peptides....";
    for(int k = 0; k < peptides_->length(); k++)
      cout << "=" << (*peptides_)[k] << "=";
    cout << endl;
  }
  if(parse_all_) {
    *wt = 1.1; // placeholder
    return True;
  }
  if(peptides_ == NULL)
    return False;
  for(int k = 0; k < peptides_->length(); k++)
    if(! strcmp((*peptides_)[k]->peptide_.c_str(), pep)) {
      if(verbose)
	cout << "returning true" << endl;
      *wt = (*peptides_)[k]->weight_;

      return True;
    }
  if(verbose) {
    cout << "-" << pep << "-" << (*peptides_)[0]->peptide_ << "-" << endl;
    cout << strlen(pep) << " vs " << (*peptides_)[0]->peptide_.length() << endl;
  }

  return False;
}


Boolean QuanticGroupPeptideParser::possiblePeptideListMember(const char* data) {
  if (parse_all_)
    return True;
  if(peptides_ == NULL)
    return False;
  for(int k = 0; k < peptides_->length(); k++)
    if(strstr(data, (*peptides_)[k]->peptide_.c_str()) != NULL)
      return True;

  return False;
}

void QuanticGroupPeptideParser::parseTuple(vector<double>* result, std::string input) {
 std::istringstream inp(input); 
 if (result) result->clear();
 while (inp) {
    std::string dat;
    if (!std::getline (inp, dat, ',')) break;

    result->push_back(atof(dat.c_str()));
 }

}



void QuanticGroupPeptideParser::setPepDataStruct(Array<Tag*>* tags, string spec_name, long scan, int precursor_charge) {
  Tag* next;
  int charge;

  //DDS: calculate pepArea variable
  bool use_this = false;

  if (data_ == NULL) {
    data_ = new QuanticPSMData();
  }

  data_->scan_ = scan;
  data_->chg_ = precursor_charge;
  data_->specName_ = spec_name;

  for(int k = 0; k < tags->length(); k++) {
    next = (*tags)[k];
    if(next->isStart()) {
      if(! strcmp(next->getName(), "quantic_result")) {
	parseTuple(data_->quantmeans_ , next->getAttributeValue("norm_quant_mean"));
	parseTuple(data_->quantstdvs_ , next->getAttributeValue("norm_quant_stdev"));
	parseTuple(data_->rawquants_ , next->getAttributeValue("quant"));


      }
      
    } // if start
  } // next tag

}
