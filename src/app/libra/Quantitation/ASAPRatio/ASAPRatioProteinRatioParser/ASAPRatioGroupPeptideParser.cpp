/*
Program       : ASAPRatioGroupPeptideParser
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: ASAPRatioGroupPeptideParser.cpp 8245 2020-08-21 04:55:02Z real_procopio $

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

#include "ASAPRatioGroupPeptideParser.h"

ASAPRatioGroupPeptideParser::ASAPRatioGroupPeptideParser(const char* pepxmlfile, Array<UniquePeptide*>* peptides, double minprob, double minwt, Boolean heavy2light) : Parser(NULL) {
  parse_all_ = False;
  memset(&data_,0,sizeof(data_));
  peptides_ = peptides;
  min_probability_ = minprob;
  heavy2light_ = heavy2light;

  ratio_sum_ = 0.0;
  ratio_square_sum_ = 0.0;
  ratio_num_ = 0;

  pepxmlfiles_ = new Array<const char*>;
  char* next = new char[strlen(pepxmlfile)+1];
  strcpy(next, pepxmlfile);
  pepxmlfiles_->insertAtEnd(next);
  single_input_ = True;
  ratio_ = NULL;

  delete[] next;
  init(NULL);
}


ASAPRatioGroupPeptideParser::ASAPRatioGroupPeptideParser(Array<const char*>* pepxmlfiles, Array<UniquePeptide*>* peptides, double minprob, double minwt, Boolean heavy2light) : Parser(NULL) {
  parse_all_ = False;
  memset(&data_,0,sizeof(data_));
  peptides_ = peptides;
  min_probability_ = minprob;
  heavy2light_ = heavy2light;
  min_weight_ = minwt;

  ratio_sum_ = 0.0;
  ratio_square_sum_ = 0.0;
  inv_ratio_sum_ = 0.0;
  inv_ratio_square_sum_ = 0.0;
  ratio_num_ = 0;

  pepxmlfiles_ = pepxmlfiles;
  single_input_ = False;

  ratio_ = new ASAPProteinRatio(min_probability_, min_weight_);

  init(NULL);
}


ASAPRatioGroupPeptideParser::ASAPRatioGroupPeptideParser(Array<const char*>* pepxmlfiles, double minprob, double minwt, Boolean heavy2light) : Parser(NULL) {
  parse_all_ = True;
  memset(&data_,0,sizeof(data_));
  peptides_ = NULL;
  min_probability_ = minprob;
  heavy2light_ = heavy2light;
  min_weight_ = minwt;

  ratio_sum_ = 0.0;
  ratio_square_sum_ = 0.0;
  inv_ratio_sum_ = 0.0;
  inv_ratio_square_sum_ = 0.0;
  ratio_num_ = 0;

  pepxmlfiles_ = pepxmlfiles;
  single_input_ = False;

  ratio_ = NULL;

  init(NULL);
}


ASAPRatioGroupPeptideParser::~ASAPRatioGroupPeptideParser() {
  for (map <string, vector<psm> >::const_iterator it = all_peptides_map_.begin(); it != all_peptides_map_.end(); it++) {
    for (vector<psm>::const_iterator p = it->second.begin(); p != it->second.end(); p++) {
      delete p->modinfo;
    }
  }

  if(ratio_ != NULL)
    delete ratio_;
}


void ASAPRatioGroupPeptideParser::parse(const char* xmlfile) {
  Tag* tag = NULL;

  //  int line_width = 10000;
  char *nextline = new char[line_width_];
  char* data = NULL;

  Boolean analyze = False;
  Boolean collection = False;

  double probability = -1.0;
  double weight = -1.0;

  Array<Tag*>* tags = NULL;
  int elution;
  long scan;
  int precursor_charge;
  char peptide[200];
  int index;

  char asap_sum_match[]      = "asapratio_summary";
  char msms_sum_match[]      = "msms_run_summary";
  char asap_time_match[]     = "asapratio_timestamp";
  char search_result_match[] = "spectrum_query";

  Boolean setpepstruct = False;

  //for(int k = 0; k < peptides_->length(); k++)
  //  cout << (*peptides_)[k] << " ";
  // cout << endl;
#ifdef USE_STD_MODS
  ModificationInfo* modinfo = NULL;
  char quant_labels[500];
  quant_labels[0] = 0;
  Array<Tag*>* modification_tags = NULL;
  Boolean mod_on = False;
  double lightmass = 0.0;
#endif

  int msms_run_idx=-1;

  for(int k = 0; k < pepxmlfiles_->length(); k++) {
    RACI fin((*pepxmlfiles_)[k]); // can read gzipped xml
    if(! fin) {
      cout << "error opening " << (*pepxmlfiles_)[k] << endl;
      exit(1);
    }
 
    while(fin.getline(nextline, line_width_)) {

      if(analyze ||
	 strstr(nextline, asap_sum_match)      != NULL ||
	 strstr(nextline, msms_sum_match)      != NULL ||
	 strstr(nextline, asap_time_match)     != NULL ||
	 strstr(nextline, search_result_match) != NULL ||
	 possiblePeptideListMember(nextline)) {

	data = strstr(nextline, "<");
	while(data != NULL) {
	  tag = new Tag(data);
	  if(tag != NULL) {

	    if(tag->isStart() && ! strcmp(tag->getName(), "asapratio_summary")) {
	      elution = atoi(tag->getAttributeValue("elution"));
	      strcpy(quant_labels, tag->getAttributeValue("labeled_residues"));
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
#ifdef USE_STD_MODS
	    else if(tag->isStart() && ! strcmp("modification_info", tag->getName())) {
	      if(modification_tags == NULL)
		modification_tags = new Array<Tag*>;
	      modification_tags->insertAtEnd(tag);
	      mod_on = !tag->isEnd();
	      if (!mod_on) { // tag already closed, process it now
		modinfo = new ModificationInfo(modification_tags);
	      }
	    }
	    else if(mod_on && tag->isEnd() && ! strcmp("modification_info", tag->getName())) {
	      modification_tags->insertAtEnd(tag);
	      modinfo = new ModificationInfo(modification_tags);
	      mod_on = False;
	    }
	    else if(mod_on) {
	      modification_tags->insertAtEnd(tag);
	    }
#endif

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
	      if(tag->isStart() && ! strcmp(tag->getName(), "asapratio_result")) {
		tags->insertAtEnd(tag);
		collection = True;
	      }
#ifdef USE_STD_MODS
	      else if(tag->isStart() && ! strcmp(tag->getName(), "asapratio_peptide_data")) {
		lightmass = atof(tag->getAttributeValue("light_mass"));
		tags->insertAtEnd(tag);
	      }
#endif
	      else if(tag->isEnd() && ! strcmp(tag->getName(), "asapratio_result")) {
		// process
		tags->insertAtEnd(tag);
		setPepDataStruct(tags, elution, scan, precursor_charge);
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
#ifdef USE_STD_MODS
		    pep.lightmass = lightmass;
		    pep.modinfo = modinfo;
		    pep.quant_labels = string(quant_labels);
#endif

		    all_peptides_map_[pep.sequence].push_back(pep);
		  }
		  else {
		    ratio_->enter(peptide,
#ifdef USE_STD_MODS
				  lightmass,
				  modinfo,
				  quant_labels,
#endif
				  data_,
				  index,
				  k,
				  weight,
				  probability,
				  msms_run_idx);
		  }
		}

		// reset
#ifdef USE_STD_MODS
		lightmass = 0.0;
		if(modification_tags != NULL) {
		  for(int k = 0; k < modification_tags->length(); k++)
		    if((*modification_tags)[k] != NULL)
		      delete (*modification_tags)[k];
		  delete modification_tags;
		  modification_tags = NULL;
		}
		quant_labels[0] = 0;
		modinfo = NULL;
#endif
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

#ifdef USE_STD_MODS
		if(modification_tags != NULL) {
		  for(int k = 0; k < modification_tags->length(); k++)
		    if((*modification_tags)[k] != NULL)
		      delete (*modification_tags)[k];
		  delete modification_tags;
		  modification_tags = NULL;
		}
		delete modinfo;
		modinfo = NULL;
#endif
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

#ifdef USE_STD_MODS
  delete modinfo;
#endif

  delete[] nextline;
}


void ASAPRatioGroupPeptideParser::setFilter(Tag* tag) {
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


void ASAPRatioGroupPeptideParser::clearRatio() {
  if(ratio_ != NULL)
    delete ratio_;
}


proDataStrct* ASAPRatioGroupPeptideParser::getProDataStruct(Array<UniquePeptide*>* peptides) {
  Boolean verbose = False;
  peptides_ = NULL;
  peptides_ = peptides;

  clearRatio();
  ratio_ = new ASAPProteinRatio(min_probability_, min_weight_);

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
	const char* quant_labels = all_peptides_map_.at(pep_k)[i].quant_labels.c_str();
	const pepDataStrct data = all_peptides_map_.at(pep_k)[i].data;
	double wt = (*peptides_)[k]->weight_;


	ratio_->enter(peptideseq,
#ifdef USE_STD_MODS
		      all_peptides_map_.at(pep_k)[i].lightmass,
		      modinfo,
		      quant_labels,
#endif
		      data,
		      all_peptides_map_.at(pep_k)[i].index,
		      all_peptides_map_.at(pep_k)[i].xml_index,
		      wt,
		      all_peptides_map_.at(pep_k)[i].prob, 
		      all_peptides_map_.at(pep_k)[i].msms_run_idx);
      }
    }
  }


  if(ratio_ == NULL) {
    cout << "error: null ratio" << endl;
    return NULL;
  }

  return ratio_->getProDataStruct();
}


proDataStrct* ASAPRatioGroupPeptideParser::getProDataStruct() {
  if(ratio_ == NULL) {
    cout << "error: null ratio" << endl;
    return NULL;
  }
  return ratio_->getProDataStruct();
}

RatioStruct ASAPRatioGroupPeptideParser::getRatio() {
  RatioStruct ratio;

  if(ratio_num_ == 0) {
      ratio.iNumPeptides = 0;
      ratio.dRatio = -9.9;
      ratio.dStdDev = -9.9;
  }
  else {
      ratio.dRatio = ratio_sum_ / ratio_num_;
      ratio.dh2lRatio = inv_ratio_sum_ / ratio_num_;
      if(ratio_num_ == 1) {
	ratio.dStdDev = 0.0;
	ratio.dh2lStdDev = 0.0;
      }
      else {
	ratio.dStdDev = sqrt((ratio_square_sum_ / ratio_num_) - (ratio.dRatio * ratio.dRatio));
	ratio.dh2lStdDev = sqrt((inv_ratio_square_sum_ / ratio_num_) - (ratio.dh2lRatio * ratio.dh2lRatio));
      }
      ratio.iNumPeptides = ratio_num_;
  }
  return ratio;
}

double ASAPRatioGroupPeptideParser::getRatioSum() {
  return ratio_sum_;
}

double ASAPRatioGroupPeptideParser::getRatioSquareSum() {
  return ratio_square_sum_;
}

int ASAPRatioGroupPeptideParser::getRatioNum() {
  return ratio_num_;
}


Boolean ASAPRatioGroupPeptideParser::peptideListMember(const char* pep, double* wt) {
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
    if(! strcmp((*peptides_)[k]->peptide_, pep)) {
      if(verbose)
	cout << "returning true" << endl;
      *wt = (*peptides_)[k]->weight_;

      return True;
    }
  if(verbose) {
    cout << "-" << pep << "-" << (*peptides_)[0]->peptide_ << "-" << endl;
    cout << strlen(pep) << " vs " << strlen((*peptides_)[0]->peptide_) << endl;
  }

  return False;
}


Boolean ASAPRatioGroupPeptideParser::possiblePeptideListMember(const char* data) {
  if (parse_all_)
    return True;
  if(peptides_ == NULL)
    return False;
  for(int k = 0; k < peptides_->length(); k++)
    if(strstr(data, (*peptides_)[k]->peptide_) != NULL)
      return True;

  return False;
}


void ASAPRatioGroupPeptideParser::setPepDataStruct(Array<Tag*>* tags, int elution, long scan, int precursor_charge) {
  Tag* next;
  int charge;

  //DDS: calculate pepArea variable
  bool use_this = false;
  double maxArea = 0;
  double tmpArea = 0;

  double tmp_ltime = 0;
  double tmp_ltime_wd = 0;
  double tmp_htime = 0;
  double tmp_htime_wd = 0;

  double ltime = 0;
  double ltime_wd = 0;
  double htime = 0;
  double htime_wd = 0;

  data_.scan = scan;
  data_.chrg = precursor_charge;
  data_.eltn = elution;

  for(int k = 0; k < tags->length(); k++) {
    next = (*tags)[k];
    if(next->isStart()) {
      if(! strcmp(next->getName(), "asapratio_result")) {
	data_.pepRatio[0] = atof(next->getAttributeValue("mean"));
	data_.pepRatio[1] = atof(next->getAttributeValue("error"));
	data_.pepH2LRatio[0] = atof(next->getAttributeValue("heavy2light_mean"));
	data_.pepH2LRatio[1] = atof(next->getAttributeValue("heavy2light_error"));
	// make change here
	if(data_.pepRatio[0] == -1) {
	  data_.pepRatio[0] = -2;
	  data_.pepH2LRatio[0] = -2;
	}
	else if(data_.pepRatio[0] >= 999.0) {
	  data_.pepRatio[0] = -1;
	  data_.pepH2LRatio[0] = -1;
	}

      }
      else if(! strcmp(next->getName(), "asapratio_peptide_data")) {
	data_.indx = atoi(next->getAttributeValue("status"));
	data_.cidIndx = atoi(next->getAttributeValue("cidIndex"));
	data_.msLight = atof(next->getAttributeValue("light_mass"));
	data_.msHeavy = atof(next->getAttributeValue("heavy_mass"));
	data_.areaFlag = atoi(next->getAttributeValue("area_flag"));
      }
      else if(! strcmp(next->getName(), "asapratio_contribution")) {
	charge = atoi(next->getAttributeValue("charge"));
	data_.pkRatio[charge-1] = atof(next->getAttributeValue("ratio"));
	data_.pkError[charge-1] = atof(next->getAttributeValue("error"));
	data_.pkCount[charge-1] = atoi(next->getAttributeValue("use"));
	
	if (use_this && tmpArea > maxArea) {
	    maxArea = tmpArea;
	    ltime = tmp_ltime;
	    htime = tmp_htime;
	    ltime_wd = tmp_ltime_wd;
	    htime_wd = tmp_htime_wd;
	}

	if (data_.pkCount[charge-1] == 1) {
	  use_this = true;
	  tmpArea = 0;
	  tmp_ltime = 0;
	  tmp_htime = 0;
	  tmp_ltime_wd = 0;
	  tmp_htime_wd = 0;
	}
	else {
	  use_this = false;
	  tmpArea = 0;
	  tmp_ltime = 0;
	  tmp_htime = 0;
	  tmp_ltime_wd = 0;
	  tmp_htime_wd = 0;
	}

      }
      else if(! strcmp(next->getName(), "asapratio_lc_lightpeak")) {
	int label = 0;
	data_.peaks[charge-1][label].indx = atoi(next->getAttributeValue("status"));
	data_.peaks[charge-1][label].valley[0] = atoi(next->getAttributeValue("left_valley"));
	data_.peaks[charge-1][label].valley[1] = atoi(next->getAttributeValue("right_valley"));
	data_.peaks[charge-1][label].bckgrnd = atof(next->getAttributeValue("background"));
	data_.peaks[charge-1][label].area[0] = atof(next->getAttributeValue("area"));
	data_.peaks[charge-1][label].area[1] = atof(next->getAttributeValue("area_error"));
	data_.peaks[charge-1][label].time[0] = atof(next->getAttributeValue("time"));
	data_.peaks[charge-1][label].time[1] = atof(next->getAttributeValue("time_width"));
	data_.peaks[charge-1][label].peak = atoi(next->getAttributeValue("is_heavy"));
	if (use_this) {
	  tmpArea += data_.peaks[charge-1][label].area[0];
	  tmp_ltime = data_.peaks[charge-1][label].time[0];
	  tmp_ltime_wd = data_.peaks[charge-1][label].time[1];
	}
      }
      else if(! strcmp(next->getName(), "asapratio_lc_heavypeak")) {
	int label = 1;
	data_.peaks[charge-1][label].indx = atoi(next->getAttributeValue("status"));
	data_.peaks[charge-1][label].valley[0] = atoi(next->getAttributeValue("left_valley"));
	data_.peaks[charge-1][label].valley[1] = atoi(next->getAttributeValue("right_valley"));
	data_.peaks[charge-1][label].bckgrnd = atof(next->getAttributeValue("background"));
	data_.peaks[charge-1][label].area[0] = atof(next->getAttributeValue("area"));
	data_.peaks[charge-1][label].area[1] = atof(next->getAttributeValue("area_error"));
	data_.peaks[charge-1][label].time[0] = atof(next->getAttributeValue("time"));
	data_.peaks[charge-1][label].time[1] = atof(next->getAttributeValue("time_width"));
	data_.peaks[charge-1][label].peak = atoi(next->getAttributeValue("is_heavy"));
	if (use_this) {
	  tmpArea += data_.peaks[charge-1][label].area[0];
	  tmp_htime = data_.peaks[charge-1][label].time[0];
	  tmp_htime_wd = data_.peaks[charge-1][label].time[1];
	}
      }
    } // if start
  } // next tag

  if (use_this && tmpArea > maxArea) {
    maxArea = tmpArea;
    ltime = tmp_ltime;
    htime = tmp_htime;
    ltime_wd = tmp_ltime_wd;
    htime_wd = tmp_htime_wd;
  }

  data_.pepArea = maxArea;
  data_.pepTime[0][0] = ltime;
  data_.pepTime[0][1] = ltime_wd;
  data_.pepTime[1][0] = htime;
  data_.pepTime[1][1] = htime_wd;
}
