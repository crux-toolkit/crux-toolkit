#include "RTCatalogParser.h"

/*

Program       : RTCatalog                                                       
Author        : David Shteynberg <dshteynb  AT systemsbiology.org>                                                       
Date          : 12.12.07

RTCatalogParser implementation
Copyright (C) 2007, 2017 David Shteynberg

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
#include "Parsers/Parser/TagListComparator.h" // regression test stuff - bpratt Insilicos LLC, Nov 2005
//#include "pwiz/utility/misc/random_access_compressed_ifstream.hpp" // for potentially reading pep.xml.gz
#include "gzstream.h" // for producing .gz files if indicated by output filename
#include "Util/RACI/RACI.h"

#pragma  warning(disable: 4800)

#define PROTONMASS 1.007276466879

RTCatalogParser::RTCatalogParser(double minProb, double minPTMProb, int minRT, int maxRT, 
				 bool gradientCorr, bool tables, bool xics, bool matchedions,
				 string* ignoreFile, string * acnFile, string * pepsbyrunFile,
				 string * gradPepsFile, string * iRTsFile, string* worklistFiles, int minGradPeps) {
  gradCorr_ = gradientCorr;
 
  minProb_ = minProb;
  minPTMProb_ = minPTMProb;

  run_idx_ = -1;

  minGradPeps_ = minGradPeps;

  XICs_ = xics;

  matchedions_ = matchedions;

  byrun_peps_ = NULL;
  byrun_pep_q1s_ = NULL;
  rtmix_chroms_ = NULL;
  track_peps_ = NULL;
  ignore_runs_ = NULL;
  worklist_runs_ = NULL;
  prev_rtmix_byrun_ = NULL;  
  next_rtmix_byrun_ = NULL;
  track_pepq1q3_hash_ = NULL;
  
  rtmix_peps_q1q3_ = NULL;

  rt_cat_ =  new RTCatalog(minProb, minRT, maxRT, acnFile);
  input_files_ = new Array<string*>();
  grad_peps_hash_ = new bool_hash;
  irt_peps_hash_ = new dbl_hash;
  ok_runs_ = new bool_hash;// new Array<string*>();
  //  rt_calcs_= new Array<RTCalculator*>(rt_cat_->num_runs_);
  byrun_rt_calcs_= new rtcalc_hash();


  if (ignoreFile!=NULL) {
    ignore_runs_ = new Array<string*>();

    string run_name;

    ifstream fin(ignoreFile->c_str());

    bool_hash::iterator it;

    while (1) {
      fin >> run_name; 
      if (fin.fail()) break;

      size_t pos1 =  run_name.find_last_of("/\\");
      
      run_name = run_name.substr(pos1+1);
      
      //if (run_name.find(".mzML") == string::npos) {
      //	run_name += ".mzML"; 
      //}

      ignore_runs_->insertAtEnd(new string(run_name));
    }
  }

  if (worklistFiles!=NULL) {
    worklist_runs_ = new Array<Array<string*>*>();
    
    size_t pos1=0;
    size_t pos2=0;
    string listFile;
    string run_name;
    int list_idx = -1;
    
  

    do {
      pos2 = worklistFiles->find(",", pos1);
      listFile = worklistFiles->substr(pos1, pos2-pos1);
      
      list_idx++;
      worklist_runs_->insertAtEnd(new Array<string*>());
      
      
      
      ifstream fin(listFile.c_str());
      
      while (1) {
	fin >> run_name; 
	if (fin.fail()) break;
	(*worklist_runs_)[list_idx]->insertAtEnd(new string(run_name));
      }
      
      pos1 = pos2+1;
    }
    while (pos2 != string::npos);
      
 }
  
  if (gradPepsFile!=NULL) {
    
    string grad_pep;
    
    ifstream fin(gradPepsFile->c_str());
    
    bool_hash::iterator it;
    
    while (1) {
      fin >> grad_pep; 
      if (fin.fail()) break;
      it = grad_peps_hash_->find(grad_pep);
      if (it == grad_peps_hash_->end()) {
	grad_peps_hash_->insert(make_pair(*(new string(grad_pep)), true));
      }
    }
  }
  
  
  
  if (iRTsFile!=NULL) {
    
    string grad_pep;
    double iRT;
    
    ifstream fin(iRTsFile->c_str());
    
    dbl_hash::iterator it;
    bool_hash::iterator it2;
    while (1) {
      fin >> grad_pep;
      if (fin.fail()) break;
      fin >> iRT;
      if (fin.fail()) break;
      
      it = irt_peps_hash_->find(grad_pep);
      if (it == irt_peps_hash_->end()) {
	irt_peps_hash_->insert(make_pair(*(new string(grad_pep)), iRT));
      }
      
      it2 = grad_peps_hash_->find(grad_pep);
      if (it2 == grad_peps_hash_->end()) {
	grad_peps_hash_->insert(make_pair(*(new string(grad_pep)), true));
      }
    }
  }

 if (pepsbyrunFile!=NULL) {
    
    string mod_pep;
    string run_name;
    
    ifstream fin(pepsbyrunFile->c_str());
    byrun_peps_ = new bool_hash_hash();
    bool_hash_hash::iterator it;
    bool_hash::iterator it2;
    while (1) {
      fin >> mod_pep;
      if (fin.fail()) break;
      fin >> run_name;
      if (fin.fail()) break;
      
      it = byrun_peps_->find(run_name);
      if (it == byrun_peps_->end()) {
	byrun_peps_->insert(make_pair(*(new string(run_name)), new bool_hash));
      }
      
      it2 = (*byrun_peps_)[run_name]->find(mod_pep);
      if (it2 == (*byrun_peps_)[run_name]->end()) {
	(*byrun_peps_)[run_name]->insert(make_pair(*(new string(mod_pep)), true));
      }
    }
  }
  
  tablesIn_ = tables;
}

RTCatalogParser::~RTCatalogParser() {
  //TODO implement me
}

void RTCatalogParser::setTolerances(float ppm, float dal) {
  if (ppm < 0 && dal > 0) {
    rt_cat_->setDaltonTolerance();
    XICtol_ = dal;
  }
  else if (ppm < 0) {
    ppm = 5;
    XICtol_ = ppm;
    
  }
  else {
    XICtol_ = ppm;
  }

}

void RTCatalogParser::setRTMixChromatograms(str_hash* rtmix_chroms) {
  rtmix_chroms_ = rtmix_chroms;
  if (!rtmix_peps_q1q3_) 
    rtmix_peps_q1q3_ = new dblarr_hash_hash();
  for (str_hash::iterator itr = rtmix_chroms_->begin(); itr !=rtmix_chroms_->end(); itr++) {
    rtmix_peps_q1q3_->insert(make_pair(itr->first, new dblarr_hash()));
    string* tag = new string(itr->first);
    string* file = new string(itr->second);
    processRTMixChromatograms(tag, file);
    delete tag;
    delete file;
  }
}
  
//Read Q1 Q3 pairs for the rtmixes  
void RTCatalogParser::processRTMixChromatograms(string* tag, string* file) {
  
  ifstream fin(file->c_str());
  string pep;
  double q1, q3;
  dblarr_hash::iterator itr;
  while (1) {

    fin >> pep ;

    if (fin.fail()) break;

    fin >> q1 >> q3;
    
    itr = (*rtmix_peps_q1q3_)[*tag]->find(pep);
    
    if (itr == (*rtmix_peps_q1q3_)[*tag]->end()) {
      (*rtmix_peps_q1q3_)[*tag]->insert(make_pair(*(new string(pep)), new Array<double>()));
    }
    (*(*rtmix_peps_q1q3_)[*tag])[pep]->insertAtEnd(q1);
    (*(*rtmix_peps_q1q3_)[*tag])[pep]->insertAtEnd(q3);
  }
}


void RTCatalogParser::computeRTMixRuns() {
   dblarr_hash_hash::iterator itr1;
   dblarr_hash::iterator itr2;
   size_t pos1, pos2;
   Array<string*>* runs = new Array<string*>();
   for (itr1 = rtmix_peps_q1q3_->begin(); itr1 !=  rtmix_peps_q1q3_->end(); itr1++) {
     runs->clear();
     for (int x=0; x < worklist_runs_->size(); x++) {
       for (int y=0; y < (*worklist_runs_)[x]->size(); y++) {
	 if ( (*(*worklist_runs_)[x])[y]->find(itr1->first.c_str()) != string::npos) {
	   pos1 =  (*(*worklist_runs_)[x])[y]->find_last_of("/\\");
	   //pos2 =  (*(*worklist_runs_)[x])[y]->find_last_of(".");
	   string* rn = new string((*(*worklist_runs_)[x])[y]->substr(pos1+1));
	   runs->insertAtEnd(rn);
	   if (rt_cat_->addRun(*(*(*worklist_runs_)[x])[y])) {
	     if (ok_runs_->find(*rn) == ok_runs_->end()) {
	       string * key = new string(*rn);
	       ok_runs_->insert(make_pair(*key,true));
	     }
	     run_idx_++;
	   }
	   for (itr2 = itr1->second->begin(); itr2 != itr1->second->end();  itr2++) {
	     rt_cat_->insertResult(*(*(*worklist_runs_)[x])[y], itr2->first);
	   }
	 }
       }
     }
  
     rt_cat_->trackPeptidesChromatograms(runs, itr1->second);
   }

}


void RTCatalogParser::addFile(const char* filename) {
  string* name = new string(filename);
  //cerr << "DDS DEBUG: inserting file " << filename << endl;
  input_files_->insertAtEnd(name);
}

void RTCatalogParser::parseTables(const char* c) {

  char *nextline = new char[line_width_];
  string line;
  string mod_pep = "";
  string rt_str = "";
  string spectrum_name = "";
  string exp_lbl = "";
  string charge = "";
  int pos;
  size_t pos1, pos2;
  double rt;
  string run_name;
  for(int k = 0; k < input_files_->length(); k++) {
    RACI fin((*input_files_)[k]->c_str()); // read possibly gzipped files
    if(! fin) {
      cerr << "fin: error opening " << (*input_files_)[k]->c_str() << endl;
      exit(1);
    }
    cout << "Parsing Table File: " << (*input_files_)[k]->c_str() << endl;
    
    rt_cat_->addRun(*((*input_files_)[k])); 

    run_name = *((*input_files_)[k]);
    
    //pos1 =  run_name.find_last_of("/\\");

    //    run_name = run_name.substr(pos1+1);
					   
    int pos;
    
    //TODO: this means the directory names cannot have a '.'
    //pos = run_name.find_first_of('.');

    //if (pos!=string::npos) {
    //  run_name = run_name.substr(0, pos);
    //}
    //ok_runs_->insertAtEnd(new string(run_name));
    
    if (ok_runs_->find(run_name) == ok_runs_->end()) {
      string * key = new string(run_name);
      ok_runs_->insert(make_pair(*key,true));
    }

    while(fin.getline(nextline, line_width_)) {
      line = string(nextline);
      pos = line.find_first_of('\t');
      
      if (pos!=string::npos) {
	mod_pep = line.substr(0, pos);
	line = line.substr(pos+1);
      }
      
      pos = line.find_first_of('\t');

      if (pos==string::npos) continue;

      while (pos!=string::npos) {
	rt_str = line.substr(0, pos);
	rt = atof(rt_str.c_str());
	line = line.substr(pos+1);
	pos = line.find_first_of('\t');
	rt_cat_->insertResult(run_name, spectrum_name, 1., (Array<double>*)NULL, mod_pep, mod_pep, -1., rt, -1, -1, 1, exp_lbl, charge);
      }
      rt_str = line.substr(0);
      rt = atof(rt_str.c_str());
      rt_cat_->insertResult(run_name, spectrum_name, 1., (Array<double>*)NULL, mod_pep, mod_pep, -1., rt, -1, -1, 1,exp_lbl, charge);

    }
  }
}

void RTCatalogParser::useXICs(bool set) {
  XICs_ = set;
}

void RTCatalogParser::parse(const char* c) {
  Array<Tag*>* tags = NULL;
  Tag* tag = NULL;

  char *nextline = new char[line_width_];
  char* data = NULL;

  size_t pos1, pos2;
  string spectrum_name = "";
  string exp_lbl = "";
  string charge = "";
  double prob = 0;
  double iprob = 0;
  Array<double>* allntt_prob=NULL ;
  double calcnmass = -1;
  double rt = -1;
  double ion_injection_time = -1;
  string pep_seq = "";
  string mod_pep = "";
  bool get_pep = false;

  bool skip = false;
  bool top_hit = false;
  bool ignore_run = false;
  int k=0;
  //int run_idx=0;

  int maxntt = 0;
  double prec_intens = -1;
  double collision_eng = -1;
  double matchedions_frac = -1;

  TPP_HASHMAP_T<unsigned int, double>*  modPTMprobs = new TPP_HASHMAP_T<unsigned int, double>();


  if (XICs_)
    parsePeptideQ1s(c);

  // TODO  double allntt_prob[3] = {-100, -100, -100};
  for(k = 0; k < input_files_->length(); k++) {
    RACI fin((*input_files_)[k]->c_str()); // read possibly gzipped files
    if(! fin) {
      cerr << "fin: error opening " << (*input_files_)[k]->c_str() << endl;
      exit(1);
    }
    cout << "Parsing PepXML File: " << (*input_files_)[k]->c_str() << endl;
    string run_name_full;
    string run_name;
    while(fin.getline(nextline, line_width_)) {
      data = strstr(nextline, "<");
      while(data != NULL) {
	tag = new Tag(data);
	
	if (tag != NULL) {

	  if (tag->isStart() && 
	      ! strcmp(tag->getName(), "msms_run_summary")) {
	    //(*input_files_)[k]->assign(tag->getAttributeValue("summary_xml"));
	    //TODO: This has to be more robust
	    run_name = tag->getAttributeValue("base_name");

	    run_name_full = tag->getAttributeValue("raw_data");
	    
	    pos1 =  run_name.find_last_of("/\\");
	    
	    run_name = run_name.substr(pos1+1);
	    
	    run_name += run_name_full;
	    run_name_full = run_name;

	    //if (run_name.find(".mzML") == string::npos) {
	    //  run_name += ".mzML"; 
	    //}
	    //if (run_name_full.find(".mzML") == string::npos) {
	    //  run_name_full += ".mzML"; 
	    //}
	    
	    ignore_run = false;
	    if (ignore_runs_ != NULL || worklist_runs_!= NULL) {
	      if (ignore_runs_ != NULL) {
		for (int t = 0 ; t < ignore_runs_->size(); t++) {
		  if (run_name_full == *(*ignore_runs_)[t] || run_name.find((*ignore_runs_)[t]->c_str()) == 0) {
		    ignore_run = true;
		    break;
		  }
		}
	      }
	      bool present = false;
	      if (worklist_runs_!= NULL) {

		for (int t = 0 ; t < worklist_runs_->size(); t++) {
		  for (int s = 0 ; s < (*worklist_runs_)[t]->size(); s++) {
		    size_t pos1, pos2;
		    string wl_name="";
		    if (run_name.find((*(*worklist_runs_)[t])[s]->c_str()) != string::npos) {
		      present = true;
		      break;
		    }
		    else if ( (pos1 =  (*(*worklist_runs_)[t])[s]->find_last_of("/\\")) != string::npos) {
		      if ( (pos2 =  (*(*worklist_runs_)[t])[s]->find_last_of(".")) != string::npos) {
			wl_name = (*(*worklist_runs_)[t])[s]->substr(pos1+1, pos2-pos1-1);
		      }
		      if (!wl_name.empty()) {
			if (run_name.find(wl_name.c_str()) == 0) {
			  present = true;
			  break;
			}
		      }
		    }
		    
		  }
		}
	      }
	      else {
		present = true;
	      }
	      if (!ignore_run && present) {
		if (ok_runs_->find(run_name) == ok_runs_->end()) {
		  string * key = new string(run_name);
		  ok_runs_->insert(make_pair(*key,true));		
		}
		
		//ok_runs_->insertAtEnd(new string(run_name));
	      }
	      else {
		if (ignore_runs_&&!ignore_run) {
		  ignore_runs_->insertAtEnd(new string(run_name));
		}
		cout << "Ignoring Run: " << run_name << endl;
		ignore_run = true;
	      }

	    }	  
	    else {
	      ignore_run = false;
	      if (ok_runs_->find(run_name) == ok_runs_->end()) {
		string * key = new string(run_name);
		ok_runs_->insert(make_pair(*key,true));
	      }
	      //ok_runs_->insertAtEnd(new string(run_name));
	    }
	    
	    if (!ignore_run) {
	      if (rt_cat_->addRun(run_name_full)) {
		run_idx_++;
	      }	
	    }
	  }
	  else if (!ignore_run && tag->isEnd() && 
	      ! strcmp(tag->getName(), "msms_run_summary")) {
	   
	    ignore_run = false;
	  }
	  else if (tag->isStart() && 
		   ! strcmp(tag->getName(), "spectrum_query")) {
	    top_hit = true;
	    maxntt = 0;
	    get_pep = true;
	    pep_seq = "";
	    mod_pep = "";
	    spectrum_name = "";
	    modPTMprobs->clear();
	    ion_injection_time = -1;
	    if (tag->getAttributeValue("ion_injection_time") != NULL) 
	      collision_eng = atof(tag->getAttributeValue("ion_injection_time"));

	    rt = -1;
	    if (tag->getAttributeValue("retention_time_sec") != NULL) 
	      rt = atof(tag->getAttributeValue("retention_time_sec"))/60;

	    prec_intens = -1;
	    if (tag->getAttributeValue("precursor_intensity") != NULL) 
	      prec_intens = atof(tag->getAttributeValue("precursor_intensity"));

	    collision_eng = -1;
	    if (tag->getAttributeValue("collision_energy") != NULL) 
	      collision_eng = atof(tag->getAttributeValue("collision_energy"));

	    // store the spectrum name for later
	    int len = strlen(tag->getAttributeValue("spectrum"));
	    
	    spectrum_name += tag->getAttributeValue("spectrum");

	    //if (spectrum_name.find("IR-AR1SCXF10-Z100-T1.11135.11135.2",0) != string::npos)
	    //  cerr << "DEBUG: this one !!!" << endl;
	      
	    if (tag->getAttributeValue("experiment_label") != NULL) {
	      len = strlen(tag->getAttributeValue("experiment_label"));
	    }
	    else {
	      len = 0;
	    }
	    exp_lbl = "";
	    if (len > 0)
	      exp_lbl += tag->getAttributeValue("experiment_label");
	    //char* tmp = strrchr(spectrum_name, '.');
	    //tmp = '\0';
	    charge = tag->getAttributeValue("assumed_charge");
	    
	  }  
	  else if (tag->isStart() && 
		   ! strcmp(tag->getName(), "search_hit") && 
		   atoi(tag->getAttributeValue("hit_rank")) > 1) {
	    top_hit = false;
	  }
	  else if (top_hit && tag->isStart() && 
		   ! strcmp(tag->getName(), "alternative_protein")) {
	    if (maxntt < atoi(tag->getAttributeValue("num_tol_term")))
	      maxntt = atoi(tag->getAttributeValue("num_tol_term"));
	      
	  }
	  else if (tag->isStart() && 
		   ! strcmp(tag->getName(), "search_hit") && 
		   atoi(tag->getAttributeValue("hit_rank")) == 1) {
	    top_hit = true;
	    pep_seq = "";
	    prob = 0;
	    iprob = 0;
	    pep_seq += tag->getAttributeValue("peptide");
	    calcnmass = atof(tag->getAttributeValue("calc_neutral_pep_mass"));
	    maxntt = atoi(tag->getAttributeValue("num_tol_term"));
	    //get_pep = false;

	    if (tag->getAttributeValue("num_matched_ions") && tag->getAttributeValue("tot_num_ions")) {
	      matchedions_frac = atof(tag->getAttributeValue("num_matched_ions")) / 
		atof(tag->getAttributeValue("tot_num_ions"));
	      
	      if (ion_injection_time > 0) {
		matchedions_frac /= ion_injection_time;
	      }
	    }

	    if (isnan(matchedions_frac) ) {
	      matchedions_frac  = -1;
	    }
	  }

	  else if (tag->isStart() && get_pep &&
		   ! strcmp(tag->getName(), "modification_info")) {      
	    get_pep = false;
	    mod_pep = "";
	    mod_pep += tag->getAttributeValue("modified_peptide");
	  }

	  else if (tag->isStart() && 
		   ! strcmp(tag->getName(), "mod_aminoacid_mass")) { 
	    modPTMprobs->insert(make_pair(atoi(tag->getAttributeValue("position")), -1));

	  }
	  else if (tag->isStart() && 
		   ! strcmp(tag->getName(), "mod_nterm_mass")) { 
	    modPTMprobs->insert(make_pair(0, -1));

	  }
	  else if (tag->isStart() &&
		   ! strcmp(tag->getName(), "mod_cterm_mass")) { 
	    modPTMprobs->insert(make_pair(pep_seq.length()+1, -1));

	  }

	  else if (tag->isStart() && 
		   ! strcmp(tag->getName(), "mod_aminoacid_probability")) { 
	    unsigned int pos = atoi(tag->getAttributeValue("position"));
	    double ptmprob = atof(tag->getAttributeValue("probability"));
	    if (modPTMprobs->find(pos) != modPTMprobs->end()) {
	      if (ptmprob > (*modPTMprobs)[pos]) {
		(*modPTMprobs)[pos] = ptmprob;
	      }
	      
	    }


	  }

	  else if (tag->isStart() && 
		   ! strcmp(tag->getName(), "peptideprophet_result")) {

	    // got the spectrum name and probability
	    prob = atof( tag->getAttributeValue("probability") );
	    if (mod_pep == "") {
	      mod_pep = pep_seq;
	    }
	    allntt_prob = new Array<double>(3);
	    //TODO: This parsing is fragile, perhaps move to boost regex
 	    const char* nttprob =  tag->getAttributeValue("all_ntt_prob");
	    char* buf = new char[strlen(nttprob)];
	    strcpy(buf, nttprob+1);
	    char* c = strchr(buf, ',');
	    *c = '\0';
	    (*allntt_prob)[0] = atof(buf);
	    int len = strlen(buf)+2;
	    strcpy(buf, nttprob+len);
	    c = strchr(buf, ',');
	    *c = '\0';
	    (*allntt_prob)[1] = atof(buf);
	    strcpy(buf, strrchr(nttprob,',')+1);
	    c = strchr(buf, ')');
	    *c = '\0';
	    (*allntt_prob)[2] = atof(buf);
	    delete [] buf;

	    if ( tag->getAttributeValue("analysis") == NULL || ( strcmp(tag->getAttributeValue("analysis"), "none") && strcmp(tag->getAttributeValue("analysis"), "incomplete")) ) {
	      skip = false;//  rt_cat_->insertResult(k, spectrum_name, prob, allntt_prob, pep_seq, mod_pep, calcnmass, exp_lbl, charge);
	    }
	    else {
	      skip = true;
	    }


	  }
	  else if (tag->isStart() && 
		   ! strcmp(tag->getName(), "interprophet_result")) {

	    // got the spectrum name and probability
	    iprob = atof( tag->getAttributeValue("probability") );
	    if (mod_pep == "") {
	      mod_pep = pep_seq;
	    }
	    allntt_prob = new Array<double>(3);
	    //TODO: This parsing is fragile, perhaps move to boost regex
 	    const char* nttprob =  tag->getAttributeValue("all_ntt_prob");
	    char* buf = new char[strlen(nttprob)];
	    strcpy(buf, nttprob+1);
	    char* c = strchr(buf, ',');
	    *c = '\0';
	    (*allntt_prob)[0] = atof(buf);
	    int len = strlen(buf)+2;
	    strcpy(buf, nttprob+len);
	    c = strchr(buf, ',');
	    *c = '\0';
	    (*allntt_prob)[1] = atof(buf);
	    strcpy(buf, strrchr(nttprob,',')+1);
	    c = strchr(buf, ')');
	    *c = '\0';
	    (*allntt_prob)[2] = atof(buf);
	    delete [] buf;

	    if ( tag->getAttributeValue("analysis") == NULL ||
		 ( strcmp(tag->getAttributeValue("analysis"), "none") && 
		   strcmp(tag->getAttributeValue("analysis"), "incomplete")) ) {
	      skip = false;//  rt_cat_->insertResult(k, spectrum_name, prob, allntt_prob, pep_seq, mod_pep, calcnmass, exp_lbl, charge);
	    }
	    else {
	      skip = true;
	    }

	  }
	  else if (tag->isEnd() && 
		   //		   ! strcmp(tag->getName(), "spectrum_query")) {
		   		   ! strcmp(tag->getName(), "search_hit")) {

	    //if (pep_seq.find("DLFSVLK",0) != string::npos)
	    //  cerr << "DEBUG: this one !!! " << endl;

	    if (byrun_peps_) {
	      skip = true;
	      
	      if (byrun_peps_->find(run_name) != byrun_peps_->end() &&
		  (*byrun_peps_)[run_name]->find(mod_pep) != (*byrun_peps_)[run_name]->end())
		skip = false;

	      if (byrun_peps_->find(run_name.substr(0,run_name.find_last_of("."))) != byrun_peps_->end() &&
		  (*byrun_peps_)[run_name.substr(0,run_name.find_last_of("."))]->find(mod_pep) != (*byrun_peps_)[run_name.substr(0,run_name.find_last_of("."))]->end())
		skip = false;
	    }

	    for (TPP_HASHMAP_T<unsigned int, double>::iterator itr = modPTMprobs->begin();
		 itr !=  modPTMprobs->end(); itr++) {
	      if (itr->second > -1 && itr->second < minPTMProb_) {
		skip = true;
		break;
	      }

	    }
	  
	    if (iprob < minProb_ || prob < minProb_) {
	      skip = true;
	    }

	    if ( !ignore_run && !skip ) {
	      prec_intens = prec_intens+1;
	      
	      if (!matchedions_) {
		matchedions_frac = -1;
	      }
	      prob = (*allntt_prob)[maxntt];
	      rt_cat_->insertResult(run_name_full, spectrum_name, prob,
				    allntt_prob, pep_seq, mod_pep, calcnmass,
				    rt, prec_intens, collision_eng,  matchedions_frac,
				    exp_lbl, charge);
	    }
	    rt = -1;
	    skip = false;
	    top_hit = false;
	    maxntt = 0;
	    exp_lbl = "";
	    spectrum_name = "";
	    charge = "";
	    get_pep = false;
	    pep_seq = "";
	    mod_pep = "";
	    prob = 0;
	     // TODO  allntt_prob[0] = -100; allntt_prob[1] = -100; allntt_prob[2] = -100;
	  }
		 	
	}
	delete tag;
	data = strstr(data+1, "<");      
      }
    }
    fin.close();

  }
  delete [] nextline;

}

void RTCatalogParser::parsePeptideQ1s(const char* c) {
  Array<Tag*>* tags = NULL;
  Tag* tag = NULL;

  char *nextline = new char[line_width_];
  char* data = NULL;

  size_t pos1, pos2;
  string spectrum_name = "";
  string exp_lbl = "";
  string charge = "";
  double prob = 0;
  Array<double>* allntt_prob=NULL ;
  double calcnmass = -1;
  double rt = -1;
  string pep_seq = "";
  string mod_pep = "";
  bool get_pep = false;

  bool skip = false;
  bool ignore_run = false;
  int k=0;
  int run_idx=0;

  double prec_intens = -1;
  double prec_mass = -1;
  double collision_eng = -1;

  if (!byrun_pep_q1s_) {
    byrun_pep_q1s_ = new dblarr_hash_hash();
  }

  // TODO  double allntt_prob[3] = {-100, -100, -100};
  for(k = 0; k < input_files_->length(); k++) {
    RACI fin((*input_files_)[k]->c_str()); // read possibly gzipped files
    if(! fin) {
      cerr << "fin: error opening " << (*input_files_)[k]->c_str() << endl;
      exit(1);
    }
    cout << "Extracting Q1 values for identified peptides in PepXML File: " << (*input_files_)[k]->c_str() << endl;
    string run_name_full;
    while(fin.getline(nextline, line_width_)) {
      data = strstr(nextline, "<");
      while(data != NULL) {
	tag = new Tag(data);
	
	if (tag != NULL) {

	  if (tag->isStart() && 
	      ! strcmp(tag->getName(), "msms_run_summary")) {
	    //(*input_files_)[k]->assign(tag->getAttributeValue("summary_xml"));
	    //TODO: This has to be more robust
	    string run_name = tag->getAttributeValue("base_name");

	    run_name_full = tag->getAttributeValue("raw_data");
	    
	    
	    pos1 =  run_name.find_last_of("/\\");
	    
	    run_name = run_name.substr(pos1+1);

	    run_name += run_name_full;
	    run_name_full = run_name;
	    
	    //if (run_name.find(".mzML") == string::npos) {
	    //  run_name += ".mzML"; 
	    //}
	    //if (run_name_full.find(".mzML") == string::npos) {
	    //  run_name_full += ".mzML"; 
	    //}

	    ignore_run = false;
	    if (ignore_runs_ != NULL || worklist_runs_!= NULL) {
	      if (ignore_runs_ != NULL) {
		for (int t = 0 ; t < ignore_runs_->size(); t++) {
		  if (run_name == *(*ignore_runs_)[t] || run_name.find((*ignore_runs_)[t]->c_str()) == 0) {
		    ignore_run = true;
		    break;
		  }
		}
	      }
	      bool present = false;
	      if (worklist_runs_!= NULL) {

		for (int t = 0 ; t < worklist_runs_->size(); t++) {
		  for (int s = 0 ; s < (*worklist_runs_)[t]->size(); s++) {
		    size_t pos1, pos2;
		    string wl_name="";
		    if (run_name.find((*(*worklist_runs_)[t])[s]->c_str()) != string::npos) {
		      present = true;
		      break;
		    }
		    else if ( (pos1 =  (*(*worklist_runs_)[t])[s]->find_last_of("/\\")) != string::npos) {
		      if ( (pos2 =  (*(*worklist_runs_)[t])[s]->find_last_of(".")) != string::npos) {
			wl_name = (*(*worklist_runs_)[t])[s]->substr(pos1+1, pos2-pos1-1);
		      }
		      if (!wl_name.empty()) {
			if (run_name.find(wl_name.c_str()) != string::npos) {
			  present = true;
			  break;
			}
		      }
		    }
		    
		  }
		}
	      }
	      else {
		present = true;
	      }
	      if (!ignore_run && present) {
		if (ok_runs_->find(run_name) == ok_runs_->end()) {
		  string * key = new string(run_name);
		  ok_runs_->insert(make_pair(*key,true));		
		}
		
		//ok_runs_->insertAtEnd(new string(run_name));
	      }
	      else {
		if (ignore_runs_&&!ignore_run) {
		  ignore_runs_->insertAtEnd(new string(run_name));
		}
		cout << "Ignoring Run: " << run_name << endl;
		ignore_run = true;
	      }

	    }	  
	    else {
	      ignore_run = false;
	      if (ok_runs_->find(run_name) == ok_runs_->end()) {
		string * key = new string(run_name);
		ok_runs_->insert(make_pair(*key,true));
	      }
	      //ok_runs_->insertAtEnd(new string(run_name));
	    }
	    
	    if (!ignore_run) {
	      //if (rt_cat_->addRun(run_name_full)) {
		run_idx++;
	      //}	
	    }
	  }
	  else if (!ignore_run && tag->isEnd() && 
	      ! strcmp(tag->getName(), "msms_run_summary")) {
	   
	    ignore_run = false;
	  }
	  else if (tag->isStart() && 
		   ! strcmp(tag->getName(), "spectrum_query")) {
	    get_pep = true;
	    pep_seq = "";
	    mod_pep = "";
	    spectrum_name = "";
	    rt = -1;
	    if (tag->getAttributeValue("retention_time_sec") != NULL) 
	      rt = atof(tag->getAttributeValue("retention_time_sec"))/60;

	    prec_intens = -1;
	    if (tag->getAttributeValue("precursor_intensity") != NULL) 
	      prec_intens = atof(tag->getAttributeValue("precursor_intensity"));
	    
	    prec_mass = -1;
	    if (tag->getAttributeValue("precursor_neutral_mass") != NULL) 
	      prec_mass = atof(tag->getAttributeValue("precursor_neutral_mass"));


	    collision_eng = -1;
	    if (tag->getAttributeValue("collision_energy") != NULL) 
	      collision_eng = atof(tag->getAttributeValue("collision_energy"));

	    // store the spectrum name for later
	    int len = strlen(tag->getAttributeValue("spectrum"));
	    
	    spectrum_name += tag->getAttributeValue("spectrum");
	    
	    if (tag->getAttributeValue("experiment_label") != NULL) {
	      len = strlen(tag->getAttributeValue("experiment_label"));
	    }
	    else {
	      len = 0;
	    }
	    exp_lbl = "";
	    if (len > 0)
	      exp_lbl += tag->getAttributeValue("experiment_label");
	    //char* tmp = strrchr(spectrum_name, '.');
	    //tmp = '\0';
	    charge = tag->getAttributeValue("assumed_charge");
	    
	  }
	  else if (tag->isStart() && 
		   ! strcmp(tag->getName(), "search_hit") && 
		   atoi(tag->getAttributeValue("hit_rank")) == 1) {
	    pep_seq = "";
	    pep_seq += tag->getAttributeValue("peptide");
	    calcnmass = atof(tag->getAttributeValue("calc_neutral_pep_mass"));
	    //get_pep = false;
	    
	  }

	  else if (tag->isStart() && get_pep &&
		   ! strcmp(tag->getName(), "modification_info")) {      
	    get_pep = false;
	    mod_pep = "";
	    mod_pep += tag->getAttributeValue("modified_peptide");
	  }

	  else if (tag->isStart() && 
		   ! strcmp(tag->getName(), "peptideprophet_result")) {

	    // got the spectrum name and probability
	    prob = atof( tag->getAttributeValue("probability") );
	    if (mod_pep == "") {
	      mod_pep = pep_seq;
	    }
	    allntt_prob = new Array<double>(3);
	    //TODO: This parsing is fragile, perhaps move to boost regex
 	    const char* nttprob =  tag->getAttributeValue("all_ntt_prob");
	    char* buf = new char[strlen(nttprob)];
	    strcpy(buf, nttprob+1);
	    char* c = strchr(buf, ',');
	    *c = '\0';
	    (*allntt_prob)[0] = atof(buf);
	    int len = strlen(buf)+2;
	    strcpy(buf, nttprob+len);
	    c = strchr(buf, ',');
	    *c = '\0';
	    (*allntt_prob)[1] = atof(buf);
	    strcpy(buf, strrchr(nttprob,',')+1);
	    c = strchr(buf, ')');
	    *c = '\0';
	    (*allntt_prob)[2] = atof(buf);
	    delete [] buf;

	    if ( tag->getAttributeValue("analysis") == NULL ||
		 ( strcmp(tag->getAttributeValue("analysis"), "none") 
		   && strcmp(tag->getAttributeValue("analysis"), "incomplete")) ) {
	      skip = false;//  rt_cat_->insertResult(k, spectrum_name, prob, allntt_prob, pep_seq, mod_pep, calcnmass, exp_lbl, charge);
	    }
	    else {
	      skip = true;
	    }


	  }
	  else if (tag->isStart() && 
		   ! strcmp(tag->getName(), "interprophet_result")) {

	    // got the spectrum name and probability
	    prob = atof( tag->getAttributeValue("probability") );
	    if (mod_pep == "") {
	      mod_pep = pep_seq;
	    }
	    allntt_prob = new Array<double>(3);
	    //TODO: This parsing is fragile, perhaps move to boost regex
 	    const char* nttprob =  tag->getAttributeValue("all_ntt_prob");
	    char* buf = new char[strlen(nttprob)];
	    strcpy(buf, nttprob+1);
	    char* c = strchr(buf, ',');
	    *c = '\0';
	    (*allntt_prob)[0] = atof(buf);
	    int len = strlen(buf)+2;
	    strcpy(buf, nttprob+len);
	    c = strchr(buf, ',');
	    *c = '\0';
	    (*allntt_prob)[1] = atof(buf);
	    strcpy(buf, strrchr(nttprob,',')+1);
	    c = strchr(buf, ')');
	    *c = '\0';
	    (*allntt_prob)[2] = atof(buf);
	    delete [] buf;

	    if ( tag->getAttributeValue("analysis") == NULL || 
		 ( strcmp(tag->getAttributeValue("analysis"), "none") && 
		   strcmp(tag->getAttributeValue("analysis"), "incomplete")) ) {
	      skip = false;//  rt_cat_->insertResult(k, spectrum_name, prob, allntt_prob, pep_seq, mod_pep, calcnmass, exp_lbl, charge);
	    }
	    else {
	      skip = true;
	    }

	  }
	  else if (tag->isEnd() && 
		   ! strcmp(tag->getName(), "spectrum_query")) {
	    
	    if (!rt_cat_->rejectResult(prob,rt) &&  !ignore_run && !skip ) {
	      string* rn = new string(run_name_full);  
	      
	      size_t pos1, pos2;
	      //TODO: set minimum probability option
	      if (byrun_pep_q1s_->find(*rn) == byrun_pep_q1s_->end()) {
	
		pos1 =  rn->find_last_of("/\\");
		//pos2 =  rn->find_last_of(".");
		*rn =  rn->substr(pos1+1);//, pos2-pos1-1);
		
		
		
		//if (run_name_full.find(".mzML") == string::npos) {
		//  *rn += ".mzML"; 		  
		//}

		byrun_pep_q1s_->insert(make_pair(*rn, new dblarr_hash)); 
	      }

	      if ((*byrun_pep_q1s_)[*rn]->find(mod_pep) == (*byrun_pep_q1s_)[*rn]->end()) { 
		     (*byrun_pep_q1s_)[*rn]->insert(make_pair(mod_pep, new Array<double>));  
	      }	 
	      (*(*byrun_pep_q1s_)[*rn])[mod_pep]->insertAtEnd((prec_mass + atof(charge.c_str())*PROTONMASS) / atof(charge.c_str()));
	    }
	    rt = -1;
	    skip = false;
	    exp_lbl = "";
	    spectrum_name = "";
	    charge = "";
	    get_pep = false;
	    pep_seq = "";
	    mod_pep = "";
	    prob = 0;
	     // TODO  allntt_prob[0] = -100; allntt_prob[1] = -100; allntt_prob[2] = -100;
	  }
		 	
	}
	delete tag;
	data = strstr(data+1, "<");      
      }
    }
    fin.close();

  }
  delete [] nextline;

}

void  RTCatalogParser::setMods(bool pv, bool cys) {
  rt_cat_->setMods(pv,cys);
}

void RTCatalogParser::writeRTCatalog(const char* file) {
 ifstream fin(file);
  if (fin.good()) {
    cerr << "ERROR: cannot write output to existing file: " << file << endl;
    exit(1);
  }
  fin.close();
  
  ofstream fout(file);
  if (! fout) {
    cerr << "ERROR: cannot write output to file: " << file << endl;
    exit(1);
  }
  cout << "Writing Catalog to file: " << file << endl;

  if (gradCorr_) { 
    //    correctRTsByRun();
    string res = string(file) + string(".gradCorr.html");
    ofstream gout(res.c_str());
    gradientCorrectReport(gout);
    gout.close();
  }
  rt_cat_->calcRTStatsCombined(fout);
  if (!gradCorr_) {
    rt_cat_->calcRTStatsByRun();
  }
  rt_cat_->writeRunRTStats(file);
}


void RTCatalogParser::gradientCorrectReport(ostream& out) {
 out << "<html>" << endl;
  out << "<head>" << endl;
  out << "<title>  </title>" << endl;
  out << "<link type='text/css' rel='stylesheet' href='/tpp/html/css/ex.css?3.2'/>" << endl;
  out << "<script type='text/javascript' src='/tpp/html/js/protovis.js '></script>" << endl;
  out << "</head>" << endl;
  out << "<body>" << endl;
  out << "<center><div id='id'>" << endl;
  out << "<script type='text/javascript+protovis'>" << endl;
  
  // Sizing and scales.

  out << "var w = 800," << endl;
  out << " h = 600;" << endl;
  //  out << " x = pv.Scale.linear(0, 400).range(0, w)," << endl;
  //  out << " y = pv.Scale.linear(0, 4000).range(0, h);" << endl;

  out << "var run_grads = [" << endl;
  for (int k = 0; k < rt_cat_->num_runs_ ; k++) {
    out << "[{run_index: " << k ;
    out <<", run: '" << (*rt_cat_->run_names_)[k]->substr((*rt_cat_->run_names_)[k]->find_last_of("/\\")+1) << "', ";
    out << "slope: " << (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[k]]->getRunSlope() << ", ";
    out << "intercept: " << (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[k]]->getRunInt() << ", ";
    out << "rsq: " << (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[k]]->getRunRSQ() << ", ";
    out << "npeps: " << (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[k]]->getUsedRTCount() << "}]";
    
    if (k < rt_cat_->num_runs_-1) {
      out << ", ";
    }
    out << endl;
  }
  out << "];" << endl;

  out << "var runs = [ " << endl;
  for (int k = 0; k < rt_cat_->num_runs_ ; k++) {
    out << "[" << endl;
    for (int i = 0; i < (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[k]]->getAllRTCount(); i++) {
      out << "{pep: '" << (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[k]]->getModPeptide(i)->c_str() 
	  << "', obs_rt: " << (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[k]]->getObsRT(i) << ", "
	  << "calc_rt: " << (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[k]]->getPredRT(i) << ", "
	  << "lreg_rt: " << ((*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[k]]->getObsRT(i)*(*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[k]]->getRunSlope()+(*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[k]]->getRunInt()) << ", "
	  << "outlier: " << (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[k]]->isOutlier(i) << "}";
      
      if (i < (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[k]]->getAllRTCount()-1) {
	out << ", ";
      }
      out << endl;
    }
    out << "]";
    if (k < rt_cat_->num_runs_-1) {
	out << ", ";
    }
    out << endl;
  }
  out << "];" << endl; 


  out << "var size = 180;" << endl;
  out << "var pad = 40;" << endl;
  out << "var w = 6;" << endl;
  out << "var maxRT = 120" << endl;//<< rt_cat_->acn_gradient_->getAcn(5000) << ";" << endl;
  out << "var maxY = 120;" << endl;


  out << "var run_idx_=-1;" << endl;

  out << "var vis = new pv.Panel()" << endl;
  out << ".width((size+pad)*w)" << endl;
  out << ".height((size+pad)*(1+Math.floor(runs.length/w)))" << endl;
  out << ".left(10)" << endl;
  out << ".top(5)" << endl;
  out << ".events('all');" << endl;

  out << "var plot = vis.add(pv.Panel)" << endl;
  out << "  .data(runs)" << endl;
  out << "  .top(function() { return Math.floor(this.index / w) * (size + pad) + pad / 2; } )" << endl;
  out << "  .left(function() { return this.index % w * (size + pad) + pad / 2; } )" << endl;
  out << "  .height(size)" << endl;
  out << "  .width(size);" << endl;

  out << "var plot_area = vis.add(pv.Panel)" << endl;
  out << "  .data(run_grads)" << endl;
  out << "  .top(function() { return Math.floor(this.index / w) * (size + pad) + pad / 2; } )" << endl;
  out << "  .left(function() { return this.index % w * (size + pad) + pad / 2; } )" << endl;
  out << "  .height(size)" << endl;
  out << "  .width(size);" << endl;

  out << "//  .event('mouseover',  function() { run_idx_ = this.index; vis.render(); console.log(run_idx_);})" << endl;
  out << "//  .event('mouseout',  function() { run_idx_ = -1; vis.render(); console.log(run_idx_);});" << endl;

  out << "var cell = plot.add(pv.Panel)" << endl;
  out << "  .data(function(c) { return c; } );" << endl;

  out << "var area = plot_area.add(pv.Panel)" << endl;
  out << " .data(function(c) { return c; } );" << endl;

  out << "var scale = pv.range(0, 120, 10);" << endl; // rt_cat_->acn_gradient_->getAcn(4850) << "," << rt_cat_->acn_gradient_->getAcn(4850)/10 << ");" << endl;
  
  out << "var yscale = pv.range(0, 120, 10);" << endl;

  out << "//console.log( run_grads);" << endl;

  out << "var scatr = cell.add(pv.Dot)" << endl;
  out << "  .left(function(d) { return d.obs_rt*size/maxRT; } )" << endl;
  if (irt_peps_hash_->size()==0) {
    out << "  .bottom(function(d) { return d.calc_rt*size/maxRT ; } )" << endl;
  }
  else {
    out << "  .bottom(function(d) { return d.calc_rt*size/maxY ; } )" << endl;
  }
  out << "  .shape(function(d) { if (d.outlier == 1) { return 'cross'; } else { return 'circle' ;} } )" << endl;
  out << "  .strokeStyle(function(d) { if (d.outlier == 1) { return  pv.rgb(128,0,0, 0.8) ; } else { return  pv.rgb(128,128,128, 0.8); ;} } )" << endl;
  out << "  .title(function(d) { return d.pep; } );" << endl;

  out << "var name;" << endl;
  out << "var slope;" << endl;
  out << "var intercept;" << endl;
  out << "var rsq;" << endl;

  out << "var lreg = area.add(pv.Line)" << endl;
  out << "  .data(function(d) {  rsq = d.rsq; name = d.run; slope = d.slope; intercept = d.intercept;  " << endl;
  out << "		       return scale.map(function(x) {  while (slope*x+intercept < 0 && slope > 0) { x++; }  while (slope*x+intercept > 120" //rt_cat_->acn_gradient_->getAcn(4800)
      << " && slope > 0) { x--; }  while (slope*x+intercept < 0 && slope < 0) { x--; } return { x: x, y: d.slope*x+d.intercept };  ; } ) ; " << endl;
  out << "		    } )" << endl;
  if (irt_peps_hash_->size()==0) {
    out << "  .bottom(function(d) { return d.y*size/maxRT;; } )" << endl;  
  }
  else {
    out << "  .bottom(function(d) { return d.y*size/maxY;; } )" << endl;   
  }

  out << "  .left(function(d) { return d.x*size/maxRT; ;} )" << endl;
  out << "  .title(function(d) { return name + ', ' + 'y = ' + slope + '*x ' + ' + ' + intercept + ', r_sq = ' + rsq; } );;" << endl;

  out << "var xax = area.add(pv.Rule)" << endl;
  out << ".data(scale)" << endl;
  out << ".left(function(d) { return d*size/maxRT ; } )" << endl;
  out << ".strokeStyle(function(d) { return d ? 'silver' : 'black' ; } )" << endl;
  out << ".anchor('bottom')" << endl;
  out << ".add(pv.Label)" << endl;
  out << ".text(function(d) { return d; })" << endl;
  out << ".font('7px sans-serif');" << endl;

  out << "var xlab = area.anchor('bottom')" << endl;
  out << ".top(size + pad / 2)  " << endl;
  out << ".add(pv.Label)" << endl;
  out << ".text('Observed Median RT').font('7px sans-serif');" << endl;

  out << "var yax = area.add(pv.Rule)" << endl;
  if (irt_peps_hash_->size()==0) {
    out << ".data(scale)" << endl;
    out << ".bottom(function(d) { return d*size/maxRT ;  } )" << endl;
  }
  else {
    out << ".data(yscale)" << endl;
    out << ".bottom(function(d) { return d*size/maxY ;  } )" << endl;
  }

  out << ".strokeStyle(function(d) { return d ? 'silver' : 'black' ; } )" << endl;
  out << ".anchor('left')" << endl;
  out << ".add(pv.Label)" << endl;
  out << ".text(function(d) { return d; })" << endl;

  out << ".font('7px sans-serif');" << endl;

  out << "var ylab = area.anchor('left').add(pv.Label)" << endl;
  out << ".left(-20)" << endl;
  out << ".text('RTCalc').font('7px sans-serif')" << endl;
  out << ".textAngle(-Math.PI/2);" << endl;

  out << "var run_labl = area.anchor('top')" << endl;
  out << "   .add(pv.Label)" << endl;
  out << "   .top(-10)" << endl;
  out << "   .text(function(d) { return  d.run; } )" << endl;
  out << "   .font('7px sans-serif');" << endl;

  out << "vis.render();" << endl;

}

void RTCatalogParser::sortRunNames() {
  rt_cat_->sortRunNames();

}


void RTCatalogParser::trackPeptidesReport(string* file) {
  string res2 = *file + ".rtcat.REPORT.html";
  ofstream fout2(res2.c_str());
  rt_cat_->trackPeptidesReportPV(fout2, track_peps_, irt_peps_hash_->size()>0);
  fout2.close();
}

void RTCatalogParser::chromPeptidesReport(string* file) {
  string res2 = *file + ".rtcat.CHROM.html";
  Array<string*> oks;
  //  for (bool_hash::iterator it = ok_runs_->begin(); it != ok_runs_->end(); it++) {
  for (int i=0; i< rt_cat_->run_names_->size(); i++) {
    string* rn = new string(*(*rt_cat_->run_names_)[i]);
    oks.insertAtEnd(rn);
  }
  ofstream fout2(res2.c_str());
  rt_cat_->trackPeptidesChromatograms(fout2, &oks, track_pepq1q3_hash_);
  fout2.close();
}


void RTCatalogParser::trackPeptidesAcrossRuns(string* chromFile, string* trackFile) {
  Array<string*> oks;
  if (chromFile != NULL) {
    ifstream finC(chromFile->c_str());
  
    
    
    //string res1 = *file + ".rtcat.CHROMS.html";
    //  string res2 = *file + ".rtcat.REPORT.html";
    
    if (track_peps_ != NULL) {
      track_peps_->clear();
    }
    else { 
      track_peps_ = new Array<string*>();
    }
    
    
    if( track_pepq1q3_hash_ != NULL) {
      track_pepq1q3_hash_->clear();
    }
    else {
      track_pepq1q3_hash_ = new dblarr_hash();
    }
    
    string nextpep;
    string nextline;
    string q1q3;
    int pos;
    int cpos;
    double q1;
    double q3;
    dblarr_hash::iterator itr;
    while (1) {
      //fin >> nextpep;
      getline(finC, nextline);
      if (finC.fail()) break;
      
      pos = nextline.find_first_of('\t');
      if (pos!=string::npos) {
	nextpep = nextline.substr(0, pos);
	
	while (pos!=string::npos) {
	  
	  nextline = nextline.substr(pos+1);
	  pos = nextline.find_first_of('\t');
	  
	  if (pos!=string::npos) {
	    q1q3 = nextline.substr(0, pos);
	  }
	  else {
	    q1q3 = nextline;
	  }
	  
	  cpos = q1q3.find_first_of(',');
	  
	  q1 = atof(q1q3.substr(0, cpos).c_str());
	  q3 = atof(q1q3.substr(cpos+1).c_str());
	  
	  itr = track_pepq1q3_hash_->find(nextpep);
	  
	  if (itr == track_pepq1q3_hash_->end()) {
	    track_pepq1q3_hash_->insert(make_pair(*(new string(nextpep)), new Array<double>()));
	    itr = track_pepq1q3_hash_->find(nextpep);
	  }
	  
	  itr->second->insertAtEnd(q1);	itr->second->insertAtEnd(q3);
	  
	}
	
	
      }
      else {
	nextpep = nextline;
      }
      
      track_peps_->insertAtEnd(new string(nextpep));
      
      
    }

    for (dblarr_hash_hash::iterator itr1 = rtmix_peps_q1q3_->begin(); itr1 != rtmix_peps_q1q3_->end(); itr1++) {
      for (dblarr_hash::iterator itr2 = itr1->second->begin(); itr2 != itr1->second->end();  itr2++) {
	nextpep = itr2->first;

	itr = track_pepq1q3_hash_->find(nextpep);
	if (itr == track_pepq1q3_hash_->end()) {
	  track_pepq1q3_hash_->insert(make_pair(*(new string(nextpep)), new Array<double>()));
	  itr = track_pepq1q3_hash_->find(nextpep);
	}
	for (int x = 0; x < itr2->second->size(); x+=2) {
	  q1 = (*itr2->second)[x];	  q3 = (*itr2->second)[x+1];
	  itr->second->insertAtEnd(q1);	itr->second->insertAtEnd(q3);
	}
	track_peps_->insertAtEnd(new string(nextpep));
      }
    }
    
    for (bool_hash::iterator it = ok_runs_->begin(); it != ok_runs_->end(); it++) {
      string* okr = new string(it->first);
      oks.insertAtEnd(okr);
    }
    
    rt_cat_->trackPeptidesChromatograms(&oks, track_pepq1q3_hash_);

  }
  
  if (trackFile != NULL) {
    ifstream finT(trackFile->c_str());
    
    if (track_peps_ != NULL) {
      track_peps_->clear();
    }
    else { 
      track_peps_ = new Array<string*>();
    }
    
    string nextpep;
    string nextline;
    string q1q3;
    int pos;
    int cpos;
    double q1;
    double q3;

    while (1) {
      //fin >> nextpep;
      getline(finT, nextline);
      if (finT.fail()) break;
      
      pos = nextline.find_first_of('\t');
      if (pos!=string::npos) {
	nextpep = nextline.substr(0, pos);
	
	while (pos!=string::npos) {
	  
	  nextline = nextline.substr(pos+1);
	  pos = nextline.find_first_of('\t');
	  
	  if (pos!=string::npos) {
	    q1q3 = nextline.substr(0, pos);
	  }
	  else {
	    q1q3 = nextline;
	  }
	  
	  cpos = q1q3.find_first_of(',');
	  
	  q1 = atof(q1q3.substr(0, cpos).c_str());
	  q3 = atof(q1q3.substr(cpos+1).c_str());
	  
	  
	}
	
	
      }
      else {
	nextpep = nextline;
      }
      
      
      track_peps_->insertAtEnd(new string(nextpep));
      
      
    }
  }


  
  if (XICs_) {
    cerr << "INFO: Tracking XICS ... patience please ... " << endl;
    for (bool_hash::iterator it = ok_runs_->begin(); it != ok_runs_->end(); it++) {
      string* okr = new string(it->first);
      oks.insertAtEnd(okr);
    }

    rt_cat_->calcRTStatsByRun();
    if (irt_peps_hash_->size() == 0) {
      rt_cat_->calcRTStatsCombined(cerr);
    }
    
    rt_cat_->trackPeptidesXICs(&oks, byrun_pep_q1s_, XICtol_); //Last param is maxPPM set it to really narrow 1
  }
    

}


// void RTCatalogParser::trackPeptides(string* chromFile, string* trackFile) {
//   Array<string*> oks;
//   if (chromFile != NULL) {
//     ifstream finC(chromFile->c_str());
  
    
    
//     //string res1 = *file + ".rtcat.CHROMS.html";
//     //  string res2 = *file + ".rtcat.REPORT.html";
    
//     if (track_peps_ != NULL) {
//       track_peps_->clear();
//     }
//     else { 
//       track_peps_ = new Array<string*>();
//     }
    
    
//     if( track_pepq1q3_hash_ != NULL) {
//       track_pepq1q3_hash_->clear();
//     }
//     else {
//       track_pepq1q3_hash_ = new dblarr_hash();
//     }
    
//     string nextpep;
//     string nextline;
//     string q1q3;
//     int pos;
//     int cpos;
//     double q1;
//     double q3;
//     dblarr_hash::iterator itr;
//     while (1) {
//       //fin >> nextpep;
//       getline(finC, nextline);
//       if (finC.fail()) break;
      
//       pos = nextline.find_first_of('\t');
//       if (pos!=string::npos) {
// 	nextpep = nextline.substr(0, pos);
	
// 	while (pos!=string::npos) {
	  
// 	  nextline = nextline.substr(pos+1);
// 	  pos = nextline.find_first_of('\t');
	  
// 	  if (pos!=string::npos) {
// 	    q1q3 = nextline.substr(0, pos);
// 	  }
// 	  else {
// 	    q1q3 = nextline;
// 	  }
	  
// 	  cpos = q1q3.find_first_of(',');
	  
// 	  q1 = atof(q1q3.substr(0, cpos).c_str());
// 	  q3 = atof(q1q3.substr(cpos+1).c_str());
	  
// 	  itr = track_pepq1q3_hash_->find(nextpep);
	  
// 	  if (itr == track_pepq1q3_hash_->end()) {
// 	    track_pepq1q3_hash_->insert(make_pair(*(new string(nextpep)), new Array<double>()));
// 	    itr = track_pepq1q3_hash_->find(nextpep);
// 	  }
	  
// 	  itr->second->insertAtEnd(q1);	itr->second->insertAtEnd(q3);
	  
// 	}
	
	
//       }
//       else {
// 	nextpep = nextline;
//       }
      
//       track_peps_->insertAtEnd(new string(nextpep));
      
      
//     }

//     for (dblarr_hash_hash::iterator itr1 = rtmix_peps_q1q3_->begin(); itr1 != rtmix_peps_q1q3_->end(); itr1++) {
//       for (dblarr_hash::iterator itr2 = itr1->second->begin(); itr2 != itr1->second->end();  itr2++) {
// 	nextpep = itr2->first;

// 	itr = track_pepq1q3_hash_->find(nextpep);
// 	if (itr == track_pepq1q3_hash_->end()) {
// 	  track_pepq1q3_hash_->insert(make_pair(*(new string(nextpep)), new Array<double>()));
// 	  itr = track_pepq1q3_hash_->find(nextpep);
// 	}
// 	for (int x = 0; x < itr2->second->size(); x+=2) {
// 	  q1 = (*itr2->second)[x];	  q3 = (*itr2->second)[x+1];
// 	  itr->second->insertAtEnd(q1);	itr->second->insertAtEnd(q3);
// 	}
// 	track_peps_->insertAtEnd(new string(nextpep));
//       }
//     }
    
//     for (bool_hash::iterator it = ok_runs_->begin(); it != ok_runs_->end(); it++) {
//       string* okr = new string(it->first);
//       oks.insertAtEnd(okr);
//     }
    
//     rt_cat_->trackPeptidesChromatograms(&oks, track_pepq1q3_hash_);

//   }
  
//   if (trackFile != NULL) {
//     ifstream finT(trackFile->c_str());
    
//     if (track_peps_ != NULL) {
//       track_peps_->clear();
//     }
//     else { 
//       track_peps_ = new Array<string*>();
//     }
    
//     string nextpep;
//     string nextline;
//     string q1q3;
//     int pos;
//     int cpos;
//     double q1;
//     double q3;

//     while (1) {
//       //fin >> nextpep;
//       getline(finT, nextline);
//       if (finT.fail()) break;
      
//       pos = nextline.find_first_of('\t');
//       if (pos!=string::npos) {
// 	nextpep = nextline.substr(0, pos);
	
// 	while (pos!=string::npos) {
	  
// 	  nextline = nextline.substr(pos+1);
// 	  pos = nextline.find_first_of('\t');
	  
// 	  if (pos!=string::npos) {
// 	    q1q3 = nextline.substr(0, pos);
// 	  }
// 	  else {
// 	    q1q3 = nextline;
// 	  }
	  
// 	  cpos = q1q3.find_first_of(',');
	  
// 	  q1 = atof(q1q3.substr(0, cpos).c_str());
// 	  q3 = atof(q1q3.substr(cpos+1).c_str());
	  
	  
// 	}
	
	
//       }
//       else {
// 	nextpep = nextline;
//       }
      
      
//       track_peps_->insertAtEnd(new string(nextpep));
      
      
//     }
//   }


// }



void RTCatalogParser::correctRTsByRun() {
  cout << "Computing Gradient Corrections" << endl;
  rt_cat_->calcRTStatsByRun();
  if (irt_peps_hash_->size() == 0) {
    rt_cat_->calcRTStatsCombined(cerr);
  }
  for (int i=0; i<rt_cat_->num_runs_; i++) {
    byrun_rt_calcs_->insert(make_pair(*(*rt_cat_->run_names_)[i], new RTCalculator((*rt_cat_->run_names_)[i])));
    //cerr << "DEBUG: run_name = " << *(*rt_cat_->run_names_)[i] << endl;
    for (rtinfo_hash::iterator it  = (*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]]->begin();  
	 it != (*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]]->end(); it++) {

      bool_hash::iterator it2 = grad_peps_hash_->find(it->first);

      if ( grad_peps_hash_->size() == 0 || it2 != grad_peps_hash_->end()) {
	(*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->addPeptide_RT((const char*)it->first.c_str(),
								      (const char*)it->first.c_str(), 
								      -1, 
								      it->second->med_, 
								      (*rt_cat_->run_names_)[i]
								      );
      }

    }
    double rsq = -1;
    //(*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->read_RTcoeff();
    if (irt_peps_hash_->size()==0) {
      (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->batchRTCatalog(*rt_cat_);
      rsq = (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->calc_GradientCorrection();
    }
    else {
      (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->batch_iRT(*irt_peps_hash_);
      //rsq = (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->calc_GradientCorrection(0.05, -60);
      //rsq = (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->calc_GradientCorrection(0.04, 0);
      rsq = (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->calc_GradientCorrection();
    }
    //(*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->batchRTCalc(rt_cat_->acn_gradient_);
    //(*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->batchRTCalc();
    
    //(*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->batchSSRCalcHP();
    // if ( (*rt_cat_->run_names_)[i]->find("NG_DP_Reto_1sty-NG_reto_1668_13b_10_429") != string::npos) {
    //  cerr << "DDS: DEBUG" <<endl;
    //}
    
  //  double rsq = (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->calc_GradientCorrection();
    double slope = (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->getRunSlope();
    double intercept = (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->getRunInt();
    //double rsq = (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->calc_GradientOffset();
    if ( (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->getAllRTCount() == 1) {
      rsq = 1;
    }
    
    if (rsq < 0.3 || (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->getUsedRTCount() < minGradPeps_) {
      cout << "WARNING: Gradient Correction Skipped for File: " << (*rt_cat_->run_names_)[i]->c_str() << ", Rsq = " << rsq << endl;
      slope = 0.;
      intercept = -1e9; 
    }
    (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->setRunSlope(slope);
    (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->setRunInt(intercept);

	
    cout << "Gradient Correction for File: " << (*rt_cat_->run_names_)[i]->c_str()  << endl;
    cout << "\tRTCalc(rt_Obs) = A * rt_Obs + B, where A = " <<slope << ", B = " <<intercept << ", Rsq = " << rsq << endl;
    
    
    //(*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->InvertLine();
    
    for (dblarr_hash::iterator dit = (*(*rt_cat_->byrun_peprts_hash_)[*(*rt_cat_->run_names_)[i]]).begin(); dit !=  (*(*rt_cat_->byrun_peprts_hash_)[*(*rt_cat_->run_names_)[i]]).end(); dit++) {
      //	  double RTCalc = (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->calc_PepRT((string& )dit->first);
      
      //if (dit->first=="HWYITTGPVREK") {
      //	cerr << "DDS: DEBUG" << endl;
      //}
      
      //correct the points
      for (int j = 0; j < dit->second->size(); j++) {
	double corrected =slope * (*dit->second)[j] +intercept;
	(*(*(*rt_cat_->byrun_peprts_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first])[j] = corrected;
	
      }

      double med = (*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).med_;
      (*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).origMed_ = (*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).med_;
      (*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).origMean_ = (*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).mean_;
      //Correct RTmed for peptides by run
      (*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).med_ = slope *    
	(*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).med_  +intercept;      
      (*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).min_ = slope *    
	(*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).min_  +intercept;      
      (*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).mean_ = slope *    
	(*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).mean_  +intercept;      
      (*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).q25_ = slope *    
	(*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).q25_  +intercept;      
      (*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).q75_ = slope *    
	(*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).q75_  +intercept;      
      (*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).ComputeSIQR();
    }
 
 
        
  }

  
}

void RTCatalogParser::calcRTsByRun() {
   cout << "Computing BYRUN RT Statistics" << endl;
  rt_cat_->calcRTStatsByRun();
}

void RTCatalogParser::correctRTsByAdjacentRun() {
  str_hash::iterator  prevrun_itr, nextrun_itr;
  rt_cat_->calcRTStatsByRun();
  rt_cat_->calcRTStatsCombined(cerr);
  string* mappedrun = NULL;
  const string* mixTag = NULL;
  size_t pos1, pos2;
  string run_name = "";

  for (int i=0; i<rt_cat_->num_runs_; i++) {
    byrun_rt_calcs_->insert(make_pair(*(*rt_cat_->run_names_)[i], new RTCalculator((*rt_cat_->run_names_)[i])));
    //(*rt_cat_->byrun_peprtinfo_hash_)[i];

    
    //pos1 =  (*rt_cat_->run_names_)[i]->find_last_of("/\\");
    //pos2 =  (*rt_cat_->run_names_)[i]->find_last_of(".");
    run_name =  *(*rt_cat_->run_names_)[i];
	
    

    //run_itr = run_map->find(*(*rt_cat_->run_names_)[i]);

    prevrun_itr = prev_rtmix_byrun_->find(run_name);
    nextrun_itr = next_rtmix_byrun_->find(run_name);

    if (prevrun_itr != prev_rtmix_byrun_->end()) { //Prefer previous RTmix run
      mappedrun = &prevrun_itr->second;
    }
    else if (nextrun_itr != next_rtmix_byrun_->end()) { //Defer to next RTmix run
      mappedrun = &nextrun_itr->second;
    }
    else {
      mappedrun = (*rt_cat_->run_names_)[i];
    }

    

    for (str_hash::iterator itr = rtmix_chroms_->begin(); itr !=rtmix_chroms_->end(); itr++) {
      if (mappedrun->find(itr->first) != string::npos) {
	mixTag = &itr->first;
      }
    }

    for (rtinfo_hash::iterator it  = (*rt_cat_->byrun_peprtinfo_hash_)[*mappedrun]->begin();  
	 it != (*rt_cat_->byrun_peprtinfo_hash_)[*mappedrun]->end(); it++) {
      

      if (grad_peps_hash_->size() > 0) {
	bool_hash::iterator it2 = grad_peps_hash_->find(it->first);
	
	if (irt_peps_hash_->size()==0) {
	  if (it2 != grad_peps_hash_->end()) {
	    (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->addPeptide_RT((const char*)it->first.c_str(),
									  (const char*)it->first.c_str(), 
									  -1, it->second->med_, (*rt_cat_->run_names_)[i]);
	  }
	}
	else {
	  dbl_hash::iterator it3 = irt_peps_hash_->find(it->first);
	  if (it2 != grad_peps_hash_->end()&& it3 != irt_peps_hash_->end()) {
	    (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->addPeptide_RT((const char*)it->first.c_str(),
									  (const char*)it->first.c_str(), 
									  -1, 
									it3->second, 
									  (*rt_cat_->run_names_)[i]
									  );
	    
	  }
	}
      }
      else {
	//get them from the RTMix
	dblarr_hash::iterator it2 = (*rtmix_peps_q1q3_)[*mixTag]->find(it->first);
	if (it2 != (*rtmix_peps_q1q3_)[*mixTag]->end()) {
	  (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->addPeptide_RT((const char*)it->first.c_str(),
									(const char*)it->first.c_str(), 
									  -1, it->second->med_, (*rt_cat_->run_names_)[i]);
	  }

      }

    }

    (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->read_RTcoeff();

    //    (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->batchRTCalc(rt_cat_->acn_gradient_);
    //    (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->batchRTCalc();
    
    (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->batchRTCatalog(*rt_cat_);

    //(*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->batchSSRCalcHP();
    double rsq = -1;
    if (irt_peps_hash_->size()!=0) {
      rsq = (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->calc_GradientCorrection(0.04, 0);
    }
    else {
      rsq = (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->calc_GradientCorrection();
    }
    double slope = (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->getRunSlope();
    double intercept = (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->getRunInt();
    //double rsq = (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->calc_GradientOffset();
    if ( (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->getAllRTCount() == 1) {
      rsq = 1;
    }
    
    if (rsq < 0.45) {// || (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->getUsedRTCount() < minGradPeps_) {
      cout << "WARNING: Gradient Correction Skipped for File: " << (*rt_cat_->run_names_)[i]->c_str() << ", Rsq = " << rsq << endl;
      slope = 0.;
      intercept = -1e9; 
    }
    (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->setRunSlope(slope);
    (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->setRunInt(intercept);

	
    cout << "Gradient Correction for File: " << (*rt_cat_->run_names_)[i]->c_str()  << endl;
    cout << "\tRTCalc(rt_Obs) = A * rt_Obs + B, where A = " <<slope << ", B = " <<intercept << ", Rsq = " << rsq << endl;
    
    
    //(*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->InvertLine();
    
    for (dblarr_hash::iterator dit = (*(*rt_cat_->byrun_peprts_hash_)[*(*rt_cat_->run_names_)[i]]).begin(); dit !=  (*(*rt_cat_->byrun_peprts_hash_)[*(*rt_cat_->run_names_)[i]]).end(); dit++) {
      //	  double RTCalc = (*byrun_rt_calcs_)[*(*rt_cat_->run_names_)[i]]->calc_PepRT((string& )dit->first);
      
      //if (dit->first=="HWYITTGPVREK") {
      //	cerr << "DDS: DEBUG" << endl;
      //}
      
      //correct the points
      for (int j = 0; j < dit->second->size(); j++) {
	double corrected =slope * (*dit->second)[j] +intercept;
	(*(*(*rt_cat_->byrun_peprts_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first])[j] = corrected;
	
      }

      //Correct RTmed for peptides by run
      if (*(*rt_cat_->run_names_)[i] != *mappedrun) {
        (*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).med_ = slope *   (*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).med_  +intercept;      
        (*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).min_ = slope *   (*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).min_  +intercept;      
        (*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).mean_ = slope *   (*(*(*rt_cat_->byrun_peprtinfo_hash_)[*(*rt_cat_->run_names_)[i]])[dit->first]).mean_  +intercept;      
      }
    }
 
 
        
  }

  
}



void RTCatalogParser::setUpRTMixRuns() {
  

    string nextline;
    int pos;
    int j;
    int allRun_i=0;
    string prev_rtmix = "";
    string next_rtmix = "";
    string wl_name="";
    int i = 0;

    size_t pos1, pos2;



    prev_rtmix_byrun_ = new str_hash();    
    next_rtmix_byrun_ = new str_hash();    
    Array<string*>* allRuns = new Array<string*>();
    
    str_hash::iterator it;
    dblarr_hash_hash::iterator itr1;

    //Update worklist name with runname from pepXML file for non-mixes
    for (int x= 0; x < worklist_runs_->size() ; x++) {
      
      for (int y=0; y < (*worklist_runs_)[x]->size(); y++) {
	pos1 =  (*(*worklist_runs_)[x])[y]->find_last_of("/\\");
	pos2 =  (*(*worklist_runs_)[x])[y]->find_last_of(".");
	wl_name = (*(*worklist_runs_)[x])[y]->substr(pos1+1, pos2-pos1-1);
	for (int z=0; z < rt_cat_->run_names_->size(); z++) {
	  if ( (*rt_cat_->run_names_)[z]->find(wl_name) != string::npos) {
	    (*(*worklist_runs_)[x])[y] = (*rt_cat_->run_names_)[z];
	    break;
	  }	  
	}	
      }
    }

  
    for (int x= 0; x < worklist_runs_->size() ; x++) {
      prev_rtmix = "";
      allRuns->clear();
      for (int y=0; y < (*worklist_runs_)[x]->size(); y++) {
	nextline = *(*(*worklist_runs_)[x])[y];
	// pos1 =  (*(*worklist_runs_)[x])[y]->find_last_of("/\\");
	//pos2 =  (*(*worklist_runs_)[x])[y]->find_last_of(".");
	//wl_name = (*(*worklist_runs_)[x])[y]->substr(pos1+1, pos2-pos1-1);
	wl_name = *(*(*worklist_runs_)[x])[y];
	
	allRuns->insertAtEnd(new string(wl_name));
	//rt_cat_->addRun(nextline);
	
	for (itr1 = rtmix_peps_q1q3_->begin(); itr1 !=  rtmix_peps_q1q3_->end(); itr1++) {
	  pos = nextline.find(itr1->first.c_str());
	  
	  if (pos!=string::npos) {
	    prev_rtmix = wl_name;
	    break;
	  }
	}
	
	if (pos==string::npos && !prev_rtmix.empty()) {
	  prev_rtmix_byrun_->insert(make_pair(*(*allRuns)[y], prev_rtmix));
	}      
	
      }
    


      next_rtmix = "";
      for (int y= (*worklist_runs_)[x]->size()-1 ; y >= 0; y--) {
	nextline = *(*(*worklist_runs_)[x])[y];
	// pos1 =  (*(*worklist_runs_)[x])[y]->find_last_of("/\\");
	//pos2 =  (*(*worklist_runs_)[x])[y]->find_last_of(".");
	//wl_name = (*(*worklist_runs_)[x])[y]->substr(pos1+1, pos2-pos1-1);
	wl_name = *(*(*worklist_runs_)[x])[y];
	
	//allRuns->insertAtEnd(new string(wl_name));
	//rt_cat_->addRun(nextline);
	
	for (itr1 = rtmix_peps_q1q3_->begin(); itr1 !=  rtmix_peps_q1q3_->end(); itr1++) {
	  pos = nextline.find(itr1->first.c_str());
	  
	  if (pos!=string::npos) {
	    next_rtmix = wl_name;
	    break;
	  }
	}
	
	if (pos==string::npos && !next_rtmix.empty()) {
	  next_rtmix_byrun_->insert(make_pair(*(*allRuns)[y], next_rtmix));
	}      
	
      }
    }
      
      
      //    prev_rtmix = next_rtmix;
      //     next_rtmix = "";
    
//     for (i = i; i < allRuns->size(); i++) {
      
//       it = prev_rtmix_byrun_->find(*(*allRuns)[i]);
      
//       if (it == prev_rtmix_byrun_->end() && !prev_rtmix.empty()) {
// 	prev_rtmix_byrun_->insert(make_pair(*(*allRuns)[i], prev_rtmix));
	
//       }
      
//       it = next_rtmix_byrun_->find(*(*allRuns)[i]);
      
//       if (it == next_rtmix_byrun_->end() && !next_rtmix.empty()) {
// 	next_rtmix_byrun_->insert(make_pair(*(*allRuns)[i], next_rtmix));
	
//       }
	  
//     }
	
    //Ignore Runs not in the WorkList
    bool present = false;
    for (j=0; j<rt_cat_->num_runs_; j++) {
      present = false;
      for (i=0; i< allRuns->size(); i++) {
	if ((*rt_cat_->run_names_)[j]->find((*allRuns)[i]->c_str()) != string::npos) {
	  present = true;
	  break;
	}
      }
      if (!present) {
	if (ignore_runs_==NULL) 
	  ignore_runs_ = new Array<string*>();
	ignore_runs_->insertAtEnd((*rt_cat_->run_names_)[j]);
      }
    }
    
      
//     int thisRun=-1;
//     int prevRTmix=-1;
//     int nextRTmix=-1;
//     for (i=0; i< allRuns->size(); i++) {
//       thisRun = -1;
//       prevRTmix = -1;
//       nextRTmix = -1;
//       for (j=0; j<rt_cat_->num_runs_; j++) {
// 	pos = (*rt_cat_->run_names_)[j]->find((*allRuns)[i]->c_str());
      
// 	if (pos != string::npos) {
// 	  thisRun = j;
// 	}
	
// 	it = next_rtmix_byrun_->find(*(*allRuns)[i]);
	
// 	if (it != next_rtmix_byrun_->end()) {

// 	  pos = (*rt_cat_->run_names_)[j]->find((*next_rtmix_byrun_)[*(*allRuns)[i]].c_str());
	  
// 	  if (pos != string::npos) {
// 	    nextRTmix = j;
// 	  }
// 	}

// 	it = prev_rtmix_byrun_->find(*(*allRuns)[i]);
	
// 	if (it != prev_rtmix_byrun_->end()) {
// 	  pos = (*rt_cat_->run_names_)[j]->find((*prev_rtmix_byrun_)[*(*allRuns)[i]].c_str());
	  
// 	  if (pos != string::npos) {
// 	    prevRTmix = j;
// 	  }
// 	}

//       }


//       if (thisRun != -1) {
// 	  //	  prev_rtmix_byrun_->remove(make_pair(*(*allRuns)[i], prev_rtmix));
	
// 	//pos1 =  (*rt_cat_->run_names_)[thisRun]->find_last_of("/\\");
// 	//pos2 =  (*rt_cat_->run_names_)[thisRun]->find_last_of(".");
// 	//wl_name =  (*rt_cat_->run_names_)[thisRun]->substr(pos1+1, pos2-pos1-1);
	
// 	wl_name =  *(*rt_cat_->run_names_)[thisRun];

// 	if (prevRTmix  != -1) {
// 	  //	  prev_rtmix_byrun_->insert(make_pair(*(*rt_cat_->run_names_)[thisRun], *(*rt_cat_->run_names_)[prevRTmix]));
	  
	

// 	  //pos1 =  (*rt_cat_->run_names_)[prevRTmix]->find_last_of("/\\");
// 	  //pos2 =  (*rt_cat_->run_names_)[prevRTmix]->find_last_of(".");
// 	  //prev_rtmix =  (*rt_cat_->run_names_)[prevRTmix]->substr(pos1+1, pos2-pos1-1);
// 	  prev_rtmix =  *(*rt_cat_->run_names_)[prevRTmix];


// 	  prev_rtmix_byrun_->insert(make_pair(wl_name, prev_rtmix));
// 	}
	
// 	else if (nextRTmix  == -1)  {
// 	  //prev_rtmix_byrun_->insert(make_pair(*(*rt_cat_->run_names_)[thisRun], *(*rt_cat_->run_names_)[thisRun]));
// 	  prev_rtmix_byrun_->insert(make_pair(wl_name, wl_name));
// 	}	
// 	else {
// 	  //prev_rtmix_byrun_->insert(make_pair(*(*rt_cat_->run_names_)[thisRun], *(*rt_cat_->run_names_)[thisRun]));
// 	  next_rtmix =  *(*rt_cat_->run_names_)[nextRTmix];	 
// 	  prev_rtmix_byrun_->insert(make_pair(wl_name, next_rtmix));
// 	}	

// 	if (nextRTmix  != -1) {

// 	  //pos1 =  (*rt_cat_->run_names_)[nextRTmix]->find_last_of("/\\");
// 	  //pos2 =  (*rt_cat_->run_names_)[nextRTmix]->find_last_of(".");
// 	  //next_rtmix =  (*rt_cat_->run_names_)[nextRTmix]->substr(pos1+1, pos2-pos1-1);

// 	  next_rtmix =  *(*rt_cat_->run_names_)[nextRTmix];

// 	  //next_rtmix_byrun_->insert(make_pair(*(*rt_cat_->run_names_)[thisRun], *(*rt_cat_->run_names_)[nextRTmix]));

// 	  next_rtmix_byrun_->insert(make_pair(wl_name, next_rtmix));
// 	}
// 	else if (prevRTmix  == -1)  {
// 	  //next_rtmix_byrun_->insert(make_pair(*(*rt_cat_->run_names_)[thisRun], *(*rt_cat_->run_names_)[thisRun]));
// 	  next_rtmix_byrun_->insert(make_pair(wl_name, wl_name));
// 	}
// 	else {
// 	  //prev_rtmix_byrun_->insert(make_pair(*(*rt_cat_->run_names_)[thisRun], *(*rt_cat_->run_names_)[thisRun]));
// 	  prev_rtmix =  *(*rt_cat_->run_names_)[prevRTmix];	 
// 	  next_rtmix_byrun_->insert(make_pair(wl_name, prev_rtmix));
// 	}		

//       }
//       delete (*allRuns)[i];
//     }
  
    

}
