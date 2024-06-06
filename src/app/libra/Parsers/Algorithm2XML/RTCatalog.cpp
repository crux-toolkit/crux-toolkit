/*

Program       : RTCatalog
Author        : David Shteynberg <dshteynb  AT systemsbiology.org>                                                       
Date          : 09.29.2010



Copyright (C) 2010 David Shteynberg

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
\
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

#include "RTCatalog.h"
using namespace std;
//using namespace pwiz::msdata;
using namespace mzParser;

bool compare_mz (  MZIntensityPair mz1,  MZIntensityPair mz2) {
  return (mz1.mz < mz2.mz);
}

bool compare_runs( const string* elem1, const string * elem2 ){
  size_t pos1, pos2;
  // pos1 = elem1->find_first_of("_-.");
  //pos2 = elem2->find_first_of("_-.");

  //  if (pos1 != string::npos && pos2 !=  string::npos) {
  //  min = pos1 < pos2 ? pos1 : pos2;
  //}
  //else {
  //  min = elem1->length() <  elem2->length() ? elem1->length() :  elem2->length();
  //}

  //return strncmp(elem1->c_str(), elem2->c_str(), min) < 0; // longest first
    
    string s1 ;
    string s2 ;
    
    pos1 = elem1->find_last_of("/\\");
    if (pos1 == string::npos) {
      pos1 = elem1->length();
    }
    
    pos2 = elem2->find_last_of("/\\");
    if (pos2 == string::npos) {
      pos2 = elem2->length();
    }
    
    s1 = elem1->substr(pos1);
    s2 = elem2->substr(pos2);

    return std::lexicographical_compare(s1.c_str(), s1.c_str()+pos1, 
					s2.c_str(), s2.c_str()+pos2);

}


RTCatalog::RTCatalog(double minProb, double minRT, double maxRT, string * acnFile) {
  dalTOL_ = false;
  minProb_ = minProb;
  //  peprts_hash_ = new Array<dblarr_hash*>();
  //pepintens_hash_ = new Array<dblarr_hash*>();

  byrun_peprts_hash_ = new dblarr_hash_hash();
  byrun_pep_maxIntens_ = new dbl_hash_hash();
  byrun_pep_maxMatchedIons_ = new dbl_hash_hash();
  byrun_pepintens_hash_ = new dblarr_hash_hash();

  
  byrun_pepMatchedIons_hash_ = new dblarr_hash_hash();

  
  run_names_ = new Array<string*>();

  run_files_ = new str_hash();

  pepruns_hash_ = new strparr_hash();

  ismix_run_hash_ = new bool_hash();

  //  rt_calc_ = new RTCalculator();

  peprtinfo_hash_  = new rtinfo_hash();

  mod_table_ = new mod_hash();

  //  byrun_peprtinfo_hash_  = new Array<rtinfo_hash*>();
  byrun_peprtinfo_hash_  = new rtinfo_hash_hash();

  acn_gradient_ = NULL;
  
  if (acnFile!=NULL) {
    acn_gradient_ = new GradientProgram(acnFile->c_str());
  }
  else {
    acn_gradient_ = new GradientProgram();
  }

  minRT_ = minRT;
  maxRT_ =  maxRT;
  num_runs_ = 0;
  pv_mods_ = false;
  cys_cam_ = false;
}

RTCatalog::RTCatalog(const char* file) {
  num_runs_ = -1;
  dalTOL_ = false;
  pv_mods_ = false;
  cys_cam_ = false;
  peprtinfo_hash_  = new rtinfo_hash();
  mod_table_ = new mod_hash();
  run_names_ = new Array<string*>();
  
  rtinfo_hash::iterator it;

  string pep;
  int n;
  double med, q25, q75, siqr, mean, min, stdev;
  ifstream fin(file);

  for (int i=0; i<9; i++) {
    string tmp;
    fin >> tmp;
    if (fin.fail()) break;
  }

  while (1) {
    fin >> pep >> q25 >> med >> q75 >> siqr >> mean >> stdev >> min >> n; 
    if (fin.fail()) break;
    
    it = peprtinfo_hash_->find(pep);
    if (it == peprtinfo_hash_->end()) {
      peprtinfo_hash_->insert(make_pair(pep, new RTInfo(n, q25, med, q75, mean, stdev, min)));
    }
    else {
      cerr << "WARNING: Catalog contains multiple entries for peptide " << pep << endl;
    }

  }

}

int RTCatalog::getRTCount(string& pep) {
  int rtn = -1;
  rtinfo_hash::iterator it= peprtinfo_hash_->find(pep);
  if (it != peprtinfo_hash_->end()) { 
    rtn = it->second->n_;
  }
  return rtn;
  
}

void  RTCatalog::setMods(bool pv, bool cys) {
  if (pv) {
    initMods();
  }
  pv_mods_ = pv;
  cys_cam_ = cys;
}

string RTCatalog::getPeptideString(const string& in_pep) {
  string out_pep = "";
  string token;
  if (!pv_mods_ && !cys_cam_) {
    return in_pep;
  }
  size_t pos = 0;
  size_t cl_bracket = string::npos;
  size_t op_bracket = string::npos;
  for (pos = 0; pos < in_pep.length() ;) {
    op_bracket = in_pep.find('[', pos);
    
    if (op_bracket == pos + 1) {
      cl_bracket = in_pep.find(']', pos);
    }
    else {
      cl_bracket = string::npos;
    }

    if (cl_bracket != string::npos) {
      token = in_pep.substr(pos, cl_bracket-pos+1);
    }
    else {
      token = in_pep.substr(pos, 1);
    }
    
    pos += token.length();
    if (token.length() == 1  && 
	token == "C" && cys_cam_) {
      token = "C[CAM]";
    }
    if (mod_table_->find(token) ==  mod_table_->end()) {
      out_pep += token;
    }
    else {
      out_pep += (*mod_table_)[token];
    }

  }
  return out_pep;
}

RTInfo* RTCatalog::getRTInfo(string& pep) {
  RTInfo* rtn = NULL;
  rtinfo_hash::iterator it= peprtinfo_hash_->find(pep);
  if (it != peprtinfo_hash_->end()) { 
    rtn = it->second;
  }
  return rtn;
  
}

int RTCatalog::getNumRuns() {
  return num_runs_;
}

void RTCatalog::initMods() {
  mod_table_->insert(make_pair("C[160]", "C[CAM]"));
  mod_table_->insert(make_pair("M[147]", "M[Oxi]"));
  mod_table_->insert(make_pair("K[136]", "K[+08]"));
  mod_table_->insert(make_pair("R[166]", "R[+10]"));
  mod_table_->insert(make_pair("E[111]", "E[PGE]"));
  mod_table_->insert(make_pair("Q[111]", "Q[PGQ]"));
  mod_table_->insert(make_pair("C[143]", "C[PCm]"));
  mod_table_->insert(make_pair("n[43]", "[1Ac]-"));
  mod_table_->insert(make_pair("S[167]", "S[Pho]"));
  mod_table_->insert(make_pair("T[181]", "T[Pho]"));
  mod_table_->insert(make_pair("Y[243]", "Y[Pho]"));
 

}

RTInfo* RTCatalog::getRTInfo(string& pep, string& run) {
  RTInfo* rtn = NULL;

  rtinfo_hash_hash::iterator it1= byrun_peprtinfo_hash_->find(run);

  if (it1 == byrun_peprtinfo_hash_->end()) {
    return NULL;
  }
  
  rtinfo_hash::iterator it2= (*byrun_peprtinfo_hash_)[run]->find(pep);

  if (it2 == (*byrun_peprtinfo_hash_)[run]->end() ) { 
    return NULL;
  }
  rtn = it2->second;
  return rtn;
  
}

double RTCatalog::getRTMedian(string& pep) const {
  double rtn = -1;
  rtinfo_hash::iterator it= peprtinfo_hash_->find(pep);
  if (it != peprtinfo_hash_->end()) { 
    rtn = it->second->med_;
  }
  return rtn;
}

double RTCatalog::getRTMin(string& pep) const {
  double rtn = -1;
  rtinfo_hash::iterator it= peprtinfo_hash_->find(pep);
  if (it != peprtinfo_hash_->end()) { 
    rtn = it->second->min_;
  }
  return rtn;
}

double RTCatalog::getRTSIQR(string& pep) const{
  double rtn = -1;
  rtinfo_hash::iterator it= peprtinfo_hash_->find(pep);
  if (it != peprtinfo_hash_->end()) { 
    rtn = it->second->siqr_;
  }
  return rtn;
}

double RTCatalog::getRTMean(string& pep) const {
  double rtn = -1;
  rtinfo_hash::iterator it= peprtinfo_hash_->find(pep);
  if (it != peprtinfo_hash_->end()) { 
    rtn = it->second->mean_;
  }
  return rtn;
}

double RTCatalog::getRTStdev(string& pep) const {
  double rtn = -1;
  rtinfo_hash::iterator it= peprtinfo_hash_->find(pep);
  if (it != peprtinfo_hash_->end()) { 
    rtn = it->second->stdev_;
  }
  return rtn;

}

RTCatalog::~RTCatalog() {
  //TODO: write me!
  dblarr_hash_hash::iterator it1;
  dblarr_hash::iterator it2;


  for (it1 = byrun_peprts_hash_->begin(); it1 != byrun_peprts_hash_->end(); it1++) {
    for (it2 = it1->second->begin(); it2 != it1->second->end(); it2++) {
      delete it2->second;
    }
    delete it1->second;
  }
  
  
}


void RTCatalog::insertResult(const string& run, const string& modpep) {
  string* key = new string(modpep);

  size_t pos1, pos2;
  //pos1 =  run.find_last_of("/\\");
  //pos2 =  run.find_last_of(".");
  string rn = run;//.substr(pos1+1);

  dblarr_hash::iterator it = (*(*byrun_peprts_hash_)[rn]).find(*key);
  
  if (it == (*(*byrun_peprts_hash_)[rn]).end() ) {
    (*byrun_pep_maxIntens_)[rn]->insert(make_pair(*key, 0));
    (*byrun_pep_maxMatchedIons_)[rn]->insert(make_pair(*key, 0));
    (*byrun_peprts_hash_)[rn]->insert(make_pair(*key, new Array<double>()));

    (*byrun_pepMatchedIons_hash_)[rn]->insert(make_pair(*key, new Array<double>()));
    (*byrun_pepintens_hash_)[rn]->insert(make_pair(*key, new Array<double>()));
  }
 
  strparr_hash::iterator s_it = (*pepruns_hash_).find(*key);
  
  if (s_it == (*pepruns_hash_).end() ) {
    pepruns_hash_->insert(make_pair(*key, new Array<string*>()));
    (*(*pepruns_hash_)[*key]).insertAtEnd(new string(rn));
  }
  else {
    if ( *(*(*pepruns_hash_)[*key])[(*(*pepruns_hash_)[*key]).size()-1] != rn) {
      (*(*pepruns_hash_)[*key]).insertAtEnd(new string(rn));
    }
    
  }
}

bool RTCatalog::rejectResult(double prob, double rt) {
  return (prob < minProb_ || rt < minRT_ || rt > maxRT_);
}

void RTCatalog::insertResult(string& run, string& spectrum, double prob, 
			     Array<double>* allntt_prob, string& pepseq, 
			     string& modpep, double calcnmass, double rt, 
			     double prec_intens, double collision_eng, 
			     double matchedions_frac, string& exp_lbl, 
string& charge) {

  //  if (pepseq.find("DLFSVLK",0) != string::npos)
  // cerr << "DEBUG: this one !!! " << endl;

   if (prob < minProb_ || rt < minRT_ || rt > maxRT_) 
    return;

   //  rt = acn_gradient_->getAcn(rt);

   string* key = new string(modpep);

  dblarr_hash::iterator it = (*(*byrun_peprts_hash_)[run]).find(*key);

  if (it == (*(*byrun_peprts_hash_)[run]).end() ) {
    (*byrun_peprts_hash_)[run]->insert(make_pair(*key, new Array<double>()));
    (*byrun_pepintens_hash_)[run]->insert(make_pair(*key, new Array<double>()));
    (*byrun_pepMatchedIons_hash_)[run]->insert(make_pair(*key, new Array<double>()));
  }


  if ((*byrun_pep_maxIntens_)[run]->find(*key) != (*byrun_pep_maxIntens_)[run]->end()) {
    if ( (*(*byrun_pep_maxIntens_)[run])[*key] < prec_intens) {
      (*(*byrun_pep_maxIntens_)[run])[*key] = prec_intens;
    }
  }
  else {
    (*byrun_pep_maxIntens_)[run]->insert(make_pair(*key, prec_intens));
  }
  

  //TODO:  DDS: substitute here with reading it from the XIC os 
  (*(*(*byrun_peprts_hash_)[run])[*key]).insertAtEnd(acn_gradient_->getAcn(rt));
  (*(*(*byrun_pepintens_hash_)[run])[*key]).insertAtEnd(prec_intens);

  string *spkey = new string(spectrum);
  
  
  if ((*byrun_pep_maxMatchedIons_)[run]->find(*key) != (*byrun_pep_maxMatchedIons_)[run]->end()) {
    if ( (*(*byrun_pep_maxMatchedIons_)[run])[*key] < matchedions_frac) {
      (*(*byrun_pep_maxMatchedIons_)[run])[*key] = matchedions_frac;
    }
  }
  else {
    (*byrun_pep_maxMatchedIons_)[run]->insert(make_pair(*key, matchedions_frac));
  }
  
   
   
   
  (*(*(*byrun_pepMatchedIons_hash_)[run])[*key]).insertAtEnd(matchedions_frac);
  
  
  
  strparr_hash::iterator s_it = (*pepruns_hash_).find(*key);
  
  if (s_it == (*pepruns_hash_).end() ) {
    pepruns_hash_->insert(make_pair(*key, new Array<string*>()));
    (*(*pepruns_hash_)[*key]).insertAtEnd(new string(run));
  }
  else {
    if ( *(*(*pepruns_hash_)[*key])[(*(*pepruns_hash_)[*key]).size()-1] != run) {
      (*(*pepruns_hash_)[*key]).insertAtEnd(new string(run));
    }
    
  }

  //  else {
    //assume runs are entered in order
  //   if ( (*(*pepruns_hash_)[*key])[(*(*pepruns_hash_)[*key]).size()-1] < (*run_names_)[run_idx] ) {
  //    (*(*pepruns_hash_)[*key]).insertAtEnd((*run_names_)[run_idx]);
  //  }
    
  //}
  
  
}

void  RTCatalog::trackPeptidesChromatograms(Array<string*>* runs, dblarr_hash* peps_q1q3) {
  using namespace std;

  int i=0;
  int j=0;
  int k=0;

  MSDataFile* msd;
  ChromatogramListPtr sl;
  ChromatogramPtr s;
  string q1, q3;
  size_t pos;
  bool got_match = false;
  bool got_run_match = false;
  dblarr_hash::iterator itr;

  int last_run = -1;

  dblarr_hash* bychrome_rtsecs = new dblarr_hash(); 
  dblarr_hash* bychrome_intens = new dblarr_hash(); ;
  dbl_hash* bychrome_maxIntens = new dbl_hash();
  dbl_hash* bychrome_peakRT = new dbl_hash();


  Array<double>* runchrome_rtsecs = new Array<double>(); 
  Array<double>* runchrome_intens = new Array<double>(); ;



  //  runchrome_intWtdStdev_ = new dbl_hash();
  //  runchrome_intWtdMean_hash_ = new dbl_hash();

  string pepQ1Q3;
  string pep;
  string run_pepQ1Q3;
  
  double w = 0;
  double num = 0;
  double denom = 0;
  double mean = -1;
  double stdev = -1;
  double run_maxIntens = -1;

  double run_peakRT = -1;
  double run_sumIntens = 0;
  double all_sumIntens = 0;

  double all_maxIntens = 0;

  double all_num = 0;

  double all_denom = 0;

  double all_sumsq = 0;

  double sumsq = 0;

  for (k=0; k<runs->size(); k++) {
    string rn = *(*runs)[i];
    for (dblarr_hash::iterator it = (*(*byrun_peprts_hash_)[rn]).begin(); it != (*(*byrun_peprts_hash_)[rn]).end(); it++) {
      
      (*(*(*byrun_peprts_hash_)[rn])[it->first]).clear();
      (*(*(*byrun_pepintens_hash_)[rn])[it->first]).clear();


    }
  }

  for (i=0; i<runs->size(); i++) {
    got_run_match = false;
    //TODO: this means the directory names cannot have a '.'
    string run_name = (*run_files_)[*(*runs)[i]];
    pos = run_name.find_first_of('.');

    if (pos!=string::npos) {
      run_name = run_name.substr(0, pos);
    }
    
    try {
      msd = new MSDataFile(run_name);//+".mzML");
    }
    catch (...) {
      cerr << "WARNING: Unable to open file: " << (*run_files_)[*(*runs)[i]] << endl;
      continue;
    }
    
    

    if (!msd->run.chromatogramListPtr->get())
      throw runtime_error("[trackPeptidesChromatograms] Null chromatogramListPtr.");

    
    sl = msd->run.chromatogramListPtr;
    
    cout << "Parsing Data File: " << run_name <<  endl;
   
    for (j=1; j<(int)sl->size(); j++) {

      pepQ1Q3 = "";
      q1 = "";
      q3 = "";

      s = sl->chromatogram(j, true);

      q1 = s->id;
      pos = q1.find("Q1=");
      if (pos!=string::npos) {
	q1 = q1.substr(pos);
	pos = q1.find_first_of(" \t");
	if (pos!=string::npos) {
	  q1 = q1.substr(3, pos-3);
	}
      }
      q3 = s->id;
      pos = q3.find("Q3=");
      if (pos!=string::npos) {
	q3 = q3.substr(pos);
	pos = q3.find_first_of(" \t");
	if (pos!=string::npos) {
	  q3 = q3.substr(3, pos-3);
	}
      }

      itr = peps_q1q3->begin();

      while(itr != peps_q1q3->end()) {
	for (k=0; k<itr->second->size(); k+= 2) {
	  if (fabs((*itr->second)[k] - atof(q1.c_str())) < 0.05 && fabs((*itr->second)[k+1] - atof(q3.c_str())) < 0.05 ) {
	    got_match = true;
	    pos = q1.find(".");
	    if (pos!=string::npos) {
	      q1 = q1.substr(0, pos);
	    }
	    pos = q3.find(".");
	    if (pos!=string::npos) {
	      q3 = q3.substr(0, pos);
	    }
	    pepQ1Q3 += itr->first+q1+q3;
	    pep = itr->first;

	    break;
	  }
	}
	if (got_match) break;
	itr++;
      }
      //      run_pepQ1Q3 = (*run_names_)[i] + pepQ1Q3;
      got_run_match = got_match || got_run_match ;
      
    
      if (got_match) {



	run_maxIntens = -1;
	run_peakRT = -1;
	run_sumIntens = 0;
	runchrome_rtsecs->clear();
	runchrome_intens->clear();
	vector<TimeIntensityPair> pairs;
	s->getTimeIntensityPairs(pairs);
	for (size_t ix=0; ix < pairs.size(); ix++) {
	  if (pairs[ix].intensity > 0 && pairs[ix].time*60 <= maxRT_ && 
	      pairs[ix].time*60 >= minRT_) {


	    if ((*bychrome_maxIntens).find(pepQ1Q3.c_str()) != (*bychrome_maxIntens).end()) {
	      if ((*bychrome_maxIntens)[pepQ1Q3.c_str()] <  pairs[ix].intensity) {
		(*bychrome_maxIntens)[pepQ1Q3.c_str()] =  pairs[ix].intensity;
		(*bychrome_peakRT)[pepQ1Q3.c_str()] =  pairs[ix].time;
		
	      }
	    }
	    else {
	      (*bychrome_maxIntens)[pepQ1Q3.c_str()] =  pairs[ix].intensity;
	      (*bychrome_peakRT)[pepQ1Q3.c_str()] =  pairs[ix].time;
	      (*bychrome_rtsecs)[pepQ1Q3.c_str()] = new Array<double>();
	      (*bychrome_intens)[pepQ1Q3.c_str()] = new Array<double>();
	    }

	    run_sumIntens += pairs[ix].intensity;



	    if (run_maxIntens <  pairs[ix].intensity) {
	     run_maxIntens =  pairs[ix].intensity;
	     run_peakRT = pairs[ix].time;
	    }

	    if (all_maxIntens < run_maxIntens) {
	      all_maxIntens = run_maxIntens;
	    }
	 
	    runchrome_rtsecs->insertAtEnd(pairs[ix].time);
	    runchrome_intens->insertAtEnd(pairs[ix].intensity);

	    if ((*(*byrun_peprts_hash_)[*(*runs)[i]]).find(pep) != (*(*byrun_peprts_hash_)[*(*runs)[i]]).end()) {

	      double RTsec = pairs[ix].time*60;
	      (*(*byrun_peprts_hash_)[*(*runs)[i]])[pep]->insertAtEnd(acn_gradient_->getAcn(RTsec));
	      (*(*byrun_pepintens_hash_)[*(*runs)[i]])[pep]->insertAtEnd(pairs[ix].intensity);
	    }

	    

	  }
	}
	
	if ((*(*byrun_pep_maxIntens_)[*(*runs)[i]]).find(pep) !=(*(*byrun_pep_maxIntens_)[*(*runs)[i]]).end() ) {
	  if ((*(*byrun_pep_maxIntens_)[*(*runs)[i]])[pep] < run_maxIntens) {
	    (*(*byrun_pep_maxIntens_)[*(*runs)[i]])[pep] = run_maxIntens;
	  }
	}

      
	if (i != last_run) {
	  last_run = i;
	}

	if (bychrome_rtsecs->find(pepQ1Q3.c_str()) != bychrome_rtsecs->end()) {
	  (*bychrome_rtsecs)[pepQ1Q3.c_str()]->insertAtEnd(run_peakRT);
	  (*bychrome_intens)[pepQ1Q3.c_str()]->insertAtEnd(run_maxIntens);
	}
      
	sumsq = 0;
	num = 0;
	denom = 0;
	mean = -1;
	stdev = -1;
	
	//Compute wtdMean for each chromatoram in each run
	num = 0; denom = 0; 
	all_sumIntens += run_sumIntens;
	for (int idx=0; idx < runchrome_rtsecs->size(); idx++) {
	  w =  (*runchrome_intens)[idx] / run_maxIntens;
	  num += (*runchrome_rtsecs)[idx] * w; 
	  sumsq += (*runchrome_rtsecs)[idx] * (*runchrome_rtsecs)[idx] * w ;
	  denom += w;
	  
	}
	
	if (denom > 0) {
	  mean = num / denom;
	  stdev = sqrt(sumsq/denom - mean*mean);
	}
	
	//	all_num += run_maxIntens * run_peakRT;
	//all_sumsq += run_maxIntens * run_peakRT * run_peakRT;
	//all_denom += run_maxIntens;

	
      }
      got_match = false;
    }


    delete msd;
  }



//   sumsq = 0;
//   num = 0;
//   denom = 0;
//   mean = -1;
//   stdev = -1;

//   //compute intensity weighted (normalized against maxInt) RT mean
//   for (itr = bychrome_rtsecs->begin(); itr != bychrome_rtsecs->end(); itr++) {
//     num = 0; denom = 0;
    
//     for (int idx=0; idx < itr->second->size(); idx++) {
//       w =  (*(*bychrome_intens)[itr->first])[idx] / (*bychrome_maxIntens)[itr->first];
//       num += (*(*bychrome_rtsecs)[itr->first])[idx] * w; 
//       sumsq += (*(*bychrome_rtsecs)[itr->first])[idx] * (*(*bychrome_rtsecs)[itr->first])[idx] * w ;
//       denom += w;
//     }
//     if (denom > 0) {
//       mean = num / denom;
//       stdev = sqrt(sumsq/denom - mean*mean);
//     }


//   }



}

void RTCatalog::progress(int tic, int step, int &tot) {
  
  if (!step || tic % step == 0) {
    tot++;
    if (!step && tot % 10 == 0) {
      cerr << tot;
    }
    else if (tot % 10 == 0) {     
      cerr << tot << "%";
    }
    else {
      cerr << ".";
    }
  }
 
 

}



void  RTCatalog::trackPeptidesXICs(Array<string*>* runs, dblarr_hash_hash* byrun_peps_q1s, double maxPPM) {
  using namespace std;

  int i=0;
  int j=0;
  int k=0;

  MSDataFile* msd;
  SpectrumListPtr sl;// = *msd.run.spectrumListPtr;
  SpectrumPtr s;

  //  ChromatogramListPtr sl;
  //ChromatogramPtr s;
  string q1, q3;
  size_t pos;
  bool got_match = false;
  bool got_run_match = false;
  dblarr_hash::iterator itr;

  int last_run = -1;

  dblarr_hash* bychrome_rtsecs = new dblarr_hash(); 
  dblarr_hash* bychrome_intens = new dblarr_hash(); ;
  dbl_hash* bychrome_maxIntens = new dbl_hash();
  dbl_hash* bychrome_peakRT = new dbl_hash();


  Array<double>* runchrome_rtsecs = new Array<double>(); 
  Array<double>* runchrome_intens = new Array<double>(); ;



  //  runchrome_intWtdStdev_ = new dbl_hash();
  //  runchrome_intWtdMean_hash_ = new dbl_hash();

  string pepQ1Q3;
  string pep;
  string run;
  string run_pepQ1Q3;
  
  double w = 0;
  double num = 0;
  double denom = 0;
  double mean = -1;
  double stdev = -1;
  double run_maxIntens = -1;

  double run_peakRT = -1;
  double run_sumIntens = 0;
  double all_sumIntens = 0;

  double all_maxIntens = 0;

  double all_num = 0;

  double all_denom = 0;

  double all_sumsq = 0;

  double sumsq = 0;

  for (k=0; k<runs->size(); k++) {
    string rn = *(*runs)[k];
    for (dblarr_hash::iterator it = (*(*byrun_peprts_hash_)[rn]).begin(); it != (*(*byrun_peprts_hash_)[rn]).end(); it++) {
      
      (*(*(*byrun_peprts_hash_)[rn])[it->first]).clear();
      (*(*(*byrun_pepintens_hash_)[rn])[it->first]).clear();


   }
  }

  for (i=0; i<runs->size(); i++) {
    got_run_match = false;
    //TODO: this means the directory names cannot have a '.'
    string rn = (*run_files_)[*(*runs)[i]];
    pos =  rn.find_last_of("/\\");
		//pos2 =  rn->find_last_of(".");
    rn =  rn.substr(pos+1);//, pos2-pos1-1);

    dblarr_hash_hash::iterator itr1 = byrun_peps_q1s->find(rn);


    if (itr1 == byrun_peps_q1s->end()) {
      continue;
    }

  
    
    try {
      msd = new MSDataFile((*run_files_)[*(*runs)[i]]);
    }
    catch (...) {
      cerr << "WARNING: Unable to open file: " << (*run_files_)[*(*runs)[i]] << endl;
      continue;
    }
    
   
    if (!msd->run.spectrumListPtr->get())
      throw runtime_error("[trackPeptidesXICs] Null spectrumListPtr.");

    
    sl = msd->run.spectrumListPtr;
    
    cerr << "Parsing Data File: " << (*run_files_)[*(*runs)[i]]   << endl;

    int step = (int)sl->size() / 100;
    int tot = 0;


    for (j=1; j<(int)sl->size(); j++) {
      progress(j-1, step, tot);
      s = sl->spectrum(j, true);

      SpectrumInfo info;
      
      info.update(*s, true);
      
      if (info.msLevel != 1) {
	continue;
      }


      double rt = info.retentionTime;

   



      double intens = 0;

      std::vector< MZIntensityPair >  data =  info.data;	
      
      std::sort(data.begin(), data.end(), compare_mz);

      run = itr1->first; 
      size_t id = 0;
      size_t last_id = 0;
      size_t left = 0;
      size_t right = info.dataSize ? info.dataSize-1 : 0;
      double look_mz = -1;

      if (right <= left) {
	continue;
      }
      //      for (size_t id = 0; id < info.dataSize; id++) {A

      for (dblarr_hash::iterator itr2 = itr1->second->begin() ; itr2 != itr1->second->end(); itr2++) {

	for (int q1_x = 0;  q1_x < itr2->second->size(); q1_x++) {

	  intens = 0;
	  pep = itr2->first;

	  //  double prevRT = (*(*(*byrun_peprts_hash_)[run])[pep])
	  //(*(*(*byrun_peprts_hash_)[run])[pep]).clear();
	  //(*(*(*byrun_pepintens_hash_)[run])[pep]).clear();

	  double prevRT = (*(*byrun_peprtinfo_hash_)[run])[pep]->origMed_;
	  
	  if (fabs(rt - prevRT ) > 180) {
	    continue;
	  }




	  look_mz = (*itr2->second)[q1_x];
	  id = 0;
	  left = 0;
	  right  = info.dataSize ? info.dataSize-1 : 0;

	  id = (right+left)/2;
	  
	  double got_mz = data[id].mz;
	  double ppm; 
	  
	  //Binary Search
	  while (1) {

	     got_mz = data[id].mz;
	     ppm = 1e6*(got_mz-look_mz)/look_mz;

	    if (right <= left) {
	      if (fabs(ppm) < maxPPM && data[id].intensity > intens) {
		intens = data[id].intensity;
	      }
	      break;
	    }
	    else if (fabs(ppm) < maxPPM && data[id].intensity > intens) {
	      intens = data[id].intensity;
	      break;
	    }
	    else if (look_mz < data[id].mz) {
	      right = id;
	    }
	    else {
	      left = id;
	    }

	    last_id = id;
	    id = (right+left)/2;
	    
	    if (last_id == id) {
	      id = right;
	      right = left;
	    }
	    
	  }


      
	
	  if (intens > 0) {
	    (*(*byrun_peprts_hash_)[run])[pep]->insertAtEnd(acn_gradient_->getAcn(rt));
	    (*(*byrun_pepintens_hash_)[run])[pep]->insertAtEnd(intens);
	  }
	  //  if (pep == "VFQQYAGTEVK" && intens > 0) {
	  //  cerr << "DDS: DEBUG here!!!" << endl;
	  //}

	  if ((*(*byrun_pep_maxIntens_)[run])[pep] < intens ) {
	    if (pep == "VFQQYAGTEVK") {
	      cerr << "DDS: DEBUG here!!!" << endl;
	    }
	    (*(*byrun_pep_maxIntens_)[run])[pep] = intens;
	  }

	}
      }
    
// 	for (dblarr_hash::iterator itr2 = itr1->second->begin() ; itr2 != itr1->second->end(); itr2++) {
// 	  intens = 0;
// 	  pep = itr2->first;
// 	  look_mz = (*itr2->second)[q1_x];
// 	  //if ((*(*byrun_peprts_hash_)[run]).find(pep) == (*(*byrun_peprts_hash_)[run]).end()) {
// 	  //  continue;
// 	  //}
// 	  for (int q1_x = 0;  q1_x < itr2->second->size(); q1_x++) {
// 	    if (fabs((*itr2->second)[q1_x]-data[id].mz) < 0.05 && data[id].intensity > intens) {
// 	      intens = data[id].intensity;
// 	    }
// 	  }
// 	  if (intens > 0) {
// 	    (*(*byrun_peprts_hash_)[run])[pep]->insertAtEnd(rt);
// 	    (*(*byrun_pepintens_hash_)[run])[pep]->insertAtEnd(intens);
// 	  }

// 	  if ((*(*byrun_pep_maxIntens_)[run])[pep] < intens ) {
// 	    (*(*byrun_pep_maxIntens_)[run])[pep] = intens;
// 	  }
// 	}
//       }
    }

    delete msd;
    cerr << " done" << endl;
  }
      
 




}



void  RTCatalog::trackPeptidesChromatograms(ostream& out, Array<string*>* runs, dblarr_hash* peps_q1q3) {
  using namespace std;

  int i=0;
  int j=0;
  int k=0;
  MSDataFile* msd;
  ChromatogramListPtr sl;
  ChromatogramPtr s;
  string q1, q3;
  size_t pos;
  string pep;
  bool got_match = false;
  bool got_run_match = false;
  dblarr_hash::iterator itr;

out << "<html>" << endl;
out << "<head> " << endl;
out << "<title>  </title> " << endl;
out << "<link type='text/css' rel='stylesheet' href='http://regis-web/tpp-dshteynb/html/ex.css?3.2'/> " << endl;
out << "<script type='text/javascript' src='http://regis-web/tpp-dshteynb/html/protovis-d3.2.js'></script> " << endl;
out << "</head> " << endl;
out << "<body> " << endl;
out << "<center><div id='id'> " << endl;
out << "<script type='text/javascript'> " << endl;


  out << "var runs = [" << endl;
  for (k = 0; k < runs->size() ; k++) {
    
    //    for (dblarr_hash::iterator it = (*(*byrun_peprts_hash_)[rn]).begin(); it != (*(*byrun_peprts_hash_)[rn]).end(); it++) {
      
    //  (*(*(*byrun_peprts_hash_)[rn])[it->first]).clear();
    //  (*(*(*byrun_pepintens_hash_)[rn])[it->first]).clear();


    //}    

    out << "{run_index: " << k ;
    out <<", run: '" << (*runs)[k]->substr((*runs)[k]->find_last_of("/\\")+1) << "'}";
    if (k < runs->size()-1) {
      out << ", ";
    }
    out << endl;
  }
  out << "];" << endl;

  out << "var chroms = [" << endl;
  
  int last_run = -1;

  dblarr_hash* bychrome_rtsecs = new dblarr_hash(); 
  dblarr_hash* bychrome_intens = new dblarr_hash(); ;
  dbl_hash* bychrome_maxIntens = new dbl_hash();
  dbl_hash* bychrome_peakRT = new dbl_hash();


  Array<double>* runchrome_rtsecs = new Array<double>(); 
  Array<double>* runchrome_intens = new Array<double>(); ;



  //  runchrome_intWtdStdev_ = new dbl_hash();
  //  runchrome_intWtdMean_hash_ = new dbl_hash();

  string pepQ1Q3;
  string run_pepQ1Q3;
  
  double w = 0;
  double num = 0;
  double denom = 0;
  double mean = -1;
  double stdev = -1;
  double run_maxIntens = -1;

  double run_peakRT = -1;
  double run_sumIntens = 0;
  double all_sumIntens = 0;

  double all_maxIntens = 0;

  double all_num = 0;

  double all_denom = 0;

  double all_sumsq = 0;

  double sumsq = 0;
  for (i=0; i<runs->size(); i++) {
    got_run_match = false;
    string run_name = (*run_files_)[*(*runs)[i]];
    pos = run_name.find_first_of('.');

    if (pos!=string::npos) {
      run_name = run_name.substr(0, pos);
    }
    
    try {
      msd = new MSDataFile(run_name); //+".mzML");
    }
    catch (...) {
      cerr << "WARNING: Unable to open file: " << (*run_files_)[*(*runs)[i]] << endl;
      continue;
    }
    

    if (!msd->run.chromatogramListPtr->get())
      throw runtime_error("[trackPeptidesChromatograms] Null chromatogramListPtr.");

    sl = msd->run.chromatogramListPtr;

    cout << "Parsing Data File: " << (*run_files_)[*(*runs)[i]] << endl;
   
    for (j=1; j<(int)sl->size(); j++) {

      pepQ1Q3 = "";
      q1 = "";
      q3 = "";

      s = sl->chromatogram(j, true);

      q1 = s->id;
      pos = q1.find("Q1=");
      if (pos!=string::npos) {
	q1 = q1.substr(pos);
	pos = q1.find_first_of(" \t");
	if (pos!=string::npos) {
	  q1 = q1.substr(3, pos-3);
	}
      }
      q3 = s->id;
      pos = q3.find("Q3=");
      if (pos!=string::npos) {
	q3 = q3.substr(pos);
	pos = q3.find_first_of(" \t");
	if (pos!=string::npos) {
	  q3 = q3.substr(3, pos-3);
	}
      }

      itr = peps_q1q3->begin();

      while(itr != peps_q1q3->end()) {
	for (k=0; k<itr->second->size(); k+= 2) {
	  if (fabs((*itr->second)[k] - atof(q1.c_str())) < 0.05 && fabs((*itr->second)[k+1] - atof(q3.c_str())) < 0.05 ) {
	    got_match = true;
	    pos = q1.find(".");
	    if (pos!=string::npos) {
	      q1 = q1.substr(0, pos);
	    }
	    pos = q3.find(".");
	    if (pos!=string::npos) {
	      q3 = q3.substr(0, pos);
	    }
	    pepQ1Q3 += itr->first+q1+q3;
	    pep = itr->first;

	    break;
	  }
	}
	if (got_match) break;
	itr++;
      }
      //      run_pepQ1Q3 = (*run_names_)[i] + pepQ1Q3;
      got_run_match = got_match || got_run_match ;
      
    
      if (got_match) {

	if (i != last_run) {
	  if (last_run < 0) {
	    out << "[" << endl;
	  }
	  else {
	    out << "], [" << endl;
	  }
	}
	else {
	  out << ", " << endl;
	}
	
	out << "{name: '" << s->id << "', run_index: " << i ;
	out << ", chrome: [" ;
	run_maxIntens = -1;
	run_peakRT = -1;
	run_sumIntens = 0;
	runchrome_rtsecs->clear();
	runchrome_intens->clear();
	vector<TimeIntensityPair> pairs;
	s->getTimeIntensityPairs(pairs);
	for (size_t ix=0; ix < pairs.size(); ix++) {
	  if (pairs[ix].intensity > 0 && 
		pairs[ix].time*60 <= maxRT_ && 
		pairs[ix].time*60 >= minRT_) {
	    out << "{time: " << pairs[ix].time;
	    out << ", intensity: " << pairs[ix].intensity << "}";

	    if ((*bychrome_maxIntens).find(pepQ1Q3.c_str()) != (*bychrome_maxIntens).end()) {
	      if ((*bychrome_maxIntens)[pepQ1Q3.c_str()] <  pairs[ix].intensity) {
		(*bychrome_maxIntens)[pepQ1Q3.c_str()] =  pairs[ix].intensity;
		(*bychrome_peakRT)[pepQ1Q3.c_str()] =  pairs[ix].time;

	      }
	    }
	    else {
	      (*bychrome_maxIntens)[pepQ1Q3.c_str()] =  pairs[ix].intensity;
	      (*bychrome_peakRT)[pepQ1Q3.c_str()] =  pairs[ix].time;
	      (*bychrome_rtsecs)[pepQ1Q3.c_str()] = new Array<double>();
	      (*bychrome_intens)[pepQ1Q3.c_str()] = new Array<double>();
	    }

	    (*bychrome_rtsecs)[pepQ1Q3.c_str()]->insertAtEnd(pairs[ix].time);
	    (*bychrome_intens)[pepQ1Q3.c_str()]->insertAtEnd(pairs[ix].intensity);

	    run_sumIntens = pairs[ix].intensity;



	    if (run_maxIntens <  pairs[ix].intensity) {
	     run_maxIntens =  pairs[ix].intensity;
	     run_peakRT = pairs[ix].time;
	    }

	    if (all_maxIntens < run_maxIntens) {
	      all_maxIntens = run_maxIntens;
	    }
	 
	    runchrome_rtsecs->insertAtEnd(pairs[ix].time);
	    runchrome_intens->insertAtEnd(pairs[ix].intensity);

	    
 	    if ((*(*byrun_peprts_hash_)[*(*runs)[i]]).find(pep) != (*(*byrun_peprts_hash_)[*(*runs)[i]]).end() ) {
	      double RTsec = pairs[ix].time*60;
 	      (*(*byrun_peprts_hash_)[*(*runs)[i]])[pep]->insertAtEnd(acn_gradient_->getAcn(RTsec));
 	      (*(*byrun_pepintens_hash_)[*(*runs)[i]])[pep]->insertAtEnd(pairs[ix].intensity);
 	    }

	    

	    if (ix <pairs.size()-1) {
	      out << ", " << endl;
	    }	  
	  }
	}
	out << "]" << endl;

	if (i != last_run) {
	  last_run = i;
	}

	if ((*bychrome_rtsecs)[pepQ1Q3.c_str()] != NULL ) {
	  (*bychrome_rtsecs)[pepQ1Q3.c_str()]->insertAtEnd(run_peakRT);
	  (*bychrome_intens)[pepQ1Q3.c_str()]->insertAtEnd(run_maxIntens);
	}
      
	sumsq = 0;
	num = 0;
	denom = 0;
	mean = -1;
	stdev = -1;
	
	//Compute wtdMean for each chromatoram in each run
	num = 0; denom = 0; 
	all_sumIntens += run_sumIntens;
	for (int idx=0; idx < runchrome_rtsecs->size(); idx++) {
	  w =  (*runchrome_intens)[idx] / run_maxIntens;
	  num += (*runchrome_rtsecs)[idx] * w; 
	  sumsq += (*runchrome_rtsecs)[idx] * (*runchrome_rtsecs)[idx] * w ;
	  denom += w;
	}
	if (denom > 0) {
	  mean = num / denom;
	  stdev = sqrt(sumsq/denom - mean*mean);
	}
	out << ", mean: " << mean << ", stdev: " << stdev  << ", peakRT: " << run_peakRT << ", peakIntensity: " << run_maxIntens << endl;

	//	all_num += run_maxIntens * run_peakRT;
	//all_sumsq += run_maxIntens * run_peakRT * run_peakRT;
	//all_denom += run_maxIntens;

	
	out << "}";
      }
      got_match = false;
    }


    delete msd;
  }
  out << "]];" << endl;

  out << "var wtdMeans = [" << endl;
  sumsq = 0;
  num = 0;
  denom = 0;
  mean = -1;
  stdev = -1;

  //compute intensity weighted (normalized against maxInt) RT mean
  for (itr = bychrome_rtsecs->begin(); itr != bychrome_rtsecs->end(); ) {
    num = 0; denom = 0;
    if (itr->second != NULL) {
      for (int idx=0; idx < itr->second->size(); idx++) {
	w =  (*(*bychrome_intens)[itr->first])[idx] / (*bychrome_maxIntens)[itr->first];
	num += (*(*bychrome_rtsecs)[itr->first])[idx] * w; 
	sumsq += (*(*bychrome_rtsecs)[itr->first])[idx] * (*(*bychrome_rtsecs)[itr->first])[idx] * w ;
	denom += w;
      }
      if (denom > 0) {
	mean = num / denom;
	stdev = sqrt(sumsq/denom - mean*mean);
      }
      out << "{ name: '" << itr->first << "', mean: " << mean << ", stdev: " << stdev << "}";
      itr++;
      if (itr != bychrome_rtsecs->end()) {
	out << ", " << endl;
      }
    }
    else {
      itr++;
    }

  }
  out << "];" << endl;
  

  

out << "  var size = 600;" << endl;
out << "  var pad = 40;" << endl;
out << "  var w = 1;" << endl;
out << "  var maxTIME = 5000;" << endl;
out << "  var maxINT = 15000;" << endl;
out << "  var vis = new pv.Panel()" << endl;
out << "    .width((size+pad)*w)" << endl;
out << "    .height((size+pad)*(1+Math.floor(runs.length/w)))" << endl;
out << "    .left(10)" << endl;
out << "    .top(5)" << endl;
out << "    .events('all');" << endl;

 
out << " var plot = vis.add(pv.Panel)" << endl;
out << "    .data(chroms)" << endl;
out << "    .top(function() { return Math.floor(this.index / w) * (size + pad) + pad / 2; } )" << endl;
out << "    .left(function() { return this.index % w * (size + pad) + pad / 2; } )" << endl;
out << "    .height(size)" << endl;
out << "    .width(size);" << endl;
out << " var run_plot = plot.add(pv.Panel)" << endl;
out << "   .data(function() {" << endl;
out << "	   maxINT = pv.max(pv.values(chroms[this.parent.index]).map(function(e) { " << endl;
out << "								      return e.chrome; " << endl;
out << "								    } ) , " << endl;
out << "			   function(f) { " << endl;
out << "			     return pv.max(f , function(g) { return g.intensity; } ); " << endl;
out << "			   } );" << endl;
out << "	   return chroms[this.parent.index];" << endl;
out << "	 }" << endl;
out << "	 );" << endl;
out << " var xscale = pv.range(1000, maxTIME, maxTIME/10);" << endl;
out << " var xax = plot.add(pv.Rule)" << endl;
out << "    .data(xscale)" << endl;
out << "    .left(function(d) { return d*size/maxTIME ; } )" << endl;
out << "    .strokeStyle(function(d) { return d ? '#DDD' : 'black' ; } )" << endl;
out << "    .anchor('bottom')" << endl;
out << "    .add(pv.Label)" << endl;
out << "    .text(function(d) { return d; })" << endl;
out << "    .font('7px sans-serif');" << endl;
out << "  var yax = run_plot.add(pv.Rule)" << endl;
out << "    .data(function(c) { " << endl;
out << "                         return pv.range(0, maxINT, maxINT/7) ; " << endl;
out << "		       })" << endl;
out << "    .bottom(function(d) { return d*size/maxINT ;  } )" << endl;
out << "    .strokeStyle(function(d) { return d ? '#DDD' : 'black' ; } )" << endl;
out << "    .anchor('left')" << endl;
out << "    .add(pv.Label)" << endl;
out << "    .text(function(d) { " << endl;
out << "                         return Math.round(d); " << endl;
out << "                      }" << endl;
out << "          )" << endl;
out << "    .font('7px sans-serif');" << endl;
out << "  var name = '';" << endl;
out << "  var line1 = run_plot.add(pv.Line)" << endl;
out << "   .data(function(c) {" << endl;
out << "                        name = c.name;" << endl;
out << "			return c.chrome;" << endl;
out << "	               }" << endl;
out << "          )" << endl;
out << "    .lineWidth(1)" << endl;
out << "    .title(function(d) {" << endl;
out << "                          return name; " << endl;
out << "                       }" << endl;
out << "           )" << endl;
out << "    .left(function(d) { " << endl;
out << "                          return d.time*60*size/maxTIME; " << endl;
out << "                      } )" << endl;
out << "    .bottom(function(d) { " << endl;
out << "                          return d.intensity*size/maxINT ; " << endl;
out << "                      } );" << endl;
out << "  var run_labl = plot.anchor('top')" << endl;
out << "    .add(pv.Label)" << endl;
out << "    .top(-15)" << endl;
out << "    .text(function(d) { return  runs[d[0].run_index].run; } )" << endl;
out << "    .font('12px sans-serif');" << endl;
out << "  var ylabl = plot.anchor('left')" << endl;
out << "    .add(pv.Label)" << endl;
out << "    .left(-33)" << endl;
out << "    .text('Intensity')" << endl;
out << "    .font('7px sans-serif');" << endl;
out << "  var ylabl = plot.anchor('bottom')" << endl;
out << "    .add(pv.Label)" << endl;
out << "    .top(size+pad/2)" << endl;
out << "    .text('Retention Time (seconds)')" << endl;
out << "vis.render();" << endl;
out << "</script>" << endl;
out << "</body>" << endl;
out << "</html>" << endl;
  
}



// void  RTCatalog::trackPeptideXICs(ostream& out, Array<string*>* runs, dblarr_hash* peps_q1q3) {
//   using namespace std;
//   using namespace pwiz::msdata;

//   int i=0;
//   int j=0;
//   int k=0;
//   pwiz::msdata::MSDataFile* msd;
//   ChromatogramListPtr sl;
//   ChromatogramPtr s;
//   string q1, q3;
//   size_t pos;
//   string pep;
//   bool got_match = false;
//   bool got_run_match = false;
//   dblarr_hash::iterator itr;

// out << "<html>" << endl;
// out << "<head> " << endl;
// out << "<title>  </title> " << endl;
// out << "<link type='text/css' rel='stylesheet' href='http://regis-web/tpp-dshteynb/html/ex.css?3.2'/> " << endl;
// out << "<script type='text/javascript' src='http://regis-web/tpp-dshteynb/html/protovis-d3.2.js'></script> " << endl;
// out << "</head> " << endl;
// out << "<body> " << endl;
// out << "<center><div id='id'> " << endl;
// out << "<script type='text/javascript'> " << endl;


//   out << "var runs = [" << endl;
//   for (k = 0; k < runs->size() ; k++) {
    
//     //    for (dblarr_hash::iterator it = (*(*byrun_peprts_hash_)[rn]).begin(); it != (*(*byrun_peprts_hash_)[rn]).end(); it++) {
      
//     //  (*(*(*byrun_peprts_hash_)[rn])[it->first]).clear();
//     //  (*(*(*byrun_pepintens_hash_)[rn])[it->first]).clear();


//     //}    

//     out << "{run_index: " << k ;
//     out <<", run: '" << (*runs)[k]->substr((*runs)[k]->find_last_of("/\\")+1) << "'}";
//     if (k < runs->size()-1) {
//       out << ", ";
//     }
//     out << endl;
//   }
//   out << "];" << endl;

//   out << "var chroms = [" << endl;
  
//   int last_run = -1;

//   dblarr_hash* bychrome_rtsecs = new dblarr_hash(); 
//   dblarr_hash* bychrome_intens = new dblarr_hash(); ;
//   dbl_hash* bychrome_maxIntens = new dbl_hash();
//   dbl_hash* bychrome_peakRT = new dbl_hash();


//   Array<double>* runchrome_rtsecs = new Array<double>(); 
//   Array<double>* runchrome_intens = new Array<double>(); ;



//   //  runchrome_intWtdStdev_ = new dbl_hash();
//   //  runchrome_intWtdMean_hash_ = new dbl_hash();

//   string pepQ1Q3;
//   string run_pepQ1Q3;
  
//   double w = 0;
//   double num = 0;
//   double denom = 0;
//   double mean = -1;
//   double stdev = -1;
//   double run_maxIntens = -1;

//   double run_peakRT = -1;
//   double run_sumIntens = 0;
//   double all_sumIntens = 0;

//   double all_maxIntens = 0;

//   double all_num = 0;

//   double all_denom = 0;

//   double all_sumsq = 0;

//   double sumsq = 0;
//   for (i=0; i<runs->size(); i++) {
//     got_run_match = false;
//     string run_name = (*run_files_)[*(*runs)[i]];
//     pos = run_name.find_first_of('.');

//     if (pos!=string::npos) {
//       run_name = run_name.substr(0, pos);
//     }
    
//     try {
//       msd = new pwiz::msdata::MSDataFile(run_name+".mzML");
//     }
//     catch (...) {
//       cerr << "WARNING: Unable to open file: " << (*run_files_)[*(*runs)[i]] << endl;
//       continue;
//     }
    

//     if (!msd->run.chromatogramListPtr.get())
//       throw runtime_error("[trackPeptidesChromatograms] Null chromatogramListPtr.");

//     sl = msd->run.chromatogramListPtr;

//     cout << "Parsing mzML File: " << (*run_files_)[*(*runs)[i]] << endl;
   
//     for (j=1; j<(int)sl->size(); j++) {

//       pepQ1Q3 = "";
//       q1 = "";
//       q3 = "";

//       s = sl->chromatogram(j, true);

//       q1 = s->id;
//       pos = q1.find("Q1=");
//       if (pos!=string::npos) {
// 	q1 = q1.substr(pos);
// 	pos = q1.find_first_of(" \t");
// 	if (pos!=string::npos) {
// 	  q1 = q1.substr(3, pos-3);
// 	}
//       }
//       q3 = s->id;
//       pos = q3.find("Q3=");
//       if (pos!=string::npos) {
// 	q3 = q3.substr(pos);
// 	pos = q3.find_first_of(" \t");
// 	if (pos!=string::npos) {
// 	  q3 = q3.substr(3, pos-3);
// 	}
//       }

//       itr = peps_q1q3->begin();

//       while(itr != peps_q1q3->end()) {
// 	for (k=0; k<itr->second->size(); k+= 2) {
// 	  if (fabs((*itr->second)[k] - atof(q1.c_str())) < 0.05 && fabs((*itr->second)[k+1] - atof(q3.c_str())) < 0.05 ) {
// 	    got_match = true;
// 	    pos = q1.find(".");
// 	    if (pos!=string::npos) {
// 	      q1 = q1.substr(0, pos);
// 	    }
// 	    pos = q3.find(".");
// 	    if (pos!=string::npos) {
// 	      q3 = q3.substr(0, pos);
// 	    }
// 	    pepQ1Q3 += itr->first+q1+q3;
// 	    pep = itr->first;

// 	    break;
// 	  }
// 	}
// 	if (got_match) break;
// 	itr++;
//       }
//       //      run_pepQ1Q3 = (*run_names_)[i] + pepQ1Q3;
//       got_run_match = got_match || got_run_match ;
      
    
//       if (got_match) {

// 	if (i != last_run) {
// 	  if (last_run < 0) {
// 	    out << "[" << endl;
// 	  }
// 	  else {
// 	    out << "], [" << endl;
// 	  }
// 	}
// 	else {
// 	  out << ", " << endl;
// 	}
	
// 	out << "{name: '" << s->id << "', run_index: " << i ;
// 	out << ", chrome: [" ;
// 	run_maxIntens = -1;
// 	run_peakRT = -1;
// 	run_sumIntens = 0;
// 	runchrome_rtsecs->clear();
// 	runchrome_intens->clear();
// 	vector<pwiz::msdata::TimeIntensityPair> pairs;
// 	s->getTimeIntensityPairs(pairs);
// 	for (size_t ix=0; ix < pairs.size(); ix++) {
// 	  if (pairs[ix].intensity > 0 && 
// 		pairs[ix].time*60 <= maxRT_ && 
// 		pairs[ix].time*60 >= minRT_) {
// 	    out << "{time: " << pairs[ix].time;
// 	    out << ", intensity: " << pairs[ix].intensity << "}";

// 	    if ((*bychrome_maxIntens).find(pepQ1Q3.c_str()) != (*bychrome_maxIntens).end()) {
// 	      if ((*bychrome_maxIntens)[pepQ1Q3.c_str()] <  pairs[ix].intensity) {
// 		(*bychrome_maxIntens)[pepQ1Q3.c_str()] =  pairs[ix].intensity;
// 		(*bychrome_peakRT)[pepQ1Q3.c_str()] =  pairs[ix].time;

// 	      }
// 	    }
// 	    else {
// 	      (*bychrome_maxIntens)[pepQ1Q3.c_str()] =  pairs[ix].intensity;
// 	      (*bychrome_peakRT)[pepQ1Q3.c_str()] =  pairs[ix].time;
// 	      (*bychrome_rtsecs)[pepQ1Q3.c_str()] = new Array<double>();
// 	      (*bychrome_intens)[pepQ1Q3.c_str()] = new Array<double>();
// 	    }

// 	    (*bychrome_rtsecs)[pepQ1Q3.c_str()]->insertAtEnd(pairs[ix].time);
// 	    (*bychrome_intens)[pepQ1Q3.c_str()]->insertAtEnd(pairs[ix].intensity);

// 	    run_sumIntens = pairs[ix].intensity;



// 	    if (run_maxIntens <  pairs[ix].intensity) {
// 	     run_maxIntens =  pairs[ix].intensity;
// 	     run_peakRT = pairs[ix].time;
// 	    }

// 	    if (all_maxIntens < run_maxIntens) {
// 	      all_maxIntens = run_maxIntens;
// 	    }
	 
// 	    runchrome_rtsecs->insertAtEnd(pairs[ix].time);
// 	    runchrome_intens->insertAtEnd(pairs[ix].intensity);

	    
//  	    if ((*(*byrun_peprts_hash_)[*(*runs)[i]]).find(pep) != (*(*byrun_peprts_hash_)[*(*runs)[i]]).end() ) {
// 	      double RTsec = pairs[ix].time*60;
//  	      (*(*byrun_peprts_hash_)[*(*runs)[i]])[pep]->insertAtEnd(RTsec);
//  	      (*(*byrun_pepintens_hash_)[*(*runs)[i]])[pep]->insertAtEnd(pairs[ix].intensity);
//  	    }

	    

// 	    if (ix <pairs.size()-1) {
// 	      out << ", " << endl;
// 	    }	  
// 	  }
// 	}
// 	out << "]" << endl;

// 	if (i != last_run) {
// 	  last_run = i;
// 	}

// 	if ((*bychrome_rtsecs)[pepQ1Q3.c_str()] != NULL ) {
// 	  (*bychrome_rtsecs)[pepQ1Q3.c_str()]->insertAtEnd(run_peakRT);
// 	  (*bychrome_intens)[pepQ1Q3.c_str()]->insertAtEnd(run_maxIntens);
// 	}
      
// 	sumsq = 0;
// 	num = 0;
// 	denom = 0;
// 	mean = -1;
// 	stdev = -1;
	
// 	//Compute wtdMean for each chromatoram in each run
// 	num = 0; denom = 0; 
// 	all_sumIntens += run_sumIntens;
// 	for (int idx=0; idx < runchrome_rtsecs->size(); idx++) {
// 	  w =  (*runchrome_intens)[idx] / run_maxIntens;
// 	  num += (*runchrome_rtsecs)[idx] * w; 
// 	  sumsq += (*runchrome_rtsecs)[idx] * (*runchrome_rtsecs)[idx] * w ;
// 	  denom += w;
// 	}
// 	if (denom > 0) {
// 	  mean = num / denom;
// 	  stdev = sqrt(sumsq/denom - mean*mean);
// 	}
// 	out << ", mean: " << mean << ", stdev: " << stdev  << ", peakRT: " << run_peakRT << ", peakIntensity: " << run_maxIntens << endl;

// 	//	all_num += run_maxIntens * run_peakRT;
// 	//all_sumsq += run_maxIntens * run_peakRT * run_peakRT;
// 	//all_denom += run_maxIntens;

	
// 	out << "}";
//       }
//       got_match = false;
//     }


//     delete msd;
//   }
//   out << "]];" << endl;

//   out << "var wtdMeans = [" << endl;
//   sumsq = 0;
//   num = 0;
//   denom = 0;
//   mean = -1;
//   stdev = -1;

//   //compute intensity weighted (normalized against maxInt) RT mean
//   for (itr = bychrome_rtsecs->begin(); itr != bychrome_rtsecs->end(); ) {
//     num = 0; denom = 0;
//     if (itr->second != NULL) {
//       for (int idx=0; idx < itr->second->size(); idx++) {
// 	w =  (*(*bychrome_intens)[itr->first])[idx] / (*bychrome_maxIntens)[itr->first];
// 	num += (*(*bychrome_rtsecs)[itr->first])[idx] * w; 
// 	sumsq += (*(*bychrome_rtsecs)[itr->first])[idx] * (*(*bychrome_rtsecs)[itr->first])[idx] * w ;
// 	denom += w;
//       }
//       if (denom > 0) {
// 	mean = num / denom;
// 	stdev = sqrt(sumsq/denom - mean*mean);
//       }
//       out << "{ name: '" << itr->first << "', mean: " << mean << ", stdev: " << stdev << "}";
//       itr++;
//       if (itr != bychrome_rtsecs->end()) {
// 	out << ", " << endl;
//       }
//     }
//     else {
//       itr++;
//     }

//   }
//   out << "];" << endl;
  

  

// out << "  var size = 600;" << endl;
// out << "  var pad = 40;" << endl;
// out << "  var w = 1;" << endl;
// out << "  var maxTIME = 5000;" << endl;
// out << "  var maxINT = 15000;" << endl;
// out << "  var vis = new pv.Panel()" << endl;
// out << "    .width((size+pad)*w)" << endl;
// out << "    .height((size+pad)*(1+Math.floor(runs.length/w)))" << endl;
// out << "    .left(10)" << endl;
// out << "    .top(5)" << endl;
// out << "    .events('all');" << endl;

 
// out << " var plot = vis.add(pv.Panel)" << endl;
// out << "    .data(chroms)" << endl;
// out << "    .top(function() { return Math.floor(this.index / w) * (size + pad) + pad / 2; } )" << endl;
// out << "    .left(function() { return this.index % w * (size + pad) + pad / 2; } )" << endl;
// out << "    .height(size)" << endl;
// out << "    .width(size);" << endl;
// out << " var run_plot = plot.add(pv.Panel)" << endl;
// out << "   .data(function() {" << endl;
// out << "	   maxINT = pv.max(pv.values(chroms[this.parent.index]).map(function(e) { " << endl;
// out << "								      return e.chrome; " << endl;
// out << "								    } ) , " << endl;
// out << "			   function(f) { " << endl;
// out << "			     return pv.max(f , function(g) { return g.intensity; } ); " << endl;
// out << "			   } );" << endl;
// out << "	   return chroms[this.parent.index];" << endl;
// out << "	 }" << endl;
// out << "	 );" << endl;
// out << " var xscale = pv.range(1000, maxTIME, maxTIME/10);" << endl;
// out << " var xax = plot.add(pv.Rule)" << endl;
// out << "    .data(xscale)" << endl;
// out << "    .left(function(d) { return d*size/maxTIME ; } )" << endl;
// out << "    .strokeStyle(function(d) { return d ? '#DDD' : 'black' ; } )" << endl;
// out << "    .anchor('bottom')" << endl;
// out << "    .add(pv.Label)" << endl;
// out << "    .text(function(d) { return d; })" << endl;
// out << "    .font('7px sans-serif');" << endl;
// out << "  var yax = run_plot.add(pv.Rule)" << endl;
// out << "    .data(function(c) { " << endl;
// out << "                         return pv.range(0, maxINT, maxINT/7) ; " << endl;
// out << "		       })" << endl;
// out << "    .bottom(function(d) { return d*size/maxINT ;  } )" << endl;
// out << "    .strokeStyle(function(d) { return d ? '#DDD' : 'black' ; } )" << endl;
// out << "    .anchor('left')" << endl;
// out << "    .add(pv.Label)" << endl;
// out << "    .text(function(d) { " << endl;
// out << "                         return Math.round(d); " << endl;
// out << "                      }" << endl;
// out << "          )" << endl;
// out << "    .font('7px sans-serif');" << endl;
// out << "  var name = '';" << endl;
// out << "  var line1 = run_plot.add(pv.Line)" << endl;
// out << "   .data(function(c) {" << endl;
// out << "                        name = c.name;" << endl;
// out << "			return c.chrome;" << endl;
// out << "	               }" << endl;
// out << "          )" << endl;
// out << "    .lineWidth(1)" << endl;
// out << "    .title(function(d) {" << endl;
// out << "                          return name; " << endl;
// out << "                       }" << endl;
// out << "           )" << endl;
// out << "    .left(function(d) { " << endl;
// out << "                          return d.time*60*size/maxTIME; " << endl;
// out << "                      } )" << endl;
// out << "    .bottom(function(d) { " << endl;
// out << "                          return d.intensity*size/maxINT ; " << endl;
// out << "                      } );" << endl;
// out << "  var run_labl = plot.anchor('top')" << endl;
// out << "    .add(pv.Label)" << endl;
// out << "    .top(-15)" << endl;
// out << "    .text(function(d) { return  runs[d[0].run_index].run; } )" << endl;
// out << "    .font('12px sans-serif');" << endl;
// out << "  var ylabl = plot.anchor('left')" << endl;
// out << "    .add(pv.Label)" << endl;
// out << "    .left(-33)" << endl;
// out << "    .text('Intensity')" << endl;
// out << "    .font('7px sans-serif');" << endl;
// out << "  var ylabl = plot.anchor('bottom')" << endl;
// out << "    .add(pv.Label)" << endl;
// out << "    .top(size+pad/2)" << endl;
// out << "    .text('Retention Time (seconds)')" << endl;
// out << "vis.render();" << endl;
// out << "</script>" << endl;
// out << "</body>" << endl;
// out << "</html>" << endl;
  
// }

void  RTCatalog::trackPeptidesReport(ostream& out, Array<string*>* peps) {
  out << "Peptide" << "\t";
  
  for (int k = 0; k < num_runs_ ; k++) {
    out << (*run_names_)[k]->c_str() << "_RTMedian";

    if (k < num_runs_-1) {
      out << "\t";
    }
   
  }
  out << endl;
  rtinfo_hash::iterator it;

  for (int i = 0; i < peps->size() ; i++) {
    out << (*peps)[i]->c_str() << "\t";
    for (int k = 0; k < num_runs_ ; k++) {
      it = (*byrun_peprtinfo_hash_)[*(*run_names_)[k]]->find(*(*peps)[i]);
      
      if (it == (*byrun_peprtinfo_hash_)[*(*run_names_)[k]]->end()) {
	out << "N/A";
      }
      else {
	out << it->second->med_;
      }

      if (k < num_runs_-1) {
	out << "\t";
      }
      
    }
    out << endl;
    
  }
     
}

void RTCatalog::sortRunNames() {
  std::vector<string*>* tmp = new std::vector<string*>(num_runs_);
  for (int k = 0; k < num_runs_ ; k++) {
    (*tmp)[k] = (*run_names_)[k];
  }

  std::sort(tmp->begin(), tmp->end(), compare_runs);

  for (int k = 0; k < num_runs_ ; k++) {
    (*run_names_)[k] = (*tmp)[k];
  }

  
  delete tmp;

}

void  RTCatalog::trackPeptidesReportPV(ostream& out, Array<string*>* peps, bool iRT) {
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
  out << " h = 600," << endl;
  out << " x = pv.Scale.linear(0, 400).range(0, w)," << endl;
  if (!iRT) {
    out << " y = pv.Scale.linear(0, 4000).range(0, h);" << endl;
  }
  else {
    out << " y = pv.Scale.linear(0, 120).range(0, h);" << endl;
  }

  //Data 
  
  out << "var runs = [";
  for (int k = 0; k < num_runs_ ; k++) {
    out << "\"" << (*run_names_)[k]->substr((*run_names_)[k]->find_last_of("/\\")+1) << "\"";

    if (k < num_runs_-1) {
      out << ", ";
    }
    
  }
  out << "];" << endl;

  out << "var stats = [\"RTmed\",\"RTmean\",\"RTsiqr\",\"RTstdev\",\"RTmin\", \"RTnobs\"];" << endl ;

  out << "var control_peptides = [ " << "\n";  
  rtinfo_hash::iterator it;

  for (int i = 0; i < peps->size() ; i++) {
    for (int k = 0; k < num_runs_ ; k++) {
      it = (*byrun_peprtinfo_hash_)[*(*run_names_)[k]]->find(*(*peps)[i]);
      
      if (it == (*byrun_peprtinfo_hash_)[*(*run_names_)[k]]->end()) {
	//	out << "N/A";
      }
      else {
	out << "{seq:\""<< (*peps)[i]->c_str() << "\", ";
	out << "run:\"" << (*run_names_)[k]->substr((*run_names_)[k]->find_last_of("/\\")+1) << "\", ";
	out << "RTmed: "   << it->second->med_ << ", ";
	out << "RTmean: "  << it->second->mean_ << ", ";
	out << "RTsiqr: "  << it->second->siqr_ << ", ";
	out << "RTstdev: " << it->second->stdev_ << ", ";
	out << "RTmin: "   << it->second->min_ << ", ";
	out << "RTnobs:"   << it->second->n_ << "}";
	if (k < num_runs_-1 || i <  peps->size() - 1 ) {
	  out << ",";
	}
	out  << endl;
      }
 
    }

  }
  out << "];" << endl;



  // Plot
out << "var vis = new pv.Panel()" << endl;
out << "    .width(w+200).height(h+50).left(50).bottom(50)" << endl;
out << "    .events('all');" << endl;
out << "//    .event('mousemove', pv.Behavior.point());" << endl;
out << "" << endl;
out << "var Xax   = vis.add(pv.Rule) " << endl;
out << ".data([ {  a  : 0, x : 1 },{ a : 50, x : 2 },{ a : 100, x : 3 },{ a : 150, x : 4 },{ a : 200, x : 5 },{ a : 250, x : 6 },{ a : 300, x : 7 },{ a : 350, x : 8 },{ a : 400, x : 9 } ])" << endl;
out << ".left(function(d) pv.Scale.linear(0, 400).range(0,w)(d.a))" << endl;
out << ".strokeStyle('rgba(128,128,128,.2)');" << endl;
out << "" << endl;
out << "Xax.anchor('bottom')" << endl;
out << ".add(pv.Label) " << endl;
out << ".text(function(d) (d.a));" << endl;
out << "" << endl;
out << "var Yax   = vis.add(pv.Rule) " << endl;

 if (!iRT) {
   out << ".data([ { y : 0, x : 1 },{ y : 500, x : 2 },{ y : 1000, x : 3 },{ y : 1500, x : 4 },{ y : 2000, x : 5 },{ y : 2500, x : 6 },{ y : 3000, x : 7 },{ y : 3500, x : 8 },{ y : 4000, x : 9 } ])" << endl;
   out << ".bottom(function(d) pv.Scale.linear(0, 4000).range(0,h)(d.y))" << endl;
 }
 else {
   out << ".data([ { y : 0, x : 1 },{ y : 15, x : 2 },{ y : 30, x : 3 },{ y : 45, x : 4 },{ y : 60, x : 5 },{ y : 75, x : 6 },{ y : 90, x : 7 },{ y : 105, x : 8 },{ y : 120, x : 9 } ])" << endl;
   out << ".bottom(function(d) pv.Scale.linear(0, 120).range(0,h)(d.y))" << endl;
 }
out << ".strokeStyle('rgba(128,128,128,.2)');" << endl;
out << "" << endl;
out << "Yax.anchor('left')" << endl;
out << ".add(pv.Label) " << endl;
out << ".text(function(d) (d.y));" << endl;
out << "" << endl;
out << "vis.add(pv.Label)" << endl;
out << ".bottom(h/2)" << endl;
out << ".left(-40)" << endl;
out << ".text('Median RT')" << endl;
out << ".textAngle(-Math.PI/2);" << endl;
out << "" << endl;
out << "vis.add(pv.Label)" << endl;
out << ".bottom(-25)" << endl;
out << ".left(w/2)" << endl;
out << ".text('MS Run Index');" << endl;
out << "" << endl;
out << "" << endl;
out << "" << endl;
out << "var activePep = '';" << endl;
out << "var activeText = '';" << endl;
out << "" << endl;
out << "var dot_active = -1;" << endl;
out << "" << endl;
out << "var medNobs = pv.median(control_peptides.map(function(d) { return d.RTnobs+1;}));" << endl;
out << "//console.log(medNobs);" << endl;
out << "" << endl;
out << "var maxNobs = pv.max(control_peptides.map(function(d) { return d.RTnobs+1;}));" << endl;
out << "//console.log(pv.Scale.linear(0, 1, 1).range(pv.rgb(0,255,255, 0.8), pv.rgb(255,0,0, 0.8)));" << endl;
out << "" << endl;
out << "var dot_obs = '';" << endl;
out << "" << endl;
out << "" << endl;
out << "" << endl;
out << "" << endl;
out << "var vis62 = vis" << endl;
out << ".add(pv.Bar) " << endl;
out << ".data(control_peptides)" << endl;
out << "//    .def('active', -1)" << endl;
out << ".fillStyle(function() { if (dot_active == this.index)  { return pv.rgb(128,0,0, 0.8); } else { return null ; } } )" << endl;
out << ".strokeStyle(function(d) { tmp = Math.floor(255*d.RTnobs+1/medNobs); return pv.rgb(tmp,255-tmp,255-tmp, 0.8) ; })" << endl;
out << "//.radius(function(d) { if (d.RTnobs+1 == 1) { return 2; } else { return d.RTsiqr/(4000/600) < 2 ? 2 : d.RTsiqr/(4000/600)  ;} })" << endl;
out << "//.size(function(d) { if (d.RTnobs+1 == 1) { return 5; } else { return 1.57*Math.pow(d.RTsiqr/(4000/600),2);} })" << endl;
out << "//.shape(function(d) { if (d.RTnobs+1 == 1) { return 'cross'; } else { return 'circle' ;} } )" << endl;
out << ".height(function(d) { return y(2*d.RTsiqr); } )" << endl;
out << ".width(x(1))" << endl;
out << ".left( pv.Scale.ordinal(control_peptides.map(function(d) { return d.run;}))" << endl;
out << "   .split(0, w)" << endl;
out << "   .by(function(d) { return d.run;} ))" << endl;
out << ".bottom(function(d) { return y(d.RTmed-d.RTsiqr); })    " << endl;
out << ".event('mouseover', function(d) { dot_obs = d.RTnobs+1; activeText=d.run+', SIQR='+d.RTsiqr+', NOBS='+d.RTnobs; activePep=d.seq; dot_active = this.index; this.parent.render();} ) " << endl;
out << "//.event('point', function(d) { dot_obs = d.RTnobs; activeText=d.run+', SIQR='+d.RTsiqr+', NOBS='+d.RTnobs; activePep=d.seq; dot_active = this.index; this.parent.render();} ) " << endl;
out << "//.event('unpoint', function(d) { dot_obs = ''; activeText = ''; activePep=''; dot_active = -1; vis.render();  } )" << endl;
out << ".event('mouseout',  function(d) { dot_obs = ''; activeText = ''; activePep=''; dot_active = -1; vis.render();  } )" << endl;
out << ".title( function(d) { return d.run; } )" << endl;
out << "//.anchor('right')" << endl;
out << "//.add(pv.Label)" << endl;
out << "//.visible(function() {return dot_active == this.index ;} )" << endl;
out << "//.text(function(d) { return d.run;} )" << endl;
out << "//.textAngle(-0.3);" << endl;
out << "" << endl;
out << "" << endl;
out << "var vis63 = vis" << endl;
out << ".add(pv.Dot) " << endl;
out << ".data(control_peptides)" << endl;
out << "//    .def('active', -1)" << endl;
out << ".fillStyle(function() { if (dot_active == this.index)  { return pv.rgb(128,0,0, 0.8); } else { return null ; } } )" << endl;
out << ".strokeStyle(function(d) { tmp = Math.floor(255*d.RTnobs+1/medNobs); return pv.rgb(tmp,255-tmp,255-tmp, 0.8) ; })" << endl;
out << "//.radius(function(d) { if (d.RTnobs+1 == 1) { return 2; } else { return d.RTsiqr/(4000/600) < 2 ? 2 : d.RTsiqr/(4000/600)  ;} })" << endl;
out << "//.size(function(d) { if (d.RTnobs+1 == 1) { return 5; } else { return 1.57*Math.pow(d.RTsiqr/(4000/600),2);} })" << endl;
out << ".shape(function(d) { if (d.RTnobs+1 == 1) { return 'cross'; } else { return 'tick' ;} } )" << endl;
out << ".angle(function(d) { if (d.RTnobs+1 >= 1)  return 1.57; }  )" << endl;
out << "//.height(function(d) { return y(2*d.RTsiqr); } )" << endl;
out << ".size(function(d) {if (d.RTnobs+1 == 1) { return x(5); } else { return x(1); }})" << endl;
out << ".left( pv.Scale.ordinal(control_peptides.map(function(d) { return d.run;}))" << endl;
out << "   .split(0, w)" << endl;
out << "   .by(function(d) { return d.run;} ))" << endl;
out << ".bottom(function(d) { return y(d.RTmed); })    " << endl;
out << ".event('mouseover', function(d) { dot_obs = d.RTnobs; activeText=d.run+', SIQR='+d.RTsiqr+', NOBS='+d.RTnobs; activePep=d.seq; dot_active = this.index; this.parent.render();} ) " << endl;
out << "//.event('point', function(d) { dot_obs = d.RTnobs; activeText=d.run+', SIQR='+d.RTsiqr+', NOBS='+d.RTnobs; activePep=d.seq; dot_active = this.index; this.parent.render();} ) " << endl;
out << "//.event('unpoint', function(d) { dot_obs = ''; activeText = ''; activePep=''; dot_active = -1; vis.render();  } )" << endl;
out << ".event('mouseout',  function(d) { dot_obs = ''; activeText = ''; activePep=''; dot_active = -1; vis.render();  } )" << endl;
out << ".title( function(d) { return d.run; } )" << endl;
out << "//.anchor('right')" << endl;
out << "//.add(pv.Label)" << endl;
out << "//.visible(function() {return dot_active == this.index ;} )" << endl;
out << "//.text(function(d) { return d.run;} )" << endl;
out << "//.textAngle(-0.3);" << endl;
out << "" << endl;
out << "" << endl;
out << "var tiplabel = vis.add(pv.Dot)" << endl;
out << ".data(function() { return [activeText]; })" << endl;
out << ".visible(function(){ return (dot_active != -1) ;} )" << endl;
out << ".fillStyle(function() { return pv.rgb(128,0,0, 0.8);} )" << endl;
out << ".strokeStyle(function(d) { tmp = Math.floor(255*dot_obs/medNobs); return pv.rgb(tmp,255-tmp,255-tmp, 0.8) ; })" << endl;
out << ".shape(function() { if ( dot_obs == 1) { return 'cross'; } else if (dot_obs > 1) { return 'bar' ;} else { return null; } } )" << endl;
out << ".bottom(30)" << endl;
out << ".left(20)" << endl;
out << ".add(pv.Label)" << endl;
out << ".anchor('right')" << endl;
out << ".data(function() { return [activeText]; } )" << endl;
out << ".text(function(d) { return d; } )" << endl;
out << ".textAngle(0);" << endl;
out << "" << endl;
out << "var peplabel = vis.add(pv.Label)" << endl;
out << ".data(function() { return [activePep]; })" << endl;
out << ".bottom(50)" << endl;
out << ".left(20)" << endl;
out << ".text(function(d) { return d; } )" << endl;
out << ".textAngle(0);" << endl;
out << "" << endl;
out << "" << endl;
out << "" << endl;
out << "var legend = vis.add(pv.Dot)" << endl;
out << ".shape('cross')" << endl;
out << ".strokeStyle(function() { tmp = Math.floor(255/medNobs); return pv.rgb(tmp,255-tmp,255-tmp, 0.8) ; })" << endl;
out << ".bottom(y(3900))" << endl;
out << ".left(x(20))" << endl;
out << ".anchor('right')" << endl;
out << ".add(pv.Label)" << endl;
out << ".left(x(30))" << endl;
out << ".text('NOBS = 1');" << endl;
out << "" << endl;
out << "var legend = vis.add(pv.Dot)" << endl;
out << ".shape('dot')" << endl;
out << ".strokeStyle(function() { tmp = Math.floor(255/medNobs); return pv.rgb(tmp,255-tmp,255-tmp, 0.8) ; })" << endl;
out << ".bottom(y(3700))" << endl;
out << ".left(x(20))" << endl;
out << ".anchor('right')" << endl;
out << ".add(pv.Label)" << endl;
out << ".left(x(30))" << endl;
out << ".text('NOBS > 1, Diameter = RT Inter-Quartile Range');" << endl;
out << "" << endl;
out << "" << endl;
out << "var nobs_scale=-1;" << endl;
out << "" << endl;
out << "" << endl;
out << "var nobsColors = vis.add(pv.Bar)" << endl;
out << ".data(pv.range(0, 1, 1/medNobs))" << endl;
out << ".bottom(y(3400))" << endl;
out << ".left(function() { return x(20) + this.index * x(5);})" << endl;
out << ".width(x(5))" << endl;
out << ".height(y(100))" << endl;
out << ".event('mouseover', function() { nobs_scale=this.index; vis.render(); } ) " << endl;
out << ".event('mouseout', function() { nobs_scale=-1; vis.render();} )" << endl;
out << "// .fillStyle(pv.Scale.linear(0, .5, 1).range('red', 'yellow', 'green'));" << endl;
out << ".fillStyle(pv.Scale.linear(0, 1, medNobs).range(pv.rgb(0,255,255, 0.8), pv.rgb(255,0,0, 0.8)))" << endl;
out << ".add(pv.Label)" << endl;
out << ".bottom(y(3400))" << endl;
out << ".visible(function() { return nobs_scale == this.index ;} )" << endl;
out << ".text(function() { if (nobs_scale == this.index) { tmp = this.index+1; } else { tmp = null; } return tmp; });" << endl;
out << "" << endl;
out << "" << endl;
out << "var nobsScaleLabelTitle= vis.add(pv.Label)" << endl;
out << ".bottom(y(3300))" << endl;
out << ".left(x(20))" << endl;
out << ".text('Number of Observations (NOBS) Scale');" << endl;
out << "" << endl;
out << "var nobsScaleLabelMin= vis.add(pv.Label)" << endl;
out << ".bottom(y(3400))" << endl;
out << ".left(x(20)-10)" << endl;
out << ".text('1');" << endl;
out << "" << endl;
out << "var nobsScaleLabelMax= vis.add(pv.Label)" << endl;
out << ".bottom(y(3400))" << endl;
out << ".left(x(20)+medNobs*x(5))" << endl;
out << ".text('> Median NOBS ('+medNobs+')');" << endl;





out << "vis.root.render();" << endl;
out << "</script></div></center>" << endl;
out << "</body></html>" << endl;


}

bool RTCatalog::addRun(string& run_name) {
  


  string* rn = new string(run_name);
  string* rf = new string(run_name);
  

  size_t pos1, pos2;
  pos1 =  rn->find_last_of("/\\");
  //pos2 =  rn->find_last_of(".");
  *rn =  rn->substr(pos1+1);//, pos2-pos1-1);

  if (disco_) {    
    pos2 = rn->find_last_of("_");
    
    *rn = rn->substr(0, pos2);
  }

  //if (run_name.find(".mzML") == string::npos) {
  //  *rn += ".mzML"; 
  //  *rf += ".mzML"; 
  //}


  if (byrun_peprtinfo_hash_->find(*rn) != byrun_peprtinfo_hash_->end()) {
    run_name = *rn;
    delete rf;
    delete rn;
    return false;
  }


  
  run_names_->insertAtEnd(rn);
  run_files_->insert(make_pair(*rn, *rf));
  byrun_peprtinfo_hash_->insert(make_pair(*rn, new rtinfo_hash()));
  byrun_peprts_hash_->insert(make_pair(*rn, new dblarr_hash()));
  byrun_pep_maxIntens_->insert(make_pair(*rn, new dbl_hash()));
  byrun_pep_maxMatchedIons_->insert(make_pair(*rn, new dbl_hash()));
  byrun_pepintens_hash_->insert(make_pair(*rn, new dblarr_hash()));
  byrun_pepMatchedIons_hash_->insert(make_pair(*rn, new dblarr_hash()));

  //ismix_run_hash_->insert(make_pair(*rn, mix_run));
  //  peprts_hash_->insertAtEnd(new dblarr_hash());
  //pepintens_hash_->insertAtEnd(new dblarr_hash());

  num_runs_++;
  run_name = *rn;
  return true;

}




#ifndef __LGPL__
void RTCatalog::calcRTStatsCombined(ostream& out) {
  gsl_vector* r;
  
  unsigned long i, j, k;
  unsigned long size = 0;
  peprtinfo_hash_->clear();
  //  out << "Peptide" << "\t"
  //<< "RTmed" << "\t"
  //<< "RTsiqr" << "\t"
  //<< "RTmean" << "\t"
  //<< "RTstdev" << "\t"
  //<< "RTmin" << "\t"
  //<< "RTnobs" << "\n";

  out << "PEPTIDE\tQ25\tMEDIAN\tQ75\tSIQR\tMEAN\tSTDEV\tMIN\tNOBS\n";

  for (strparr_hash::iterator s_it = (*pepruns_hash_).begin(); s_it != (*pepruns_hash_).end(); s_it++) {
  
    size = 0;
    //compute number of datapoints
    for (i=0; i<s_it->second->size(); i++) {
      for (j=0; j<(*(*(*byrun_pepintens_hash_)[*(*s_it->second)[i]])[s_it->first]).size(); j++) {
	//if ((*byrun_peprtinfo_hash_)[*(*s_it->second)[i]]->find(s_it->first) != (*byrun_peprtinfo_hash_)[*(*s_it->second)[i]]->end())  {
	if ((*(*(*byrun_peprts_hash_)[*(*s_it->second)[i]])[s_it->first])[j] > -1e6) {
	  size ++;
	}
      }
    }

    if (size == 0) {
      continue;
    }

 
    double maxIntens = -1;
    double maxMatchedIons = -1;
    size = 0;
    //recompute number of datapoints
    for (i=0; i<s_it->second->size(); i++) {
      maxIntens = (*(*byrun_pep_maxIntens_)[*(*s_it->second)[i]])[s_it->first];
      maxMatchedIons = (*(*byrun_pep_maxMatchedIons_)[*(*s_it->second)[i]])[s_it->first];

      if (maxMatchedIons > 0) {
	for (j=0; j<(*(*(*byrun_pepMatchedIons_hash_)[*(*s_it->second)[i]])[s_it->first]).size(); j++) {
	  if (fabs((*(*(*byrun_pepMatchedIons_hash_)[*(*s_it->second)[i]])[s_it->first])[j]-maxMatchedIons) < 0.01 && (*(*(*byrun_peprts_hash_)[*(*s_it->second)[i]])[s_it->first])[j] > -1e6 ) {
	    size++;
	  }
	}

      }
      else {
	for (j=0; j<(*(*(*byrun_pepintens_hash_)[*(*s_it->second)[i]])[s_it->first]).size(); j++) {
	  if ((maxIntens <= 0 || (*(*(*byrun_pepintens_hash_)[*(*s_it->second)[i]])[s_it->first])[j] / maxIntens >= 0.5) && 
	      (*(*(*byrun_peprts_hash_)[*(*s_it->second)[i]])[s_it->first])[j] > -1e6 ) {
	    size++;
	  }
	}
      }
    }


    double med = 0;
    double q25 = 0;
    double q75 = 0;
    double mean = 0;
    double min = 0;

    double sd = 0;
    double siqr = 0;

    //    gsl_vector_free(r);
  
    if (size > 0) {
      //populate vector with RTs
      r = gsl_vector_calloc(size);
      k = 0;
      for (i=0; i<s_it->second->size(); i++) {
	maxIntens = (*(*byrun_pep_maxIntens_)[*(*s_it->second)[i]])[s_it->first];
	maxMatchedIons = (*(*byrun_pep_maxMatchedIons_)[*(*s_it->second)[i]])[s_it->first];

	if (maxMatchedIons > 0) {
	  for (j=0; j<(*(*(*byrun_peprts_hash_)[*(*s_it->second)[i]])[s_it->first]).size(); j++) {
	    //Must be strong intensity
	    if (fabs((*(*(*byrun_pepMatchedIons_hash_)[*(*s_it->second)[i]])[s_it->first])[j] - maxMatchedIons) <0.01 && (*(*(*byrun_peprts_hash_)[*(*s_it->second)[i]])[s_it->first])[j] > -1e6 ) {
	      //gsl_vector_set(r, k, acn_gradient_->getAcn((*(*byrun_peprtinfo_hash_)[*(*s_it->second)[i]])[s_it->first]->med_));
	      //gsl_vector_set(r, k, acn_gradient_->getAcn((*(*(*byrun_peprts_hash_)[*(*s_it->second)[i]])[s_it->first])[j]));
	      gsl_vector_set(r, k, (*(*(*byrun_peprts_hash_)[*(*s_it->second)[i]])[s_it->first])[j]);
	      k++;
	    }
	    
	  }
	}
	else {
	  if (maxIntens <= 0) {
	    cerr << "WARNING: Couldn't find max Intensity for run=" << *(*s_it->second)[i] << ", pep=" << s_it->first << endl;
	  }
	  for (j=0; j<(*(*(*byrun_peprts_hash_)[*(*s_it->second)[i]])[s_it->first]).size(); j++) {
	    //Must be strong intensity
	    if ((maxIntens <= 0 ||(*(*(*byrun_pepintens_hash_)[*(*s_it->second)[i]])[s_it->first])[j] / maxIntens >= 0.5) && 
		(*(*(*byrun_peprts_hash_)[*(*s_it->second)[i]])[s_it->first])[j] > -1e6 ) {
	      //gsl_vector_set(r, k, acn_gradient_->getAcn((*(*byrun_peprtinfo_hash_)[*(*s_it->second)[i]])[s_it->first]->med_));
	      //gsl_vector_set(r, k, acn_gradient_->getAcn((*(*(*byrun_peprts_hash_)[*(*s_it->second)[i]])[s_it->first])[j]));
	      gsl_vector_set(r, k, (*(*(*byrun_peprts_hash_)[*(*s_it->second)[i]])[s_it->first])[j]);
	      k++;
	    }
	  }
	}
      }
      
      
      gsl_sort_vector(r);
       med = gsl_stats_median_from_sorted_data(gsl_vector_ptr(r, 0), 1, size);
       q25 = gsl_stats_quantile_from_sorted_data(gsl_vector_ptr(r, 0), 1, size, 0.25);
       q75 = gsl_stats_quantile_from_sorted_data(gsl_vector_ptr(r, 0), 1, size, 0.75);
       mean = gsl_stats_mean(gsl_vector_ptr(r, 0), 1, size);
       min = gsl_stats_min(gsl_vector_ptr(r, 0), 1, size);
      
       sd = 0;
       siqr = (q75 - q25)/2;
      
      if (size > 1) sd = gsl_stats_sd(gsl_vector_ptr(r, 0), 1, size);
      gsl_vector_free(r);
    }

    out << getPeptideString(s_it->first) << "\t"
        << q25 << "\t"
	<< med << "\t"
        << q75 << "\t"
	<< siqr << "\t"
	<< mean << "\t"
	<< sd << "\t"
	<< min << "\t"
	<< size << "\n";

    peprtinfo_hash_->insert(make_pair(s_it->first, new RTInfo(size, q25, med, q75, mean, sd, min)));
    
  
  }

}

void RTCatalog::writeRunRTStats(const char* file) {
  string name = string(file)+"_BYRUN";
  ofstream fout(name.c_str());

  fout << "PEPTIDE\tMEDIAN\t(UNCORR_MED)\tSIQR\tMEAN\t(UNCORR_MEAN)\tSTDEV\tNOBS\n";
  for (int k = 0; k < num_runs_ ; k++) {
    fout << *(*run_names_)[k] << "\n";

    
    for (rtinfo_hash::iterator it = (*byrun_peprtinfo_hash_)[*(*run_names_)[k]]->begin();
	 it != (*byrun_peprtinfo_hash_)[*(*run_names_)[k]]->end(); it++) {
      fout << getPeptideString(it->first) << "\t"
	   << it->second->med_ << "\t"
	   << "(" << it->second->origMed_ << ")" << "\t"
	   << it->second->siqr_ << "\t"
	   << it->second->mean_ << "\t"
	   << "(" << it->second->origMean_ << ")" << "\t"
	   << it->second->stdev_ << "\t"
	   << it->second->n_ << "\n";
    }


  }
  fout.close();
  
}

void RTCatalog::calcRTStatsByRun() {
  gsl_vector* r;
  
  int j, k, l;
  int size = 0;

  for (k = 0; k < num_runs_ ; k++) {
    (*byrun_peprtinfo_hash_)[*(*run_names_)[k]]->clear();
    for (dblarr_hash::iterator it = (*(*byrun_peprts_hash_)[*(*run_names_)[k]]).begin(); it != (*(*byrun_peprts_hash_)[*(*run_names_)[k]]).end(); it++) {
      size = 0;
      for (j=0; j< (*(*(*byrun_pepintens_hash_)[*(*run_names_)[k]])[it->first]).size(); j++) {
	if ((*(*(*byrun_peprts_hash_)[*(*run_names_)[k]])[it->first])[j] > 0 ) {
	  size++;
	}
      }
      //      size = (*(*(*byrun_peprts_hash_)[*(*run_names_)[k]])[it->first]).size();

      if (size == 0) {
	(*byrun_peprtinfo_hash_)[*(*run_names_)[k]]->insert(make_pair(it->first, new RTInfo(0, -1, -1, -1, -1, -1, -1)));    
	continue;
      }
      
      //compute precursor Intensity Stats
      r = gsl_vector_calloc(size);
      l = 0;
      for (j=0; j<size; j++) {
	if ((*(*(*byrun_peprts_hash_)[*(*run_names_)[k]])[it->first])[j] > 0 ) {
	  gsl_vector_set(r, l, (*(*(*byrun_pepintens_hash_)[*(*run_names_)[k]])[it->first])[j]);
	  l++;
	}
      }
      

      double maxIntens = gsl_stats_max(gsl_vector_ptr(r, 0), 1, size);
      
      size = 0;
      //recompute number of datapoints
      for (j=0; j< (*(*(*byrun_pepintens_hash_)[*(*run_names_)[k]])[it->first]).size(); j++) {
	if ((maxIntens <= 0 || (*(*(*byrun_pepintens_hash_)[*(*run_names_)[k]])[it->first])[j] / maxIntens >= 0.5) &&
	      (*(*(*byrun_peprts_hash_)[*(*run_names_)[k]])[it->first])[j] > -1e6) {
	  size++;
	}
      }


      //populate vector with RTs
      gsl_vector_free(r);
      double med = 0;
      double q25 = 0;
      double q75 = 0;
      double mean = 0;
      double min = 0;
      
      double sd = 0;
      double siqr = 0;


      if (size > 0) {
	r = gsl_vector_calloc(size);
	l = 0;
	for (j=0; j< (*(*(*byrun_peprts_hash_)[*(*run_names_)[k]])[it->first]).size(); j++) {
	  if ((maxIntens <= 0 ||  (*(*(*byrun_pepintens_hash_)[*(*run_names_)[k]])[it->first])[j] / maxIntens >= 0.5) && 
	      (*(*(*byrun_peprts_hash_)[*(*run_names_)[k]])[it->first])[j] > -1e6  ) {
	    gsl_vector_set(r, l, (*(*(*byrun_peprts_hash_)[*(*run_names_)[k]])[it->first])[j]);
	    l++;
	  }
	}
	
	
	gsl_sort_vector(r);
	med = gsl_stats_median_from_sorted_data(gsl_vector_ptr(r, 0), 1, size);
	q25 = gsl_stats_quantile_from_sorted_data(gsl_vector_ptr(r, 0), 1, size, 0.25);
	q75 = gsl_stats_quantile_from_sorted_data(gsl_vector_ptr(r, 0), 1, size, 0.75);
	mean = gsl_stats_mean(gsl_vector_ptr(r, 0), 1, size);
	min = gsl_stats_min(gsl_vector_ptr(r, 0), 1, size);
	
	sd = 0;
	//	siqr = (q75 - q25)/2;
	
	if (size > 1) sd = gsl_stats_sd(gsl_vector_ptr(r, 0), 1, size);

	gsl_vector_free(r);
      }
      //cout << s_it->first << "\t"
      //	  << med << "\t"
      //	  << siqr << "\t"
      //	  << mean << "\t"
      //	  << sd << "\t"
      //	  << size << "\n";
      
      (*byrun_peprtinfo_hash_)[*(*run_names_)[k]]->insert(make_pair(it->first, new RTInfo(size, q25, med, q75, mean, sd, min)));      
    }
  }
}

#endif
