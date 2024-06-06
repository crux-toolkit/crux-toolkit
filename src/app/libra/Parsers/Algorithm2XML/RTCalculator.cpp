/* ***********************************************************************

Program       : RTCalculator                                                       
Author        : David Shteynberg <dshteynb  AT systemsbiology.org>                                                       
Date          : 09.29.2010

Primary data object holding all mixture distributions for each precursor ion charge

Copyright (C) 2010 David Shteynberg

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

*************************************************************************  */

#include <time.h>
#include <assert.h>
#include <map>
#include <fstream>
#include <sstream>

#include "RTCalculator.h"


#define RTMAX 5600

map<string, double>* RTCalculator::SSR1Coeff_ = NULL;
map<string, double>* RTCalculator::SSR1NTermCoeff_ = NULL;

map<string, double>* RTCalculator::HelixProp_ = NULL;

map<string, double>* RTCalculator::HydroPath_ = NULL;

map<string, double>* RTCalculator::EludeIndex_ = NULL;

map<string, double>* RTCalculator::ClustProp_ = NULL;


RTCalculator::RTCalculator() { 
  rt_catalog_ = NULL;
  c_ = NULL;
  ready_ = false;
  run_name_ = new string("noname") ;
  run_RT_sum_ = 0;
  run_RT_count_ = 0;
  RTs_ = new Array<double>();
  used_ = new Array<bool>();
  learned_RTs_ = new Array<double>();
  learn_probs_ = new Array<double>();
  learn_ntts_ = new Array<int>();
  matrix_ = NULL;
  learn_matrix_ = NULL;
  min_RT_=-1;
  adjusted_ = false;
  initSSR1Coeff();
  initHelixProp();
  r_sq_ = 0;
  run_slope_ = 0;
  ann_ = NULL;
  run_intercept_ = 0;
  best_modpep_rts_ = new dblarr_hash();
  best_modpep_runs_ = new intarr_hash();
  best_modpep_ids_ = new intarr_hash();
  pIcalc_ = new pICalculator();
  peps_ = new Array<char*>();

  max_len_ = 49;

  used_count_ = 0;
  max_rt_ = 0;
  min_rt_ = 0;
  
  verbose_ = 0;

  ann_file_ = "";

#ifndef __LGPL__
  rt_fit_ = NULL;
  //rt_fit_ = gsl_multifit_linear_alloc(1, 1); 
#endif
  

}
RTCalculator::RTCalculator(string* run) {
  //DDS: RT model
  //DDS: TODO Move this stuff into the RT model class???
  c_ = NULL;
  rt_catalog_ = NULL;
  ready_ = false;
  run_name_ = run;
  run_RT_sum_ = 0;
  run_RT_count_ = 0;
  RTs_ = new Array<double>();
  used_ = new Array<bool>();
  learned_RTs_ = new Array<double>();
  learn_probs_ = new Array<double>();
  learn_ntts_ = new Array<int>();
  pIcalc_ = new pICalculator();
  peps_ = new Array<char*>();

  max_len_ = 49;

  max_rt_ = 0;
  min_rt_ = 0;
  
  used_count_ = 0;
  
  verbose_ = 0;

  ann_file_ = "";
   
  matrix_ = NULL;
  learn_matrix_ = NULL;
  r_sq_ = 0;
  run_slope_ = 0;
  run_intercept_ = 0;
#ifndef __LGPL__
  rt_fit_ = NULL;
  //  rt_fit_ = gsl_multifit_linear_alloc(1, 1); 
#endif
  adjusted_ = false;
  min_RT_=-1;
  initSSR1Coeff();
  initHelixProp();
  best_modpep_rts_ = new dblarr_hash();
  best_modpep_runs_ = new intarr_hash();
  best_modpep_ids_ = new intarr_hash();
  ann_ = NULL;
}

RTCalculator::~RTCalculator() {
  delete RTs_;
  delete used_;
  delete learned_RTs_;
  delete learn_probs_;
  delete learn_ntts_;
  delete best_modpep_rts_;
  delete best_modpep_runs_;
  delete best_modpep_ids_;
  delete pIcalc_;
  delete peps_;
  
  if (matrix_) delete (matrix_);
  if (c_) gsl_vector_free(c_);

  if (ann_) fann_destroy(ann_);
  ann_ = NULL;

}

double RTCalculator::Peptide_RT(char* seq, int scan) {
  double rtn = getSSRCalcHP(seq);
  // double rtn = getHydroEHP(seq);
  //double rtn = getHydroGHP(seq);
  RTs_->insertAtEnd(rtn);
  
  char* new_pep = new char[strlen(seq)+1];
  strcpy(new_pep, seq);
  //new_pep[strlen(seq)]='\0';
  //DDS: RT model
  scans_.insertAtEnd(scan);
  if (min_RT_ <= 0 || scan < min_RT_) {
    min_RT_ = scan;
  }
  run_RT_sum_ += rtn;
  run_RT_count_++;

  return rtn;
}

void RTCalculator::adjustRT() {
  if (! adjusted_)
    for (int i=0; i<(int)rts_.size(); i++) {
      rts_[i] = rts_[i] - min_RT_;
      
    }

  adjusted_ = true;

}

double RTCalculator::getPredRT(int i) { 
  if (RTs_ && i >= 0 && i < RTs_->size())
    return (*RTs_)[i]; 
  
  return -1;
}

double RTCalculator::getObsRT(int i) { 
  if (i >= 0 && i < (int)rts_.size())
    return rts_[i]; 
  
  return -1;
} 

const string* RTCalculator::getModPeptide(int i) { 
  if (i >= 0 && i < (int)modified_peps_.size())
    return &modified_peps_[i]; 
  
  return NULL;
}

int RTCalculator::isOutlier(int i) { 
  if (i >= 0 && i < used_->size())
    return (*used_)[i] ? 0 : 1; 
  
  return -1;
}


int RTCalculator::getPredRTCount() { 
  if (RTs_)
    return (*RTs_).size(); 
  
  return -1;
}


//void  RTCalculator::linearTxRT(Array<RTCalculator*>* run_calcs) {
//   if (! adjusted_ ) {
    
//     for (int i=0; i<(int)rts_.size(); i++) {
//       rts_[i] = (*run_calcs)[getRunNameatIdx(i)]->getRunSlope() *  (*RTs_)[i] + (*run_calcs)[getRunNameatIdx(i)]->getRunInt();
//     }
//     RTs_->clear();
//     run_RT_count_=0;
//     // nt i = 0;

// //     for (int d=0; d<(*run_calcs).size(); d++) {
// //       for (int dd=0; dd<(*run_calcs)[d]->getPredRTCount(); dd++) {
// // 	rts_[i++] = (*run_calcs)[d]->getRunSlope() *  (*run_calcs)[d]->getPredRT(dd) + (*run_calcs)[d]->getRunInt();
// //       }
// //     }



//     for (  dblarr_hash::iterator q = (*best_modpep_rts_).begin(); q != (*best_modpep_rts_).end(); q++) {
//       for (int j = 0; j < (*(*q).second).size(); j++) {
// 	(*(*best_modpep_rts_)[(*q).first])[j] = (*run_calcs)[(*(*best_modpep_runs_)[(*q).first])[j]]->getRunSlope() * 	(*(*best_modpep_rts_)[(*q).first])[j] + (*run_calcs)[(*(*best_modpep_runs_)[(*q).first])[j]]->getRunInt();
//       }
//     }
    

//   }
//   adjusted_ = true;
//}

void RTCalculator::distillModifiedPeptides(Array<double>* probs, Array<int> * ntts) {
  //scans_.clear();
  //modified_peps_.clear();
  //learn_rts_.clear();
  //RTs_->clear();
  //scans_.clear();
  //probs_->clear();
  //ntts_->clear();

  for (  dblarr_hash::iterator q = (*best_modpep_rts_).begin(); q != (*best_modpep_rts_).end(); q++) {
    cerr << (*q).first << "\t";

    double maxPr = 0;
    int ntt = 0;
    intarr_hash::iterator id = best_modpep_ids_->find((*q).first);
    for (int j = 0; j < (*(*id).second).size(); j++) {
      if (j == 0 || (*probs)[(*(*id).second)[j]] > maxPr) {
	maxPr = (*probs)[(*(*id).second)[j]];
	ntt = (*ntts)[(*(*id).second)[j]];
      }
    }


    gsl_vector* r = gsl_vector_calloc((*(*q).second).size());
    for (int j = 0; j < (*(*q).second).size(); j++) {
      if ((*probs)[(*(*id).second)[j]] == maxPr) {
	gsl_vector_set(r, j, (*(*q).second)[j]);
	cerr << (*(*q).second)[j] << "\t";
      }
      
    }
    
    gsl_sort_vector(r);
    double med = gsl_stats_median_from_sorted_data(gsl_vector_ptr(r, 0), 1, (*(*q).second).size());
    gsl_vector_free(r);
    cerr << med << endl;;
    //modified_peps_.push_back((*q).first);
    learn_peps_.push_back((*q).first);
    learn_rts_.push_back(med);
    
    

    learn_probs_->insertAtEnd(maxPr);
    learn_ntts_->insertAtEnd(ntt);
  }
  
   
  
}


void RTCalculator::addPeptide_RT(const char* seq, const char* modified_pep, int scan, double rt, string* run_name) {
  //double rtn = getSSRCalcHP(seq);
  // double rtn = getHydroEHP(seq);
  //double rtn = getHydroGHP(seq);
  //RTs_->insertAtEnd(rtn);
  string pep_seq = string(seq);
  
  string unmod_pep = stripPeptide(pep_seq);
  
  
  char* new_pep = new char[strlen(unmod_pep.c_str())+1];
  strcpy(new_pep, unmod_pep.c_str());
  //new_pep[strlen(seq)]='\0';
  peps_->insertAtEnd(new_pep);
  //DDS: RT model
  scans_.insertAtEnd(scan);
  
  modified_peps_.push_back(modified_pep);
  
  double tmp;

  // HENRY: Also score experimental RT value (0.0 means unknown, in which case
  // the scan number will be used as a pseudo-RT)
  if (rt < 0.0001) {
    tmp = scan;
  } else {
    tmp = rt;
  }
  rts_.push_back(tmp);

  used_->insertAtEnd(true);
  used_count_++;

  
  if (min_RT_ <= 0 || tmp < min_RT_) {
    min_RT_ = tmp;
  }

  run_ids_.push_back(run_name);
  //return rtn;
}

double RTCalculator::getSSRCalcHP(char* peptide) {
  double out;
  char* szBuf = new char[2000];

  sprintf(szBuf, "SSRCalc=%s", getBinPath());
  putenv(szBuf);

  sprintf(szBuf, "%sSSRCalc3.pl --alg 3.0 --seq %s --output tsv | cut -f3", getBinPath(), peptide);

  FILE* pipe;
  if ((pipe = tpplib_popen(szBuf, "r")) == NULL) {
    cout << "error calling SSRCalc3.pl" << endl;
    exit(1);
  }
  else {
    char *res = fgets(szBuf, 2000, pipe);
    pclose(pipe);
    szBuf[strlen(szBuf)-1] = 0;
    out = atof(szBuf);
  }
  delete [] szBuf;

  return out;

}


// void RTCalculator::batchSRCalcHP() {
//   //  RTs_->reserve(peps_->size());

//   RTs_->reserve(rts_.size());
//   for (int dd = 0; dd < rts_.size(); dd++) {

//     string pep =  modified_peps_[dd];

//     double predictedRT = calcSSR1(pep);
  
//     //    run_RT_count_++;
//     run_RT_sum_ += predictedRT;
    
//     RTs_->insertAtEnd(predictedRT);
//   }
//   run_RT_count_ = RTs_->size();

//   learned_RTs_->reserve(learn_rts_.size());
//   for (int dd = 0; dd < learn_peps_.size(); dd++) {

//     string pep =  learn_peps_[dd];

//     double predictedRT = calcSSR1(pep);
  
//     //    run_RT_count_++;
//     learned_RT_sum_ += predictedRT;
    
//     learned_RTs_->insertAtEnd(predictedRT);
//   }
//   learned_RT_count_ = learned_RTs_->size();
  
// }
  
void RTCalculator::batchRTCalc() {
  //  RTs_->reserve(peps_->size());
  RTs_->clear();
  RTs_->reserve((int)rts_.size());
  run_RT_sum_ = 0;
  for (int dd = 0; dd < (int)rts_.size(); dd++) {

    string pep =  modified_peps_[dd];

    double predictedRT = 0;

    if (ann_) {
      predictedRT = calcANN_PepRT2(pep);
    }
    else {
      predictedRT = calc_PepRT(pep);
    }
  
    //    run_RT_count_++;
    run_RT_sum_ += predictedRT;
    
    RTs_->insertAtEnd(predictedRT);
  }
  run_RT_count_ = RTs_->size();
  
}


void RTCalculator::batchRTCalc(GradientProgram* acnProg) {
  //  RTs_->reserve(peps_->size());
  RTs_->clear();
  RTs_->reserve((int)rts_.size());
  run_RT_sum_ = 0;
  for (int dd = 0; dd < (int)rts_.size(); dd++) {

    string pep =  modified_peps_[dd];

    double predictedRT = acnProg->getAcn((calc_PepRT(pep)));
  
    //    run_RT_count_++;
    run_RT_sum_ += predictedRT;
    
    RTs_->insertAtEnd(predictedRT);
  }
  run_RT_count_ = RTs_->size();
  
}


void RTCalculator::batchRTCalc(const char* catfile) {
  RTs_->clear();
  RTs_->reserve(peps_->size());
  run_RT_sum_ = 0;
  if (rt_catalog_ == NULL) {
    rt_catalog_ = new RTCatalog(catfile);
    
      for (int dd = 0; dd < (int)rts_.size(); dd++) {
	string pep =  modified_peps_[dd];

	RTInfo* info = rt_catalog_->getRTInfo(pep);
	double nextRT = rt_catalog_->getRTMedian(pep);
	
	if (!info) {

	  if (ann_)
	    nextRT = calcANN_PepRT2(pep);
	  else 
	    nextRT = calc_PepRT(pep);
	}
	RTs_->insertAtEnd(nextRT);
	
      }
  }
  run_RT_count_ = RTs_->size();
}
  
void RTCalculator::read_RTCatalog(const char* catfile) {
  if (rt_catalog_ == NULL) {
    rt_catalog_ = new RTCatalog(catfile);
  }
}

void RTCalculator::batchRTCatalog(const RTCatalog& catalog) {
  RTs_->clear();
  RTs_->reserve((int)rts_.size());
  run_RT_sum_ = 0;
  for (int dd = 0; dd < (int)rts_.size(); dd++) {
    
    string pep =  modified_peps_[dd];
    
    double predictedRT = catalog.getRTMedian(pep);

    if (predictedRT < 0) {
      cerr << "ERROR: Peptide " << pep << " not found in Preliminary Catalog." << endl;
    }
    
    //    run_RT_count_++;
    run_RT_sum_ += predictedRT;
    
    RTs_->insertAtEnd(predictedRT);
  }
  run_RT_count_ = RTs_->size();
}


void RTCalculator::batch_iRT(dbl_hash& irt_hash) {
  RTs_->clear();
  RTs_->reserve((int)rts_.size());
  run_RT_sum_ = 0;
  for (int dd = 0; dd < (int)rts_.size(); dd++) {
    
    string pep =  modified_peps_[dd];
    
    double predictedRT = -9999999;
    dbl_hash::iterator it = irt_hash.find(pep);
    if (it != irt_hash.end()) {
      predictedRT = (*it).second;
    }

    if (predictedRT < -1000000) {
      cerr << "ERROR: Peptide " << pep << " not found in Preliminary Catalog." << endl;
    }
    else {
    //    run_RT_count_++;
      run_RT_sum_ += predictedRT;
    
      RTs_->insertAtEnd(predictedRT);
    }
  }
  run_RT_count_ = RTs_->size();
}


double RTCalculator::get_PepRTCatalog(string& pep) {
  RTInfo* rtinfo;
  rtinfo_hash::iterator itr = rt_catalog_->peprtinfo_hash_->find(pep);
  if ( itr != rt_catalog_->peprtinfo_hash_->end()) {
    rtinfo = itr->second;
    //if (rtinfo->siqr_ > 90) {
    //  cerr << "Skipping peptide " << pep << " which has SIQR of " << rtinfo->siqr_ << endl;
    //  return -2;
    //}
    //if (rtinfo->n_  < 3) {
    //  cerr << "Skipping peptide " << pep << " which has Num Obs of " << rtinfo->n_ << endl;
    //  return -2;
    //}
    //    return rtinfo->min_;  
    return rtinfo->med_;  
    
  }
  return -1000000000000.0;
  
}


void RTCalculator::batchSSRCalcHP() {
  double out;
  char pep_file_str[1028];
  time_t now; time(&now);
  struct tm now_time = *localtime(& now);
  srandom(now_time.tm_mday + now_time.tm_hour + now_time.tm_min + now_time.tm_sec);
  sprintf(pep_file_str, "%s.SSRCalc.input.peps.%ld.tmp", getWebserverTmpPath()?getWebserverTmpPath():"./",(long)random());
  FILE* pep_file;
  if ((pep_file = fopen(pep_file_str, "w")) == NULL) {
    cout << "error calling writing tmp file" << endl;
    exit(1);
  }
  else {
    for (int i=0; i<peps_->size(); i++) {
      fprintf(pep_file, "%s\n", (*peps_)[i]);
    }
  }
  fclose(pep_file);
  char* szBuf = new char[2000];
  run_RT_sum_  = 0;
  run_RT_count_ = 0;

  sprintf(szBuf, "SSRCalc=%s", getBinPath());
  putenv(szBuf);

  sprintf(szBuf, "%sSSRCalc3.pl --alg 3.0 --source %s --output tsv | cut -f3", getBinPath(), pep_file_str);
  FILE* pipe;
  if ((pipe = tpplib_popen(szBuf, "r")) == NULL) {
    cout << "error calling SSRCalc3.pl" << endl;
    exit(1);
  }
  else {
    while (fgets(szBuf, 2000, pipe) != NULL) {
      szBuf[strlen(szBuf)-1] = 0;
      out = atof(szBuf);
      RTs_->insertAtEnd(out);
      learned_RTs_->insertAtEnd(out);
      run_RT_sum_ += out;
      run_RT_count_++;
    }
  }
  pclose(pipe);
  sprintf(szBuf, "rm -f %s", pep_file_str);
  int ret = tpplib_system(szBuf);
  delete [] szBuf;
}
  

double RTCalculator::getHydroEHP(char* peptide) {
  double out;
  char* szBuf = new char[2000];
  sprintf(szBuf, "%s -e 'use Hydro; print Hydro::calcHydrophobicityE( sequence => \"%s\" );'", getPerlBinary(), peptide);
  FILE* pipe;
  if ((pipe = tpplib_popen(szBuf, "r")) == NULL) {
    cout << "error calling Hydro.pm" << endl;
    exit(1);
  }
  else {
    char *res = fgets(szBuf, 2000, pipe);
    pclose(pipe);
    szBuf[strlen(szBuf)-1] = 0;
    out = atof(szBuf);
  }
  delete [] szBuf;
  return out;

}

double RTCalculator::getHydroGHP(char* peptide) {
  double out;
  char* szBuf = new char[2000];
  sprintf(szBuf, "%s -e 'use Hydro; print Hydro::calcHydrophobicityG( sequence => \"%s\" );'", getPerlBinary(), peptide);
  FILE* pipe;
  if ((pipe = tpplib_popen(szBuf, "r")) == NULL) {
    cout << "error calling Hydro.pm" << endl;
    exit(1);
  }
  else {
    char *res = fgets(szBuf, 2000, pipe);
    pclose(pipe);
    szBuf[strlen(szBuf)-1] = 0;
    out = atof(szBuf);
  }
  delete [] szBuf;

  return out;

}


void RTCalculator::recalc_RTstats(Array<double>* probs) {
  double numer = 0;
  double denom = 0;
  double tot_RTs = 0;
  double prob;
  run_RT_used_count_ = 0;
  // cout << "DDS: probsSize=" << probs->size() << " runRTNum=" << run_RT_count_<< endl;
  //assert(probs->size() == run_RT_count_);
  int i;
  for (i=0; i<run_RT_count_; i++) {
    if ((*probs)[i] >= 0) {
      run_RT_used_count_ ++;
      prob = (*probs)[i]*(*probs)[i];
      tot_RTs += (*RTs_)[i];
      numer += prob*(*RTs_)[i];
      denom += prob;
    }
  }

  if (denom == 0) {
    numer = tot_RTs;
    denom = run_RT_count_;
    
    run_RT_mean_ = numer / denom;
    run_RT_stddev_ = 0;

    for (i=0; i<run_RT_count_; i++) {
      run_RT_stddev_ += pow((run_RT_mean_ - (*RTs_)[i]), 2);
    }
  }
  else {
    run_RT_mean_ = numer / denom;
    run_RT_stddev_ = 0;
    
    for (i=0; i<run_RT_count_; i++) {
      if ((*probs)[i] >= 0) {
	prob = (*probs)[i]*(*probs)[i];
	run_RT_stddev_ += prob * pow((run_RT_mean_ - (*RTs_)[i]), 2);
      }
    }
  }
  run_RT_stddev_ /= denom;

  run_RT_stddev_ = pow(run_RT_stddev_, 0.5);  

}

void RTCalculator::write_RTann(string& file) {
  if (file.empty()) 
    file = getConfPath() + (string)"RTtrain";

  ann_file_ = file + ".fann";

  fann_save(ann_, ann_file_.c_str());
  
}



void RTCalculator::write_RTstats(ostream& out) {
  // out <<  run_idx_ << "\t" << run_slope_ << "\t" << run_intercept_ << "\t" << r_sq_ << "\t" << run_RT_used_count_<< endl;
  
  out << endl;
  out << "Rsq = " << r_sq_ << " slope = " << run_slope_ << " intercept = " << run_intercept_ << endl;
  out << "Coefficients:" << endl;
  
  for (map<string, double>::iterator q = coeff_.begin(); q != coeff_.end(); q++) {
    out << "  " << q->first << " : " << q->second << endl;
  }

}

void RTCalculator::write_RTcoeff(ostream& out) {
  string param; 
  for (map<string, double>::iterator q = coeff_.begin(); q != coeff_.end(); q++) {
    out << q->first << " " << q->second << endl;
  }

}

void RTCalculator::write_paramIndex(ostream& out) {
  string param; 
  long pc=0;
  out << "ANN_FILE" << endl;
  out << ann_file_ << endl;
  out << "MAX_RT" << endl;
  out << max_rt_ << endl;
  out << "AAINDEX" << endl;
  for (map<string, int>::iterator q = aa_index_.begin(); q != aa_index_.end(); q++) {
    out << q->first << " " << q->second << endl;
  }
  
  out << "PARAMINDEX" << endl;
  for (map<string, int>::iterator q = param_index_.begin(); q != param_index_.end(); q++) {
    out << q->first << " " << q->second << endl;
    if (!pc || pc < q->second+1) {
      pc = q->second+1;
    }
    
  }
  out << "PARAMCOUNT" << endl << pc << endl;

}


void RTCalculator::read_RTcoeff() {
  string coeffile =  getConfPath() + (string)"RTCalc.coeff";
  
  ifstream fin(coeffile.c_str());
  if(! fin) {
    cerr << "cannot read coefficients file " << coeffile << endl;
    exit(1);
  }
  
  read_RTcoeff(fin); 
}

void RTCalculator::read_RTcoeff(const char* coeffile) {
  ifstream fin(coeffile);
  if(! fin) {
    cerr << "cannot read coefficients file " << coeffile << endl;
    exit(1);
  }
  
  read_RTcoeff(fin); 
}



void RTCalculator::read_RTann() {
  string file =  getConfPath() + (string)"RTtrain.fann";
  ann_file_ = file;
  ann_ = fann_create_from_file(ann_file_.c_str());
}

void RTCalculator::read_RTann(string& file) {
  ann_file_ = file;
  ann_ = fann_create_from_file(ann_file_.c_str()); 
}

void RTCalculator::read_paramIndex(istream& in) {
  string param;
  int val = 0;
  bool aaindex = false;

  paramcount_ = 0;
  max_rt_ = 0;
  min_rt_ = 0;
  aacount_ = 0;
  while (1) {
    in >> param;
    
    if (param == "ANN_FILE") {
      in >> ann_file_;
      read_RTann(ann_file_);
      continue;
    }
    else if (param == "MAX_RT") {
      in >> max_rt_;
      continue;
    }
    else if (param == "AAINDEX") {
      aaindex = true;
      in >> param;
    }
    else if (param == "PARAMCOUNT") {
      in >> paramcount_;
    }
    else if (param == "PARAMINDEX") {
      aaindex = false;
      in >> param;
    }

    in >> val;
    if (in.fail()) {
      break;
    }
    if (!aaindex) {
      param_index_[param] = val;
      //paramcount_++;
    }
    else {
      aa_index_[param] = val;
      aacount_++;
    }
  }
  //paramcount_ += max_len_ * aacount_;
}

void RTCalculator::read_RTcoeff(istream& in) {
  string param;
  double val = 0;
  paramcount_ = 0;
  while (1) {
    in >> param >> val;
    if (in.fail()) {
      break;
    }
    coeff_[param] = val;
    param_index_[param] = paramcount_++;
  }
  
  c_ = gsl_vector_calloc(paramcount_);

  for (map<string, int>::iterator pp = param_index_.begin(); pp != param_index_.end(); pp++) {
    gsl_vector_set(c_, pp->second, coeff_[pp->first]);
  }
  

}

void RTCalculator::train_RTcoeff(istream& in) {
  string pep; 
  double rt;

  learn_peps_.clear();
  learn_rts_.clear();


  
  while (1) {
    in >> pep >> rt ;
    unsigned int len = (unsigned int)countAA(pep);

    if (len > max_len_) {
      max_len_ =  len;
    }

    if (rt > max_rt_) {
      max_rt_ =  rt;
    }

    if (rt < min_rt_ || min_rt_ <= 0) {
      min_rt_ =  rt;
    }
    if (in.fail()) {
      break;
    }


    dblarr_hash::iterator q = best_modpep_rts_->find(string(pep));
    
    if (q == best_modpep_rts_->end() ) {
      best_modpep_rts_->insert(make_pair(pep, new Array<double>()));
    }
    (*(*best_modpep_rts_)[pep]).insertAtEnd(rt);

    //learn_peps_.push_back(pep);
    //learn_rts_.push_back(rt);

  }



  for (  dblarr_hash::iterator q = (*best_modpep_rts_).begin(); q != (*best_modpep_rts_).end(); q++) {

    //if ((*(*q).second).size() < 2) {
    //  continue;
    //}

    gsl_vector* r = gsl_vector_calloc((*(*q).second).size());
    for (int j = 0; j < (*(*q).second).size(); j++) {
      gsl_vector_set(r, j, (*(*q).second)[j]);
    }
      
  
    gsl_sort_vector(r);
    double med = gsl_stats_median_from_sorted_data(gsl_vector_ptr(r, 0), 1, (*(*q).second).size());
    double mn = gsl_stats_mean(gsl_vector_ptr(r, 0), 1, (*(*q).second).size());
    double sd = gsl_stats_sd(gsl_vector_ptr(r, 0), 1, (*(*q).second).size());
    gsl_vector_free(r);

 

//     r = gsl_vector_calloc((*(*q).second).size());
//     for (int j = 0; j < (*(*q).second).size(); j++) {
//       if (fabs((*(*q).second)[j]-mn)/sd < 2) {
// 	gsl_vector_set(r, j, (*(*q).second)[j]);
//       }
//     }
//     gsl_sort_vector(r);

//     med = gsl_stats_median_from_sorted_data(gsl_vector_ptr(r, 0), 1, (*(*q).second).size());
//     mn = gsl_stats_mean(gsl_vector_ptr(r, 0), 1, (*(*q).second).size());
//     sd = gsl_stats_sd(gsl_vector_ptr(r, 0), 1, (*(*q).second).size());
//     gsl_vector_free(r);


    // if ((*(*q).second).size() > 1 && sd > 30 ) {
    //  continue;
    //}
 

    learn_peps_.push_back((*q).first);
    learn_rts_.push_back(med);


  }


}



void RTCalculator::train_RTCatalog(const char* catfile) {
  string pep; 
  double rt;
  RTInfo* rtinfo;
  learn_peps_.clear();
  learn_rts_.clear();

  rt_catalog_ = new RTCatalog(catfile);

  for (rtinfo_hash::iterator itr = rt_catalog_->peprtinfo_hash_->begin(); itr != rt_catalog_->peprtinfo_hash_->end(); itr++) {
    
    pep = itr->first;
    rtinfo = itr->second;
    rt = rtinfo->med_;


    unsigned int len = (unsigned int)countAA(pep);

    if (len > max_len_) {
      max_len_ =  len;
    }

    if (rt > max_rt_) {
      max_rt_ =  rt;
    }
    if (rt < min_rt_) {
      min_rt_ =  rt;
    }
    if (rtinfo->siqr_ > 1) {
      cerr << "Skipping peptide " << pep << " which has SIQR of " << rtinfo->siqr_ << endl;
      continue;
    }
    if (rtinfo->n_  < 5) {
      cerr << "Skipping peptide " << pep << " which has Num Obs of " << rtinfo->n_ << endl;
      continue;
    }
    
    


    dblarr_hash::iterator q = best_modpep_rts_->find(string(pep));
    
    if (q == best_modpep_rts_->end() ) {
      best_modpep_rts_->insert(make_pair(pep, new Array<double>()));
    }
    (*(*best_modpep_rts_)[pep]).insertAtEnd(rt);



  }



  for (  dblarr_hash::iterator q = (*best_modpep_rts_).begin(); q != (*best_modpep_rts_).end(); q++) {

    //if ((*(*q).second).size() < 2) {
    //  continue;
    //}

    gsl_vector* r = gsl_vector_calloc((*(*q).second).size());
    for (int j = 0; j < (*(*q).second).size(); j++) {
      gsl_vector_set(r, j, (*(*q).second)[j]);
    }
      
  
    gsl_sort_vector(r);
    double med = gsl_stats_median_from_sorted_data(gsl_vector_ptr(r, 0), 1, (*(*q).second).size());
    double mn = gsl_stats_mean(gsl_vector_ptr(r, 0), 1, (*(*q).second).size());
    double sd = gsl_stats_sd(gsl_vector_ptr(r, 0), 1, (*(*q).second).size());
    gsl_vector_free(r);

 

//     r = gsl_vector_calloc((*(*q).second).size());
//     for (int j = 0; j < (*(*q).second).size(); j++) {
//       if (fabs((*(*q).second)[j]-mn)/sd < 2) {
// 	gsl_vector_set(r, j, (*(*q).second)[j]);
//       }
//     }
//     gsl_sort_vector(r);

//     med = gsl_stats_median_from_sorted_data(gsl_vector_ptr(r, 0), 1, (*(*q).second).size());
//     mn = gsl_stats_mean(gsl_vector_ptr(r, 0), 1, (*(*q).second).size());
//     sd = gsl_stats_sd(gsl_vector_ptr(r, 0), 1, (*(*q).second).size());
//     gsl_vector_free(r);


    // if ((*(*q).second).size() > 1 && sd > 30 ) {
    //  continue;
    //}
 

    learn_peps_.push_back((*q).first);
    learn_rts_.push_back(med);


  }


}


double RTCalculator::calcANN_PepRT(string& pep) {
  fann_type *input = new fann_type[paramcount_];
  for (int ii=0; ii<paramcount_; ii++) {
    input[ii]=0;
  }
  int paramcount = paramcount_;
  int totAA = countAA(pep);
  string::size_type pos = (string::size_type)0;
  double score;
  map<string, int>::iterator found;
  
  // count numbers of each type of amino acids
  
  int numAA = 0;
  string aatag, nextaa;
  while (pos != string::npos || numAA == totAA-1) {
      
      if (numAA == 0) {
	aatag = nextAAToken(pep, pos, pos);
      }
      else {
	aatag = nextaa;
      }
      
      score = 1.;
      if (numAA < totAA-1) {
	nextaa =  nextAAToken(pep, pos, pos);
	
      }
      found = param_index_.find(aatag);
      if (found != param_index_.end()) {
	input[found->second]+=score;
      }     
      
      if (aatag.find_first_of("nc") == string::npos && (numAA < 3 || numAA >= totAA-3)) {
	 score = 1;
	if (numAA == 0 ) {
	  aatag = "termN1_" + aatag;

	}
	if (numAA == totAA-1) { 
	  aatag = "termC1_" + aatag;


	}
	if (numAA == 1) {
	  aatag = "termN2_" + aatag;
	}
	if ( numAA == totAA-2)  {
	  aatag = "termC2_" + aatag;

	}
	
	if (numAA == 2 ) {
	  aatag = "termN3_" + aatag;
	}
	if (numAA == totAA-3)  {
	  aatag = "termC3_" + aatag;
	}
	
	found = param_index_.find(aatag);
	if (found != param_index_.end()) {
	  input[found->second]+=score;
	}     
      }
      
      
      numAA++;
      
    }
    
    
    // found = param_index_.find("1");
    //if (found != param_index_.end()) {
    //  input[found->second]=1.;
    //}   
    
    found = param_index_.find("LEN");
    if (found != param_index_.end()) {
      input[found->second]=(double)numAA;
    }   
    
    found = param_index_.find("pI");
    if (found != param_index_.end()) {
      //input[found->second]=log(pIcalc_->ModifiedPeptide_pI(pep.c_str(),  0, pep.length()));  
      input[found->second]=isoelectric(pep);
      
    }   
    
    
    found = param_index_.find("helicity1");
    if (found != param_index_.end()) {
      input[found->second]=helicity1(pep);
    }   
    
    found = param_index_.find("clusterness");
    if (found != param_index_.end()) {
      input[found->second]=clusterness(pep);
    }   
    
    found = param_index_.find("helicity2");
    if (found != param_index_.end()) {
      input[found->second]=helicity2(pep);
    }   
    
    found = param_index_.find("helectric");
    if (found != param_index_.end()) {
      input[found->second]=helectric(pep);
    }   

    found = param_index_.find("hydroMoment");
    if (found != param_index_.end()) {
      input[found->second]=hydroMoment(pep);
    }   

    found = param_index_.find("EludeAmphipathicityHelixMin");
    if (found != param_index_.end()) {
      EludeAmphipathicityHelix(pep, input);
    }   
    
    found = param_index_.find("EludeHydrophobicMoment100Min");
    if (found != param_index_.end()) {
      EludeHydrophobicMoments(pep, input);
    }   

    found = param_index_.find("EludeIndexPartialSum5Max");
    if (found != param_index_.end()) {
      EludeIndexPartialSum(pep, input);
    }   

    found = param_index_.find("EludeIndexSum");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexSum(pep);
    }   
    found = param_index_.find("EludeIndexAvg");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexAvg(pep);
    }   
    found = param_index_.find("EludeIndexN");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexN(pep);
    }   
    found = param_index_.find("EludeIndexC");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexC(pep);
    }   
    found = param_index_.find("EludeIndexNearestNeighborPos");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexNearestNeighborPos(pep);
    }   
    found = param_index_.find("EludeIndexNearestNeighborNeg");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexNearestNeighborNeg(pep);
    } 
    found = param_index_.find("EludeIndexSumSqrDiff");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexSumSqrDiff(pep);
    }   

 
    fann_type *calc_out;

    calc_out = fann_run(ann_, input);

    double predictedRT = calc_out[0];
    
	delete[] input;
    return predictedRT*max_rt_;

}
double RTCalculator::calcANN_PepRT2(string& pep) {
  fann_type *input = new fann_type[paramcount_];
  for (int ii=0; ii<paramcount_; ii++) {
    input[ii]=0;
  }
  int paramcount = paramcount_;
  int totAA = countAA(pep);
  string::size_type pos = (string::size_type)0;
  double score;
  map<string, int>::iterator found;
  
  // count numbers of each type of amino acids
  
  int numAA = 0;
  string aatag, nextaa;
  while (pos != string::npos || numAA == totAA-1) {
            
      if (numAA == 0) {
	aatag = nextAAToken(pep, pos, pos);
      }
      else {
	aatag = nextaa;
      }
      
      score = 1.;
      if (numAA < totAA-1) {
	nextaa =  nextAAToken(pep, pos, pos);
	
      }
      found = aa_index_.find(aatag);
      int pep_pos = numAA;
      if (numAA < totAA / 2) {
	pep_pos = numAA;
      }
      else {
	pep_pos = max_len_ - (totAA - numAA);
      }
      if (found != aa_index_.end()) {
		int aa_i = found->second;
	for (int ii=0; ii<aacount_; ii++) {
	
	  if (ii == aa_i) {
	    input[pep_pos*aacount_+ii]= 1;
	  }
	  else {
	    input[pep_pos*aacount_+ii]= 0;
	  }
	    
	}     
      }     
      
      
      
      numAA++;
      
    }
    
    
    found = param_index_.find("LEN");
    if (found != param_index_.end()) {
      input[found->second]=(double)numAA;
    }   
    
    found = param_index_.find("pI");
    if (found != param_index_.end()) {
      input[found->second]=pIcalc_->ModifiedPeptide_pI(pep.c_str(),  0, (int)pep.length());  
      //      input[found->second]=isoelectric(pep);
      
    }   
    
    
    found = param_index_.find("helicity1");
    if (found != param_index_.end()) {
      input[found->second]=helicity1(pep);
    }   
    
    found = param_index_.find("clusterness");
    if (found != param_index_.end()) {
      input[found->second]=clusterness(pep);
    }   
    
    found = param_index_.find("helicity2");
    if (found != param_index_.end()) {
      input[found->second]=helicity2(pep);
    }   
    
    found = param_index_.find("helectric");
    if (found != param_index_.end()) {
      input[found->second]=helectric(pep);
    }   

    found = param_index_.find("hydroMoment");
    if (found != param_index_.end()) {
      input[found->second]=hydroMoment(pep);
    }   


    found = param_index_.find("EludeAmphipathicityHelixMin");
    if (found != param_index_.end()) {
      EludeAmphipathicityHelix(pep, input);
    }   
    
    found = param_index_.find("EludeHydrophobicMoment100Min");
    if (found != param_index_.end()) {
      EludeHydrophobicMoments(pep, input);
    }   
    
    found = param_index_.find("EludeIndexPartialSum5Max");
    if (found != param_index_.end()) {
      EludeIndexPartialSum(pep, input);
    }   

    found = param_index_.find("EludeIndexSum");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexSum(pep);
    }   
    found = param_index_.find("EludeIndexAvg");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexAvg(pep);
    }   
    found = param_index_.find("EludeIndexN");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexN(pep);
    }   
    found = param_index_.find("EludeIndexC");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexC(pep);
    }   
    found = param_index_.find("EludeIndexNearestNeighborPos");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexNearestNeighborPos(pep);
    }   
    found = param_index_.find("EludeIndexNearestNeighborNeg");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexNearestNeighborNeg(pep);
    } 
    found = param_index_.find("EludeIndexSumSqrDiff");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexSumSqrDiff(pep);
    }   

    fann_type *calc_out;

    calc_out = fann_run(ann_, input);


    double predictedRT = calc_out[0];
	delete[] input;
    
    return predictedRT*max_rt_;

}


double RTCalculator::calc_PepRT(string& pep) {

#ifdef __LGPL__
  return -1.;
#else
  
  int paramcount = paramcount_;
  
  vector<double> x(paramcount, 0.0);
  gsl_matrix *cov = gsl_matrix_calloc(paramcount, paramcount);
  gsl_vector* X =  gsl_vector_calloc(paramcount);
  int totAA = countAA(pep);
  // count numbers of each type of amino acids
  string::size_type pos = 0;
  int numAA = 0;
  double score;
  string aatag, nextaa;
  map<string, int>::iterator found;
  while (pos != string::npos|| numAA == totAA-1) {
    
    if (numAA == 0) {
      aatag = nextAAToken(pep, pos, pos);
    }
    else {
      aatag = nextaa;
    }

    score = 1.;
    if (numAA < totAA-1) {
      nextaa =  nextAAToken(pep, pos, pos);
      //       string aabond  = aatag + "-" + nextaa;
      //       map<string, int>::iterator found1 = param_index_.find(aabond);
      //       if (found1 != param_index_.end()) {
      // 	x[found1->second]++;
      //       }
    }
    found = param_index_.find(aatag);
    if (found != param_index_.end()) {
      x[found->second]+=score;
    }   
    if (aatag.find_first_of("nc") == string::npos && ( numAA < 3 || numAA >= totAA-3) ) {
      aatag = "term_" + aatag;
      if (numAA == 0 || numAA == totAA-1) {
	score = 0.5;
      }
      if (numAA == 1 || numAA == totAA-2)  {
	score = 0.3;
      }
      
      if (numAA == 2 || numAA == totAA-3)  {
	score = 0.1;
      }
      
      found = param_index_.find(aatag);
      if (found != param_index_.end()) {
	x[found->second]+=score;
      }   
	

    }
 

    numAA++;

    
  }

  found = param_index_.find("1");
  if (found != param_index_.end()) {
    x[found->second]=1.;
  }   
  
  found = param_index_.find("ln(LEN)");
  if (found != param_index_.end()) {
    x[found->second]=log((double)numAA);
  }   

  found = param_index_.find("pI");
  if (found != param_index_.end()) {
    //x[found->second]=log(pIcalc_->ModifiedPeptide_pI(pep.c_str(),  0, pep.length() ));  
    x[found->second]=isoelectric(pep);

  }   


  found = param_index_.find("helicity1");
  if (found != param_index_.end()) {
    x[found->second]=helicity1(pep);
  }   

  found = param_index_.find("clusterness");
  if (found != param_index_.end()) {
    x[found->second]=clusterness(pep);
  }   

  found = param_index_.find("helicity2");
  if (found != param_index_.end()) {
    x[found->second]=helicity2(pep);
  }   

  found = param_index_.find("helectric");
  if (found != param_index_.end()) {
    x[found->second]=helectric(pep);
  }   

  found = param_index_.find("hydroMoment");
  if (found != param_index_.end()) {
    x[found->second]=hydroMoment(pep);
  }   
    
  found = param_index_.find("EludeAmphipathicityHelixMin");
  if (found != param_index_.end()) {
    EludeAmphipathicityHelix(pep, &x);
  }   
  
  found = param_index_.find("EludeHydrophobicMoment100Min");
  if (found != param_index_.end()) {
    EludeHydrophobicMoments(pep, &x);
  }   

  found = param_index_.find("EludeIndexPartialSum5Max");
  if (found != param_index_.end()) {
    EludeIndexPartialSum(pep, &x);
  }   
    
  //  found = param_index_.find("smallness");
  //if (found != param_index_.end()) {
  //  x[found->second]=smallness(pep);
  //}   
  
  




  // x[3] = computeHelixScore(pep);

  // x[2] = pIcalc_->ModifiedPeptide_pI(pep.c_str(),  NULL);  

  //x[3] = helicity1(pep);
  //  x[4] = clusterness(pep);
    
  //  x[4] = helicity2(pep);
  //x[5] = helectric(pep);


  for (int j = 0; j < paramcount; j++) {
    gsl_vector_set(X, j, x[j]);    
  }  

  
  double predictedRT;
  double predictedDev;

  int  err = gsl_multifit_linear_est(X, c_, cov, &predictedRT, &predictedDev);
  
  gsl_vector_free(X);
  gsl_matrix_free(cov);
  return predictedRT;
#endif
}

string RTCalculator::getRunName() {
  return *run_name_;
}

string RTCalculator::getRunNameatIdx(int i) {
  return *run_ids_[i];
}

double RTCalculator::getRunRTMean() {
  return run_RT_mean_;
}

double RTCalculator::getRunRTStdDev() {
  return run_RT_stddev_;
}


void RTCalculator::InvertLine() {
  double tmp_slope = 1 / run_slope_;
  double tmp_int = -1 * run_intercept_ / run_slope_;
  
  run_slope_ = tmp_slope;
  run_intercept_ = tmp_int;

}

Boolean RTCalculator::recalc_RTgsl(Array<double>* probs, double min_prob, Array<int>* ntts, int min_ntt) {

  if (ready_) return True;

  if (r_sq_ < 0.5) {
    cerr << "WARNING: Not enough correlation in multifit linear RT model " << *run_name_ << ". RT Model has been disabled." << endl;
    return False;
  }
  ready_ = true;
  return True;

}

Boolean RTCalculator::recalc_RTgsl(double min_prob, int min_ntt) {

  if (ready_) return True;

  if (r_sq_ < 0.5) {
    cerr << "WARNING: Not enough correlation in multifit linear RT model " << *run_name_ << ". RT Model has been disabled." << endl;
    return False;
  }
  ready_ = true;
  return True;

}
Boolean RTCalculator::recalc_RTstats(Array<double>* probs, double min_prob, Array<int>* ntts, int min_ntt) {
  return recalc_RTstats( probs, min_prob, ntts,  min_ntt, RTMAX);
}

Boolean RTCalculator::recalc_RTstats(Array<double>* probs, double min_prob, Array<int>* ntts, int min_ntt, double rtMax) {
  double RTnumer = 0;
  double RTdenom = 0;
  double tot_RTs = 0;
  double all_RTs = 0;
  double SCANnumer = 0;
  double SCANdenom = 0;
  double tot_SCANs = 0;
  double all_SCANs = 0;
  double prob;  
  double RT_mean = 0;
  double RT_stddev = 0;
  double SCAN_mean = 0;
  double SCAN_stddev = 0;
  
  //if (ready_) return True;

  run_RT_used_count_ = 0;
  //   cout << "DDS: probsSize=" << probs->size() << " runRTNum=" << run_RT_count_<< endl;
  assert(probs->size() == run_RT_count_);
  int i;

  int_hash* pep_rt_ix = new int_hash();
  intarr_hash* pep_ixs = new intarr_hash();

  intarr_hash::iterator ia_it ;
  
  string pep;
  for (i=0; i<run_RT_count_; i++) {
    if ((*probs)[i] >= min_prob && (*ntts)[i] >= min_ntt && (*RTs_)[i] > 0 && rts_[i] <= rtMax) {
      pep = modified_peps_[i];
      ia_it =  pep_ixs->find(pep);
      if (ia_it == pep_ixs->end()) {
	pep_ixs->insert(make_pair(pep, new Array<int>()));
      }
      (*pep_ixs)[pep]->insertAtEnd(i);
    }

  }
  
#ifndef __LGPL__
  gsl_vector* r;
  int size = 0;
  int j=0;
  for (ia_it = pep_ixs->begin(); ia_it != pep_ixs->end(); ia_it++) {
    pep = ia_it->first;
    size = ia_it->second->size();
    r = gsl_vector_calloc(size);

    for (j=0; j < size; j++) {
      gsl_vector_set(r, j, rts_[(*(*pep_ixs)[pep])[j]]);
    }
    gsl_sort_vector(r);
    double med = gsl_stats_median_from_sorted_data(gsl_vector_ptr(r, 0), 1, size);

    double min = gsl_stats_min(gsl_vector_ptr(r, 0), 1, size);

    for (j=0; j < size; j++) {
      med = min;
      if (fabs(med - rts_[(*(*pep_ixs)[pep])[j]]) < 0.0001) {
	(*pep_rt_ix)[pep] = (*(*pep_ixs)[pep])[j];
      }
    }
    
    gsl_vector_free(r);

  }
  



#endif 
  int_hash::iterator i_it;

  for (i_it=pep_rt_ix->begin(); i_it!=pep_rt_ix->end(); i_it++) {
    i = i_it->second;
    if ((*probs)[i] >= 0 && (*RTs_)[i] >= 0) {
      prob = (*probs)[i];
      all_RTs += (*RTs_)[i];
      all_SCANs += rts_[i];
      //      all_SCANs += scans_[i];
      if (prob >= min_prob && (*ntts)[i] >= min_ntt) {
	tot_RTs += (*RTs_)[i];
	RTnumer += prob*(*RTs_)[i];
	RTdenom += prob;
	SCANnumer += prob*rts_[i];
	//	SCANnumer += prob*scans_[i];
	SCANdenom += prob;
	run_RT_used_count_ ++;
	if (verbose_)
	  cout << modified_peps_[i] << "\t" << (*probs)[i] << "\t" << rts_[i] << "\t" << (*RTs_)[i] << endl;

      }
    }
  }

  if (RTdenom < 5) {
    cerr << "WARNING: Not enough high probability IDs in run index " << *run_name_ << " to generate RT model. RT Model has been disabled." << endl;
    //return False;
  }
  else {
    RT_mean = RTnumer / RTdenom;
    RT_stddev = 0;
    SCAN_mean = SCANnumer / SCANdenom;
    SCAN_stddev = 0;
    RTdenom = 0;
    SCANdenom = 0;
    for (i_it=pep_rt_ix->begin(); i_it!=pep_rt_ix->end(); i_it++) {
      i = i_it->second;
      if ((*probs)[i] >= 0 && (*RTs_)[i] >= 0) {
	prob = (*probs)[i];
	if (prob >= min_prob && (*ntts)[i] >= min_ntt ) {
	  RT_stddev += prob * pow((RT_mean - (*RTs_)[i]), 2);
	  SCAN_stddev += prob * pow((SCAN_mean - rts_[i]), 2);
	  //	  SCAN_stddev += prob * pow((SCAN_mean - scans_[i]), 2);
	  RTdenom += prob;
	  SCANdenom += prob;
	}
      }
    }
  }

  RT_stddev /= RTdenom;

  RT_stddev = pow(RT_stddev, 0.5); 
  
  run_RT_mean_ = RT_mean;
  run_RT_stddev_ = RT_stddev;

  SCAN_stddev /= SCANdenom;

  SCAN_stddev = pow(SCAN_stddev, 0.5); 
  
  run_SCAN_mean_ = SCAN_mean;
  run_SCAN_stddev_ = SCAN_stddev;

  double SSxx = 0;
  double SSxy = 0;
  double SSyy = 0;

  
  
  for (i_it=pep_rt_ix->begin(); i_it!=pep_rt_ix->end(); i_it++) {
      i = i_it->second;
    if ((*probs)[i] >= 0 && (*RTs_)[i] >= 0) {
      if ((*probs)[i] >= min_prob && (*ntts)[i] >= min_ntt) {
	SSxx += pow((rts_[i] - run_SCAN_mean_), 2);
	//	SSyy += pow((scans_[i] - run_SCAN_mean_), 2);
	SSyy += pow(((*RTs_)[i] - run_RT_mean_), 2);
	SSxy += (rts_[i] - run_SCAN_mean_)*((*RTs_)[i] - run_RT_mean_);
      }
    }
  }  


  delete pep_rt_ix;
  delete pep_ixs;

  run_slope_ = SSxy / SSxx;
  run_intercept_ = run_SCAN_mean_ - run_slope_ * run_RT_mean_;
  r_sq_ = (SSxy*SSxy) / (SSxx*SSyy);
  if (r_sq_ < 0.5) {
    cerr << "WARNING: Not enough correlation r_sq=" << r_sq_ << " between theoretical RT values and scan numbers in run index " << *run_name_ << ". RT Model has been disabled." << endl;
    //    return False;
  }

  cerr << "Run Index: " << *run_name_ << ", slope=" << run_slope_ << ", intercept=" << run_intercept_ << ", r_sq=" << r_sq_ << "\n";
  ready_ = true;
  return True;
}




//DDS: RT model
void RTCalculator::calc_RTstats() {
  double numer = 0;
  double denom = 0;
  run_RT_used_count_ = run_RT_count_;
  int i;
  for (i=0; i<run_RT_count_; i++) {
    numer += (*RTs_)[i];
    denom ++;
  }
  
  run_RT_mean_ = numer/denom;
  run_RT_stddev_ = 0;
 
  for (i=0; i<run_RT_count_; i++) {
    run_RT_stddev_ += pow((run_RT_mean_ - (*RTs_)[i]), 2);
  }

  run_RT_stddev_ /= denom;
  
  run_RT_stddev_ = pow(run_RT_stddev_, 0.5);  
}

bool RTCalculator::linearRegressRT(Array<double>* probs, double min_prob, Array<int>* ntts, int min_ntt) {

  //  min_prob = 0;
#ifdef __LGPL__
  // cannot do fitting, so just calculate "predictedRT" by SSRCalc 1.0 algorithm
  // (published in MCP)
  return (false);



#else
  
  cerr << "Performing linear regression of retention times..." << endl;
  
  // The sequence vs RT relationship is learned by linear regression, rather
  // than predicted by SSRCalc or Hydro. This is more robust and generally applicable.
  // Note that RT predictors like SSRCalc or Hydro makes assumptions about the column used,
  // components and gradient of mobile phase, etc, and cannot deal with modified amino acids.
  
  int numDataPoints = (int)rts_.size();

  int numGood = 0;
  
  for (int d = 0; d < numDataPoints; d++) {
    if ((*probs)[d] >= min_prob && (*ntts)[d] >= min_ntt) {
      numGood++;
      if (min_RT_ <= 0 || rts_[d] < min_RT_) {
	min_RT_ = rts_[d];
      }

    }
  }

  adjustRT();

  cerr << "Fitting model to " << numGood << " good identifications with P>=" << min_prob << " and NTT>=" << min_ntt << "." << endl;
  
  if (numGood < 80) {
    cerr << "Not enough good identifications to fit retention time model." << endl;
    cerr << "Reverting to SSRCalc 1.0 predictions." << endl;
    return (false);
  }
    
  //  map<string, int> param_index;
  vector<string> params;
  int paramcount = 0;

  // add the constant
  param_index_["1"] = paramcount++;
  params.push_back("1");
  // add the log of the peptide length as one of the modeled variables
  param_index_["ln(LEN)"] = paramcount++;
  params.push_back("ln(LEN)");
  string pep;
  string aatag;
  // find out what amino acids and modified amino acids are present in the identifications
  for (int dd = 0; dd < numDataPoints; dd++) {
     pep =  modified_peps_[dd];
    string::size_type pos = 0;
    int numAA = 0;
    while (pos != string::npos) {
       aatag = nextAAToken(pep, pos, pos);
      if (param_index_.find(aatag) == param_index_.end()) {
	param_index_[aatag] = paramcount;
	paramcount++;
	params.push_back(aatag);
      }
      if (numAA < 3) {
	// create another variable to model the disproportionate influence of 
	// amino acids close to the N-terminus
	aatag = "n_" + aatag;
	if (param_index_.find(aatag) == param_index_.end()) {
	  param_index_[aatag] = paramcount;
	  paramcount++;
	  params.push_back(aatag);
	}	
      }
      numAA++;
    }
  }
   
  
  if (!matrix_) {
    
    matrix_ = new vector<vector<double> >; 

    // fill in the matrix -- this will be the input to the linear regression algorithm
    for (int d = 0; d <  (int)modified_peps_.size(); d++) {
    
      vector<double> x(paramcount + 1, 0.0);
      
      string pep =  modified_peps_[d];
      
      // count numbers of each type of amino acids
      string::size_type pos = 0;
      int numAA = 0;
      while (pos != string::npos) {
	string aatag = nextAAToken(pep, pos, pos);
        map<string, int>::iterator found = param_index_.find(aatag);
        if (found != param_index_.end()) {
	  x[found->second]++;
        }       
	
	if (numAA < 3) {
	  aatag = "n_" + aatag;
	  map<string, int>::iterator found1 = param_index_.find(aatag);
	  if (found1 != param_index_.end()) {
	    if (numAA == 0) {
	      x[found1->second] += 0.5; // approx weighting based on SSRCalc 1
	    } else if (numAA == 1) {
	      x[found1->second] += 0.3;
	    } else if (numAA == 2) {
	      x[found1->second] += 0.1;
	      
	    }
	  }  
	}
		
	numAA++;
	
      }
        
      x[0] = 1.0;
      x[1] = log((double)numAA);
      //      for (int z = 2; z < paramcount; z++) {
      //	if (params[z][0] >= 'A' && params[z][0] <= 'Z') {
      //	  x[z] /= ((double)(numAA));
      //	}
      // }
      // the very last column is the retention time itself (can also use scan number)
      x[paramcount] = rts_[d];
    
      matrix_->push_back(x);
      
    }
  }




  // using GNU Scientific Library (GSL) to do fitting.
  // note that GSL is GPL, so this feature must be disabled if LGPL is desired

  
  if (rt_fit_ != NULL) gsl_multifit_linear_free(rt_fit_);
  rt_fit_ = gsl_multifit_linear_alloc(numGood, paramcount); 
  
  gsl_vector* y = gsl_vector_calloc(numGood);
  gsl_vector_set_zero(y);
  gsl_vector* w = gsl_vector_calloc(numGood);
  gsl_vector_set_zero(w);
  gsl_matrix* X = gsl_matrix_calloc(numGood, paramcount); 
  gsl_matrix_set_zero(X);
  int count = 0;
  for(int d = 0; d < (int)matrix_->size(); d++){
    if ((*probs)[d] >= min_prob && (*ntts)[d] >= min_ntt) {
      gsl_vector_set(y, count, (*matrix_)[d][paramcount]);
      gsl_vector_set(w, count, (*probs)[d]);
      for(int j = 0; j < paramcount; j++) {
	//	if ( j == 27  &&  (*matrix_)[d][j] > 0 ) {
	// cerr << "DDS: DBG" << endl;
	//}
	  gsl_matrix_set(X, count, j, (*matrix_)[d][j]);
      }
      count++;
    }
  }
  
  /* DEBUG code to print out matrix

  cerr << endl << "==============" << endl;
  for (int row = 0; row < count; row++) {
    cerr << gsl_vector_get(y, row);
      for (int col = 0; col < 41; col++) {
      cerr << gsl_matrix_get(X, row, col) << ' ';
    }
    cerr << endl;
  }
  */
  
  gsl_vector *c = gsl_vector_calloc(paramcount);
  gsl_vector_set_zero(c);
  gsl_matrix *cov = gsl_matrix_calloc(paramcount, paramcount);
  gsl_matrix_set_zero(cov);
  double chisq;
  double tolerance = 0.01;
  size_t rank = 0;

  int err = gsl_multifit_wlinear_svd(X, w, y, tolerance, &rank, c, cov, &chisq, rt_fit_);
  int newNumGood = 0;
  
  //DDS: Outlier removal here
  
  gsl_vector *dis =  gsl_vector_calloc(numGood);
  gsl_vector_set_zero(dis);
  int ix = 0;
  double dis_mean, dis_sd;

  
  //DDS: Find distance mean and stddev

  for (int d = 0; d < (int)matrix_->size(); d++) {

    if ((*probs)[d] >= min_prob && (*ntts)[d] >= min_ntt) {
      double predictedRT;
      double predictedDev;
      gsl_vector *x = gsl_vector_calloc(paramcount);
      gsl_vector_set_zero(x);
      for (int j = 0; j < paramcount; j++) {
	gsl_vector_set(x, j, (*matrix_)[d][j]);
      }
      err = gsl_multifit_linear_est(x, c, cov, &predictedRT, &predictedDev);

      gsl_vector_set(dis, ix++, (rts_)[d] - predictedRT);

      gsl_vector_free(x);
      
    }    
  }
  
  dis_mean = gsl_stats_wmean(gsl_vector_ptr(w, 0), 1,gsl_vector_ptr(dis, 0), 1, numGood);
  run_RT_mean_ = dis_mean;
  dis_sd = gsl_stats_wsd(gsl_vector_ptr(w, 0), 1,gsl_vector_ptr(dis, 0), 1, numGood);
  run_RT_stddev_ = dis_sd;


  count = 0;
  
  for (int d = 0; d < (int)matrix_->size(); d++) {

    if ((*probs)[d] >= min_prob && (*ntts)[d] >= min_ntt) {
      double predictedRT;
      double predictedDev;
      gsl_vector *x = gsl_vector_calloc(paramcount);
      
      for (int j = 0; j < paramcount; j++) {
	gsl_vector_set(x, j, (*matrix_)[d][j]);
      }
      err = gsl_multifit_linear_est(x, c, cov, &predictedRT, &predictedDev);
      
      if (fabs(((rts_)[d] - predictedRT - dis_mean)/dis_sd) < 2) {
	count++;	
      }
      gsl_vector_free(x);
    }    
  }
  
  gsl_vector_free(y);
  gsl_vector_free(w);
  gsl_matrix_free(X);
  
  numGood = count;
  run_RT_used_count_ = count;
  count = 0;
  y = gsl_vector_calloc(numGood);
  gsl_vector_set_zero(y);
  w = gsl_vector_calloc(numGood);
  gsl_vector_set_zero(y);
  X = gsl_matrix_calloc(numGood, paramcount); 
  gsl_vector_set_zero(y);
  
  
  for (int d = 0; d < (int)matrix_->size(); d++) {
    
    if ((*probs)[d] >= min_prob && (*ntts)[d] >= min_ntt) {
      double predictedRT;
      double predictedDev;
      gsl_vector *x = gsl_vector_calloc(paramcount);
      
      for (int j = 0; j < paramcount; j++) {
	gsl_vector_set(x, j, (*matrix_)[d][j]);
      }
      err = gsl_multifit_linear_est(x, c, cov, &predictedRT, &predictedDev);
      
      if (fabs(((rts_)[d] - predictedRT - dis_mean)/dis_sd) < 2) {
	gsl_vector_set(y, count, (*matrix_)[d][paramcount]);
	gsl_vector_set(w, count, (*probs)[d]);
	for (int j = 0; j < paramcount; j++) {
	  //	  if ( j == 27  &&  (*matrix_)[d][j] > 0 ) {
	  // cerr << "DDS: DBG" << endl;
	  //}
	  gsl_matrix_set(X, count, j, (*matrix_)[d][j]);
	}

	count++;
      }
      
      gsl_vector_free(x);
    }    
  }

  gsl_vector_free(dis);

  cerr << "After outlier removal, fitting model to " << numGood << " good identifications with P>=" << min_prob << " and NTT>=" << min_ntt << "." << endl;
  
  if (numGood < 80) {
    cerr << "Not enough good identifications to fit retention time model." << endl;
    cerr << "Reverting to SSRCalc 1.0 predictions." << endl;
    return (false);
  }
  
  //  cerr << "RANK=" << rank << endl;
  //  cerr << "DONE FITTING" << endl;
  if (rt_fit_ != NULL)  gsl_multifit_linear_free(rt_fit_);
  rt_fit_ = gsl_multifit_linear_alloc(numGood, paramcount); 
  //  err = gsl_multifit_wlinear_svd(X, w, y,  tolerance, &rank, c, cov, &chisq, rt_fit_);

  err = gsl_multifit_linear_svd(X, y,  tolerance, &rank, c, cov, &chisq, rt_fit_);

 for (map<string, int>::iterator pp = param_index_.begin(); pp != param_index_.end(); pp++) {
    coeff_[pp->first] = gsl_vector_get(c, pp->second);
  }
  
  double tss = (count - 1) * gsl_stats_wvariance(gsl_vector_ptr(w, 0), 1, gsl_vector_ptr(y, 0), 1, count);
  r_sq_ = 1 - (chisq / tss);
    
  RTs_->reserve(numDataPoints);
  // now calculate the "predicted" RTs from this regressed model
  for (int d = 0; d < numDataPoints; d++) {
    double predictedRT;
    double predictedDev;
    gsl_vector *x = gsl_vector_calloc(paramcount);
    for (int j = 0; j < paramcount; j++) {
      gsl_vector_set(x, j, (*matrix_)[d][j]);
    }
      
    err = gsl_multifit_linear_est(x, c, cov, &predictedRT, &predictedDev);
  
    run_RT_count_++;
    run_RT_sum_ += predictedRT;
    
    RTs_->insertAtEnd(predictedRT);
  
    gsl_vector_free(x);
  }
    
  gsl_vector_free(c);
  gsl_matrix_free(cov);
  gsl_vector_free(y);
  gsl_matrix_free(X);
// gsl_multifit_linear_free(work);
  
//  plotRegressionFits(probs, min_prob, ntts, min_ntt);
  
  if (r_sq_ < 0.6) {
    // too much scatter to trust this?
    return (false);
  }
#endif  


  return (true);
}

double RTCalculator::getRTScore(double RT,  double rt) {
  return getRTScore(RT, 0, rt);

}

double RTCalculator::getRTScore(double RT, int scan, double rt) {
  //TODO: DDS hack for now should divide by average stdev over all runs
  //if (run_RT_stddev_ <= 0.05) run_RT_stddev_ = 0.5;
  
 // double rtn = (rt - RT)/1000.0;

  //#ifdef __LGPL__
  // double rtn = (rt - (run_slope_ * RT + run_intercept_));
  //#else
    
  if (isnan(run_RT_stddev_) || run_RT_stddev_ <= 0) {
    return 0;
  }

  double rtn = (rt - RT) / run_RT_stddev_;


  //  double rtn = (scan - (run_slope_ * RT + run_intercept_))/1000;
  //  cout << "peptide seq: " << pep << " RT: " << RT << " Score: " << rtn <<endl;

  //   if (rtn > 6000) {
  // rtn = 6000;
  // }
   
  // if (rtn < -3000) {
  //   rtn = -3000;
  // }

  if (rtn < -2) {
    rtn = -2;
  }
  if (rtn > 6) {
    rtn = 6;
  }
  //#endif

  return rtn;
  
}


double RTCalculator::getRTScore(int index) {
  //TODO: DDS hack for now should divide by average stdev over all runs
  //if (run_RT_stddev_ <= 0.05) run_RT_stddev_ = 0.5;
  
 // double rtn = (rt - RT)/1000.0;

  //#ifdef __LGPL__
  // double rtn = (rt - (run_slope_ * RT + run_intercept_));
  //#else

  double rtn;
  string pep = modified_peps_[index];
  double rt = rts_[index];
  double RTmn = rt_catalog_->getRTMean(pep);
  double RTsd = rt_catalog_->getRTStdev(pep);
  int RTnum =  rt_catalog_->getRTCount(pep);

  

 
  if (RTnum < 2 || RTsd <= 0) {
    // TODO: ???
    // is this the right value 
    // for the peptides that aren't in the catalog
    rtn = 5;
     
 }
  else {
    rtn = (rt - RTmn) / RTsd;
  }

  rtn = fabs(rtn);

  //  double rtn = (scan - (run_slope_ * RT + run_intercept_))/1000;
  //  cout << "peptide seq: " << pep << " RT: " << RT << " Score: " << rtn <<endl;

  //   if (rtn > 6000) {
  // rtn = 6000;
  // }
   
  // if (rtn < -3000) {
  //   rtn = -3000;
  // }

  // if (rtn < -5) {
  //  rtn = -5;
  //}

  if (rtn > 5) {
    rtn = 5;
  }
  //#endif

  return rtn;
  
}


double RTCalculator::EXP10(double value)
{
   return( pow(10.0,value) );
}

void RTCalculator::plotRegressionFits(Array<double>* probs, double min_prob, Array<int>* ntts, int min_ntt) {
  
  stringstream fileNameSS;
  fileNameSS << "RTFit_" << *run_name_;
  string specFileName(fileNameSS.str() + ".spec");
  string gpFileName(fileNameSS.str() + ".gp");
  
  ofstream fout(specFileName.c_str());
  ofstream gpout(gpFileName.c_str());
  
  double xmax = -1000000000000.0;
  double xmin = 1000000000000.0;
  double ymax = -1000000000000.0;
  double ymin = 1000000000000.0;
  for (int d = 0; d < probs->size(); d++) {
    
    if ((*RTs_)[d] < xmin) xmin = (*RTs_)[d];
    if ((*RTs_)[d] > xmax) xmax = (*RTs_)[d];
    if (rts_[d] < ymin) ymin = rts_[d];
    if (rts_[d] > ymax) ymax = rts_[d];
    
  if ((*probs)[d] >= min_prob && (*ntts)[d] >= min_ntt) {
    fout << scans_[d] << '\t' << (*RTs_)[d] << '\t' << rts_[d] << '\t' << -1.0 << '\t' << -1.0 << endl;
  } else if ((*probs)[d] >= min_prob) {
    fout << scans_[d] << '\t' << (*RTs_)[d] << '\t' << -1.0 << '\t' << rts_[d] << '\t' << -1.0 << endl;
  } else {
    fout << scans_[d] << '\t' << (*RTs_)[d] << '\t' << -1.0 << '\t' << -1.0 << '\t' << rts_[d] << endl;
  }
    
  }
  gpout << "set terminal png" << endl;
  gpout << "set output \"" << fileNameSS.str() + ".png" << "\"" << endl;
  gpout << "set nokey" << endl;
  gpout << "set border 1" << endl;
  gpout << "set xtics border nomirror" << endl;
  gpout << "set ytics border nomirror" << endl;
  gpout << "set xzeroaxis linetype -1 linewidth 1.0" << endl;
  gpout << "set yzeroaxis linetype -1 linewidth 1.0" << endl;
  gpout << "set mxtics" << endl;
  gpout << "set mytics" << endl;
  
  gpout << "set origin 0.0,0.0" << endl;
  gpout << "set yrange [0.0:" << ymax << "]" << endl;
  gpout << "set xlabel \"" << "Predicted RT" << endl;
  gpout << "set ylabel \"" << "Actual RT" << endl;
  
  gpout.precision(3);
  gpout << "set label \"" << "R2=" << r_sq_ << " (" << r_sq_ << ")\" at " << xmin << ", " << ymax << " front" << endl; 
  
  
  if (!coeff_.empty()) {
  double move = (ymax - ymin) / 42;
  double shift = move;
  for (map<string, double>::iterator q = coeff_.begin(); q != coeff_.end(); q++) {
    gpout.precision(1);
    gpout << "set label \"" << q->first << ":" << fixed << q->second << "\" at ";
    gpout.precision(3);
    gpout << xmin << ", " << ymax - shift << " front" << endl;
    shift += move;  
  }
  }
  
  gpout << "plot \"" << specFileName << "\" using 2:5 with dots lc -1";
  gpout << ", \"" << specFileName << "\" using 2:4 with points pt 1 lc 4";
  gpout << ", \"" << specFileName << "\" using 2:3 with points pt 1 lc 1";
  
  if (coeff_.empty()) {
    gpout << ", " << run_intercept_ << " + " <<  run_slope_  << " * x with line lc 1" << endl;
  } else {
    gpout << ", x with line lc 1" << endl;
  }
    
  string gnuCmd("gnuplot " + gpFileName);
  int return_result = tpplib_system(gnuCmd.c_str());
  
  string rmCmd("rm " + specFileName + " " + gpFileName);
  
}

int RTCalculator::countAA(string& pep) {
   int totAA = 0;
   string::size_type pos = 0;     
   while (pos != string::npos) {
     string aatag = nextAAToken(pep, pos, pos);
     totAA++;
   } 
   
  return totAA;

}




// nextAAToken - simply extracts the next AA token, starting to look from position 'start' in string s. 
// An AA token starts at an AA and ends before the next AA. That is, C*, CJ, C[339], M(O), KJJJJJJ are 
// all AA tokens. When function returns, 'end' will point to the next AA or string::npos if there's nothing left
string RTCalculator::nextAAToken(string& s, string::size_type start, string::size_type& end) {

  string allAA("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
  allAA += "nc"; // add the n- and c-termini
  
  if (allAA.find(s[start], 0) == string::npos) {
    start = s.find_first_of(allAA, start);
    
    if (start == string::npos) {
      // no more AA
      end = string::npos;
      return ("");
    }
  } 
  
  
  end = s.find_first_of(allAA + "[", start + 1);
  
  if(end == string::npos) {
    
    return (s.substr(start));

  } else if (s[end] == '[') {

    string::size_type closeBracket = s.find(']', end + 1);
    if (closeBracket == string::npos) {
      // [ not closed, just stop at first AA after the [
      end = s.find_first_of(allAA, end + 1); 
    } else {
      end = closeBracket + 1;
    }
  }
  
  if (end == string::npos || end >= (int)s.length()) {
    end = string::npos;
    return (s.substr(start));
  }

  return (s.substr(start, end - start));
}

string RTCalculator::getAATokenIndex(string& s, int i) {

  string::size_type pos = 0;
  int c = 0;
  string rtn="";
  
  while (pos != string::npos && c++ <= i) {
    rtn = nextAAToken(s, pos, pos);
  }
  
  return rtn;
  
}


int  RTCalculator::getPeptideLength(string& s) {

  string::size_type pos = 0;
  int rtn=-1;
  
  while (pos != string::npos) {
    nextAAToken(s, pos, pos);
    rtn++;
  }
  
  return rtn;

}

string  RTCalculator::stripPeptide(string& s) {

  string::size_type pos = 0;
  string rtn="";
  int len = getPeptideLength(s);
  
  for (int i=0; i<len; i++) {
    rtn += getAATokenIndex(s,i).substr(0,1);

  }
  
  return rtn;

}

string  RTCalculator::reversePeptide(string& s) {

  string::size_type pos = 0;
  string rtn="";
  int len = getPeptideLength(s);
  
  for (int i=len-1; i>-1; i--) {
    rtn += getAATokenIndex(s,i);

  }

  
  return rtn;

}


double RTCalculator::calcSSR1(string& pep) {
 
  string::size_type pos = 0;
  int numAA = 0;
  double H = 0.0;

  while (pos != string::npos) {
    string aatag = nextAAToken(pep, pos, pos);
    
    if (SSR1Coeff_->find(aatag) != SSR1Coeff_->end()) {
      H += (*SSR1Coeff_)[aatag];
    }

    if (numAA < 3) {
      if (SSR1NTermCoeff_->find(aatag) != SSR1NTermCoeff_->end()) {
	double adjustment = (*SSR1NTermCoeff_)[aatag];
	if (numAA == 0) {
	  H += 0.42 * adjustment;
	} else if (numAA == 1) {
	  H += 0.22 * adjustment;
	} else {
	  H += 0.05 * adjustment;
	}
      }
    }

    numAA++;
  }

  double KL = 1.0;
  if (numAA < 10) {
    KL = 1 - 0.027 * (10 - numAA);
  } else if (numAA > 20) {
    KL = 1 - 0.014 * (numAA - 20);
  }

  if (H >= 38.0) {
    H = H - 0.3 * (H - 38.0);
  }

  return (H);

}

void RTCalculator::initSSR1Coeff() {

  if (SSR1Coeff_ && SSR1NTermCoeff_) {
    return;
  }

  SSR1Coeff_ = new map<string, double>;
  SSR1NTermCoeff_ = new map<string, double>;

  (*SSR1Coeff_)["W"] = 11.0;
  (*SSR1Coeff_)["F"] = 10.5;
  (*SSR1Coeff_)["L"] = 9.6;
  (*SSR1Coeff_)["I"] = 8.4;

  (*SSR1Coeff_)["M"] = 5.8;
  (*SSR1Coeff_)["M[147]"] = 5.8; // assume M[147] behaves as M

  (*SSR1Coeff_)["V"] = 5.0;
  (*SSR1Coeff_)["Y"] = 4.0;
  (*SSR1Coeff_)["A"] = 0.8;
  (*SSR1Coeff_)["T"] = 0.4;
  (*SSR1Coeff_)["P"] = 0.2;
  (*SSR1Coeff_)["E"] = 0.0;
  (*SSR1Coeff_)["D"] = -0.5;

  (*SSR1Coeff_)["C"] = -0.8; // assume C behaves as C[161]
  (*SSR1Coeff_)["C[160]"] = -0.8; // assume C[160] behaves as C[161]
  (*SSR1Coeff_)["C[161]"] = -0.8; // original SSRCalc actually uses C[161]

  (*SSR1Coeff_)["S"] = -0.8;
  (*SSR1Coeff_)["Q"] = -0.9;
  (*SSR1Coeff_)["G"] = -0.9;
  (*SSR1Coeff_)["N"] = -1.2;
  (*SSR1Coeff_)["R"] = -1.3;
  (*SSR1Coeff_)["H"] = -1.3;
  (*SSR1Coeff_)["K"] = -1.9;

  ////

  (*SSR1NTermCoeff_)["W"] = -4.0;
  (*SSR1NTermCoeff_)["F"] = -7.0;
  (*SSR1NTermCoeff_)["L"] = -9.0;
  (*SSR1NTermCoeff_)["I"] = -8.0;

  (*SSR1NTermCoeff_)["M"] = -5.5;
  (*SSR1NTermCoeff_)["M[147]"] = -5.5; // assume M[147] behaves as M

  (*SSR1NTermCoeff_)["V"] = -5.5;
  (*SSR1NTermCoeff_)["Y"] = -3.0;
  (*SSR1NTermCoeff_)["A"] = -1.5;
  (*SSR1NTermCoeff_)["T"] = 5.0;
  (*SSR1NTermCoeff_)["P"] = 4.0;
  (*SSR1NTermCoeff_)["E"] = 7.0;
  (*SSR1NTermCoeff_)["D"] = 9.0;

  (*SSR1NTermCoeff_)["C"] = 4.0; // assume C behaves as C[161]
  (*SSR1NTermCoeff_)["C[160]"] = 4.0; // assume C[160] behaves as C[161]
  (*SSR1NTermCoeff_)["C[161]"] = 4.0; // original SSRCalc actually uses C[161]

  (*SSR1NTermCoeff_)["S"] = 5.0;
  (*SSR1NTermCoeff_)["Q"] = 1.0;
  (*SSR1NTermCoeff_)["G"] = 5.0;
  (*SSR1NTermCoeff_)["N"] = 5.0;
  (*SSR1NTermCoeff_)["R"] = 8.0;
  (*SSR1NTermCoeff_)["H"] = 4.0;
  (*SSR1NTermCoeff_)["K"] = 4.6;

}

void RTCalculator::initHelixProp() {

  if (HelixProp_) {
    return;
  }

  HelixProp_ = new map<string, double>;


  (*HelixProp_)["W"] = -0.49;
  (*HelixProp_)["F"] = -0.54;
  (*HelixProp_)["L"] = -0.21;
  (*HelixProp_)["I"] = -0.41;

  (*HelixProp_)["M"] = -0.24;
  (*HelixProp_)["M[147]"] = -0.24; // assume M[147] behaves as M

  (*HelixProp_)["V"] = -0.61;
  (*HelixProp_)["Y"] = -0.53;
  (*HelixProp_)["A"] = 0;
  (*HelixProp_)["T"] = -0.66;
  (*HelixProp_)["P"] = -3.16;
  (*HelixProp_)["E"] = -0.40;
  (*HelixProp_)["D"] = -0.69;

  (*HelixProp_)["C"] = -0.68; // assume C behaves as C[161]
  (*HelixProp_)["C[160]"] = -0.68; // assume C[160] behaves as C[161]
  (*HelixProp_)["C[161]"] = -0.68; // original SSRCalc actually uses C[161]

  (*HelixProp_)["S"] = -0.5;
  (*HelixProp_)["Q"] = -0.39;
  (*HelixProp_)["G"] = -1;
  (*HelixProp_)["N"] = -0.65;
  (*HelixProp_)["R"] = -0.21;
  (*HelixProp_)["H"] = -0.56;
  (*HelixProp_)["K"] = -0.26;


  (*HelixProp_)["W_bsc"] = 1.6;
  (*HelixProp_)["F_bsc"] = 0.5;
  (*HelixProp_)["L_bsc"] = 1.6;
  (*HelixProp_)["I_bsc"] = 3.5;

  (*HelixProp_)["M_bsc"] = 1.8;
  (*HelixProp_)["M[147]_bsc"] = 1.8; // assume M[147] behaves as M

  (*HelixProp_)["V_bsc"] = 1.4;
  (*HelixProp_)["Y_bsc"] = 0.2;

  (*HelixProp_)["C_bsc"] = 0; // assume C behaves as C[161]
  (*HelixProp_)["C[160]_bsc"] = 0; // assume C[160] behaves as C[161]
  (*HelixProp_)["C[161]_bsc"] = 0; // original SSRCalc actually uses C[161]
  
  (*HelixProp_)["P_bsc"] = 0;

  (*HelixProp_)["A_bsc"] = 1;
  (*HelixProp_)["E_bsc"] = 0;
  (*HelixProp_)["Z_bsc"] = 0;
  (*HelixProp_)["T_bsc"] = 0;
  (*HelixProp_)["B_bsc"] = 0;
  (*HelixProp_)["D_bsc"] = 0;

  (*HelixProp_)["Q_bsc"] = 0;
  (*HelixProp_)["S_bsc"] = 0;
  (*HelixProp_)["G_bsc"] = 0;
  (*HelixProp_)["R_bsc"] = 0;
  (*HelixProp_)["N_bsc"] = 0;
  (*HelixProp_)["H_bsc"] = 0;
  (*HelixProp_)["K_bsc"] = 0;
  (*HelixProp_)["X_bsc"] = 0;


  (*HelixProp_)["W_cmu"] = 1.;
  (*HelixProp_)["F_cmu"] = 1.;
  (*HelixProp_)["L_cmu"] = 1.6;
  (*HelixProp_)["I_cmu"] = 1.4;

  (*HelixProp_)["M_cmu"] = 1.;
  (*HelixProp_)["M[147]_cmu"] = 1.; // assume M[147] behaves as M

  (*HelixProp_)["V_cmu"] = 1.2;
  (*HelixProp_)["Y_cmu"] = 1.;

  (*HelixProp_)["C_cmu"] = 1.; // assume C behaves as C[161]
  (*HelixProp_)["C[160]_cmu"] = 1.; // assume C[160] behaves as C[161]
  (*HelixProp_)["C[161]_cmu"] = 1.; // original SSRCalc actually uses C[161]
  
  (*HelixProp_)["P_cmu"] = 0.3;

  (*HelixProp_)["A_cmu"] = 1.2;
  (*HelixProp_)["E_cmu"] = 1.1;
  (*HelixProp_)["Z_cmu"] = 1.1;
  (*HelixProp_)["T_cmu"] = 1.;
  (*HelixProp_)["B_cmu"] = 1.1;
  (*HelixProp_)["D_cmu"] = 1.1;

 

  (*HelixProp_)["Q_cmu"] = 1.;
  (*HelixProp_)["S_cmu"] = 1.;
  (*HelixProp_)["G_cmu"] = 0.3;
  (*HelixProp_)["R_cmu"] = 1.;
  (*HelixProp_)["N_cmu"] = 0.4;
  (*HelixProp_)["H_cmu"] = 0.6;
  (*HelixProp_)["K_cmu"] = 1.;
  (*HelixProp_)["X_cmu"] = 1.;

(*HelixProp_)["XXUX"] =  0.8;
(*HelixProp_)["XZOX"] =  0.8;
(*HelixProp_)["XUXX"] =  0.8;

(*HelixProp_)["XXOX"] =  0.7;
(*HelixProp_)["XOXX"] =  0.7;
(*HelixProp_)["XZUX"] =  0.7;
(*HelixProp_)["XXOZ"] =  0.7;
(*HelixProp_)["ZXOX"] =  0.7;
(*HelixProp_)["XOZZ"] =  0.7;
(*HelixProp_)["ZOXX"] =  0.7;
(*HelixProp_)["ZOZX"] =  0.7;
(*HelixProp_)["ZUXX"] =  0.7;

(*HelixProp_)["ZXUX"] =  0.5;
(*HelixProp_)["XOZX"] =  0.5;
(*HelixProp_)["XZOZ"] =  0.5;
(*HelixProp_)["XUZX"] =  0.5;

(*HelixProp_)["ZZOX"] =  0.2;
(*HelixProp_)["ZXOZ"] =  0.2;
(*HelixProp_)["ZOXZ"] =  0.2;
(*HelixProp_)["XOXZ"] =  0.2;

(*HelixProp_)["ZZUZ"] =  0.2;
(*HelixProp_)["XUXZ"] =  0.2;
(*HelixProp_)["ZUXZ"] =  0.2;
(*HelixProp_)["XZUZ"] =  0.2;
(*HelixProp_)["XUZZ"] =  0.2;
(*HelixProp_)["ZXUZ"] =  0.2;
(*HelixProp_)["ZOZZ"] =  0.2;
(*HelixProp_)["ZZOZ"] =  0.2;
(*HelixProp_)["ZZUX"] =  0.2;
(*HelixProp_)["ZUZX"] =  0.2;
(*HelixProp_)["XXUZ"] =  0.2;
(*HelixProp_)["ZUZZ"] =  0.2;

(*HelixProp_)["XXOXX"] =  3.75;
(*HelixProp_)["XXOXZ"] =  3.75;
(*HelixProp_)["XXOZX"] =  3.75;
(*HelixProp_)["XZOXX"] =  3.75;
(*HelixProp_)["ZXOXX"] =  3.75;
(*HelixProp_)["XXOZZ"] =  2.7;
(*HelixProp_)["XZOXZ"] =  2.7;
(*HelixProp_)["XZOZX"] =  2.7;
(*HelixProp_)["ZXOXZ"] =  2.7;
(*HelixProp_)["ZXOZX"] =  2.7;
(*HelixProp_)["ZZOXX"] =  2.7;
(*HelixProp_)["ZXOZZ"] =  1.3;
(*HelixProp_)["XZOZZ"] =  1.3;
(*HelixProp_)["ZZOXZ"] =  1.3;
(*HelixProp_)["ZZOZX"] =  1.3;
(*HelixProp_)["ZZOZZ"] =  1.3;
(*HelixProp_)["XXUXX"] =  3.75;
(*HelixProp_)["XXUXZ"] =  3.75;
(*HelixProp_)["XXUZX"] =  3.75;
(*HelixProp_)["XZUXX"] =  3.75;
(*HelixProp_)["ZXUXX"] =  3.75;
(*HelixProp_)["XXUZZ"] =  1.1;
(*HelixProp_)["XZUXZ"] =  1.1;
(*HelixProp_)["XZUZX"] =  1.1;
(*HelixProp_)["ZXUZX"] =  1.1;
(*HelixProp_)["ZXUXZ"] =  1.1;
(*HelixProp_)["ZZUXX"] =  1.1;

(*HelixProp_)["XZUZZ"] =  1.3;
(*HelixProp_)["ZXUZZ"] =  1.3;
(*HelixProp_)["ZZUXZ"] =  1.3;
(*HelixProp_)["ZZUZX"] =  1.3;
(*HelixProp_)["ZZUZZ"] =  1.3;

(*HelixProp_)["XXOOX"] =  1.25;
(*HelixProp_)["ZXOOX"] =  1.25;
(*HelixProp_)["XZOOX"] =  1.25;
(*HelixProp_)["XOOXX"] =  1.25;
(*HelixProp_)["XOOXZ"] =  1.25;
(*HelixProp_)["XOOZX"] =  1.25;

(*HelixProp_)["XXOOZ"] =  1.25;
(*HelixProp_)["ZXOOZ"] =  1.25;
(*HelixProp_)["XZOOZ"] =  1.25;
(*HelixProp_)["ZZOOX"] =  1.25;
(*HelixProp_)["ZZOOZ"] =  1.25;
(*HelixProp_)["ZOOXX"] =  1.25;
(*HelixProp_)["ZOOXZ"] =  1.25;
(*HelixProp_)["ZOOZX"] =  1.25;
(*HelixProp_)["XOOZZ"] =  1.25;
(*HelixProp_)["ZOOZZ"] =  1.25;

(*HelixProp_)["XXOUX"] =  1.25;
(*HelixProp_)["ZXOUX"] =  1.25;
(*HelixProp_)["XXUOX"] =  1.25;
(*HelixProp_)["ZXUOX"] =  1.25;
(*HelixProp_)["XOUXX"] =  1.25;
(*HelixProp_)["XOUXZ"] =  1.25;
(*HelixProp_)["XUOXX"] =  1.25;
(*HelixProp_)["XUOXZ"] =  1.25;

(*HelixProp_)["XXOUZ"] =  0.75;
(*HelixProp_)["ZXOUZ"] =  0.75;
(*HelixProp_)["XZOUX"] =  0.75;
(*HelixProp_)["XZOUZ"] =  0.75;
(*HelixProp_)["ZZOUX"] =  0.75;
(*HelixProp_)["ZZOUZ"] =  0.75;
(*HelixProp_)["XXUOZ"] =  0.75;
(*HelixProp_)["ZXUOZ"] =  0.75;
(*HelixProp_)["XZUOX"] =  0.75;
(*HelixProp_)["XZUOZ"] =  0.75;
(*HelixProp_)["ZZUOX"] =  0.75;
(*HelixProp_)["ZZUOZ"] =  0.75;
(*HelixProp_)["ZOUXX"] =  0.75;
(*HelixProp_)["ZOUXZ"] =  0.75;
(*HelixProp_)["XOUZX"] =  0.75;
(*HelixProp_)["ZOUZX"] =  0.75;
(*HelixProp_)["XOUZZ"] =  0.75;
(*HelixProp_)["ZOUZZ"] =  0.75;
(*HelixProp_)["ZUOXX"] =  0.75;
(*HelixProp_)["ZUOXZ"] =  0.75;
(*HelixProp_)["XUOZX"] =  0.75;
(*HelixProp_)["ZUOZX"] =  0.75;
(*HelixProp_)["XUOZZ"] =  0.75;
(*HelixProp_)["ZUOZZ"] =  0.75;
(*HelixProp_)["XUUXX"] =  1.25;
(*HelixProp_)["XXUUX"] =  1.25;

(*HelixProp_)["XXUUZ"] =  0.6;
(*HelixProp_)["ZXUUX"] =  0.6;
(*HelixProp_)["ZXUUZ"] =  0.6;
(*HelixProp_)["XZUUX"] =  0.6;
(*HelixProp_)["XZUUZ"] =  0.6;
(*HelixProp_)["ZZUUX"] =  0.6;
(*HelixProp_)["ZZUUZ"] =  0.6;
(*HelixProp_)["ZUUXX"] =  0.6;
(*HelixProp_)["XUUXZ"] =  0.6;
(*HelixProp_)["ZUUXZ"] =  0.6;
(*HelixProp_)["XUUZX"] =  0.6;
(*HelixProp_)["ZUUZX"] =  0.6;
(*HelixProp_)["XUUZZ"] =  0.6;
(*HelixProp_)["ZUUZZ"] =  0.6;

(*HelixProp_)["XXOOXX"] =  3;
(*HelixProp_)["XXOOXZ"] =  3;
(*HelixProp_)["ZXOOXX"] =  3;
(*HelixProp_)["ZXOOXZ"] =  3;
(*HelixProp_)["XXOUXX"] =  3;
(*HelixProp_)["XXOUXZ"] =  3;
(*HelixProp_)["XXUOXX"] =  3;
(*HelixProp_)["XXUOXZ"] =  3;
(*HelixProp_)["ZXUOXX"] =  3;
(*HelixProp_)["ZXOUXX"] =  3;

(*HelixProp_)["XXOOZX"] =  1.6;
(*HelixProp_)["XXOOZZ"] =  1.6;
(*HelixProp_)["XZOOXX"] =  1.6;
(*HelixProp_)["XZOOXZ"] =  1.6;
(*HelixProp_)["XZOOZX"] =  1.6;
(*HelixProp_)["XZOOZZ"] =  1.6;
(*HelixProp_)["ZXOOZX"] =  1.6;
(*HelixProp_)["ZXOOZZ"] =  1.6;
(*HelixProp_)["ZZOOXX"] =  1.6;
(*HelixProp_)["ZZOOXZ"] =  1.6;
(*HelixProp_)["ZXOUXZ"] =  1.6;
(*HelixProp_)["XZUOXX"] =  1.6;
(*HelixProp_)["ZXUOXZ"] =  1.6;

(*HelixProp_)["ZZOOZX"] =  1.5;
(*HelixProp_)["ZZOOZZ"] =  1.5;
(*HelixProp_)["XXOUZX"] =  1.5;
(*HelixProp_)["XXOUZZ"] =  1.5;
(*HelixProp_)["XZOUXX"] =  1.5;
(*HelixProp_)["XZOUXZ"] =  1.5;
(*HelixProp_)["ZXOUZX"] =  1.5;
(*HelixProp_)["ZXOUZZ"] =  1.5;
(*HelixProp_)["ZZOUXX"] =  1.5;
(*HelixProp_)["ZZOUXZ"] =  1.5;
(*HelixProp_)["XXUOZX"] =  1.5;
(*HelixProp_)["XXUOZZ"] =  1.5;
(*HelixProp_)["XZUOXZ"] =  1.5;
(*HelixProp_)["ZXUOZX"] =  1.5;
(*HelixProp_)["ZXUOZZ"] =  1.5;
(*HelixProp_)["ZZUOXX"] =  1.5;
(*HelixProp_)["ZZUOXZ"] =  1.5;

(*HelixProp_)["ZZUOZX"] =  1.25;
(*HelixProp_)["ZZUOZZ"] =  1.25;
(*HelixProp_)["ZZOUZX"] =  1.25;
(*HelixProp_)["ZZOUZZ"] =  1.25;
(*HelixProp_)["XZOUZX"] =  1.25;
(*HelixProp_)["XZOUZZ"] =  1.25;
(*HelixProp_)["XZUOZX"] =  1.25;
(*HelixProp_)["XZUOZZ"] =  1.25;
(*HelixProp_)["XXUUXX"] =  1.25;
(*HelixProp_)["XXUUXZ"] =  1.25;
(*HelixProp_)["ZXUUXX"] =  1.25;

(*HelixProp_)["XXUUZX"] =  1.25;
(*HelixProp_)["XXUUZZ"] =  1.25;
(*HelixProp_)["XZUUXX"] =  1.25;
(*HelixProp_)["XZUUXZ"] =  1.25;
(*HelixProp_)["XZUUZX"] =  0.75;
(*HelixProp_)["XZUUZZ"] =  0.75;
(*HelixProp_)["ZXUUXZ"] =  1.25;
(*HelixProp_)["ZXUUZX"] =  1.25;
(*HelixProp_)["ZXUUZZ"] =  1.25;
(*HelixProp_)["ZZUUXX"] =  1.25;
(*HelixProp_)["ZZUUXZ"] =  1.25;
(*HelixProp_)["ZZUUZX"] =  0.75;
(*HelixProp_)["ZZUUZZ"] =  0.75;


  if (ClustProp_) {
    return;
  }

  ClustProp_ = new map<string, double>;
 (*ClustProp_)["11"] =  0.3;
 (*ClustProp_)["15"] =  0.4;
 (*ClustProp_)["51"] =  0.4;
 (*ClustProp_)["55"] =  1.3;
 (*ClustProp_)["111"] =  0.5;
 (*ClustProp_)["115"] =  0.7;
 (*ClustProp_)["151"] =  0.7;
 (*ClustProp_)["155"] =  2.1;
 (*ClustProp_)["511"] =  0.7;
 (*ClustProp_)["515"] =  2.1;
 (*ClustProp_)["551"] =  2.1;
 (*ClustProp_)["555"] =  2.8;
 (*ClustProp_)["1111"] =  0.7;
 (*ClustProp_)["1115"] =  0.9;
 (*ClustProp_)["1151"] =  0.9;
 (*ClustProp_)["1155"] =  2.2;
 (*ClustProp_)["1511"] =  0.9;
 (*ClustProp_)["1515"] =  2.2;
 (*ClustProp_)["1551"] =  0.9;
 (*ClustProp_)["1555"] =  3;
 (*ClustProp_)["5111"] =  0.9;
 (*ClustProp_)["5115"] =  2.2;
 (*ClustProp_)["5151"] =  2.2;
 (*ClustProp_)["5155"] =  3;
 (*ClustProp_)["5511"] =  2.2;
 (*ClustProp_)["5515"] =  3;
 (*ClustProp_)["5551"] =  3;
 (*ClustProp_)["5555"] =  3.5;
 (*ClustProp_)["11111"] =  0.9;
 (*ClustProp_)["11115"] =  1;
 (*ClustProp_)["11151"] =  1;
 (*ClustProp_)["11155"] =  2.3;
 (*ClustProp_)["11511"] =  1;
 (*ClustProp_)["11515"] =  2.3;
 (*ClustProp_)["11551"] =  2.3;
 (*ClustProp_)["11555"] =  3.1;
 (*ClustProp_)["15111"] =  1;
 (*ClustProp_)["15115"] =  2.3;
 (*ClustProp_)["15151"] =  2.3;
 (*ClustProp_)["15155"] =  3.1;
 (*ClustProp_)["15511"] =  2.3;
 (*ClustProp_)["15515"] =  3.1;
 (*ClustProp_)["15551"] =  3.1;
 (*ClustProp_)["15555"] =  3.6;
 (*ClustProp_)["51111"] =  1.0;
 (*ClustProp_)["51115"] =  2.3;
 (*ClustProp_)["51151"] =  2.3;
 (*ClustProp_)["51155"] =  3.1;
 (*ClustProp_)["51511"] =  3.6;
 (*ClustProp_)["51515"] =  2.3;
 (*ClustProp_)["51551"] =  3.1;
 (*ClustProp_)["51555"] =  3.6;
 (*ClustProp_)["55111"] =  2.3;
 (*ClustProp_)["55115"] =  3.1;
 (*ClustProp_)["55151"] =  3.1;
 (*ClustProp_)["55155"] =  3.6;
 (*ClustProp_)["55511"] =  3.1;
 (*ClustProp_)["55515"] =  3.6;
 (*ClustProp_)["55551"] =  3.6;
 (*ClustProp_)["55555"] =  4.0;
 (*ClustProp_)["111111"] =  1.1;
 (*ClustProp_)["111115"] =  1.7;
 (*ClustProp_)["111151"] =  1.7;
 (*ClustProp_)["111155"] =  2.5;
 (*ClustProp_)["111511"] =  1.7;
 (*ClustProp_)["111515"] =  2.5;
 (*ClustProp_)["111551"] =  2.5;
 (*ClustProp_)["111555"] =  3.3;
 (*ClustProp_)["115111"] =  1.7;
 (*ClustProp_)["115115"] =  2.5;
 (*ClustProp_)["115151"] =  2.5;
 (*ClustProp_)["115155"] =  3.3;
 (*ClustProp_)["115511"] =  2.5;
 (*ClustProp_)["115515"] =  3.3;
 (*ClustProp_)["115551"] =  3.3;
 (*ClustProp_)["115555"] =  3.7;
 (*ClustProp_)["151111"] =  1.7;
 (*ClustProp_)["151115"] =  2.5;
 (*ClustProp_)["151151"] =  2.5;
 (*ClustProp_)["151155"] =  3.3;
 (*ClustProp_)["151511"] =  2.5;
 (*ClustProp_)["151515"] =  3.3;
 (*ClustProp_)["151551"] =  3.3;
 (*ClustProp_)["151555"] =  3.7;
 (*ClustProp_)["155111"] =  2.5;
 (*ClustProp_)["155115"] =  3.3;
 (*ClustProp_)["155151"] =  3.3;
 (*ClustProp_)["155155"] =  3.7;
 (*ClustProp_)["155511"] =  3.3;
 (*ClustProp_)["155515"] =  3.7;
 (*ClustProp_)["155551"] =  3.7;
 (*ClustProp_)["155555"] =  4.1;
 (*ClustProp_)["511111"] =  1.7;
 (*ClustProp_)["511115"] =  2.5;
 (*ClustProp_)["511151"] =  2.5;
 (*ClustProp_)["511155"] =  3.3;
 (*ClustProp_)["511511"] =  2.5;
 (*ClustProp_)["511515"] =  3.3;
 (*ClustProp_)["511551"] =  3.3;
 (*ClustProp_)["511555"] =  3.7;
 (*ClustProp_)["515111"] =  2.5;
 (*ClustProp_)["515115"] =  3.3;
 (*ClustProp_)["515151"] =  3.3;
 (*ClustProp_)["515155"] =  3.7;
 (*ClustProp_)["515511"] =  3.3;
 (*ClustProp_)["515515"] =  3.7;
 (*ClustProp_)["515551"] =  3.7;
 (*ClustProp_)["515555"] =  4.1;
 (*ClustProp_)["551111"] =  2.5;
 (*ClustProp_)["551115"] =  3.3;
 (*ClustProp_)["551151"] =  3.3;
 (*ClustProp_)["551155"] =  3.7;
 (*ClustProp_)["551511"] =  3.3;
 (*ClustProp_)["551515"] =  3.7;
 (*ClustProp_)["551551"] =  3.7;
 (*ClustProp_)["551555"] =  4.1;
 (*ClustProp_)["555111"] =  3.3;
 (*ClustProp_)["555115"] =  3.7;
 (*ClustProp_)["555151"] =  3.7;
 (*ClustProp_)["555155"] =  4.1;
 (*ClustProp_)["555511"] =  3.7;
 (*ClustProp_)["555515"] =  4.1;
 (*ClustProp_)["555551"] =  4.1;
 (*ClustProp_)["555555"] =  4.5;


 
 if (HydroPath_) {
   return;
 }
 
 HydroPath_ = new map<string, double>;
 
 (*HydroPath_)["W"] = 0.81;
 (*HydroPath_)["F"] = 1.2;
 (*HydroPath_)["L"] = 1.1;
 (*HydroPath_)["I"] = 1.4;
 
 (*HydroPath_)["M"] = 0.64;
 (*HydroPath_)["M[147]"] = 0.64; // assume M[147] behaves as M
 
 (*HydroPath_)["V"] = 1.1;
 (*HydroPath_)["Y"] = 0.26;
 (*HydroPath_)["A"] = 0.62;
 (*HydroPath_)["T"] = -0.05;
 (*HydroPath_)["P"] = 0.12;
 (*HydroPath_)["E"] = -0.74;
 (*HydroPath_)["D"] = -0.90;
 
 (*HydroPath_)["C"] = 0.29; // assume C behaves as C[161]
 (*HydroPath_)["C[160]"] = 0.29; // assume C[160] behaves as C[161]
 (*HydroPath_)["C[161]"] = 0.29; // original SSRCalc actually uses C[161]
 
 (*HydroPath_)["S"] = -0.18;
 (*HydroPath_)["Q"] = -0.85;
 (*HydroPath_)["G"] = 0.48;
 (*HydroPath_)["N"] = -0.78;
 (*HydroPath_)["R"] = -2.5;
 (*HydroPath_)["H"] = -0.4;
 (*HydroPath_)["K"] = -1.5;

 if (EludeIndex_) {
   return;
 }

  EludeIndex_ = new map<string, double>;

  (*EludeIndex_)["W"] = 0.698326;
  (*EludeIndex_)["F"] = 0.896282;
  (*EludeIndex_)["L"] = 1.21493;
  (*EludeIndex_)["I"] = 0.793001;

  (*EludeIndex_)["M"] = 0.41043;
  (*EludeIndex_)["M[147]"] = 0.41043; // assume M[147] behaves as M

  (*EludeIndex_)["V"] = 0.790132;
  (*EludeIndex_)["Y"] = 0.38229;
  (*EludeIndex_)["A"] = 0.294636;
  (*EludeIndex_)["T"] = 0.077593;
  (*EludeIndex_)["P"] = -0.00652169;
  (*EludeIndex_)["E"] = 0.14578;
  (*EludeIndex_)["D"] = 0.109216;

  (*EludeIndex_)["C"] = -0.0960899; // assume C behaves as C[161]
  (*EludeIndex_)["C[160]"] = -0.0960899; // assume C[160] behaves as C[161]
  (*EludeIndex_)["C[161]"] = -0.0960899; // original SSRCalc actually uses C[161]

  (*EludeIndex_)["S"] = -0.142022;
  (*EludeIndex_)["Q"] = -0.0524759;
  (*EludeIndex_)["G"] = -0.148832;
  (*EludeIndex_)["N"] = -0.079866;
  (*EludeIndex_)["R"] = -0.0468878;
  (*EludeIndex_)["H"] = -1.04939;
  (*EludeIndex_)["K"] = -0.521984;


}



bool RTCalculator::linearRegressRT(double min_prob, int min_ntt) {

  //  min_prob = 0;
#ifdef __LGPL__
  // cannot do fitting, so just calculate "predictedRT" by SSRCalc 1.0 algorithm
  // (published in MCP)
  return (false);



#else
  
  cerr << "Performing linear regression of retention times..." << endl;
  
  // The sequence vs RT relationship is learned by linear regression, rather
  // than predicted by SSRCalc or Hydro. This is more robust and generally applicable.
  // Note that RT predictors like SSRCalc or Hydro makes assumptions about the column used,
  // components and gradient of mobile phase, etc, and cannot deal with modified amino acids.
  
  int numDataPoints = (int)learn_rts_.size();

  int numGood = 0;
  
  for (int d = 0; d < numDataPoints; d++) {
    if ((*learn_probs_)[d] >= min_prob && (*learn_ntts_)[d] >= min_ntt) {
      numGood++;
      //     if (min_RT_ <= 0 || learn_rts_[d] < min_RT_) {
      //	min_RT_ = learn_rts_[d];
      //}

    }
  }

  //adjustRT();

  cerr << "Fitting model to " << numGood << " good identifications with P>=" << min_prob << " and NTT>=" << min_ntt << "." << endl;
  
  if (numGood < 80) {
    cerr << "Not enough good identifications to fit retention time model." << endl;
    cerr << "Reverting to SSRCalc 1.0 predictions." << endl;
    return (false);
  }
    
  map<string, int> param_index;
  vector<string> params;
  int paramcount = 0;

  // add the constant
  param_index["1"] = paramcount++;
  params.push_back("1");
  // add the log of the peptide length as one of the modeled variables
  param_index["ln(LEN)"] = paramcount++;
  params.push_back("ln(LEN)");

  // find out what amino acids and modified amino acids are present in the identifications
  for (int dd = 0; dd < (int)modified_peps_.size(); dd++) {
    string pep =  modified_peps_[dd];
    string::size_type pos = 0;
    int numAA = 0;
    while (pos != string::npos) {
      string aatag = nextAAToken(pep, pos, pos);
      if (param_index.find(aatag) == param_index.end()) {
	param_index[aatag] = paramcount;
	paramcount++;
	params.push_back(aatag);
      }
      if (numAA < 3) {
	// create another variable to model the disproportionate influence of 
	// amino acids close to the N-terminus
	aatag = "n_" + aatag;
	if (param_index.find(aatag) == param_index.end()) {
	  param_index[aatag] = paramcount;
	  paramcount++;
	  params.push_back(aatag);
	}	
      }
      numAA++;
    }
  }
   
  
  if (!matrix_) {
    
    matrix_ = new vector<vector<double> >; 

    // fill in the matrix -- this will be the input to the linear regression algorithm
    for (int d = 0; d < (int)modified_peps_.size(); d++) {
    
      vector<double> x(paramcount + 1, 0.0);
      
      string pep =  modified_peps_[d];
      
      // count numbers of each type of amino acids
      string::size_type pos = 0;
      int numAA = 0;
      while (pos != string::npos) {
	string aatag = nextAAToken(pep, pos, pos);
        map<string, int>::iterator found = param_index.find(aatag);
        if (found != param_index.end()) {
	  x[found->second]++;
        }       
	
	if (numAA < 3) {
	  aatag = "n_" + aatag;
	  map<string, int>::iterator found1 = param_index.find(aatag);
	  if (found1 != param_index.end()) {
	    if (numAA == 0) {
	      x[found1->second] += 0.5; // approx weighting based on SSRCalc 1
	    } else if (numAA == 1) {
	      x[found1->second] += 0.3;
	    } else if (numAA == 2) {
	      x[found1->second] += 0.1;
	      
	    }
	  }  
	}
		
	numAA++;
	
      }
        
      x[0] = 1.0;
      x[1] = log((double)numAA);
      //      for (int z = 2; z < paramcount; z++) {
      //	if (params[z][0] >= 'A' && params[z][0] <= 'Z') {
      //	  x[z] /= ((double)(numAA));
      //	}
      // }
      // the very last column is the retention time itself (can also use scan number)
      x[paramcount] = rts_[d];
    
      matrix_->push_back(x);
      
    }
  }

  if (!learn_matrix_) {
    
    learn_matrix_ = new vector<vector<double> >; 

    // fill in the matrix -- this will be the input to the linear regression algorithm
    for (int d = 0; d < (int)learn_peps_.size(); d++) {
    
      vector<double> x(paramcount + 1, 0.0);
      
      string pep =  learn_peps_[d];
      
      // count numbers of each type of amino acids
      string::size_type pos = 0;
      int numAA = 0;
      while (pos != string::npos) {
	string aatag = nextAAToken(pep, pos, pos);
        map<string, int>::iterator found = param_index.find(aatag);
        if (found != param_index.end()) {
	  x[found->second]++;
        }       
	
	if (numAA < 3) {
	  aatag = "n_" + aatag;
	  map<string, int>::iterator found1 = param_index.find(aatag);
	  if (found1 != param_index.end()) {
	    if (numAA == 0) {
	      x[found1->second] += 0.5; // approx weighting based on SSRCalc 1
	    } else if (numAA == 1) {
	      x[found1->second] += 0.3;
	    } else if (numAA == 2) {
	      x[found1->second] += 0.1;
	      
	    }
	  }  
	}
		
	numAA++;
	
      }
        
      x[0] = 1.0;
      x[1] = log((double)numAA);
      //      for (int z = 2; z < paramcount; z++) {
      //	if (params[z][0] >= 'A' && params[z][0] <= 'Z') {
      //	  x[z] /= ((double)(numAA));
      //	}
      // }
      // the very last column is the retention time itself (can also use scan number)
      x[paramcount] = learn_rts_[d];
    
      learn_matrix_->push_back(x);
      
    }
  }


  // using GNU Scientific Library (GSL) to do fitting.
  // note that GSL is GPL, so this feature must be disabled if LGPL is desired
  if (rt_fit_ != NULL ) gsl_multifit_linear_free(rt_fit_);
  rt_fit_ = gsl_multifit_linear_alloc(numGood, paramcount); 
  
  gsl_vector* y = gsl_vector_calloc(numGood);
  gsl_vector* w = gsl_vector_calloc(numGood);
  gsl_matrix* X = gsl_matrix_calloc(numGood, paramcount); 
  int count = 0;
  for(int d = 0; d < (int)learn_matrix_->size(); d++){
    if ((*learn_probs_)[d] >= min_prob && (*learn_ntts_)[d] >= min_ntt) {
      gsl_vector_set(y, count, (*learn_matrix_)[d][paramcount]);
      gsl_vector_set(w, count, (*learn_probs_)[d]);
      for(int j = 0; j < paramcount; j++) {
	gsl_matrix_set(X, count, j, (*learn_matrix_)[d][j]);
      }
      count++;
    }
  }
  
  /* DEBUG code to print out matrix

  cerr << endl << "==============" << endl;
  for (int row = 0; row < count; row++) {
    cerr << gsl_vector_get(y, row);
      for (int col = 0; col < 41; col++) {
      cerr << gsl_matrix_get(X, row, col) << ' ';
    }
    cerr << endl;
  }
  */
  
  gsl_vector *c = gsl_vector_calloc(paramcount);
  gsl_vector_set_zero(c);
  gsl_matrix *cov = gsl_matrix_calloc(paramcount, paramcount);
  gsl_matrix_set_zero(cov);
  double chisq;
  double tolerance = 0.01;
  size_t rank = 0;

  int err = gsl_multifit_wlinear_svd(X, w, y, tolerance, &rank, c, cov, &chisq, rt_fit_);
  int newNumGood = 0;
  
  //DDS: Outlier removal here
  
  gsl_vector *dis =  gsl_vector_calloc(numGood);
  int ix = 0;
  double dis_mean, dis_sd;

  
  //DDS: Find distance mean and stddev

  for (int d = 0; d < (int)learn_matrix_->size(); d++) {

    if ((*learn_probs_)[d] >= min_prob && (*learn_ntts_)[d] >= min_ntt) {
      double predictedRT;
      double predictedDev;
      gsl_vector *x = gsl_vector_calloc(paramcount);

      for (int j = 0; j < paramcount; j++) {
	gsl_vector_set(x, j, (*learn_matrix_)[d][j]);
      }
      err = gsl_multifit_linear_est(x, c, cov, &predictedRT, &predictedDev);

      gsl_vector_set(dis, ix++, (learn_rts_)[d] - predictedRT);

      gsl_vector_free(x);
      
    }    
  }
  
  dis_mean = gsl_stats_wmean(gsl_vector_ptr(w, 0), 1,gsl_vector_ptr(dis, 0), 1, numGood);
  run_RT_mean_ = dis_mean;
  dis_sd = gsl_stats_wsd(gsl_vector_ptr(w, 0), 1,gsl_vector_ptr(dis, 0), 1, numGood);
  run_RT_stddev_ = dis_sd;


  count = 0;
  
  for (int d = 0; d < (int)learn_matrix_->size(); d++) {

    if ((*learn_probs_)[d] >= min_prob && (*learn_ntts_)[d] >= min_ntt) {
      double predictedRT;
      double predictedDev;
      gsl_vector *x = gsl_vector_calloc(paramcount);

      for (int j = 0; j < paramcount; j++) {
	gsl_vector_set(x, j, (*learn_matrix_)[d][j]);
      }
      err = gsl_multifit_linear_est(x, c, cov, &predictedRT, &predictedDev);

      if (fabs(((learn_rts_)[d] - predictedRT - dis_mean)/dis_sd) < 2) {
	gsl_vector_set(y, count, (*learn_matrix_)[d][paramcount]);
	gsl_vector_set(w, count, (*learn_probs_)[d]);
	for(int j = 0; j < paramcount; j++) {
	  gsl_matrix_set(X, count, j, (*learn_matrix_)[d][j]);
	}
	count++;
	
      }
      gsl_vector_free(x);
    }    
  }
  
  gsl_vector_free(y);
  gsl_vector_free(w);
  gsl_matrix_free(X);
  
  numGood = count;
  run_RT_used_count_ = count;
  count = 0;
  y = gsl_vector_calloc(numGood);
  w = gsl_vector_calloc(numGood);
  X = gsl_matrix_calloc(numGood, paramcount); 
   
  for (int d = 0; d < (int)learn_matrix_->size(); d++) {

    if ((*learn_probs_)[d] >= min_prob && (*learn_ntts_)[d] >= min_ntt) {
      double predictedRT;
      double predictedDev;
      gsl_vector *x = gsl_vector_calloc(paramcount);

      for (int j = 0; j < paramcount; j++) {
	gsl_vector_set(x, j, (*learn_matrix_)[d][j]);
      }
      err = gsl_multifit_linear_est(x, c, cov, &predictedRT, &predictedDev);

      if (fabs(((learn_rts_)[d] - predictedRT - dis_mean)/dis_sd) < 2) {
	gsl_vector_set(y, count, (*learn_matrix_)[d][paramcount]);
	gsl_vector_set(w, count, (*learn_probs_)[d]);
	for (int j = 0; j < paramcount; j++) {
	  gsl_matrix_set(X, count, j, (*learn_matrix_)[d][j]);
	}
	count++;
      }
      
      gsl_vector_free(x);
    }    
  }

  gsl_vector_free(dis);

  cerr << "After outlier removal, fitting model to " << numGood << " good identifications with P>=" << min_prob << " and NTT>=" << min_ntt << "." << endl;
  
  if (numGood < 80) {
    cerr << "Not enough good identifications to fit retention time model." << endl;
    cerr << "Reverting to SSRCalc 1.0 predictions." << endl;
    return (false);
  }
  
  //  cerr << "RANK=" << rank << endl;
  //  cerr << "DONE FITTING" << endl;
  if (rt_fit_ != NULL )gsl_multifit_linear_free(rt_fit_);
  rt_fit_ = gsl_multifit_linear_alloc(numGood, paramcount); 
  err = gsl_multifit_wlinear_svd(X, w, y,  tolerance, &rank, c, cov, &chisq, rt_fit_);

 for (map<string, int>::iterator pp = param_index.begin(); pp != param_index.end(); pp++) {
    coeff_[pp->first] = gsl_vector_get(c, pp->second);
  }
  
  double tss = (count - 1) * gsl_stats_wvariance(gsl_vector_ptr(w, 0), 1, gsl_vector_ptr(y, 0), 1, count);
  r_sq_ = 1 - (chisq / tss);
    
  RTs_->reserve((int)rts_.size());
  run_RT_count_=0;
  run_RT_sum_=0;
  learned_RT_count_=0;
  learned_RT_sum_=0;
  // now calculate the "predicted" RTs from this regressed model
  for (int d = 0; d < (int)rts_.size(); d++) {
    double predictedRT;
    double predictedDev;
    gsl_vector *x = gsl_vector_calloc(paramcount);
    for (int j = 0; j < paramcount; j++) {
      gsl_vector_set(x, j, (*matrix_)[d][j]);
    }
      
    err = gsl_multifit_linear_est(x, c, cov, &predictedRT, &predictedDev);
  
    run_RT_count_++;
    run_RT_sum_ += predictedRT;
    
    RTs_->insertAtEnd(predictedRT);
  
    gsl_vector_free(x);
  }

  // now calculate the "predicted" RTs from this regressed model
  for (int d = 0; d < (int)learn_rts_.size(); d++) {
    double predictedRT;
    double predictedDev;
    gsl_vector *x = gsl_vector_calloc(paramcount);
    for (int j = 0; j < paramcount; j++) {
      gsl_vector_set(x, j, (*matrix_)[d][j]);
    }
      
    err = gsl_multifit_linear_est(x, c, cov, &predictedRT, &predictedDev);
  
    learned_RT_count_++;
    learned_RT_sum_ += predictedRT;
    
    learned_RTs_->insertAtEnd(predictedRT);

  
    gsl_vector_free(x);
  }
    
  gsl_vector_free(c);
  gsl_matrix_free(cov);
  gsl_vector_free(y);
  gsl_matrix_free(X);
// gsl_multifit_linear_free(work);
  
//  plotRegressionFits(probs_, min_prob, ntts_, min_ntt);
  
  if (r_sq_ < 0.6) {
    // too much scatter to trust this?
    return (false);
  }
#endif  


  return (true);
}


Boolean RTCalculator::recalc_RTstats( double min_prob, int min_ntt) {
  double RTnumer = 0;
  double RTdenom = 0;
  double tot_RTs = 0;
  double all_RTs = 0;
  double SCANnumer = 0;
  double SCANdenom = 0;
  double tot_SCANs = 0;
  double all_SCANs = 0;
  double prob;  
  double RT_mean = 0;
  double RT_stddev = 0;
  double SCAN_mean = 0;
  double SCAN_stddev = 0;
  
  if (ready_) return True;

  run_RT_used_count_ = 0;
  //   cout << "DDS: probs_Size=" << learn_probs_->size() << " runRTNum=" << learned_RT_count_<< endl;
  assert(learn_probs_->size() == learned_RT_count_);
  int i;
  for (i=0; i<learned_RT_count_; i++) {
    if ((*learn_probs_)[i] >= 0) {
      prob = (*learn_probs_)[i];
      all_RTs += (*learned_RTs_)[i];
      all_SCANs += learn_rts_[i];
      //      all_SCANs += scans_[i];
      if (prob >= min_prob && (*learn_ntts_)[i] >= min_ntt) {
	tot_RTs += (*learned_RTs_)[i];
	RTnumer += prob*(*learned_RTs_)[i];
	RTdenom += prob;
	SCANnumer += prob*learn_rts_[i];
	//	SCANnumer += prob*scans_[i];
	SCANdenom += prob;
	run_RT_used_count_ ++;
      }
    }
  }

  if (RTdenom < 5) {
    cerr << "WARNING: Not enough high probability IDs in run index " << *run_name_ << " to generate RT model. RT Model has been disabled." << endl;
    return False;
  }
  else {
    RT_mean = RTnumer / RTdenom;
    RT_stddev = 0;
    SCAN_mean = SCANnumer / SCANdenom;
    SCAN_stddev = 0;
    RTdenom = 0;
    SCANdenom = 0;
    for (i=0; i<learned_RT_count_; i++) {
      if ((*learn_probs_)[i] >= 0) {
	prob = (*learn_probs_)[i];
	if (prob >= min_prob && (*learn_ntts_)[i] >= min_ntt ) {
	  RT_stddev += prob * pow((RT_mean - (*learned_RTs_)[i]), 2);
	  SCAN_stddev += prob * pow((SCAN_mean - learn_rts_[i]), 2);
	  //	  SCAN_stddev += prob * pow((SCAN_mean - scans_[i]), 2);
	  RTdenom += prob;
	  SCANdenom += prob;
	}
      }
    }
  }

  RT_stddev /= RTdenom;

  RT_stddev = pow(RT_stddev, 0.5); 
  
  run_RT_mean_ = RT_mean;
  run_RT_stddev_ = RT_stddev;

  SCAN_stddev /= SCANdenom;

  SCAN_stddev = pow(SCAN_stddev, 0.5); 
  
  run_SCAN_mean_ = SCAN_mean;
  run_SCAN_stddev_ = SCAN_stddev;

  double SSxx = 0;
  double SSxy = 0;
  double SSyy = 0;

  
  
  for (i=0; i<learned_RT_count_; i++) {
    if ((*learn_probs_)[i] >= 0) {
      if ((*learn_probs_)[i] >= min_prob && (*learn_ntts_)[i] >= min_ntt) {
	SSyy += pow((learn_rts_[i] - run_SCAN_mean_), 2);
	//	SSyy += pow((scans_[i] - run_SCAN_mean_), 2);
	SSxx += pow(((*learned_RTs_)[i] - run_RT_mean_), 2);
	SSxy += (learn_rts_[i] - run_SCAN_mean_)*((*learned_RTs_)[i] - run_RT_mean_);
      }
    }
  }  
  run_slope_ = SSxy / SSxx;
  run_intercept_ = run_SCAN_mean_ - run_slope_ * run_RT_mean_;
  r_sq_ = (SSxy*SSxy) / (SSxx*SSyy);
  if (r_sq_ < 0.5) {
    cerr << "WARNING: Not enough correlation between theoretical RT values and scan numbers in run index " << *run_name_ << ". RT Model has been disabled." << endl;
    return False;
  }
  ready_ = true;
  return True;
}



bool RTCalculator::linearRegressRT() {
  double numSD = 20;

  //  min_prob = 0;
#ifdef __LGPL__
  // cannot do fitting, so just calculate "predictedRT" by SSRCalc 1.0 algorithm
  // (published in MCP)
  return (false);



#else
  
  cerr << "Performing linear regression of retention times..." << endl;
  
  // The sequence vs RT relationship is learned by linear regression, rather
  // than predicted by SSRCalc or Hydro. This is more robust and generally applicable.
  // Note that RT predictors like SSRCalc or Hydro makes assumptions about the column used,
  // components and gradient of mobile phase, etc, and cannot deal with modified amino acids.
  
  int numDataPoints = (int)learn_rts_.size();

  int numGood = numDataPoints;
  

  //adjustRT();

  cerr << "Fitting model to " << numGood << "  peptides ." << endl;
  
  if (numGood < 80) {
    cerr << "Not enough good identifications to fit retention time model." << endl;
    cerr << "Reverting to SSRCalc 1.0 predictions." << endl;
    return (false);
  }
    
  //map<string, int> param_index;
  param_index_.clear();
  vector<string> params;
  paramcount_ = 0;

  // add the constant
  param_index_["1"] = paramcount_++;
  params.push_back("1");

 
  // add the log of the peptide length as one of the modeled variables
  param_index_["ln(LEN)"] = paramcount_++;
  params.push_back("ln(LEN)");

  param_index_["pI"] = paramcount_++;
  params.push_back("pI");

  param_index_["helicity1"] = paramcount_++;
  params.push_back("helicity1");
  
  param_index_["clusterness"] = paramcount_++;
  params.push_back("clusterness");

   
  //param_index_["helicity2"] = paramcount_++;
  //params.push_back("helicity2");

  param_index_["helectric"] = paramcount_++;
  params.push_back("helectric");


  param_index_["hydroMoment"] = paramcount_++;
  params.push_back("hydroMoment");

  
  param_index_["EludeAmphipathicityHelixMin"] = paramcount_++;
  params.push_back("EludeAmphipathicityHelixMin");


  param_index_["EludeAmphipathicityHelixMax"] = paramcount_++;
  params.push_back("EludeAmphipathicityHelixMax");


  param_index_["EludeHydrophobicMoment180Min"] = paramcount_++;
  params.push_back("EludeHydrophobicMoment180Min");


  param_index_["EludeHydrophobicMoment180Max"] = paramcount_++;
  params.push_back("EludeHydrophobicMoment180Max");


  param_index_["EludeHydrophobicMoment100Min"] = paramcount_++;
  params.push_back("EludeHydrophobicMoment100Min");


  param_index_["EludeHydrophobicMoment100Max"] = paramcount_++;
  params.push_back("EludeHydrophobicMoment100Max");

//  param_index_["EludeIndexPartialSum5Max"] = paramcount_++;
//   params.push_back("EludeIndexPartialSum5Max");

//   param_index_["EludeIndexPartialSum2Max"] = paramcount_++;
//   params.push_back("EludeIndexPartialSum2Max");

//   param_index_["EludeIndexPartialSum5Min"] = paramcount_++;
//   params.push_back("EludeIndexPartialSum5Min");

//   param_index_["EludeIndexPartialSum2Min"] = paramcount_++;
//   params.push_back("EludeIndexPartialSum2Min");

  //param_index_["EludeIndexSum"] = paramcount_++;
  //params.push_back("EludeIndexSum");

  //param_index_["EludeIndexAvg"] = paramcount_++;
  //params.push_back("EludeIndexAvg");

  //param_index_["EludeIndexN"] = paramcount_++;
  //params.push_back("EludeIndexN");

  //param_index_["EludeIndexC"] = paramcount_++;
  //params.push_back("EludeIndexC");

  //param_index_["EludeIndexNearestNeighborPos"] = paramcount_++;
  //params.push_back("EludeIndexNearestNeighborPos");

  //param_index_["EludeIndexNearestNeighborNeg"] = paramcount_++;
  //params.push_back("EludeIndexNearestNeighborNeg");

  //param_index_["EludeIndexSumSqrDiff"] = paramcount_++;
  //params.push_back("EludeIndexSumSqrDiff");

  //  param_index_["smallness"] = paramcount_++;
  //params.push_back("smallness");


  string::size_type pos;
  // find out what amino acids and modified amino acids are present in the identifications
  for (int dd = 0; dd < (int)learn_peps_.size(); dd++) {
    string pep =  learn_peps_[dd];
    int totAA = countAA(pep);;
    pos = 0;
    int numAA = 0;
    string aatag, nextaa;
    while (pos != string::npos || numAA == totAA-1) {
      
      if (numAA == 0) {
	aatag = nextAAToken(pep, pos, pos);
      }
      else {
	aatag = nextaa;
      }
      

      if (numAA < totAA-1) {
 	nextaa =  nextAAToken(pep, pos, pos);
	
      }

      if (param_index_.find(aatag) == param_index_.end()) {
	param_index_[aatag] = paramcount_;
	paramcount_++;
	params.push_back(aatag);
      }

      
      if (aatag.find_first_of("nc") == string::npos && (numAA < 3 || numAA >= totAA-3)) {
	aatag = "term_" + aatag;
	//if (numAA<3) {
	//  aatag = "Nterm_" + aatag;
	//}
	//else {
	//  aatag = "Cterm_" + aatag;
	//}

	if (param_index_.find(aatag) == param_index_.end()) {
	  param_index_[aatag] = paramcount_;
	  paramcount_++;
	  params.push_back(aatag);
	}	
	
      }
//       else if (numAA >= totAA-3) {
// 	if (numAA == totAA-3) {
// 	  aatag = "Cnearterm_" + aatag;
// 	}
// 	if (numAA == totAA-2) {
// 	  aatag = "Cnearterm_" + aatag;
// 	}
// 	if (numAA == totAA-1) {
// 	  aatag = "Cterm_" + aatag;
// 	}
	
//       }
      
  


      numAA++;

	
      
    }
  }
  
  double score = 1.;
  map<string, int>::iterator found;
  if (!learn_matrix_) {
    
    learn_matrix_ = new vector<vector<double> >; 

    // fill in the matrix -- this will be the input to the linear regression algorithm
    for (int d = 0; d < (int)learn_peps_.size(); d++) {
      pos = 0;
      vector<double> x(paramcount_ + 1, 0.0);
      
      string pep =  learn_peps_[d];
      
      int totAA = countAA(pep);
      
      // count numbers of each type of amino acids

      int numAA = 0;
      string aatag, nextaa;
      while (pos != string::npos || numAA == totAA-1) {
	
	if (numAA == 0) {
	  aatag = nextAAToken(pep, pos, pos);
	}
	else {
	  aatag = nextaa;
	}
	
	score = 1.;
 	if (numAA < totAA-1) {
 	  nextaa =  nextAAToken(pep, pos, pos);

 	}
	found = param_index_.find(aatag);
        if (found != param_index_.end()) {
	  x[found->second]+=score;
        }     

	if (aatag.find_first_of("nc") == string::npos && (numAA < 3 || numAA >= totAA-3)) {
	  aatag = "term_" + aatag;
	  if (numAA == 0 || numAA == totAA-1) {
	    score = 0.5;
	  }
	  if (numAA == 1 || numAA == totAA-2)  {
	    score = 0.3;
	  }
	  
	  if (numAA == 2 || numAA == totAA-3)  {
	    score = 0.1;
	  }

// 	  if (numAA<3) {
// 	    aatag = "Nterm_" + aatag;
// 	  }
// 	  else {
// 	    aatag = "Cterm_" + aatag;
// 	  }

	  found = param_index_.find(aatag);
	  if (found != param_index_.end()) {
	    x[found->second]+=score;
	  }     
	}
// 	else if (numAA >= totAA-3) {
// 	  if (numAA == totAA-3) {
// 	    aatag = "Cnearterm_" + aatag;
// 	    score = 0.5;
// 	  }
// 	  if (numAA == totAA-2) {
// 	    aatag = "Cnearterm_" + aatag;
// 	    score = 1.;
// 	  }
// 	  if (numAA == totAA-1) {
// 	    aatag = "Cterm_" + aatag;
// 	    score = 1.;
// 	  }

// 	}
	
	numAA++;

      }
	

      found = param_index_.find("1");
      if (found != param_index_.end()) {
	x[found->second]=1.;
      }   
      
      found = param_index_.find("ln(LEN)");
      if (found != param_index_.end()) {
	x[found->second]=log((double)numAA);
      }   
      
      found = param_index_.find("pI");
      if (found != param_index_.end()) {
	//x[found->second]=log(pIcalc_->ModifiedPeptide_pI(pep.c_str(),  0, pep.length()));  
	x[found->second]=isoelectric(pep);
    
      }   
      
      
      found = param_index_.find("helicity1");
      if (found != param_index_.end()) {
	x[found->second]=helicity1(pep);
      }   
      
      found = param_index_.find("clusterness");
      if (found != param_index_.end()) {
	x[found->second]=clusterness(pep);
      }   
      
      //   found = param_index_.find("helicity2");
      //if (found != param_index_.end()) {
      //	x[found->second]=helicity2(pep);
      //}   
      
      found = param_index_.find("helectric");
      if (found != param_index_.end()) {
	x[found->second]=helectric(pep);
      }   
  
      found = param_index_.find("hydroMoment");
      if (found != param_index_.end()) {
      	x[found->second]=hydroMoment(pep);
      }  

      found = param_index_.find("EludeAmphipathicityHelixMin");
      if (found != param_index_.end()) {
	EludeAmphipathicityHelix(pep, &x);
      }   
      
      found = param_index_.find("EludeHydrophobicMoment100Min");
      if (found != param_index_.end()) {
	EludeHydrophobicMoments(pep, &x);
      }   
  
      found = param_index_.find("EludeIndexPartialSum5Max");
      if (found != param_index_.end()) {
	EludeIndexPartialSum(pep, &x);
      }   

      found = param_index_.find("EludeIndexSum");
      if (found != param_index_.end()) {
	x[found->second]=EludeIndexSum(pep);
      }   
      found = param_index_.find("EludeIndexAvg");
      if (found != param_index_.end()) {
	x[found->second]=EludeIndexAvg(pep);
      }   
      found = param_index_.find("EludeIndexN");
      if (found != param_index_.end()) {
	x[found->second]=EludeIndexN(pep);
      }   
      found = param_index_.find("EludeIndexC");
      if (found != param_index_.end()) {
	x[found->second]=EludeIndexC(pep);
      }   
      found = param_index_.find("EludeIndexNearestNeighborPos");
      if (found != param_index_.end()) {
	x[found->second]=EludeIndexNearestNeighborPos(pep);
      }   
      found = param_index_.find("EludeIndexNearestNeighborNeg");
      if (found != param_index_.end()) {
	x[found->second]=EludeIndexNearestNeighborNeg(pep);
      } 
      found = param_index_.find("EludeIndexSumSqrDiff");
      if (found != param_index_.end()) {
	x[found->second]=EludeIndexSumSqrDiff(pep);
      }   
      
 

      //      found = param_index_.find("smallness");
      //if (found != param_index_.end()) {
      //	x[found->second]=smallness(pep);
      //}   
  
        
// 	x[0] = 1.0;
// 	x[1] = log((double)numAA);
	
	
// 	//	x[3] = computeHelixScore(pep);
// 	x[2] = pIcalc_->ModifiedPeptide_pI(pep.c_str(),  NULL);
// 	x[3] = helicity1(pep);
// 	x[4] = clusterness(pep);




	//	x[4] = helicity2(pep);
	//x[5] = helectric(pep);

	//      for (int z = 2; z < paramcount_; z++) {
	//	if (params[z][0] >= 'A' && params[z][0] <= 'Z') {
	//	  x[z] /= ((double)(numAA));
	//	}
	// }
	// the very last column is the retention time itself (can also use scan number)
	x[paramcount_] = learn_rts_[d];
	
	learn_matrix_->push_back(x);
	
     
    }
  }


  // using GNU Scientific Library (GSL) to do fitting.
  // note that GSL is GPL, so this feature must be disabled if LGPL is desired
  if (rt_fit_ != NULL ) gsl_multifit_linear_free(rt_fit_);
  rt_fit_ = gsl_multifit_linear_alloc(numGood, paramcount_); 
  
  gsl_vector* y = gsl_vector_calloc(numGood);
  gsl_vector* w = gsl_vector_calloc(numGood);
  gsl_matrix* X = gsl_matrix_calloc(numGood, paramcount_); 
  int count = 0;
  for(int d = 0; d < (int)learn_matrix_->size(); d++){
      gsl_vector_set(y, count, (*learn_matrix_)[d][paramcount_]);
      //      gsl_vector_set(w, count, (*learn_probs_)[d]);
      gsl_vector_set(w, count, 1);
      for(int j = 0; j < paramcount_; j++) {
	gsl_matrix_set(X, count, j, (*learn_matrix_)[d][j]);
      }
      count++;
  }
  
  
  /* DEBUG code to print out matrix

  cerr << endl << "==============" << endl;
  for (int row = 0; row < count; row++) {
    cerr << gsl_vector_get(y, row);
      for (int col = 0; col < 41; col++) {
      cerr << gsl_matrix_get(X, row, col) << ' ';
    }
    cerr << endl;
  }
  */
  
  gsl_vector* c = gsl_vector_calloc(paramcount_);
  gsl_vector_set_zero(c);
  gsl_matrix *cov = gsl_matrix_calloc(paramcount_, paramcount_);
  gsl_matrix_set_zero(cov);
  double chisq;
  double tolerance = 0.01;
  size_t rank = 0;

  int err = gsl_multifit_wlinear_svd(X, w, y, tolerance, &rank, c, cov, &chisq, rt_fit_);
  int newNumGood = 0;
  
  //DDS: Outlier removal here
  
  gsl_vector *dis =  gsl_vector_calloc(numGood);
  int ix = 0;
  double dis_mean, dis_sd;

  
  //DDS: Find distance mean and stddev

  for (int d = 0; d < (int)learn_matrix_->size(); d++) {

      double predictedRT;
      double predictedDev;
      gsl_vector *x = gsl_vector_calloc(paramcount_);

      for (int j = 0; j < paramcount_; j++) {
	gsl_vector_set(x, j, (*learn_matrix_)[d][j]);
      }
      err = gsl_multifit_linear_est(x, c, cov, &predictedRT, &predictedDev);

      gsl_vector_set(dis, ix++, (learn_rts_)[d] - predictedRT);

      gsl_vector_free(x);
      
  }
  
  dis_mean = gsl_stats_wmean(gsl_vector_ptr(w, 0), 1,gsl_vector_ptr(dis, 0), 1, numGood);
  run_RT_mean_ = dis_mean;
  dis_sd = gsl_stats_wsd(gsl_vector_ptr(w, 0), 1,gsl_vector_ptr(dis, 0), 1, numGood);
  run_RT_stddev_ = dis_sd;


  count = 0;
  
  for (int d = 0; d < (int)learn_matrix_->size(); d++) {

      double predictedRT;
      double predictedDev;
      gsl_vector *x = gsl_vector_calloc(paramcount_);

      for (int j = 0; j < paramcount_; j++) {
	gsl_vector_set(x, j, (*learn_matrix_)[d][j]);
      }
      err = gsl_multifit_linear_est(x, c, cov, &predictedRT, &predictedDev);

      if (fabs(((learn_rts_)[d] - predictedRT - dis_mean)/dis_sd) < numSD) {
	gsl_vector_set(y, count, (*learn_matrix_)[d][paramcount_]);
	//gsl_vector_set(w, count, (*learn_probs_)[d]);
	gsl_vector_set(w, count, 1);
	for(int j = 0; j < paramcount_; j++) {
	  gsl_matrix_set(X, count, j, (*learn_matrix_)[d][j]);
	}
	count++;
	
      }
      gsl_vector_free(x);
      
  }
  
  gsl_vector_free(y);
  gsl_vector_free(w);
  gsl_matrix_free(X);
  
  numGood = count;
  run_RT_used_count_ = count;
  count = 0;
  y = gsl_vector_calloc(numGood);
  w = gsl_vector_calloc(numGood);
  X = gsl_matrix_calloc(numGood, paramcount_); 
   
  for (int d = 0; d < (int)learn_matrix_->size(); d++) {

    double predictedRT;
    double predictedDev;
    gsl_vector *x = gsl_vector_calloc(paramcount_);
    
    for (int j = 0; j < paramcount_; j++) {
      gsl_vector_set(x, j, (*learn_matrix_)[d][j]);
    }
    err = gsl_multifit_linear_est(x, c, cov, &predictedRT, &predictedDev);
    
    if (fabs(((learn_rts_)[d] - predictedRT - dis_mean)/dis_sd) < numSD) {
      gsl_vector_set(y, count, (*learn_matrix_)[d][paramcount_]);
      //      gsl_vector_set(w, count, (*learn_probs_)[d]);
      gsl_vector_set(w, count, 1);
      for (int j = 0; j < paramcount_; j++) {
	gsl_matrix_set(X, count, j, (*learn_matrix_)[d][j]);
      }
      count++;
    }
    
    gsl_vector_free(x);
  }    
  

  gsl_vector_free(dis);

  cerr << "After outlier removal, fitting model to " << numGood << " peptides." << endl;
  
  if (numGood < 80) {
    cerr << "Not enough good identifications to fit retention time model." << endl;
    cerr << "Reverting to SSRCalc 1.0 predictions." << endl;
    return (false);
  }
  
  //  cerr << "RANK=" << rank << endl;
  //  cerr << "DONE FITTING" << endl;
  if (rt_fit_ != NULL )gsl_multifit_linear_free(rt_fit_);
  rt_fit_ = gsl_multifit_linear_alloc(numGood, paramcount_); 
  err = gsl_multifit_wlinear_svd(X, w, y,  tolerance, &rank, c, cov, &chisq, rt_fit_);

 for (map<string, int>::iterator pp = param_index_.begin(); pp != param_index_.end(); pp++) {
    coeff_[pp->first] = gsl_vector_get(c, pp->second);
  }
  
  double tss = (count - 1) * gsl_stats_wvariance(gsl_vector_ptr(w, 0), 1, gsl_vector_ptr(y, 0), 1, count);
  r_sq_ = 1 - (chisq / tss);
    
  RTs_->reserve((int)rts_.size());
  run_RT_count_=0;
  run_RT_sum_=0;
  learned_RT_count_=0;
  learned_RT_sum_=0;
  // now calculate the "predicted" RTs from this regressed model
  for (int d = 0; d < (int)rts_.size(); d++) {
    double predictedRT;
    double predictedDev;
    gsl_vector *x = gsl_vector_calloc(paramcount_);
    for (int j = 0; j < paramcount_; j++) {
      gsl_vector_set(x, j, (*learn_matrix_)[d][j]);
    }
      
    err = gsl_multifit_linear_est(x, c, cov, &predictedRT, &predictedDev);
  
    run_RT_count_++;
    run_RT_sum_ += predictedRT;
    
    RTs_->insertAtEnd(predictedRT);
  
    gsl_vector_free(x);
  }

  // now calculate the "predicted" RTs from this regressed model
  for (int d = 0; d < (int)learn_rts_.size(); d++) {
    double predictedRT;
    double predictedDev;
    gsl_vector *x = gsl_vector_calloc(paramcount_);
    for (int j = 0; j < paramcount_; j++) {
      gsl_vector_set(x, j, (*learn_matrix_)[d][j]);
    }
      
    err = gsl_multifit_linear_est(x, c, cov, &predictedRT, &predictedDev);
  
    learned_RT_count_++;
    learned_RT_sum_ += predictedRT;
    
    learned_RTs_->insertAtEnd(predictedRT);

  
    gsl_vector_free(x);
  }
    
  //  gsl_vector_free(c);
  c_ = c;
  gsl_matrix_free(cov);
  gsl_vector_free(y);
  gsl_matrix_free(X);
// gsl_multifit_linear_free(work);
  
//  plotRegressionFits(probs_, min_prob, ntts_, min_ntt);
  
  if (r_sq_ < 0.6) {
    // too much scatter to trust this?
    return (false);
  }
#endif  


  return (true);
}


bool RTCalculator::learnNeuralNet() {
  
  cerr << "Performing neural net training of retention times..." << endl;
  
  // The sequence vs RT relationship is learned by linear regression, rather
  // than predicted by SSRCalc or Hydro. This is more robust and generally applicable.
  // Note that RT predictors like SSRCalc or Hydro makes assumptions about the column used,
  // components and gradient of mobile phase, etc, and cannot deal with modified amino acids.
  
  int numDataPoints = (int)learn_rts_.size();

  int numGood = numDataPoints;
  

  //adjustRT();

  cerr << "Fitting model to " << numGood << "  peptides ." << endl;
  
  if (numGood < 80) {
    cerr << "Not enough good identifications to fit retention time model." << endl;
    cerr << "Reverting to SSRCalc 1.0 predictions." << endl;
    return (false);
  }
    
  //map<string, int> param_index;
  param_index_.clear();
  vector<string> params;
  paramcount_ = 0;

  // add the constant
  //  param_index_["1"] = paramcount_++;
  //params.push_back("1");

 
  // add the log of the peptide length as one of the modeled variables
  param_index_["LEN"] = paramcount_++;
  params.push_back("LEN");

  param_index_["pI"] = paramcount_++;
  params.push_back("pI");

   param_index_["helicity1"] = paramcount_++;
   params.push_back("helicity1");
  
   param_index_["clusterness"] = paramcount_++;
   params.push_back("clusterness");

   
   //param_index_["helicity2"] = paramcount_++;
   //params.push_back("helicity2");

   param_index_["helectric"] = paramcount_++;
   params.push_back("helectric");

   param_index_["hydroMoment"] = paramcount_++;
   params.push_back("hydroMoment");

  
  param_index_["EludeAmphipathicityHelixMin"] = paramcount_++;
  params.push_back("EludeAmphipathicityHelixMin");


  param_index_["EludeAmphipathicityHelixMax"] = paramcount_++;
  params.push_back("EludeAmphipathicityHelixMax");


  param_index_["EludeHydrophobicMoment180Min"] = paramcount_++;
  params.push_back("EludeHydrophobicMoment180Min");


  param_index_["EludeHydrophobicMoment180Max"] = paramcount_++;
  params.push_back("EludeHydrophobicMoment180Max");


  param_index_["EludeHydrophobicMoment100Min"] = paramcount_++;
  params.push_back("EludeHydrophobicMoment100Min");


  param_index_["EludeHydrophobicMoment100Max"] = paramcount_++;
  params.push_back("EludeHydrophobicMoment100Max");



  //  param_index_["smallness"] = paramcount_++;
  //params.push_back("smallness");


  string::size_type pos;
  double maxRT = max_rt_;
  double minRT = min_rt_;
  // find out what amino acids and modified amino acids are present in the identifications
  for (int dd = 0; dd < (int)learn_peps_.size(); dd++) {
    string pep =  learn_peps_[dd];
    //if (learn_rts_[dd] > maxRT) {
    // maxRT = learn_rts_[dd];
    //}
    int totAA = countAA(pep);;
    pos = 0;
    int numAA = 0;
    string aatag, nextaa;
    while (pos != string::npos || numAA == totAA-1) {
      
      if (numAA == 0) {
	aatag = nextAAToken(pep, pos, pos);
      }
      else {
	aatag = nextaa;
      }
      

      if (numAA < totAA-1) {
 	nextaa =  nextAAToken(pep, pos, pos);
	
      }

      if (param_index_.find(aatag) == param_index_.end()) {
	param_index_[aatag] = paramcount_;
	paramcount_++;
	params.push_back(aatag);
      }

      
      if (aatag.find_first_of("nc") == string::npos && (numAA < 3 || numAA >= totAA-3)) {
	//aatag = "term_" + aatag;
	if (numAA == 0 ) {
	  aatag = "termN1_" + aatag;

	}
	if (numAA == totAA-1) { 
	  aatag = "termC1_" + aatag;


	}
	if (numAA == 1) {
	  aatag = "termN2_" + aatag;
	}
	if ( numAA == totAA-2)  {
	  aatag = "termC2_" + aatag;

	}
	
	if (numAA == 2 ) {
	  aatag = "termN3_" + aatag;
	}
	if (numAA == totAA-3)  {
	  aatag = "termC3_" + aatag;
	}
	//if (numAA<3) {
	//  aatag = "Nterm_" + aatag;
	//}
	//else {
	//  aatag = "Cterm_" + aatag;
	//}

	if (param_index_.find(aatag) == param_index_.end()) {
	  param_index_[aatag] = paramcount_;
	  paramcount_++;
	  params.push_back(aatag);
	}	
	
      }

  


      numAA++;

	
      
    }
  }
  
  double score = 1.;


  //Create a training data file


  map<string, int>::iterator found;

  ofstream trout("train.dat");
 
  trout << learn_peps_.size() <<  " " << paramcount_ << " 1" << endl; 


    // fill in the matrix -- this will be the input to the linear regression algorithm
  for (int d = 0; d < (int)learn_peps_.size(); d++) {
    pos = 0;
    vector<double> x(paramcount_ + 1, 0.0);
    
    string pep =  learn_peps_[d];
    
    int totAA = countAA(pep);
    
    // count numbers of each type of amino acids
    
    int numAA = 0;
    string aatag, nextaa;
    while (pos != string::npos || numAA == totAA-1) {
      
      if (numAA == 0) {
	aatag = nextAAToken(pep, pos, pos);
      }
      else {
	aatag = nextaa;
      }
      
      score = 1.;
      if (numAA < totAA-1) {
	nextaa =  nextAAToken(pep, pos, pos);
	
      }
      found = param_index_.find(aatag);
      if (found != param_index_.end()) {
	x[found->second]+=score;
      }     
      
      if (aatag.find_first_of("nc") == string::npos && (numAA < 3 || numAA >= totAA-3)) {
	 score = 1;
	if (numAA == 0 ) {
	  aatag = "termN1_" + aatag;

	}
	if (numAA == totAA-1) { 
	  aatag = "termC1_" + aatag;


	}
	if (numAA == 1) {
	  aatag = "termN2_" + aatag;
	}
	if ( numAA == totAA-2)  {
	  aatag = "termC2_" + aatag;

	}
	
	if (numAA == 2 ) {
	  aatag = "termN3_" + aatag;
	}
	if (numAA == totAA-3)  {
	  aatag = "termC3_" + aatag;
	}
	
	found = param_index_.find(aatag);
	if (found != param_index_.end()) {
	  x[found->second]+=score;
	}     
      }
      
      
      numAA++;
      
    }
    
    
    // found = param_index_.find("1");
    //if (found != param_index_.end()) {
    //  x[found->second]=1.;
    //}   
    
    found = param_index_.find("LEN");
    if (found != param_index_.end()) {
      x[found->second]=(double)numAA;
    }   
    
    found = param_index_.find("pI");
    if (found != param_index_.end()) {
      //x[found->second]=log(pIcalc_->ModifiedPeptide_pI(pep.c_str(),  0, pep.length()));  
      x[found->second]=isoelectric(pep);
      
    }   
    
    
    found = param_index_.find("helicity1");
    if (found != param_index_.end()) {
      x[found->second]=helicity1(pep);
    }   
    
    found = param_index_.find("clusterness");
    if (found != param_index_.end()) {
      x[found->second]=clusterness(pep);
    }   
    
    found = param_index_.find("helicity2");
    if (found != param_index_.end()) {
      x[found->second]=helicity2(pep);
    }   
    
    found = param_index_.find("helectric");
    if (found != param_index_.end()) {
      x[found->second]=helectric(pep);
    }   
    
    found = param_index_.find("hydroMoment");
    if (found != param_index_.end()) {
      x[found->second]=hydroMoment(pep);
    }   
  
    found = param_index_.find("EludeAmphipathicityHelixMin");
    if (found != param_index_.end()) {
      EludeAmphipathicityHelix(pep, &x);
    }   
    
    found = param_index_.find("EludeHydrophobicMoment100Min");
    if (found != param_index_.end()) {
      EludeHydrophobicMoments(pep, &x);
    }   
      
    found = param_index_.find("EludeIndexPartialSum5Max");
    if (found != param_index_.end()) {
      EludeIndexPartialSum(pep, &x);
    }   

    found = param_index_.find("EludeIndexSum");
    if (found != param_index_.end()) {
      x[found->second]=EludeIndexSum(pep);
    }   
    found = param_index_.find("EludeIndexAvg");
    if (found != param_index_.end()) {
      x[found->second]=EludeIndexAvg(pep);
    }   
    found = param_index_.find("EludeIndexN");
    if (found != param_index_.end()) {
      x[found->second]=EludeIndexN(pep);
    }   
    found = param_index_.find("EludeIndexC");
    if (found != param_index_.end()) {
      x[found->second]=EludeIndexC(pep);
    }   
    found = param_index_.find("EludeIndexNearestNeighborPos");
    if (found != param_index_.end()) {
      x[found->second]=EludeIndexNearestNeighborPos(pep);
    }   
    found = param_index_.find("EludeIndexNearestNeighborNeg");
    if (found != param_index_.end()) {
      x[found->second]=EludeIndexNearestNeighborNeg(pep);
    } 
    found = param_index_.find("EludeIndexSumSqrDiff");
    if (found != param_index_.end()) {
      x[found->second]=EludeIndexSumSqrDiff(pep);
    }   

 
    x[paramcount_] = learn_rts_[d]/(maxRT-minRT) + (minRT / (minRT-maxRT));
    
    
    for (int i = 0; i < (int)x.size()-1; i++) {
      trout <<  x[i] << " ";

    }
    trout  << "\n" <<  x[paramcount_] << endl;
  }
  trout.close();

  const unsigned int num_input = paramcount_;
  const unsigned int num_output = 1;
  const unsigned int num_layers = 3;
  unsigned int layers[num_layers] = {paramcount_, 24, 1};
  const float desired_error = (const float) 0.0005;
  const unsigned int max_epochs = 10000;
  const unsigned int epochs_between_reports = 1000;
  

  ann_ = fann_create_standard_array(num_layers, layers);
  
  fann_set_activation_function_hidden(ann_, FANN_SIGMOID_SYMMETRIC);
  fann_set_activation_function_output(ann_, FANN_SIGMOID);
  fann_set_training_algorithm(ann_, FANN_TRAIN_QUICKPROP);
  fann_train_on_file(ann_, "train.dat", max_epochs, epochs_between_reports, desired_error);
  std::string pep;
  fann_type *input = new fann_type[paramcount_];
  //test on training data
  for (int d = 0; d < (int)learn_peps_.size(); d++) {
    
    for (int ii=0; ii<paramcount_; ii++) {
      input[ii]=0;
    }
    pep =  learn_peps_[d];
    
    int totAA = countAA(pep);
    
    // count numbers of each type of amino acids
    pos = (unsigned int) 0;
    int numAA = 0;
    string aatag, nextaa;
    while (pos != string::npos || numAA == totAA-1) {
      
      if (numAA == 0) {
	aatag = nextAAToken(pep, pos, pos);
      }
      else {
	aatag = nextaa;
      }
      
      score = 1.;
      if (numAA < totAA-1) {
	nextaa =  nextAAToken(pep, pos, pos);
	
      }
      found = param_index_.find(aatag);
      if (found != param_index_.end()) {
	input[found->second]+=score;
      }     
      
      if (aatag.find_first_of("nc") == string::npos && (numAA < 3 || numAA >= totAA-3)) {
	 score = 1;
	if (numAA == 0 ) {
	  aatag = "termN1_" + aatag;

	}
	if (numAA == totAA-1) { 
	  aatag = "termC1_" + aatag;


	}
	if (numAA == 1) {
	  aatag = "termN2_" + aatag;
	}
	if ( numAA == totAA-2)  {
	  aatag = "termC2_" + aatag;

	}
	
	if (numAA == 2 ) {
	  aatag = "termN3_" + aatag;
	}
	if (numAA == totAA-3)  {
	  aatag = "termC3_" + aatag;
	}
	
	
	found = param_index_.find(aatag);
	if (found != param_index_.end()) {
	  input[found->second]+=score;
	}     
      }
      
      
      numAA++;
      
    }
    
    

    found = param_index_.find("LEN");
    if (found != param_index_.end()) {
      input[found->second]=(double)numAA;
    }   
    
    found = param_index_.find("pI");
    if (found != param_index_.end()) {
      //input[found->second]=log(pIcalc_->ModifiedPeptide_pI(pep.c_str(),  0, pep.length()));  
      input[found->second]=isoelectric(pep);
      
    }   
    
    
    found = param_index_.find("helicity1");
    if (found != param_index_.end()) {
      input[found->second]=helicity1(pep);
    }   
    
    found = param_index_.find("clusterness");
    if (found != param_index_.end()) {
      input[found->second]=clusterness(pep);
    }   
    
    found = param_index_.find("helicity2");
    if (found != param_index_.end()) {
      input[found->second]=helicity2(pep);
    }   
    
    found = param_index_.find("helectric");
    if (found != param_index_.end()) {
      input[found->second]=helectric(pep);
    }   
    
    found = param_index_.find("hydroMoment");
    if (found != param_index_.end()) {
      input[found->second]=hydroMoment(pep);
    }   
    
    found = param_index_.find("EludeAmphipathicityHelixMin");
    if (found != param_index_.end()) {
      EludeAmphipathicityHelix(pep, input);
    }   
    
    found = param_index_.find("EludeHydrophobicMoment100Min");
    if (found != param_index_.end()) {
      EludeHydrophobicMoments(pep, input);
    }   

    found = param_index_.find("EludeIndexPartialSum5Max");
    if (found != param_index_.end()) {
      EludeIndexPartialSum(pep, input);
    }   
    found = param_index_.find("EludeIndexSum");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexSum(pep);
    }   
    found = param_index_.find("EludeIndexAvg");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexAvg(pep);
    }   
    found = param_index_.find("EludeIndexN");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexN(pep);
    }   
    found = param_index_.find("EludeIndexC");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexC(pep);
    }   
    found = param_index_.find("EludeIndexNearestNeighborPos");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexNearestNeighborPos(pep);
    }   
    found = param_index_.find("EludeIndexNearestNeighborNeg");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexNearestNeighborNeg(pep);
    } 
    found = param_index_.find("EludeIndexSumSqrDiff");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexSumSqrDiff(pep);
    }   

 
    fann_type *calc_out;

    calc_out = fann_run(ann_, input);
      
    learned_RTs_->insertAtEnd((calc_out[0]-(minRT/(minRT-maxRT)))*(maxRT-minRT));   



  }
    




  //fann_save(ann_, "train.net");
  
  //fann_destroy(ann_);
  delete[] input;
  return true;
}


bool RTCalculator::learnNeuralNet2() {
  
  cerr << "Performing neural net training of retention times..." << endl;
  
  // The sequence vs RT relationship is learned by linear regression, rather
  // than predicted by SSRCalc or Hydro. This is more robust and generally applicable.
  // Note that RT predictors like SSRCalc or Hydro makes assumptions about the column used,
  // components and gradient of mobile phase, etc, and cannot deal with modified amino acids.
  
  int numDataPoints = (int)learn_rts_.size();

  int numGood = numDataPoints;
  

  //adjustRT();

  cerr << "Fitting model to " << numGood << "  peptides ." << endl;
  
  if (numGood < 80) {
    cerr << "Not enough good identifications to fit retention time model." << endl;
    cerr << "Reverting to SSRCalc 1.0 predictions." << endl;
    return (false);
  }
    
  //map<string, int> param_index;
  param_index_.clear();

  aa_index_.clear();
  vector<string> params;
  paramcount_ = 0;
  aacount_ = 0;

  //  param_index_["smallness"] = paramcount_++;
  //params.push_back("smallness");
 
  string::size_type pos;
  double maxRT = max_rt_;
  double minRT = min_rt_;
  // find out what amino acids and modified amino acids are present in the identifications
  for (int dd = 0; dd < (int)learn_peps_.size(); dd++) {
    string pep =  learn_peps_[dd];
    //    if (learn_rts_[dd] > maxRT) {
    //   maxRT = learn_rts_[dd];
    //}
    int totAA = countAA(pep);;
    pos = 0;
    int numAA = 0;
    string aatag, nextaa;
    while (pos != string::npos || numAA == totAA-1) {
      
      if (numAA == 0) {
	aatag = nextAAToken(pep, pos, pos);
      }
      else {
	aatag = nextaa;
      }
      

      if (numAA < totAA-1) {
 	nextaa =  nextAAToken(pep, pos, pos);
	
      }

      if (aa_index_.find(aatag) == aa_index_.end()) {
	aa_index_[aatag] = aacount_;
	aacount_++;
	params.push_back(aatag);
      }
      numAA++;

	
      
    }
  }

  paramcount_ += aacount_*max_len_;

  param_index_["LEN"] = paramcount_++;
  params.push_back("LEN");

  param_index_["pI"] = paramcount_++;
  params.push_back("pI");

  param_index_["helicity1"] = paramcount_++;
  params.push_back("helicity1");
  
  param_index_["clusterness"] = paramcount_++;
  params.push_back("clusterness");

   
  //param_index_["helicity2"] = paramcount_++;
  //params.push_back("helicity2");
  
  param_index_["helectric"] = paramcount_++;
  params.push_back("helectric");
 
  param_index_["hydroMoment"] = paramcount_++;
  params.push_back("hydroMoment");
  
  param_index_["EludeAmphipathicityHelixMin"] = paramcount_++;
  params.push_back("EludeAmphipathicityHelixMin");


  param_index_["EludeAmphipathicityHelixMax"] = paramcount_++;
  params.push_back("EludeAmphipathicityHelixMax");


  param_index_["EludeHydrophobicMoment180Min"] = paramcount_++;
  params.push_back("EludeHydrophobicMoment180Min");


  param_index_["EludeHydrophobicMoment180Max"] = paramcount_++;
  params.push_back("EludeHydrophobicMoment180Max");


  param_index_["EludeHydrophobicMoment100Min"] = paramcount_++;
  params.push_back("EludeHydrophobicMoment100Min");


  param_index_["EludeHydrophobicMoment100Max"] = paramcount_++;
  params.push_back("EludeHydrophobicMoment100Max");


  //param_index_["EludeIndexPartialSum5Max"] = paramcount_++;
  //params.push_back("EludeIndexPartialSum5Max");

  //param_index_["EludeIndexPartialSum2Max"] = paramcount_++;
  //params.push_back("EludeIndexPartialSum2Max");

  //param_index_["EludeIndexPartialSum5Min"] = paramcount_++;
  //params.push_back("EludeIndexPartialSum5Min");

  //param_index_["EludeIndexPartialSum2Min"] = paramcount_++;
  //params.push_back("EludeIndexPartialSum2Min");

  //param_index_["EludeIndexSum"] = paramcount_++;
  //params.push_back("EludeIndexSum");

  //param_index_["EludeIndexAvg"] = paramcount_++;
  //params.push_back("EludeIndexAvg");

  //param_index_["EludeIndexN"] = paramcount_++;
  //params.push_back("EludeIndexN");

  //param_index_["EludeIndexC"] = paramcount_++;
  //params.push_back("EludeIndexC");

  //param_index_["EludeIndexNearestNeighborPos"] = paramcount_++;
  //params.push_back("EludeIndexNearestNeighborPos");

  //param_index_["EludeIndexNearestNeighborNeg"] = paramcount_++;
  //params.push_back("EludeIndexNearestNeighborNeg");

  //param_index_["EludeIndexSumSqrDiff"] = paramcount_++;
  //params.push_back("EludeIndexSumSqrDiff");

  double score = 1.;


  //Create a training data file


  map<string, int>::iterator found;

  ofstream trout;

  trout.open("train.dat");
 

  trout << learn_peps_.size() <<  " " << paramcount_ << " 1" << endl; 

  vector<double>* x = new vector<double>(paramcount_ + 1, 0.0);
    // fill in the matrix -- this will be the input to the linear regression algorithm
  for (unsigned int d = 0; d < learn_peps_.size(); d++) {
    pos = 0;

    for (int dd = 0; dd<(int)(*x).size(); dd++) {
      (*x)[dd] = 0;
    }
    string pep =  learn_peps_[d];
    
    int totAA = countAA(pep);
    
    // count positions of each type of amino acids
    
    int numAA = 0;
    string aatag, nextaa;
    while (pos != string::npos || numAA == totAA-1) {
      
      if (numAA == 0) {
	aatag = nextAAToken(pep, pos, pos);
      }
      else {
	aatag = nextaa;
      }
      
      score = 1.;
      if (numAA < totAA-1) {
	nextaa =  nextAAToken(pep, pos, pos);
	
      }
      int pep_pos = numAA;
      if (numAA < totAA / 2) {
	pep_pos=numAA;
      }
      else {
	pep_pos=max_len_-(totAA - numAA);
      }
      found = aa_index_.find(aatag);
      
      if (found != aa_index_.end()) {
	int aa_i = found->second;
	for (int ii=0; ii<aacount_; ii++) {
	
	  if (ii == aa_i) {
	    (*x)[pep_pos*aacount_+ii]= 1;
	  }
	  else {
	    (*x)[pep_pos*aacount_+ii]= 0;
	  }
	    
	}     
      }      
      
      numAA++;
      
    }

    found = param_index_.find("LEN");
    if (found != param_index_.end()) {
      (*x)[found->second]=(double)numAA;
    }   
    
    found = param_index_.find("pI");
    if (found != param_index_.end()) {
      (*x)[found->second]=pIcalc_->ModifiedPeptide_pI(pep.c_str(),  0, getPeptideLength(pep));  
      //x[found->second]=isoelectric(pep);
      
    }   
    
    
    found = param_index_.find("helicity1");
    if (found != param_index_.end()) {
      (*x)[found->second]=helicity1(pep);
    }   
    
    found = param_index_.find("clusterness");
    if (found != param_index_.end()) {
      (*x)[found->second]=clusterness(pep);
    }   
    
    found = param_index_.find("helicity2");
    if (found != param_index_.end()) {
      (*x)[found->second]=helicity2(pep);
    }   
    
    found = param_index_.find("helectric");
    if (found != param_index_.end()) {
      (*x)[found->second]=helectric(pep);
    }   
    
    found = param_index_.find("hydroMoment");
    if (found != param_index_.end()) {
      (*x)[found->second]=hydroMoment(pep);
    }   
  
    found = param_index_.find("EludeAmphipathicityHelixMin");
    if (found != param_index_.end()) {
      EludeAmphipathicityHelix(pep, x);
    }   
    
    found = param_index_.find("EludeHydrophobicMoment100Min");
    if (found != param_index_.end()) {
      EludeHydrophobicMoments(pep, x);
    }   

    found = param_index_.find("EludeIndexPartialSum5Max");
    if (found != param_index_.end()) {
      EludeIndexPartialSum(pep, x);
    }   

    found = param_index_.find("EludeIndexSum");
    if (found != param_index_.end()) {
      (*x)[found->second]=EludeIndexSum(pep);
    }   
    found = param_index_.find("EludeIndexAvg");
    if (found != param_index_.end()) {
      (*x)[found->second]=EludeIndexAvg(pep);
    }   
    found = param_index_.find("EludeIndexN");
    if (found != param_index_.end()) {
      (*x)[found->second]=EludeIndexN(pep);
    }   
    found = param_index_.find("EludeIndexC");
    if (found != param_index_.end()) {
      (*x)[found->second]=EludeIndexC(pep);
    }   
    found = param_index_.find("EludeIndexNearestNeighborPos");
    if (found != param_index_.end()) {
      (*x)[found->second]=EludeIndexNearestNeighborPos(pep);
    }   
    found = param_index_.find("EludeIndexNearestNeighborNeg");
    if (found != param_index_.end()) {
      (*x)[found->second]=EludeIndexNearestNeighborNeg(pep);
    } 
    found = param_index_.find("EludeIndexSumSqrDiff");
    if (found != param_index_.end()) {
      (*x)[found->second]=EludeIndexSumSqrDiff(pep);
    } 

    // found = param_index_.find("1");
    //if (found != param_index_.end()) {
    //  (*x)[found->second]=1.;
    //}   
    
  
    
    
  (*x)[paramcount_] = learn_rts_[d]/(maxRT-minRT) + (minRT / (minRT-maxRT));
    
    
    for (int i = 0; i < (int)(*x).size()-1; i++) {
      trout <<  (*x)[i] << " ";

    }
    trout  << "\n" <<  (*x)[paramcount_] << endl;
  }
  trout.close();
  delete x;

  const unsigned int num_input = paramcount_;
  const unsigned int num_output = 1;
  const unsigned int num_layers = 3;
  unsigned int layers[num_layers] = {paramcount_, 24, 1};
  const float desired_error = (const float) 0.001;
  const unsigned int max_epochs = 500;
  const unsigned int epochs_between_reports = 100;
  

  ann_ = fann_create_standard_array(num_layers, layers);
  
  fann_set_activation_function_hidden(ann_, FANN_SIGMOID_SYMMETRIC);
  fann_set_activation_function_output(ann_, FANN_SIGMOID);
  
  fann_train_on_file(ann_, "train.dat", max_epochs, epochs_between_reports, desired_error);
  std::string pep;
  fann_type *input = new fann_type[paramcount_];
  //test on training data
  for (int d = 0; d < (int)learn_peps_.size(); d++) {
    
    for (int ii=0; ii<paramcount_; ii++) {
      input[ii]=0;
    }
    pep =  learn_peps_[d];
    
    int totAA = countAA(pep);
    
    // count numbers of each type of amino acids
    pos = (unsigned int) 0;
    int numAA = 0;
    string aatag, nextaa;
    while (pos != string::npos || numAA == totAA-1) {
      
      if (numAA == 0) {
	aatag = nextAAToken(pep, pos, pos);
      }
      else {
	aatag = nextaa;
      }
      
      score = 1.;
      if (numAA < totAA-1) {
	nextaa =  nextAAToken(pep, pos, pos);
	
      }
      found = aa_index_.find(aatag);
      int pep_pos = numAA;
      if (numAA < totAA / 2) {
	pep_pos = numAA;
      }
      else {
	pep_pos = max_len_ - (totAA - numAA);
      }
      if (found != aa_index_.end()) {
		int aa_i = found->second;
	for (int ii=0; ii<aacount_; ii++) {
	
	  if (ii == aa_i) {
	    input[pep_pos*aacount_+ii]= 1;
	  }
	  else {
	    input[pep_pos*aacount_+ii]= 0;
	  }
	    
	}     
      }     
      
      
      
      numAA++;
      
    }
    
    
    found = param_index_.find("LEN");
    if (found != param_index_.end()) {
      input[found->second]=(double)numAA;
    }   
    
    found = param_index_.find("pI");
    if (found != param_index_.end()) {
      input[found->second]=pIcalc_->ModifiedPeptide_pI(pep.c_str(),  0, getPeptideLength(pep));  
      //      input[found->second]=isoelectric(pep);
      
    }   
    
    
    found = param_index_.find("helicity1");
    if (found != param_index_.end()) {
      input[found->second]=helicity1(pep);
    }   
    
    found = param_index_.find("clusterness");
    if (found != param_index_.end()) {
      input[found->second]=clusterness(pep);
    }   
    
    found = param_index_.find("helicity2");
    if (found != param_index_.end()) {
      input[found->second]=helicity2(pep);
    }   
    
    found = param_index_.find("helectric");
    if (found != param_index_.end()) {
      input[found->second]=helectric(pep);
    }   

  
    found = param_index_.find("hydroMoment");
    if (found != param_index_.end()) {
      input[found->second]=hydroMoment(pep);
    }   



    found = param_index_.find("EludeAmphipathicityHelixMin");
    if (found != param_index_.end()) {
      EludeAmphipathicityHelix(pep, input);
    }   
    
    found = param_index_.find("EludeHydrophobicMoment100Min");
    if (found != param_index_.end()) {
      EludeHydrophobicMoments(pep, input);
    }   

    found = param_index_.find("EludeIndexPartialSum5Max");
    if (found != param_index_.end()) {
      EludeIndexPartialSum(pep, input);
    }   

    found = param_index_.find("EludeIndexSum");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexSum(pep);
    }   
    found = param_index_.find("EludeIndexAvg");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexAvg(pep);
    }   
    found = param_index_.find("EludeIndexN");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexN(pep);
    }   
    found = param_index_.find("EludeIndexC");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexC(pep);
    }   
    found = param_index_.find("EludeIndexNearestNeighborPos");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexNearestNeighborPos(pep);
    }   
    found = param_index_.find("EludeIndexNearestNeighborNeg");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexNearestNeighborNeg(pep);
    } 
    found = param_index_.find("EludeIndexSumSqrDiff");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexSumSqrDiff(pep);
    } 
    fann_type *calc_out;

    calc_out = fann_run(ann_, input);
      
    learned_RTs_->insertAtEnd((calc_out[0]-(minRT/(minRT-maxRT)))*(maxRT-minRT)); 



  }
    




  //fann_save(ann_, "train.net");
  
  //fann_destroy(ann_);
  delete[] input;
  return true;
}
bool RTCalculator::learnNeuralNetCV() {
  
  cerr << "Performing neural net training of CV model..." << endl;
  
  // The sequence vs RT relationship is learned by linear regression, rather
  // than predicted by SSRCalc or Hydro. This is more robust and generally applicable.
  // Note that RT predictors like SSRCalc or Hydro makes assumptions about the column used,
  // components and gradient of mobile phase, etc, and cannot deal with modified amino acids.
  
  int numDataPoints = (int)learn_rts_.size();

  int numGood = numDataPoints;
  

  //adjustRT();


  cerr << "Fitting model to " << numGood << "  peptides ." << endl;
  
  if (numGood < 80) {
    cerr << "Not enough good identifications to fit CV Neural Net model." << endl;
    exit(1);
  }
    
  //map<string, int> param_index;
  param_index_.clear();

  aa_index_.clear();
  vector<string> params;
  paramcount_ = 0;
  aacount_ = 0;

  //  param_index_["smallness"] = paramcount_++;
  //params.push_back("smallness");
 
  string::size_type pos;
  double maxRT = max_rt_;
  double minRT = min_rt_;

  // find out what amino acids and modified amino acids are present in the identifications
  for (int dd = 0; dd < (int)learn_peps_.size(); dd++) {
    string pep =  learn_peps_[dd];
    //    if (learn_rts_[dd] > maxRT) {
    //   maxRT = learn_rts_[dd];
    //}
    int totAA = countAA(pep);;
    pos = 0;
    int numAA = 0;
    string aatag, nextaa;
    while (pos != string::npos || numAA == totAA-1) {
      
      if (numAA == 0) {
	aatag = nextAAToken(pep, pos, pos);
      }
      else {
	aatag = nextaa;
      }
      

      if (numAA < totAA-1) {
 	nextaa =  nextAAToken(pep, pos, pos);
	
      }

      if (aa_index_.find(aatag) == aa_index_.end()) {
	aa_index_[aatag] = aacount_;
	aacount_++;
	params.push_back(aatag);
      }
      numAA++;

	
      
    }
  }

  paramcount_ += aacount_*max_len_;

  param_index_["LEN"] = paramcount_++;
  params.push_back("LEN");

  param_index_["pI"] = paramcount_++;
  params.push_back("pI");

  param_index_["helicity1"] = paramcount_++;
  params.push_back("helicity1");
  
  param_index_["clusterness"] = paramcount_++;
  params.push_back("clusterness");

   
  //param_index_["helicity2"] = paramcount_++;
  //params.push_back("helicity2");
  
  param_index_["helectric"] = paramcount_++;
  params.push_back("helectric");
 
  param_index_["hydroMoment"] = paramcount_++;
  params.push_back("hydroMoment");
  
  param_index_["EludeAmphipathicityHelixMin"] = paramcount_++;
  params.push_back("EludeAmphipathicityHelixMin");


  param_index_["EludeAmphipathicityHelixMax"] = paramcount_++;
  params.push_back("EludeAmphipathicityHelixMax");


  param_index_["EludeHydrophobicMoment180Min"] = paramcount_++;
  params.push_back("EludeHydrophobicMoment180Min");


  param_index_["EludeHydrophobicMoment180Max"] = paramcount_++;
  params.push_back("EludeHydrophobicMoment180Max");


  param_index_["EludeHydrophobicMoment100Min"] = paramcount_++;
  params.push_back("EludeHydrophobicMoment100Min");


  param_index_["EludeHydrophobicMoment100Max"] = paramcount_++;
  params.push_back("EludeHydrophobicMoment100Max");


  //param_index_["EludeIndexPartialSum5Max"] = paramcount_++;
  //params.push_back("EludeIndexPartialSum5Max");

  //param_index_["EludeIndexPartialSum2Max"] = paramcount_++;
  //params.push_back("EludeIndexPartialSum2Max");

  //param_index_["EludeIndexPartialSum5Min"] = paramcount_++;
  //params.push_back("EludeIndexPartialSum5Min");

  //param_index_["EludeIndexPartialSum2Min"] = paramcount_++;
  //params.push_back("EludeIndexPartialSum2Min");

  //param_index_["EludeIndexSum"] = paramcount_++;
  //params.push_back("EludeIndexSum");

  //param_index_["EludeIndexAvg"] = paramcount_++;
  //params.push_back("EludeIndexAvg");

  //param_index_["EludeIndexN"] = paramcount_++;
  //params.push_back("EludeIndexN");

  //param_index_["EludeIndexC"] = paramcount_++;
  //params.push_back("EludeIndexC");

  //param_index_["EludeIndexNearestNeighborPos"] = paramcount_++;
  //params.push_back("EludeIndexNearestNeighborPos");

  //param_index_["EludeIndexNearestNeighborNeg"] = paramcount_++;
  //params.push_back("EludeIndexNearestNeighborNeg");

  //param_index_["EludeIndexSumSqrDiff"] = paramcount_++;
  //params.push_back("EludeIndexSumSqrDiff");

  double score = 1.;


  //Create a training data file


  map<string, int>::iterator found;

  ofstream trout;

  trout.open("train.dat");
 
 const unsigned int num_input = paramcount_;
  const unsigned int num_output = 1;
  const unsigned int num_layers = 3;

  unsigned int layers[num_layers] = {paramcount_, 22, num_output};
  const float desired_error = (const float) 0.001;
  const unsigned int max_epochs = 500;
  const unsigned int epochs_between_reports = 100;
  
  trout << learn_peps_.size() <<  " " << paramcount_ << " " << num_output << endl; 

  vector<double>* x = new vector<double>(paramcount_ + 1, 0.0);
    // fill in the matrix -- this will be the input to the linear regression algorithm
  for (unsigned int d = 0; d < learn_peps_.size(); d++) {
    pos = 0;

    for (int dd = 0; dd<(int)(*x).size(); dd++) {
      (*x)[dd] = 0;
    }
    string pep =  learn_peps_[d];
    
    int totAA = countAA(pep);
    
    // count positions of each type of amino acids
    
    int numAA = 0;
    string aatag, nextaa;
    while (pos != string::npos || numAA == totAA-1) {
      
      if (numAA == 0) {
	aatag = nextAAToken(pep, pos, pos);
      }
      else {
	aatag = nextaa;
      }
      
      score = 1.;
      if (numAA < totAA-1) {
	nextaa =  nextAAToken(pep, pos, pos);
	
      }
      int pep_pos = numAA;
      if (numAA < totAA / 2) {
	pep_pos=numAA;
      }
      else {
	pep_pos=max_len_-(totAA - numAA);
      }
      found = aa_index_.find(aatag);
      
      if (found != aa_index_.end()) {
	int aa_i = found->second;
	for (int ii=0; ii<aacount_; ii++) {
	
	  if (ii == aa_i) {
	    (*x)[pep_pos*aacount_+ii]= 1;
	  }
	  else {
	    (*x)[pep_pos*aacount_+ii]= 0;
	  }
	    
	}     
      }      
      
      numAA++;
      
    }

    found = param_index_.find("LEN");
    if (found != param_index_.end()) {
      (*x)[found->second]=(double)numAA;
    }   
    
    found = param_index_.find("pI");
    if (found != param_index_.end()) {
      (*x)[found->second]=pIcalc_->ModifiedPeptide_pI(pep.c_str(),  0, getPeptideLength(pep));  
      //x[found->second]=isoelectric(pep);
      
    }   
    
    
    found = param_index_.find("helicity1");
    if (found != param_index_.end()) {
      (*x)[found->second]=helicity1(pep);
    }   
    
    found = param_index_.find("clusterness");
    if (found != param_index_.end()) {
      (*x)[found->second]=clusterness(pep);
    }   
    
    found = param_index_.find("helicity2");
    if (found != param_index_.end()) {
      (*x)[found->second]=helicity2(pep);
    }   
    
    found = param_index_.find("helectric");
    if (found != param_index_.end()) {
      (*x)[found->second]=helectric(pep);
    }   
    
    found = param_index_.find("hydroMoment");
    if (found != param_index_.end()) {
      (*x)[found->second]=hydroMoment(pep);
    }   
  
    found = param_index_.find("EludeAmphipathicityHelixMin");
    if (found != param_index_.end()) {
      EludeAmphipathicityHelix(pep, x);
    }   
    
    found = param_index_.find("EludeHydrophobicMoment100Min");
    if (found != param_index_.end()) {
      EludeHydrophobicMoments(pep, x);
    }   

    found = param_index_.find("EludeIndexPartialSum5Max");
    if (found != param_index_.end()) {
      EludeIndexPartialSum(pep, x);
    }   

    found = param_index_.find("EludeIndexSum");
    if (found != param_index_.end()) {
      (*x)[found->second]=EludeIndexSum(pep);
    }   
    found = param_index_.find("EludeIndexAvg");
    if (found != param_index_.end()) {
      (*x)[found->second]=EludeIndexAvg(pep);
    }   
    found = param_index_.find("EludeIndexN");
    if (found != param_index_.end()) {
      (*x)[found->second]=EludeIndexN(pep);
    }   
    found = param_index_.find("EludeIndexC");
    if (found != param_index_.end()) {
      (*x)[found->second]=EludeIndexC(pep);
    }   
    found = param_index_.find("EludeIndexNearestNeighborPos");
    if (found != param_index_.end()) {
      (*x)[found->second]=EludeIndexNearestNeighborPos(pep);
    }   
    found = param_index_.find("EludeIndexNearestNeighborNeg");
    if (found != param_index_.end()) {
      (*x)[found->second]=EludeIndexNearestNeighborNeg(pep);
    } 
    found = param_index_.find("EludeIndexSumSqrDiff");
    if (found != param_index_.end()) {
      (*x)[found->second]=EludeIndexSumSqrDiff(pep);
    } 

    // found = param_index_.find("1");
    //if (found != param_index_.end()) {
    //  (*x)[found->second]=1.;
    //}   
    
  
    
    
    (*x)[paramcount_] = learn_rts_[d]/(maxRT-minRT) + (minRT / (minRT-maxRT));
    
    
    for (int i = 0; i < (int)(*x).size()-1; i++) {
      trout <<  (*x)[i] << " ";

    }
    trout  << "\n" <<  (*x)[paramcount_] << endl;
  }
  trout.close();
  delete x;

  

  ann_ = fann_create_standard_array(num_layers, layers);
  
  fann_set_activation_function_hidden(ann_, FANN_SIGMOID_SYMMETRIC);
  fann_set_activation_function_output(ann_, FANN_SIGMOID);

  const unsigned int max_neurons = 100;


  const unsigned int neurons_btr = 10;

  fann_set_training_algorithm(ann_, FANN_TRAIN_RPROP);
  fann_train_on_file(ann_, "train.dat", max_epochs, epochs_between_reports, desired_error);
  std::string pep;
  fann_type *input = new fann_type[paramcount_];
  //test on training data
  for (int d = 0; d < (int)learn_peps_.size(); d++) {
    
    for (int ii=0; ii<paramcount_; ii++) {
      input[ii]=0;
    }
    pep =  learn_peps_[d];
    
    int totAA = countAA(pep);
    
    // count numbers of each type of amino acids
    pos = (unsigned int) 0;
    int numAA = 0;
    string aatag, nextaa;
    while (pos != string::npos || numAA == totAA-1) {
      
      if (numAA == 0) {
	aatag = nextAAToken(pep, pos, pos);
      }
      else {
	aatag = nextaa;
      }
      
      score = 1.;
      if (numAA < totAA-1) {
	nextaa =  nextAAToken(pep, pos, pos);
	
      }
      found = aa_index_.find(aatag);
      int pep_pos = numAA;
      if (numAA < totAA / 2) {
	pep_pos = numAA;
      }
      else {
	pep_pos = max_len_ - (totAA - numAA);
      }
      if (found != aa_index_.end()) {
		int aa_i = found->second;
	for (int ii=0; ii<aacount_; ii++) {
	
	  if (ii == aa_i) {
	    input[pep_pos*aacount_+ii]= 1;
	  }
	  else {
	    input[pep_pos*aacount_+ii]= 0;
	  }
	    
	}     
      }     
      
      
      
      numAA++;
      
    }
    
    
    found = param_index_.find("LEN");
    if (found != param_index_.end()) {
      input[found->second]=(double)numAA;
    }   
    
    found = param_index_.find("pI");
    if (found != param_index_.end()) {
      input[found->second]=pIcalc_->ModifiedPeptide_pI(pep.c_str(),  0, getPeptideLength(pep));  
      //      input[found->second]=isoelectric(pep);
      
    }   
    
    
    found = param_index_.find("helicity1");
    if (found != param_index_.end()) {
      input[found->second]=helicity1(pep);
    }   
    
    found = param_index_.find("clusterness");
    if (found != param_index_.end()) {
      input[found->second]=clusterness(pep);
    }   
    
    found = param_index_.find("helicity2");
    if (found != param_index_.end()) {
      input[found->second]=helicity2(pep);
    }   
    
    found = param_index_.find("helectric");
    if (found != param_index_.end()) {
      input[found->second]=helectric(pep);
    }   

  
    found = param_index_.find("hydroMoment");
    if (found != param_index_.end()) {
      input[found->second]=hydroMoment(pep);
    }   



    found = param_index_.find("EludeAmphipathicityHelixMin");
    if (found != param_index_.end()) {
      EludeAmphipathicityHelix(pep, input);
    }   
    
    found = param_index_.find("EludeHydrophobicMoment100Min");
    if (found != param_index_.end()) {
      EludeHydrophobicMoments(pep, input);
    }   

    found = param_index_.find("EludeIndexPartialSum5Max");
    if (found != param_index_.end()) {
      EludeIndexPartialSum(pep, input);
    }   

    found = param_index_.find("EludeIndexSum");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexSum(pep);
    }   
    found = param_index_.find("EludeIndexAvg");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexAvg(pep);
    }   
    found = param_index_.find("EludeIndexN");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexN(pep);
    }   
    found = param_index_.find("EludeIndexC");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexC(pep);
    }   
    found = param_index_.find("EludeIndexNearestNeighborPos");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexNearestNeighborPos(pep);
    }   
    found = param_index_.find("EludeIndexNearestNeighborNeg");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexNearestNeighborNeg(pep);
    } 
    found = param_index_.find("EludeIndexSumSqrDiff");
    if (found != param_index_.end()) {
      input[found->second]=EludeIndexSumSqrDiff(pep);
    } 
    fann_type *calc_out;

    calc_out = fann_run(ann_, input);
      
    learned_RTs_->insertAtEnd((calc_out[0]-(minRT/(minRT-maxRT)))*(maxRT-minRT));



  }
    




  //fann_save(ann_, "train.net");
  
  //fann_destroy(ann_);
  delete[] input;
  return true;
}

double RTCalculator::computeHelixScore(string& pep) {
  double rtn = 0;
  string::size_type pos = 0;
  string aatag;
  while (pos != string::npos) {
    aatag = nextAAToken(pep, pos, pos);
    rtn += (*HelixProp_)[aatag];
  }

  return rtn;

}

double RTCalculator::heli1TermAdj(string & hc, int start, int len) {
  int i = 0 ;
  int where = -1;
    
  for (i=0; i<(int)(hc.length()-4); i++) {
    string m = hc.substr(i, 1);

    if (m.find_first_of("OU")!=string::npos) {
      where = i;
      break;
    }
  }
  where = where + start;

  if (where<2) { return .2; }
  if (where<3) { return .25; }
  if (where<4) { return .45; }
  
  if (where>len-3) { return .2; }
  if (where>len-4) { return .75; }
  if (where>len-5) { return .65; }
  
  return 1;

}

double RTCalculator::clusterness(string& pep) {
  string cstring = "";
  string aatag = "";
  string x1 = "";
  double sk, score=0., addit;
  int occurs;
  int found;
  string::size_type pos = 0;
  map<string, double>::iterator itr;
  while (pos != string::npos) {
    aatag = nextAAToken(pep, pos, pos);
     if (aatag.find_first_of( "LIW")!=string::npos) {
       cstring += "5";
     }
     else if (aatag.find_first_of( "AMYV")!=string::npos) {
       cstring += "1";
     }
     else {
       cstring += "0";
     }
  }

  cstring = "0"+cstring+"0";

  for (itr = (*ClustProp_).begin(); itr != (*ClustProp_).end(); itr++) {
    x1 = "0"+(*itr).first+"0";
    sk = (*itr).second;
    occurs = 0;
    found = -1;
    while ((found = (int)cstring.find(x1, found+1)) != (int)string::npos) {
      occurs++;
    }
    if (occurs>0)
    { 
      addit=sk*occurs;
      score+=addit; 
    }
  }

  return score;

}

//DDS: ported from SSRCalc3

double RTCalculator::smallness(string& pep) {
  
  if (pep.length() < 20) {
    return 2.;
  }
  if (pep.length() < 15) {
    return 3.;
  }

  return 1.;



}

double RTCalculator::helicity2(string& pep) {
  string bkpep = reversePeptide(pep);
  //  reverse(bkpep.begin(), bkpep.end());
  double fwHiscore, fwGscore;
  double bkHiscore, bkGscore;
  double h2FwBk;
  double HELIX2SCALE = 0.255;

  heli2Calc(pep, fwHiscore, fwGscore);
  heli2Calc(bkpep, bkHiscore, bkGscore);

  if (bkGscore>fwGscore)	
    { 
      h2FwBk = bkHiscore; 
    }
  else 	        
    { 
      h2FwBk = fwHiscore; 
    }

  double lenMult=0;
  if (getPeptideLength(pep)>30) { 
    lenMult=1; 
  }

  double noPMult=0.75;
  if (pep.find_first_of("P")!=string::npos) { 
    noPMult=0; 
  }
	
  double h2mult = 1+ lenMult + noPMult;	

  return HELIX2SCALE * h2mult * h2FwBk;
  
}

void  RTCalculator::heli2Calc(string& pep, double& Hiscore, double& Gscore) {
  string hstring, aatag, aasubstr, aastring="", best, traps;
  string::size_type pos = (string::size_type) 0;
  string::size_type pos1 = (string::size_type) 0;
  string::size_type pos2 = (string::size_type) 0;
  int llim = 1000;
  string lc;
  string pat = "";
  int i = 0, zap = 0, subt = 0, best_pos = -1;
  string f1, f2, f3;
  Hiscore = 0;
  Gscore = 0;
  if (getPeptideLength(pep) < 11) {
    return;
  }

  while (pos != string::npos) {
    aatag = nextAAToken(pep, pos, pos);
    if (aatag.find_first_of( "WFILYMVA")!=string::npos) {
      hstring += "1";
    }
    else     if (aatag.find_first_of( "GSPCNKQHRTDE")!=string::npos) {
      hstring += "0";
    }
    aastring += aatag.substr(0,1);
  }
  pos = (string::size_type) 0;

  for (i=0; i<(int)hstring.length(); i++) {
    if (hstring.substr(i, 1) == "1") {
      lc = hstring.substr(i, llim);
      zap = 0;
      subt = 0;
      
      while ( zap  <= llim && subt < 2) {
	f1 = "";
	f2 = "";
	f3 = "";
	
	if (zap < (int)lc.length()) {
	  f1 = lc.substr(zap, 1);	  
	  if (zap > 0) 
	    f2 = lc.substr(zap-1, 1);
	  if (zap < (int)(lc.length()-1))
	    f3 = lc.substr(zap+1, 1);
	}
	
	if (f1 == "1") {
	  if (zap>0) 
	    pat = pat + "--";

	  pat += aastring.substr(zap, 1);	  
	}
	else if (f2 == "1" && f1 =="0") {
	  subt++;
	  
	  if (subt < 2) 
	    pat += "->" + aastring.substr(zap-1,1);
	}
	else if (f3 == "1" && f1 == "0") 
        { 
          subt++; 
          if (subt<2) 
	    pat += "<-" + aastring.substr(zap+1,1);  
	}
	int sum =  0;
	if (f1 == "1") {
	  sum++;
	}
	if (f2 == "1") {
	  sum++;
	}
	if (f3 == "1") {
	  sum++;
	}

        if (!sum) 
	  zap=1000; 
        
	zap=zap+3;
	
      }
      if (getPeptideLength(pat) > 4) 
      { 
        string traps=aastring;
        double skore=evalH2pattern(pat,traps,i-1,string("*"));
        if (skore>=Hiscore) { 
	  Hiscore=skore; 
	  best=pat; 
	  best_pos=i; 
	} 
      }
    }
  }

	
  if (Hiscore>0) 
  { 
    Gscore=Hiscore; 
    traps=aastring;
    Hiscore=evalH2pattern(best,traps,best_pos-1,string("+"));
  }
  else
  {
    Hiscore = 0;
    Gscore = 0;
  }

  return;

}

double  RTCalculator::evalH2pattern(string& pat, string& testsq, int posn, string etype) {
  
  string f01 = pat.substr(0,1);
  int i, OFF1 = 2, acount=1, iss;
  double mult = 1., s3;
  string fpart = "", gpart = "", far1 = "", far2 = "";
  string testAAl = "", testAAr = "";


  map<string, double>::iterator itr;
  
  double prod1 = 0;
  
  if ((itr =  (*HelixProp_).find(f01+"_bsc")) != (*HelixProp_).end()) 
    prod1 = itr->second;

 

  if (OFF1+posn+2 >= (int)testsq.length()) {
    return 0;
  }
  testAAl=testsq.substr(OFF1+posn,1); 
  testAAr=testsq.substr(OFF1+posn+2,1);
  testsq=testsq.substr(OFF1+posn+1);

  mult=connector(f01,testAAl,testAAr,string("--"), far1, far2);	
   
  prod1=prod1*mult;
	
  if (etype == "*") { 
    prod1=prod1*25; 
  } 

  if (mult<0.00001) 
  { 
    return 0; 
  }

  iss =0;
  for (i = 1; i<(int)(pat.length()-2); i+=3) {
    fpart=pat.substr(i,2);
    gpart=pat.substr(i+2,1);
    s3 = 0;
    if ((itr =  (*HelixProp_).find(gpart+"_bsc")) != (*HelixProp_).end()) 
      s3 = itr->second;
    if (fpart == "--") { 
      iss=0; far1=""; far2="";
    } 
    else if (fpart == "<-" && i<(int)(testsq.length()-1)) { 
      iss=1; far1=testsq.substr(i+1,1); far2=""; 
    } 
    else if (fpart == "->" && i<(int)(testsq.length()-3)) { 
      iss=-1; far1=""; far2=testsq.substr(i+3,1); 
    }


    if (i < (int)(testsq.length() - (1+iss)) ) {
      testAAl=testsq.substr(i+1+iss,1);
    }
    else {
      testAAl="";
    }

    if (i < (int)(testsq.length() - (3+iss)) ) {
      testAAr=testsq.substr(i+3+iss,1);
    }
    else {
      testAAr="";
    }
    
    mult=connector(gpart, testAAl, testAAr, fpart, far1, far2);

    if (etype == "*")
    {
      if (mult!=0 || acount<3) { prod1=prod1*25*s3*mult; } 
    }
		
    if (etype == "+")
    {
      prod1=prod1+s3*mult;
    }
	
    if (mult<0.00001)
    {
      break;
    }
 
		
    acount++;
		
  }
  return prod1;

  

}

double RTCalculator::connector(string & acid, string & lp, string & rp, string  ct, string& far1, string& far2) {
  double mult = 1.;
  
  if (ct == "<-") { mult *= 0.2; }
  if (ct == "->") { mult *= 0.1; }

  map<string, double>::iterator itr;
  if ((itr =  (*HelixProp_).find(lp+"_cmu")) != (*HelixProp_).end()) 
    mult *= itr->second;


  if (lp != rp)  { 
    if ((itr =  (*HelixProp_).find(rp+"_cmu")) != (*HelixProp_).end()) 
      mult *= itr->second;
    
    
  }

  if (acid.find_first_of("AYVM")!=string::npos) 
  {
    if (lp.find_first_of("PG")!=string::npos || 
	rp.find_first_of("PG")!=string::npos) { 
      mult = 0; 
    }
    if (ct == "->" || ct == "<-") { 
      mult=0; 
    }
  }
  	
  if (acid.find_first_of("LWFI")!=string::npos) 
  {
    if ((lp.find_first_of("PG")!=string::npos || 
	 rp.find_first_of("PG")!=string::npos) && (ct != "--") ) { 
      mult=0; 
    }
    if ((far1.find_first_of("PG")!=string::npos || 
	 far2.find_first_of("PG")!=string::npos) && (ct == "<-" || ct == "->") ) { 
      mult=0; 
    }
  }
  return mult;
}

double RTCalculator::isoelectric(string& pep) {
  double mass =  ResidueMass::getMass('n', true) +  ResidueMass::getMass('c', true);
  
  double pi1 = pIcalc_->ModifiedPeptide_pI(pep.c_str(),  0, getPeptideLength(pep) );
  
  string::size_type pos = (string::size_type) 0;
  string aastring = "";
  string masstr = "";
  
  string::size_type lpos = (string::size_type) 0;
  string::size_type rpos = (string::size_type) 0;
  while (pos != string::npos) {
    string aatag = nextAAToken(pep, pos, pos);
       //Use mod mass [xxx.xx]
    if ((lpos=aatag.find_first_of("[")) != string::npos && 
	(rpos=aatag.find_first_of("]")) != string::npos ) {
      
      masstr = aatag.substr(lpos+1, rpos-lpos-1); 
      mass += atof(masstr.c_str());
    }
    else {
      mass += ResidueMass::getMass((aatag.c_str())[0], true);
	
    }
  }

  mass = 1.08014 * log(mass);

  double delta = pi1-29.207+mass;

  return delta;

   
}



double RTCalculator::hydroMoment(string& pep) {

  double momt = 0;


  string::size_type pos = (string::size_type) 0;
  string aastring = "";
  const char* aa;
  double h;
  int n = 0;
  while (pos != string::npos) {
    string aatag = nextAAToken(pep, pos, pos);
    aa = aatag.substr(0,1).c_str(); 
    n++;
    h = (*HydroPath_)[aa];
    momt =  pow(h * sin(2*n*3.14159/3.6), 2) + pow(h * cos(2*n*3.14159/3.6), 2);
    
    

  }

  momt = pow(momt, 0.5);

  return momt;
}
double RTCalculator::EludeIndexSum(string& peptide) {
  double rtn = 0.;
  string::size_type pos = (string::size_type) 0;
  string aatag;
  while (pos != string::npos) {
    aatag = nextAAToken(peptide, pos, pos);
    rtn += (*EludeIndex_)[aatag];
  }
  return rtn;
  
}
double RTCalculator::EludeIndexAvg(string& peptide) {
  double len = (double)getPeptideLength(peptide);
  double sum = EludeIndexSum(peptide);
  return sum / len;
}

double RTCalculator::EludeIndexN(string& peptide) {
  return (*EludeIndex_)[getAATokenIndex(peptide, 0)];
}
double RTCalculator::EludeIndexC(string& peptide) {
  return (*EludeIndex_)[getAATokenIndex(peptide, getPeptideLength(peptide)-1)];
}
double RTCalculator::EludeIndexNearestNeighborPos(string& peptide) {
  double sum = 0.0;
  unsigned int len = getPeptideLength(peptide);
  string aatag="";
  for (unsigned int ix = 0; ix < len; ++ix) {
    aatag = getAATokenIndex(peptide, ix);
    if (aatag.find_first_of("RK")!=string::npos) {
      if (ix > 0) {
	if ((*EludeIndex_)[getAATokenIndex(peptide, ix-1)] > 0)
	  sum += (*EludeIndex_)[getAATokenIndex(peptide, ix-1)];
      }
      if (ix < len - 1) {
	if ((*EludeIndex_)[getAATokenIndex(peptide, ix+1)] > 0)
	  sum += (*EludeIndex_)[getAATokenIndex(peptide, ix+1)];

      }
    }
  }
  return sum;
}


double RTCalculator::EludeIndexSumSqrDiff(string& peptide) { 
  double rtn = 0.0;
  unsigned int len = getPeptideLength(peptide);
  double diff = 0.0;
  for (unsigned int ix = 0; ix < len-1; ++ix) {
    diff = (*EludeIndex_)[getAATokenIndex(peptide, ix+1)]-(*EludeIndex_)[getAATokenIndex(peptide, ix)];
    rtn += diff*diff;
  }
  rtn /= len;

  return rtn;

}

double RTCalculator::EludeIndexNearestNeighborNeg(string& peptide) {
  double sum = 0.0;
  unsigned int len = getPeptideLength(peptide);
  string aatag="";
  for (unsigned int ix = 0; ix < len; ++ix) {
    aatag = getAATokenIndex(peptide, ix);
    if (aatag.find_first_of("DE")!=string::npos) {
      if (ix > 0) {
	if ((*EludeIndex_)[getAATokenIndex(peptide, ix-1)] > 0)
	  sum += (*EludeIndex_)[getAATokenIndex(peptide, ix-1)];

      }
      if (ix < len - 1) {
	if ((*EludeIndex_)[getAATokenIndex(peptide, ix+1)] > 0)
	  sum += (*EludeIndex_)[getAATokenIndex(peptide, ix+1)];
      }
    }
  }
  return sum;
}

void RTCalculator::EludeIndexPartialSum(string& peptide, vector<double> * features) {
  double sum = 0.0;
  size_t windows[2] = { 5, 2 };
  int len = getPeptideLength(peptide);

  size_t window = 0;
  for (int w = 0; w < 2; w++) {
    window = windows[w];
    double sum = 0.0;
    string::const_iterator lead = peptide.begin(), lag = peptide.begin();
    
    int start = 0, stop = 0;

    for (; start !=  window; ++start) {
      sum += (*EludeIndex_)[getAATokenIndex(peptide,start)];
    }

    double minS = sum, maxS = sum;

    for (; start != len; ++start, ++stop) {
      sum += (*EludeIndex_)[getAATokenIndex(peptide,start)];
      sum -= (*EludeIndex_)[getAATokenIndex(peptide,stop)];
      minS = min(sum, minS);
      maxS = max(sum, maxS);
    }
      
    if (windows[w] == 5) {
      (*features)[param_index_["EludeIndexPartialSum5Max"]] = maxS;
      (*features)[param_index_["EludeIndexPartialSum5Min"]] = minS;
    }
    if (windows[w] == 2) {
      (*features)[param_index_["EludeIndexPartialSum2Max"]] = maxS;
      (*features)[param_index_["EludeIndexPartialSum2Min"]] = minS;
    }
  }
}

void RTCalculator::EludeIndexPartialSum(string& peptide,  fann_type * features) {
  double sum = 0.0;
  size_t windows[2] = { 5, 2 };
  int len = getPeptideLength(peptide);

  size_t window = 0;
  for (int w = 0; w < 2; w++) {
    window = windows[w];
    double sum = 0.0;
    string::const_iterator lead = peptide.begin(), lag = peptide.begin();
    
    int start = 0, stop = 0;

    for (; start !=  window; ++start) {
      sum += (*EludeIndex_)[getAATokenIndex(peptide,start)];
    }

    double minS = sum, maxS = sum;

    for (; start != len; ++start, ++stop) {
      sum += (*EludeIndex_)[getAATokenIndex(peptide,start)];
      sum -= (*EludeIndex_)[getAATokenIndex(peptide,stop)];
      minS = min(sum, minS);
      maxS = max(sum, maxS);
    }
      
    if (windows[w] == 5) {
            (features)[param_index_["EludeIndexPartialSum5Max"]] = maxS;
            (features)[param_index_["EludeIndexPartialSum5Min"]] = minS;
    }
    if (windows[w] == 2) {
      (features)[param_index_["EludeIndexPartialSum2Max"]] = maxS;
      (features)[param_index_["EludeIndexPartialSum2Min"]] = minS;
    }
  }
}


void RTCalculator::EludeAmphipathicityHelix( string& peptide, vector<double> * features) {
  double min = 0.0, max = 0.0, hWindow = 0.0, rtn = 0.0;
  string aatag, aatagm3, aatagp3, aatagm4, aatagp4;
  int n = getPeptideLength(peptide);
  double cos300, cos400;
  // value to be added to the peptides < 9
  double cst;
  cos300 = cos(300 * M_PI / 180);
  cos400 = cos(400 * M_PI / 180);
  // calculate the average of hydrophobicity of the (*EludeIndex_)
  double avgHydrophobicityIndex = 0.0;

  map<string, double>::iterator itr;
  
  for (itr = EludeIndex_->begin(); itr != EludeIndex_->end(); itr++) {
    avgHydrophobicityIndex += itr->second; 
  }
  avgHydrophobicityIndex = avgHydrophobicityIndex / EludeIndex_->size();
  cst = avgHydrophobicityIndex * (1 + 2 * cos300 + 2 * cos400);
  // if the peptide is too short, use cst
  if (n < 9) {
    (*features)[param_index_["EludeAmphipathicityHelixMin"]] = cst;
    (*features)[param_index_["EludeAmphipathicityHelixMax"]] = cst;
  } else {
    // min (maximum) - initialized with the max(min) possible value
  
    for (int i = 4; i <= (n - 5); ++i) {
      hWindow = (*EludeIndex_)[getAATokenIndex(peptide, i)] + 
	(cos300 * ((*EludeIndex_)[getAATokenIndex(peptide, i - 3)] + (*EludeIndex_)[getAATokenIndex(peptide, i + 3)])) + 
	(cos400 * ((*EludeIndex_)[getAATokenIndex(peptide, i - 4)] + (*EludeIndex_)[getAATokenIndex(peptide, i + 4)]));
      if (i == 4) {
        min = hWindow;
        max = hWindow;
      } else {
        if (hWindow < min) {
          min = hWindow;
        }
        if (hWindow > max) {
          max = hWindow;
        }
      }
    }
    (*features)[param_index_["EludeAmphipathicityHelixMin"]] = min;
    (*features)[param_index_["EludeAmphipathicityHelixMax"]] = max;
  }

}

void RTCalculator::EludeAmphipathicityHelix( string& peptide, fann_type * features) {
  double min = 0.0, max = 0.0, hWindow = 0.0, rtn = 0.0;
  int n = getPeptideLength(peptide);
  double cos300, cos400;
  // value to be added to the peptides < 9
  double cst;
  cos300 = cos(300 * M_PI / 180);
  cos400 = cos(400 * M_PI / 180);
  // calculate the average of hydrophobicity of the index
  double avgHydrophobicityIndex = 0.0;

  map<string, double>::iterator itr;
  
  for (itr = EludeIndex_->begin(); itr != EludeIndex_->end(); itr++) {
    avgHydrophobicityIndex += itr->second; 
  }
  avgHydrophobicityIndex = avgHydrophobicityIndex / EludeIndex_->size();
  cst = avgHydrophobicityIndex * (1 + 2 * cos300 + 2 * cos400);
  // if the peptide is too short, use cst
  if (n < 9) {
    (features)[param_index_["EludeAmphipathicityHelixMin"]] = cst;
    (features)[param_index_["EludeAmphipathicityHelixMax"]] = cst;
  } else {
    // min (maximum) - initialized with the max(min) possible value
    for (int i = 4; i <= (n - 5); ++i) {
      hWindow = (*EludeIndex_)[getAATokenIndex(peptide, i)] + 
	(cos300 * ((*EludeIndex_)[getAATokenIndex(peptide, i - 3)] + (*EludeIndex_)[getAATokenIndex(peptide, i + 3)])) + 
	(cos400 * ((*EludeIndex_)[getAATokenIndex(peptide, i - 4)] + (*EludeIndex_)[getAATokenIndex(peptide, i + 4)]));
      if (i == 4) {
        min = hWindow;
        max = hWindow;
      } else {
        if (hWindow < min) {
          min = hWindow;
        }
        if (hWindow > max) {
          max = hWindow;
        }
      }
    }
    (features)[param_index_["EludeAmphipathicityHelixMin"]] = min;
    (features)[param_index_["EludeAmphipathicityHelixMax"]] = max;
  }

}

void RTCalculator::EludeHydrophobicMoments( string& peptide,
					    vector<double> * features) {
  double sum1 = 0.0, sum2 = 0.0;
  int lengthPeptide = getPeptideLength(peptide);
  double minHMoment, maxHMoment, hMoment;
  // calculate the angle in radians
  int angles[2] = {100,180};
  int w = 11;
  double angle = 0;
  for (int a = 0; a < 2; a++) {
    angle = angles[a] * M_PI / 180;
    
    for (int i = 0; i < min(lengthPeptide, w); ++i) {
      sum1 += (*EludeIndex_)[getAATokenIndex(peptide,i)] * sin((i + 1) * angle);
      sum2 += (*EludeIndex_)[getAATokenIndex(peptide,i)] * cos((i + 1) * angle);
    }
    hMoment = sqrt((sum1 * sum1) + (sum2 * sum2));
    minHMoment = hMoment;
    maxHMoment = hMoment;
    if (lengthPeptide > w) {
      for (int i = 1; i <= lengthPeptide - w; ++i) {
	sum1 = 0.0;
	sum2 = 0.0;
	for (int j = 0; j < w; j++) {
	  sum1 += (*EludeIndex_)[getAATokenIndex(peptide,i + j)] * sin((j + 1) * angle);
	  sum2 += (*EludeIndex_)[getAATokenIndex(peptide, i + j)] * cos((j + 1) * angle);
	}
	hMoment = sqrt((sum1 * sum1) + (sum2 * sum2));
	if (hMoment > maxHMoment) {
	  maxHMoment = hMoment;
	}
	if (hMoment < minHMoment) {
	  minHMoment = hMoment;
	}
      }
    }
    
    if (angles[a] == 100) {
      (*features)[param_index_["EludeHydrophobicMoment100Min"]] = minHMoment;
      (*features)[param_index_["EludeHydrophobicMoment100Max"]] = maxHMoment;
    }
    if (angles[a] == 180) {
      (*features)[param_index_["EludeHydrophobicMoment180Min"]] = minHMoment;
      (*features)[param_index_["EludeHydrophobicMoment180Max"]] = maxHMoment;
    }


  }

}

void RTCalculator::EludeHydrophobicMoments( string& peptide,
					    fann_type * features) {
  double sum1 = 0.0, sum2 = 0.0;
  int lengthPeptide = getPeptideLength(peptide);
  double minHMoment, maxHMoment, hMoment;
  // calculate the angle in radians
  int angles[2] = {100,180};
  int w = 11;
  double angle = 0;
  for (int a = 0; a < 2; a++) {
    angle = angles[a] * M_PI / 180;
    
    for (int i = 0; i < min(lengthPeptide, w); ++i) {
      sum1 += (*EludeIndex_)[getAATokenIndex(peptide,i)] * sin((i + 1) * angle);
      sum2 += (*EludeIndex_)[getAATokenIndex(peptide,i)] * cos((i + 1) * angle);
    }
    hMoment = sqrt((sum1 * sum1) + (sum2 * sum2));
    minHMoment = hMoment;
    maxHMoment = hMoment;
    if (lengthPeptide > w) {
      for (int i = 1; i <= lengthPeptide - w; ++i) {
	sum1 = 0.0;
	sum2 = 0.0;
	for (int j = 0; j < w; j++) {
	  sum1 += (*EludeIndex_)[getAATokenIndex(peptide,i + j)] * sin((j + 1) * angle);
	  sum2 += (*EludeIndex_)[getAATokenIndex(peptide, i + j)] * cos((j + 1) * angle);
	}
	hMoment = sqrt((sum1 * sum1) + (sum2 * sum2));
	if (hMoment > maxHMoment) {
	  maxHMoment = hMoment;
	}
	if (hMoment < minHMoment) {
	  minHMoment = hMoment;
	}
      }
    }
    
    if (angles[a] == 100) {
      (features)[param_index_["EludeHydrophobicMoment100Min"]] = minHMoment;
      (features)[param_index_["EludeHydrophobicMoment100Max"]] = maxHMoment;
    }
    if (angles[a] == 180) {
      (features)[param_index_["EludeHydrophobicMoment180Min"]] = minHMoment;
      (features)[param_index_["EludeHydrophobicMoment180Max"]] = maxHMoment;
    }


  }

}




double RTCalculator::helectric(string& pep) {
  string::size_type pos = (string::size_type) 0;
  string aastring = "";
  while (pos != string::npos) {
    string aatag = nextAAToken(pep, pos, pos);
    aastring += aatag.substr(0,1);
  }


  if ((int)aastring.length()>14) { return 0; }
	
  string co = aastring.substr(aastring.length()-4,4);


  string mpart = "";
	
  if (co.substr(0,1).find_first_of("DE") != string::npos)
  {
		
    mpart=co.substr(1,2); 

    if (mpart.find_first_of("PGKRH") != string::npos) {
      return 0; 
    } 
	
    for (int i=0; i<(int)mpart.length(); i++) {
      if (mpart[i]=='L' || mpart[i]=='I') {
	mpart[i] = 'X';
      }
      if (mpart[i]=='A' || mpart[i]=='V' ||
	  mpart[i]=='Y' || mpart[i]=='F' ||
	  mpart[i]=='W' || mpart[i]=='M' ) {
	mpart[i] = 'Z';
      }
      if (mpart[i]=='G' || mpart[i]=='S' ||
	  mpart[i]=='P' || mpart[i]=='C' ||
	  mpart[i]=='N' || mpart[i]=='K' ||
	  mpart[i]=='Q' || mpart[i]=='H' ||
	  mpart[i]=='R' || mpart[i]=='T' ||
	  mpart[i]=='D' || mpart[i]=='E') {
	mpart[i] = 'U';
      }
      
     
    }
		
    if (mpart.find("XX") != string::npos) { return 1; }
    if (mpart.find("ZX") != string::npos) { return 0.5; }
    if (mpart.find("XZ") != string::npos) { return 0.5; }
    if (mpart.find("ZZ") != string::npos) { return 0.4; }
    if (mpart.find("XU") != string::npos) { return 0.4; }
    if (mpart.find("UX") != string::npos) { return 0.4; }
    if (mpart.find("ZU") != string::npos) { return 0.2; }
    if (mpart.find("UZ") != string::npos) { return 0.2; }
		
  }
	
  return 0;
	
	
}






double RTCalculator::helicity1(string& pep) {
  string::size_type pos = (string::size_type) 0;
  string aatag;
  string hstring = "";
  int i=0;
  double sc4 = 0.;
  double sc5 = 0.;
  double sc6 = 0.;
  double term4 = 0.;
  double term5 = 0.;
  double term6 = 0.;
  double sum = 0.;
  while (pos != string::npos) {
    aatag = nextAAToken(pep, pos, pos);
    if (aatag.find_first_of( "PHRK")!=string::npos) {
      hstring += "z";
    }
    else     if (aatag.find_first_of( "WFIL")!=string::npos) {
      hstring += "X";
    }
    else     if (aatag.find_first_of( "YMVA")!=string::npos) {
      hstring += "Z";
    }
    else     if (aatag.find_first_of( "DE")!=string::npos) {
      hstring += "O";
    }
    else     if (aatag.find_first_of( "GSPCNKQHRT")!=string::npos) {
      hstring += "U";
    }
  }
  map<string, double>::iterator itr;
  for (i=0; i<(int)hstring.length()-3; i++) {
    string sub = hstring.substr(i, 4);
    sc4 = 0;
    if ((itr =  (*HelixProp_).find(sub)) != (*HelixProp_).end()) 
      sc4 = itr->second;

    if (sc4>0) {
      term4 = heli1TermAdj(sub, i, (int)hstring.length());
    }
    if ( i<(int)(hstring.length()-4) ) {
      sub = hstring.substr(i, 5);
      sc5 = 0;
      if ((itr =  (*HelixProp_).find(sub)) != (*HelixProp_).end()) 
	sc5 = itr->second;

      if (sc5>0) {
	term5 = heli1TermAdj(sub, i, (int)hstring.length());
      }
    }

    if (i<(int)(hstring.length()-5) ) {
      sub = hstring.substr(i, 6);
      //      cerr << sub << endl;
        sc6 = 0;
      if ((itr =  (*HelixProp_).find(sub)) != (*HelixProp_).end()) 
	sc6 = itr->second;
      if (sc6>0) {
	term6 = heli1TermAdj(sub, i, (int)hstring.length());
      }
    }
    if (sc6 > 0) {
      sum += sc6 * term6;
    }
    else if (sc5 > 0) {
      sum += sc5 * term5;
    }
    else if (sc4 > 0) {
      sum += sc4 * term4;
    }
    

  }
  return sum;

}

Boolean RTCalculator::recalc_RTstats() {
  double RTnumer = 0;
  double RTdenom = 0;
  double tot_RTs = 0;
  double all_RTs = 0;
  double SCANnumer = 0;
  double SCANdenom = 0;
  double tot_SCANs = 0;
  double all_SCANs = 0;
  double RT_mean = 0;
  double RT_stddev = 0;
  double SCAN_mean = 0;
  double SCAN_stddev = 0;
  
  if (ready_) return True;

  run_RT_used_count_ = 0;
  int i;
  for (i=0; i<(int)learn_rts_.size(); i++) {
    
    all_RTs += (*learned_RTs_)[i];
    all_SCANs += learn_rts_[i];
    //      all_SCANs += scans_[i];
    tot_RTs += (*learned_RTs_)[i];
    RTnumer += (*learned_RTs_)[i];
    RTdenom += 1;
    SCANnumer +=learn_rts_[i];
    //	SCANnumer += prob*scans_[i];
    SCANdenom += 1;
    run_RT_used_count_ ++;
  }


  if (RTdenom < 5) {
    cerr << "WARNING: Not enough high probability IDs in run index " << *run_name_ << " to generate RT model. RT Model has been disabled." << endl;
    return False;
  }
  else {
    RT_mean = RTnumer / RTdenom;
    RT_stddev = 0;
    SCAN_mean = SCANnumer / SCANdenom;
    SCAN_stddev = 0;
    RTdenom = 0;
    SCANdenom = 0;
    for (i=0; i<(int)learn_rts_.size(); i++) {
      RT_stddev += pow((RT_mean - (*learned_RTs_)[i]), 2);
      SCAN_stddev += pow((SCAN_mean - learn_rts_[i]), 2);
      //	  SCAN_stddev += prob * pow((SCAN_mean - scans_[i]), 2);
      RTdenom += 1;
      SCANdenom += 1;
    }
    
  }

  RT_stddev /= RTdenom;

  RT_stddev = pow(RT_stddev, 0.5); 
  
  run_RT_mean_ = RT_mean;
  run_RT_stddev_ = RT_stddev;

  SCAN_stddev /= SCANdenom;

  SCAN_stddev = pow(SCAN_stddev, 0.5); 
  
  run_SCAN_mean_ = SCAN_mean;
  run_SCAN_stddev_ = SCAN_stddev;

  double SSxx = 0;
  double SSxy = 0;
  double SSyy = 0;

  
  
  for (i=0; i<(int)learn_rts_.size(); i++) {
  
	SSyy += pow((learn_rts_[i] - run_SCAN_mean_), 2);
	//	SSyy += pow((scans_[i] - run_SCAN_mean_), 2);
	SSxx += pow(((*learned_RTs_)[i] - run_RT_mean_), 2);
	SSxy += (learn_rts_[i] - run_SCAN_mean_)*((*learned_RTs_)[i] - run_RT_mean_);

  }  

  run_slope_ = SSxy / SSxx;
  run_intercept_ = run_SCAN_mean_ - run_slope_ * run_RT_mean_;
  r_sq_ = (SSxy*SSxy) / (SSxx*SSyy);
  if (r_sq_ < 0.5) {
    cerr << "WARNING: Not enough correlation between theoretical RT values and scan numbers in run index " << *run_name_ << ". RT Model has been disabled." << endl;
    return False;
  }
  ready_ = true;
  return True;
}


double RTCalculator::calc_GradientCorrection(Array<double>* probs, double min_prob, Array<int>* ntts, int min_ntt, double rtMax) {
  double RTnumer = 0;
  double RTdenom = 0;
  double tot_RTs = 0;
  double all_RTs = 0;
  double SCANnumer = 0;
  double SCANdenom = 0;
  double tot_SCANs = 0;
  double all_SCANs = 0;
  double RT_mean = 0;
  double RT_stddev = 0;
  double SCAN_mean = 0;
  double SCAN_stddev = 0;
  
  if (ready_) return True;

  run_RT_used_count_ = 0;
  int i;
  for (i=0; i<(int)rts_.size(); i++) {
    if ((*probs)[i] >= min_prob && (*ntts)[i] >= min_ntt && (*RTs_)[i] > 0 && rts_[i] <= rtMax) {
      all_RTs += (*RTs_)[i];
      all_SCANs += rts_[i];
      //      all_SCANs += scans_[i];
      tot_RTs += (*RTs_)[i];
      RTnumer += (*RTs_)[i];
      RTdenom += 1;
      SCANnumer +=rts_[i];
      //	SCANnumer += prob*scans_[i];
      SCANdenom += 1;
      run_RT_used_count_ ++;
    }
  }


  if (RTdenom < 2) {
    cerr << "WARNING: Not enough IDs in run index " << *run_name_ << " to generate RT Gradient Correction." << endl;
    return 0;
  }
  else {
    RT_mean = RTnumer / RTdenom;
    RT_stddev = 0;
    SCAN_mean = SCANnumer / SCANdenom;
    SCAN_stddev = 0;
    RTdenom = 0;
    SCANdenom = 0;
    for (i=0; i<(int)rts_.size(); i++) {
      RT_stddev += pow((RT_mean - (*RTs_)[i]), 2);
      SCAN_stddev += pow((SCAN_mean - rts_[i]), 2);
      //	  SCAN_stddev += prob * pow((SCAN_mean - scans_[i]), 2);
      RTdenom += 1;
      SCANdenom += 1;
    }
    
  }

  RT_stddev /= RTdenom;

  RT_stddev = pow(RT_stddev, 0.5); 
  
  run_RT_mean_ = RT_mean;
  run_RT_stddev_ = RT_stddev;

  SCAN_stddev /= SCANdenom;

  SCAN_stddev = pow(SCAN_stddev, 0.5); 
  
  run_SCAN_mean_ = SCAN_mean;
  run_SCAN_stddev_ = SCAN_stddev;

  double SSxx = 0;
  double SSxy = 0;
  double SSyy = 0;

  
  
  for (i=0; i<(int)rts_.size(); i++) {
    if ((*probs)[i] >= min_prob && (*ntts)[i] >= min_ntt && (*RTs_)[i] > 0 && rts_[i] <= rtMax) {
      SSxx += pow((rts_[i] - run_SCAN_mean_), 2);
      //	SSyy += pow((scans_[i] - run_SCAN_mean_), 2);
      SSyy += pow(((*RTs_)[i] - run_RT_mean_), 2);
      SSxy += (rts_[i] - run_SCAN_mean_)*((*RTs_)[i] - run_RT_mean_);
    }
  }  

  run_slope_ = SSxy / SSxx;
  run_intercept_ = run_RT_mean_ - run_slope_ * run_SCAN_mean_;
  r_sq_ = (SSxy*SSxy) / (SSxx*SSyy);
  
  if (isnan(r_sq_)) {
    r_sq_ = -100;
  }

  if (r_sq_ < 0.5) {
    return r_sq_;
  }
  ready_ = true;
  return r_sq_;
}

double RTCalculator::getUsedForGradientRate() {
  if (rts_.size() > 0)
    return (double)used_count_ / (double)rts_.size();
  return 0.;
}

double RTCalculator::calc_GradientCorrection(double slope , double intercept) {
  double RTnumer = 0;
  double RTdenom = 0;
  double tot_RTs = 0;
  double all_RTs = 0;
  double SCANnumer = 0;
  double SCANdenom = 0;
  double tot_SCANs = 0;
  double all_SCANs = 0;
  double RT_mean = 0;
  double RT_stddev = 0;
  double SCAN_mean = 0;
  double SCAN_stddev = 0;
  int orig_used_count = used_count_;
  if (rts_.size() == 0 ) {
    return 0;
  }

  if (rts_.size() == 1 ) {
    return calc_GradientOffset();
  }
  
  
  //  if (ready_) return True;

  run_RT_used_count_ = 0;
  int i;
  for (i=0; i<(int)rts_.size(); i++) {
    if (rts_[i] > 0) {
      all_RTs += (*RTs_)[i];
      all_SCANs += rts_[i];
      //      all_SCANs += scans_[i];
      tot_RTs += (*RTs_)[i];
      RTnumer += (*RTs_)[i];
      RTdenom += 1;
      SCANnumer +=rts_[i];
      //	SCANnumer += prob*scans_[i];
      SCANdenom += 1;
      run_RT_used_count_ ++;

    }
  }
    

  if (RTdenom < 2) {
    cerr << "WARNING: Not enough IDs in run index " << *run_name_ << " to generate RT Gradient Correction." << endl;
    return 0;
  }
  else {
    RT_mean = RTnumer / RTdenom;
    RT_stddev = 0;
    SCAN_mean = SCANnumer / SCANdenom;
    SCAN_stddev = 0;
    RTdenom = 0;
    SCANdenom = 0;
    for (i=0; i<(int)rts_.size(); i++) {
      if (rts_[i] > 0) {
	RT_stddev += pow((RT_mean - (*RTs_)[i]), 2);
	SCAN_stddev += pow((SCAN_mean - rts_[i]), 2);
	//	  SCAN_stddev += prob * pow((SCAN_mean - scans_[i]), 2);
	RTdenom += 1;
	SCANdenom += 1;
      }
    }
  }

  RT_stddev /= RTdenom;

  RT_stddev = pow(RT_stddev, 0.5); 
  
  run_RT_mean_ = RT_mean;
  run_RT_stddev_ = RT_stddev;

  SCAN_stddev /= SCANdenom;

  SCAN_stddev = pow(SCAN_stddev, 0.5); 
  
  run_SCAN_mean_ = SCAN_mean;
  run_SCAN_stddev_ = SCAN_stddev;

  double SSxx = 0;
  double SSxy = 0;
  double SSyy = 0;

  
  
  for (i=0; i<(int)rts_.size(); i++) {
    if (rts_[i] > 0) {
      (*used_)[i] = true;
      SSxx += pow((rts_[i] - run_SCAN_mean_), 2);
      //	SSyy += pow((scans_[i] - run_SCAN_mean_), 2);
      SSyy += pow(((*RTs_)[i] - run_RT_mean_), 2);
      SSxy += (rts_[i] - run_SCAN_mean_)*((*RTs_)[i] - run_RT_mean_);
    }
    else {
      (*used_)[i] = false;
      used_count_--;
    }

  }  

  slope = run_slope_ = SSxy / SSxx;
  intercept = run_intercept_ = run_RT_mean_ - run_slope_ * run_SCAN_mean_;
  r_sq_ = (SSxy*SSxy) / (SSxx*SSyy);
 
  
  bool done =  r_sq_ >= 0.9;
  int count = 0;

  //  double best_rsq = r_sq_ ;
  //int last_max_i = -1;
  while ( !done && ++count < 50) {
    //Outlier Removal
    //DDS: Find distance mean and stddev
    gsl_vector *dis =  gsl_vector_calloc(run_RT_used_count_);


    int max_i = -1;
    double distMean = 0;
    double distStdev = 0;
    int usedCount = 0;
    double max_dis = 0;
    for (i=0; i<(int)rts_.size(); i++) {
      if ((*used_)[i]) {
	distMean += fabs(( slope * rts_[i] +  intercept) -(*RTs_)[i]);
	usedCount++;
      }
      if (max_dis < fabs(( slope * rts_[i] +  intercept) -(*RTs_)[i])) {
	max_dis = fabs(( slope * rts_[i] +  intercept) -(*RTs_)[i]);
	max_i = i;
      }
    }


    if (usedCount) 
      distMean /= usedCount;

    

    for (i=0; i<(int)rts_.size(); i++) {
      if ((*used_)[i]) {
	distStdev += pow(distMean-(( slope * rts_[i] +  intercept) -(*RTs_)[i]), 2);
      }

    }

    if (usedCount) {
      distStdev /= usedCount;

      distStdev = sqrt(distStdev);
    }

    
    for (i=0; i<(int)rts_.size(); i++) {
      if ((*used_)[i] && fabs(distMean-(( slope * rts_[i] +  intercept) -(*RTs_)[i])) > 1 * distStdev) {

	if (fabs(( slope * rts_[i] +  intercept) - (*RTs_)[i]) >= 0.98*max_dis) {
	  used_count_--;
	  (*used_)[i] = false ;
	}
      }
    }

    
    
    
    
    RTnumer = 0;
    RTdenom = 0;
    tot_RTs = 0;
    all_RTs = 0;
    SCANnumer = 0;
    SCANdenom = 0;
    tot_SCANs = 0;
    all_SCANs = 0;
    RT_mean = 0;
    RT_stddev = 0;
    SCAN_mean = 0;
    SCAN_stddev = 0;
    
    
    run_RT_used_count_ = 0;
    
    for (i=0; i<(int)rts_.size(); i++) {
      if ((*used_)[i]) {
	all_RTs += (*RTs_)[i];
	all_SCANs += rts_[i];
	//      all_SCANs += scans_[i];
	tot_RTs += (*RTs_)[i];
	RTnumer += (*RTs_)[i];
	RTdenom += 1;
	SCANnumer +=rts_[i];
	//	SCANnumer += prob*scans_[i];
	SCANdenom += 1;
	run_RT_used_count_ ++;
	
      }
    }
    
    
    if (RTdenom < 2) {
      cerr << "WARNING: Not enough IDs in run index " << *run_name_ << " to generate RT Gradient Correction." << endl;
      return 0;
    }
    else {
      RT_mean = RTnumer / RTdenom;
      RT_stddev = 0;
      SCAN_mean = SCANnumer / SCANdenom;
      SCAN_stddev = 0;
      RTdenom = 0;
      SCANdenom = 0;
      for (i=0; i<(int)rts_.size(); i++) {
	if ((*used_)[i]) {
	  RT_stddev += pow((RT_mean - (*RTs_)[i]), 2);
	  SCAN_stddev += pow((SCAN_mean - rts_[i]), 2);
	  //	  SCAN_stddev += prob * pow((SCAN_mean - scans_[i]), 2);
	  RTdenom += 1;
	  SCANdenom += 1;
	}
      }
    }
    
    RT_stddev /= RTdenom;
    
    RT_stddev = pow(RT_stddev, 0.5); 
    
    run_RT_mean_ = RT_mean;
    run_RT_stddev_ = RT_stddev;
    
    SCAN_stddev /= SCANdenom;
    
    SCAN_stddev = pow(SCAN_stddev, 0.5); 
    
    run_SCAN_mean_ = SCAN_mean;
    run_SCAN_stddev_ = SCAN_stddev;
    
    SSxx = 0;
    SSxy = 0;
    SSyy = 0;
    
    
    
    for (i=0; i<(int)rts_.size(); i++) {
      if ((*used_)[i]) {
	SSxx += pow((rts_[i] - run_SCAN_mean_), 2);
	//	SSyy += pow((scans_[i] - run_SCAN_mean_), 2);
	SSyy += pow(((*RTs_)[i] - run_RT_mean_), 2);
	SSxy += (rts_[i] - run_SCAN_mean_)*((*RTs_)[i] - run_RT_mean_);
	
      }  
    }
    
    slope = run_slope_ = SSxy / SSxx;
    intercept = run_intercept_ = run_RT_mean_ - run_slope_ * run_SCAN_mean_;

    r_sq_ = (SSxy*SSxy) / (SSxx*SSyy);


    done = r_sq_ > 0.9;
  }


      
  if (isnan(r_sq_)) {
    r_sq_ = -100;
  }

  if (r_sq_ < 0.5) {
    return r_sq_;
  }
  ready_ = true;
  return r_sq_;

}

double RTCalculator::calc_GradientCorrection() {
  int tmp_used = used_count_;
  calc_GradientCorrection(1, 0);
  used_count_ = tmp_used;
  return calc_GradientCorrection(run_slope_, run_intercept_);

}


double RTCalculator::calc_GradientOffset() {
  double RTnumer = 0;
  double RTdenom = 0;
  double tot_RTs = 0;
  double all_RTs = 0;
  double SCANnumer = 0;
  double SCANdenom = 0;
  double tot_SCANs = 0;
  double all_SCANs = 0;
  double RT_mean = 0;
  double RT_stddev = 0;
  double SCAN_mean = 0;
  double SCAN_stddev = 0;
  
  if (ready_) return True;

  run_RT_used_count_ = 0;
  int i;
  for (i=0; i<(int)rts_.size(); i++) {
      all_RTs += (*RTs_)[i];
      all_SCANs += rts_[i];
      //      all_SCANs += scans_[i];
      tot_RTs += (*RTs_)[i];
      RTnumer += (*RTs_)[i];
      RTdenom += 1;
      SCANnumer +=rts_[i];
      //	SCANnumer += prob*scans_[i];
      SCANdenom += 1;
      run_RT_used_count_ ++;

  }


  if (RTdenom < 1) {
    cerr << "WARNING: Not enough IDs in run index " << *run_name_ << " to generate RT Gradient Correction." << endl;
    return 0;
  }
  else {
    RT_mean = RTnumer / RTdenom;
    RT_stddev = 0;
    SCAN_mean = SCANnumer / SCANdenom;
    SCAN_stddev = 0;
    RTdenom = 0;
    SCANdenom = 0;
    for (i=0; i<(int)rts_.size(); i++) {
      RT_stddev += pow((RT_mean - (*RTs_)[i]), 2);
      SCAN_stddev += pow((SCAN_mean - rts_[i]), 2);
      //	  SCAN_stddev += prob * pow((SCAN_mean - scans_[i]), 2);
      RTdenom += 1;
      SCANdenom += 1;
    }
    
  }

  RT_stddev /= RTdenom;

  RT_stddev = pow(RT_stddev, 0.5); 
  
  run_RT_mean_ = RT_mean;
  run_RT_stddev_ = RT_stddev;

  SCAN_stddev /= SCANdenom;

  SCAN_stddev = pow(SCAN_stddev, 0.5); 
  
  run_SCAN_mean_ = SCAN_mean;
  run_SCAN_stddev_ = SCAN_stddev;

  double SSxx = 0;
  double SSxy = 0;
  double SSyy = 0;

  
  
  for (i=0; i<(int)rts_.size(); i++) {
 
      SSxx += pow((rts_[i] - run_SCAN_mean_), 2);
      //	SSyy += pow((scans_[i] - run_SCAN_mean_), 2);
      SSyy += pow(((*RTs_)[i] - run_RT_mean_), 2);
      SSxy += (rts_[i] - run_SCAN_mean_)*((*RTs_)[i] - run_RT_mean_);
  
  }  

  run_slope_ = 1;//SSxy / SSxx;
  run_intercept_ = run_RT_mean_ - run_slope_ * run_SCAN_mean_;
  r_sq_ = (SSxy*SSxy) / (SSxx*SSyy);

  
  if (isnan(r_sq_)) {
    r_sq_ = -100;
  }
  
  if (r_sq_ < 0.5) {
    return r_sq_;
  }
  ready_ = true;
  return r_sq_;
}


void RTCalculator::gradientCorrect() {
  int i;
  for (i=0; i<(int)rts_.size(); i++) {
    rts_[i] = rts_[i] * run_slope_ + run_intercept_;
    
  }
  

}

// Borrowed from Elude

// most and least hydrophobic patches
// double* RTModel::amphipathicityHelix(const float* index,
//                                      const string& peptide,
//                                      double* features) {
//   double min = 0.0, max = 0.0, hWindow = 0.0;
//   int n = peptide.length();
//   double cos300, cos400;
//   // value to be added to the peptides < 9
//   double cst;
//   cos300 = cos(300 * M_PI / 180);
//   cos400 = cos(400 * M_PI / 180);
//   // calculate the average of hydrophobicity of the index
//   double avgHydrophobicityIndex = 0.0;
//   for (int i = 0; i < 'Z' - 'A' + 1; ++i) {
//     avgHydrophobicityIndex += index[i];
//   }
//   avgHydrophobicityIndex = avgHydrophobicityIndex / aaAlphabet.size();
//   cst = avgHydrophobicityIndex * (1 + 2 * cos300 + 2 * cos400);
//   // if the peptide is too short, use cst
//   if (n < 9) {
//     *(features++) = cst;
//     *(features++) = cst;
//   } else {
//     // min (maximum) - initialized with the max(min) possible value
//     for (int i = 4; i <= (peptide.length() - 5); ++i) {
//       hWindow = index[peptide[i] - 'A'] + (cos300 * (index[peptide[i - 3]
//           - 'A'] + index[peptide[i + 3] - 'A'])) + (cos400
//           * (index[peptide[i - 4] - 'A'] + index[peptide[i + 4] - 'A']));
//       if (i == 4) {
//         min = hWindow;
//         max = hWindow;
//       } else {
//         if (hWindow < min) {
//           min = hWindow;
//         }
//         if (hWindow > max) {
//           max = hWindow;
//         }
//       }
//     }
//     *(features++) = min;
//     *(features++) = max;
//   }
//   return features;
// }
