#ifndef _QUANTIC_H_
#define _QUANTIC_H_
/*

Program       : Quantic
Author        : David Shteynberg <dshteynb  AT systemsbiology.org>
Date          : 04.20.2018
SVN Info      : $Id$

Copyright (C) 2018 David Shteynberg

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

*/

#include "Common/tpp_hashmap.h"
#include "Search/SpectraST/PeptideUser.hpp"
#include "Search/SpectraST/SpectraSTPepXMLLibImporter.hpp"
#include "Search/SpectraST/SpectraSTPeakList.hpp"
#include "Validation/InterProphet/InterProphetParser/KDModel.h"
#include <sstream>
#include <iomanip>
#include "CMercury8.h"

#ifndef __LGPL__
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#endif

using namespace std;
using namespace hardklor;

typedef TPP_STDSTRING_HASHMAP(string*) strp_hash;

class NeutQuant {
public:
  NeutQuant() {
    loss_ = 0;
    quant_ = 0;
    mean_  = 0;
    stdev_ = 0;
  }
  NeutQuant(double l, double q, double m, double s) {
    loss_ = l;
    quant_ = q;
    mean_  = m;
    stdev_ = s;
  }
  ~NeutQuant();

  double loss_;
  double quant_;
  double mean_;
  double stdev_;
};

class ModifiedPeptide {
 public:
  ModifiedPeptide() {
    mpep_ = NULL;
    mods_ = NULL;
    nomods_ = NULL;
  }
  ~ModifiedPeptide() {
    delete mpep_;
  }
  PeptideUser* mpep_;
  vector<vector<int>*>* mods_;
  vector<vector<int>*>* nomods_;
};

class Quantic {

 public:
  Quantic(string& pep, int charge, double calc_neut_mass, 
	  cRamp* cramp, long scan,
	  vector<string>& modaas, vector<double>& shift, vector<vector<double>*>& neutlosses,
	  TPP_HASHMAP_T<char, double>* stat_mods, TPP_HASHMAP_T<char, vector<double>*>* var_mods, 
	  TPP_HASHMAP_T<char, double>* stat_prot_termods,	TPP_HASHMAP_T<char, vector<double>*>* var_prot_termods,
#ifdef MSVC
	  HANDLE* mutex,
#else
	  pthread_mutex_t* mutex,
#endif
	  bool is_nterm_pep, bool is_cterm_pep, TPP_HASHMAP_T<int, double>& pos_mod_hash);
  // bool massdiff_mode = false, bool labile_mode = false, bool direct_mode = false);



  ~Quantic();
  
  void setDiaMode(bool mode) {
    diaMode_ = mode;
    
    mercury_ = new CMercury8(NULL);
  }

  bool init();




  vector<vector<int>*>* getModCombo(vector<gsl_combination*>& compare);


  void computeMassDiffs(vector<double>* absmassdiffs);


  void setNions(string& nions) {
    for (int i=0; i<nions.size(); i++) {
      nions_.push_back(nions[i]);
    }
  }
  
  void setCions(string& cions) {
    for (int i=0; i<cions.size(); i++) {
      cions_.push_back(cions[i]);
    }
  }
 
  

  int nTermMod();
  int cTermMod();


  bool hasNTermMod() { return has_nterm_mod_; }; 
  bool hasCTermMod() { return has_cterm_mod_; }; 
  
  bool isNtermPep() { return is_nterm_pep_; }
  bool isCtermPep() { return is_cterm_pep_; }


  void setPrecision(ostringstream & stream , double& value);

  std::string numberToString(double, int);
  std::string numberToString(long);

  int NAA() { return NAA_; }
  
  PeptideUser* getPeptide() { return pep_; }
  
  void setMZTolerance(double tol) {
    mz_tol_ = tol;
  }

  void setTopPeaks( unsigned int top_peaks) {
    top_peaks_ = top_peaks;
  }


  void setVerbose(unsigned int v) {
    verbose_ = v;
  }
  
  void setPeptideCoupled(bool p) {
    pep_coupled_ = p;
  }


  void setPPMTolerance(double tol) {
    ppm_tol_ = tol;
  }

  void setMassDiff(double massdiff); 
  
  void evaluatePeptideAnnotatedTIC(map<double,double>* avoidMZ=NULL);
  
  double getAnnotatedTIC() {
    return annot_tic_ ;
  }

  double getTIC() {
    return TIC_ ;
  }
  
  string pepToMolecForm(string& pep);

  
  string getAnnotatedString() {
    return annot_str_ ;
  }

  string getTheoreticalString() {
    return theor_str_ ;
  }

  double getQuant(int type, int idx) {
    return (*neutlossquants_[type])[idx]->quant_ ;
  }

  double getQuantMean(int type, int idx) {
    return (*neutlossquants_[type])[idx]->mean_ ;
  }

  double getQuantStdev(int type, int idx) {
    return (*neutlossquants_[type])[idx]->stdev_ ;
  }

  void setMolForms(string& ptmMolForms);

  void correctIsotopeErrors(int type, CMercury8* mercury);

  void correctIsotopeErrors(vector<NeutQuant*>& input, int type, CMercury8* mercury);

  void evaluateModPep(ModifiedPeptide* mpep);

  void updateQuants(int type,
		    vector<double>& quantEvid) {
    int q_count = 0;
    for (int n=0; n < neutlosses_[type]->size(); n++) {
      if ( quantEvid[n] > 0 ) {
	q_count++;
      }
    }
    if (q_count == neutlosses_[type]->size()) {
      for (int n=0; n < neutlosses_[type]->size(); n++) {
	(*neutquants_[type])[n] +=quantEvid[n];
	(*neutlossquants_[type])[n]->quant_ +=quantEvid[n];		
      }
    }
    else {
      for (int n=0; n < neutlosses_[type]->size(); n++) {
	if ( quantEvid[n] > 0 ) {
	  (*neutquantvect_[type])[n]->pop_back();
	}
      }
    }
  }
  
 void updateQuants(int type,
		    vector<NeutQuant*>& quantEvid) {
    int q_count = 0;
    for (int n=0; n < neutlosses_[type]->size(); n++) {
      if ( quantEvid[n]->quant_ > 0 ) {
	q_count++;
      }
    }
    if (q_count == neutlosses_[type]->size()) {
      for (int n=0; n < neutlosses_[type]->size(); n++) {
	(*neutquants_[type])[n] +=quantEvid[n]->quant_;
	(*neutlossquants_[type])[n]->quant_ +=quantEvid[n]->quant_;
	(*neutquantvect_[type])[n]->push_back(quantEvid[n]->quant_);
      }
    }
    else {
      for (int n=0; n < neutlosses_[type]->size(); n++) {
	//if ( quantEvid[n]->quant_ > 0 ) {
	//(*neutquantvect_[type])[n]->push_back(quantEvid[n]->quant_);
	  //}
      }
    }
  }
  void calcQuantStats(int type) {
    vector<double> wts;
    vector<double> vals;

    double wtsum = 0;
    double wtsumsq = 0;
    double wtsqsum = 0;
    
    vector<double> sums;
    vector<double> sumsqrs;

    
    for (int n=0; n < neutlosses_[type]->size(); n++) {
      sums.push_back(0);
      sumsqrs.push_back(0);
      for (int i=0; i<(*neutquantvect_[type])[n]->size(); i++) {
	if (n == 0) {
	  wts.push_back((*(*neutquantvect_[type])[n])[i]);
	}
	else {
	  wts[i] += (*(*neutquantvect_[type])[n])[i]; //Weight is total intensity of all channels
	}		
      }
    }

    for (int n=0; n < neutlosses_[type]->size(); n++) {
      for (int i=0; i<(*neutquantvect_[type])[n]->size(); i++) {
	(*(*neutquantvect_[type])[n])[i] = (*(*neutquantvect_[type])[n])[i] / wts[i]; //Normalize
	sums[n] +=  wts[i] * (*(*neutquantvect_[type])[n])[i];

	if (n==0) {
	  wtsum += wts[i];
	  wtsqsum += wts[i]*wts[i];
	}
      }
      if (fabs(wtsum) > 0) {
	(*neutquantmeans_[type])[n] = sums[n] / wtsum;
	(*neutlossquants_[type])[n]->mean_ = sums[n] / wtsum; 

      }
      else {
	(*neutquantmeans_[type])[n] = 0; 
	(*neutlossquants_[type])[n]->mean_ = 0;
      }
    }
    
    for (int n=0; n < neutlosses_[type]->size(); n++) {
      for (int i=0; i<(*neutquantvect_[type])[n]->size(); i++) {
	if (i==0) sums[n] = 0;
	sums[n] +=  wts[i] *
	  ( (*(*neutquantvect_[type])[n])[i] - (*neutquantmeans_[type])[n] ) *
	  ( (*(*neutquantvect_[type])[n])[i] - (*neutquantmeans_[type])[n] );
	
      }
      if (fabs(wtsum*wtsum - wtsqsum) > 0) {
	(*neutquantstdvs_[type])[n] = sqrt(sums[n] * wtsum / (wtsum*wtsum - wtsqsum));

	(*neutlossquants_[type])[n]->stdev_ = sqrt(sums[n] * wtsum / (wtsum*wtsum - wtsqsum)); 

      }
      else {
	(*neutquantstdvs_[type])[n] = 0;
	(*neutlossquants_[type])[n]->stdev_ = 0;
      }
    }

    

  }

  
  map<double,double>* evaluateModPepMap(ModifiedPeptide* mpep, bool update);
  map<double,double>* evaluateModPepMap(ModifiedPeptide* mpep, bool update,  map<double,double>* avoidMZ=NULL);

  void processPeakList();
  
#ifdef MSVC
  HANDLE* mutex_;
#else
  pthread_mutex_t* mutex_;
#endif

 private:
  string theor_str_;
  string annot_str_;
  double annot_tic_;

  unsigned int top_peaks_;
  PeptideUser* pep_;
  PeptideUser* pep_unmod_;
  string* pep_str_;
  SpectraSTLibEntry* entry_;
  SpectraSTPeakList* peakList_;
  long scan_;
  int charge_;
  cRamp* cramp_;
  int NAA_;
  double mz_tol_;
  bool etd_;

  float ppm_tol_;
  bool unknown_mod_;
  
  double mass_diff_;
  bool massdiff_mode_;

  bool direct_mode_;
  bool autodirect_;
  
  vector<int> ntermod_;
  vector<int> ctermod_;

  int recur_index_; 

  int nDECOY_;

  bool labile_mode_;

  double TIC_;

  double minInt_;

  double maxInt_;

  bool pep_coupled_; //Also quantify peptide coupled neutral losses

  vector<string> modaas_;
  vector<double> shift_;

  vector<int> direct_mode_bytype_;


  vector<vector<NeutQuant*>*> neutlossquants_;
  
  vector<vector<double>*> neutlosses_;

  vector<vector<double>*> neutquants_;

  vector<vector<double>*> neutquantmeans_;
  vector<vector<double>*> neutquantstdvs_;

  vector<vector<vector<double>*>*> neutquantvect_;

  vector<strp_hash*> label_;

  vector<string> pep_prob_str_;
  
  double calc_neut_mass_;

  string pep_unmod_str_;
  
  vector<TPP_HASHMAP_T<int, double>* > pos_prob_hash_;

  vector<TPP_HASHMAP_T<int, double>* > loss_prob_hash_;

  vector<SpectraSTLibEntry*>* decoyEntries_; 
  vector<SpectraSTPeakList*>* decoyPeakLists_;

  vector<vector<double>*> lossprob_;
  vector<vector<double>*> siteprob_;
  vector<vector<double>*> sitesum_;
  vector<vector<double>*> site_MaxEvidence_;
  vector<vector<double>*> site_MinEvidence_;
  vector<vector<double>*> site_decoyMaxEvidence_;
  vector<vector<double>*> site_decoyMinEvidence_;
  vector<vector<double>*> site_decoyEvidence_;


  vector<vector<int>*> site_MaxCount_;
  vector<vector<int>*> site_nomodMaxCount_;

  vector<vector<int>*> labile_MaxCount_;
  vector<vector<int>*> labile_nomodMaxCount_;

  vector<vector<int>*> unlabile_MaxCount_;
  vector<vector<int>*> unlabile_nomodMaxCount_;


  vector<vector<double>*> labile_MaxEvidence_;
  vector<vector<double>*> unlabile_MaxEvidence_;

 vector<vector<double>*> labile_Pval_;

 vector<vector<double>*> labile_Oscore_;
 vector<vector<double>*> labile_Iscore_;

 vector<vector<double>*> labile_Mscore_;

 vector<vector<double>*> labile_ObsEvidence_;

 vector<vector<double>*> unlabile_ObsEvidence_;

 vector<vector<double>*> labile_ExpEvidence_;

 vector<vector<double>*> unlabile_ExpEvidence_;


  vector<vector<vector<double>*>*> site_nomodAllEvidence_;
  vector<vector<vector<double>*>*> site_AllEvidence_;

  vector<vector<double>*> site_Pval_;


  vector<vector<double>*> site_Iscore_;

  vector<vector<double>*> site_Oscore_;

  vector<vector<double>*> site_Mscore_;

  vector<vector<double>*> site_ObsModEvidence_;

  vector<vector<double>*> site_ObsUnModEvidence_;
  vector<vector<double>*> site_ExpModEvidence_;
  vector<vector<double>*> site_ExpUnModEvidence_;


  vector<vector<double>*> site_nomodMaxEvidence_;
  vector<vector<double>*> site_nomodMinEvidence_;

  vector<vector<ModifiedPeptide*>*> site_nomodMaxEvidPep_;
  vector<vector<ModifiedPeptide*>*> site_MaxEvidPep_;

  vector<vector<ModifiedPeptide*>*> labile_MaxEvidPep_;
  vector<vector<ModifiedPeptide*>*> unlabile_MaxEvidPep_;

  vector<vector<vector<vector<int>*>*>*> site_nomodMaxEvidMods_;
  vector<vector<vector<vector<int>*>*>*> site_MaxEvidMods_;

  vector<vector<vector<vector<int>*>*>*> labile_MaxEvidMods_;
  vector<vector<vector<vector<int>*>*>*> unlabile_MaxEvidMods_;

  vector<vector<double>*> site_nomodDecoyMaxEvidence_;
  vector<vector<double>*> site_nomodDecoyMinEvidence_;

  vector<vector<int>*> site_N_;


  vector<vector<int>*> nomodsite_;
  vector<vector<int>*> modsite_;

  vector<vector<vector<int>*>*> allmods_;
  vector<vector<vector<int>*>*> allnomods_;
  vector<vector<int>*> mods_;
  vector<vector<int>*> nomods_;
  vector<int> nmods_;


  ModifiedPeptide* mpep_;

  vector<vector<vector<int>*>*> combs_of_mods_;

  TPP_HASHMAP_T<char, double>* stat_mods_;
  TPP_HASHMAP_T<char, vector<double>*>* var_mods_;

  bool is_nterm_pep_;

  bool is_cterm_pep_;

  vector<char> nions_;
  vector<char> cions_;

  bool has_nterm_mod_;

  bool has_cterm_mod_;

  CMercury8* mercury_;

  TPP_HASHMAP_T<char, double>* stat_prot_termods_;
  TPP_HASHMAP_T<char, vector<double>*>* var_prot_termods_;

  TPP_HASHMAP_T<char, string> ptm_mol_hash_;
  
  bool diaMode_;

  int nterm_;
  
  unsigned int verbose_;
};



#endif
