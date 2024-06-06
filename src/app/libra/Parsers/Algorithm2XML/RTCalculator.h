/*

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

*/

#ifndef RTCALC_H
#define RTCALC_H

#include "Common/sysdepend.h"
#include "RTCatalog.h"
#include "GradientProgram.h"


#include <string.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include <fann.h>

#ifndef __LGPL__
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_vector.h>
#endif


#include "Common/Array.h"
#include "Common/ModificationInfo/ModificationInfo.h"
#include "Validation/MixtureDistribution/ICATMixtureDistr.h"

#include "Parsers/Algorithm2XML/pICalculator.h"

#include "Common/tpp_hashmap.h" 

typedef TPP_STDSTRING_HASHMAP(Array<double>*) dblarr_hash;
typedef TPP_STDSTRING_HASHMAP(Array<int>*) intarr_hash;

using namespace std;


class RTCalculator {
  friend class RTMixtureDistr;
  //  friend class VariableOffsetRTMixtureDistr;
  friend class KernelDensityRTMixtureDistr;

 public:

  RTCalculator();
  ~RTCalculator();
  RTCalculator(string* run);
  // seq must be stripped of all modifications
  double Peptide_RT(char* seq, int scan);
  void addPeptide_RT(const char* seq, const char* modified_pep, int scan, double rt, string* run_name);
  // can have modifications

  void EludeAmphipathicityHelix( string& peptide, vector<double> * features) ;
  void EludeHydrophobicMoments( string& peptide, vector<double> * features);
  void EludeAmphipathicityHelix( string& peptide, fann_type * features) ;
  void EludeHydrophobicMoments( string& peptide, fann_type * features);
  void EludeIndexPartialSum(string& peptide, vector<double> * features);
  void EludeIndexPartialSum(string& peptide, fann_type * features);
  double EludeIndexNearestNeighborNeg(string& peptide);
  double EludeIndexNearestNeighborPos(string& peptide);
  double EludeIndexSumSqrDiff(string& peptide);
  double EludeIndexSum(string& peptide);
  double EludeIndexAvg(string& peptide) ;
  double EludeIndexN(string& peptide) ;
  double EludeIndexC(string& peptide) ;

  double calc_GradientCorrection(Array<double>* probs, double min_prob, Array<int>* ntts, int min_ntt, double rtMax); 
  double calc_GradientCorrection();
  double calc_GradientCorrection(double slope, double intercept);
  double calc_GradientOffset();  
 
  Boolean recalc_RTgsl(double min_prob, int min_ntt) ;
  void train_RTCatalog(const char* catfile);

  void adjustRT();

  double calcSSR1(string& pep);

  //  void linearTxRT(Array<RTCalculator*>* run_calcs);

  void calc_RTstats();
  double getUsedForGradientRate();

  double calc_PepRT(string& pep);
  Boolean recalc_RTstats(Array<double>* probs, double min_prob, Array<int>* ntts, int min_ntt);

  Boolean recalc_RTstats(Array<double>* probs, double min_prob, Array<int>* ntts, int min_ntt, double rtMax);


  void write_RTcoeff(ostream& out);
  void read_RTcoeff(istream& in);

  void write_paramIndex(ostream& out);
  void read_paramIndex(istream& in);

  void read_RTcoeff();
  void read_RTcoeff(const char* coeffile);
  
  double get_PepRTCatalog(string& pep) ;

  void read_RTCatalog(const char* catfile);

  void train_RTcoeff(istream& in);

  bool learnNeuralNet();

  bool learnNeuralNet2();

  bool learnNeuralNetCV();

 Boolean recalc_RTstats( double min_prob, int min_ntt);
  Boolean recalc_RTgsl(Array<double>* probs, double min_prob, Array<int>* ntts, int min_ntt);
  //Boolean recalc_RTstats(Array<double>* probs, double min_prob); 
  double getRTScore(double RT, int scan, double rt);

  double getRTScore(double RT, double rt);
  void recalc_RTstats(Array<double>* probs);
  Boolean recalc_RTstats();
  
  void write_RTstats(ostream& out);

  void write_RTann(string& file);

  void read_RTann();
  void read_RTann(string& file);
  
  string getRunName();
  double getRunSlope() { return run_slope_; }
  double getRunInt() { return run_intercept_; }


  void setRunSlope(double sl) { run_slope_ = sl; }
  void setRunInt(double inc) { run_intercept_ = inc; }

  double getRunRSQ() { return r_sq_; }


  double getPredRT(int i) ;

  double getObsRT(int i) ;
  
  int getUsedRTCount() { return used_count_; }

  int getAllRTCount() { return (int)rts_.size(); }

  double calcANN_PepRT(string& pep);
  double calcANN_PepRT2(string& pep);

  int getPredRTCount() ;
  string getRunNameatIdx(int i);
  double getRunRTMean();
  double getRunRTStdDev();
  double getSSRCalcHP(char* peptide);
  void batchSSRCalcHP();
  void batchRTCalc();
  void batchRTCatalog(const RTCatalog& catalog);
  void batch_iRT(dbl_hash& irt_hash);
  void batchRTCalc(GradientProgram*);
  void batchRTCalc(const char* catfile);
  double getHydroEHP(char* peptide);
  double getHydroGHP(char* peptide);

  double getRTScore(int index);

  double computeHelixScore(string& pep);
  void initHelixProp();
  void InvertLine();


  bool linearRegressRT(Array<double>* probs, double min_prob, Array<int>* ntts, int min_ntt);
  bool linearRegressRT(double min_prob, int min_ntt);
  bool linearRegressRT();
  void plotRegressionFits(Array<double>* probs, double min_prob, Array<int>* ntts, int min_ntt);

  void distillModifiedPeptides(Array<double>* probs, Array<int>* ntts); 

  // DDS: Ported from perl SSRCalc3.pl
  double hydroMoment(string& pep) ;
  double helectric(string& pep) ;
  double isoelectric(string& pep) ;
  double smallness(string& pep) ;
  double helicity1(string& pep) ;
  double clusterness(string& pep) ;
  double heli1TermAdj(string & hc, int start, int len);
  double helicity2(string& pep) ;
  void heli2Calc(string& pep, double& Hiscore, double& Gscore);
  double evalH2pattern(string& pat, string& testsq, int posn, string etype);
  double connector(string & acid, string & lp, string & rp, string  ct, string& far1, string& far2) ;
  int  getPeptideLength(string& s);

  string reversePeptide(string& s);

  string stripPeptide(string& s);

  void gradientCorrect();
  
  const string* getModPeptide(int i);

  int isOutlier(int i);

 protected:
  string ann_file_;
  string* run_name_;
  double run_RT_sum_;
  double learned_RT_sum_;
  int run_RT_count_;
  int learned_RT_count_;
  double run_RT_mean_;
  double run_RT_stddev_;
  double run_SCAN_mean_;
  double run_SCAN_stddev_;
  int run_RT_used_count_;
  double run_slope_;
  double run_intercept_;
  double r_sq_;
  
  double min_RT_;
  
  bool ready_;
  bool adjusted_;

  Array<int> scans_;
  Array<double>* learn_probs_;
  Array<int>* learn_ntts_;
  Array<double>* RTs_; // predicted RT's

  Array<bool>* used_; // boolean array identifying retentions that are actually used to learn the linear model (after outlier removal)

  int used_count_;

  Array<double>* learned_RTs_; // predicted RT's

  Array<char*>* peps_;
  
  vector<string> modified_peps_;
 
  vector<string> learn_peps_;

  dblarr_hash* best_modpep_rts_;

  intarr_hash* best_modpep_runs_;

  intarr_hash* best_modpep_ids_;
 
  unsigned int max_len_;

  // HENRY - also store the experimental RT values.
  // And isn't it time to use std::vector over whatever "Array" is?
  vector<double> rts_; // experimental RT's


  vector<double> learn_rts_; // experimental RT's to learn

  gsl_vector* c_;

  vector<string*> run_ids_;
  vector<vector<double> >* matrix_;
  vector<vector<double> >* learn_matrix_;

  unsigned int paramcount_;

  int aacount_;

  double max_rt_;

  double min_rt_;

  pICalculator * pIcalc_; 

  int verbose_;
  
  map<string, int> param_index_;

  map<string, int> aa_index_;

  double EXP10(double value);
  
  map<string, double> coeff_;
  
  RTCatalog* rt_catalog_;

  struct fann *ann_;

#ifndef __LGPL__
  gsl_multifit_linear_workspace* rt_fit_;
#endif

  string nextAAToken(string& s, string::size_type start, string::size_type& end);

  string getAATokenIndex(string& s, int i);

  int countAA(string& pep);

  static map<string, double>* HelixProp_;
  static map<string, double>* HydroPath_;
  static map<string, double>* ClustProp_;
  static map<string, double>* EludeIndex_;
  static map<string, double>* SSR1Coeff_;
  static map<string, double>* SSR1NTermCoeff_;

  static void initSSR1Coeff();

  
};

#endif
