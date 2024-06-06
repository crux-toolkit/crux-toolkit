#ifndef _RTCATALOG_H_
#define _RTCATALOG_H_

/*

Program       : RTCatalog
Author        : David Shteynberg <dshteynb  AT systemsbiology.org>                                                       
Date          : 09.29.07

Primary data object holding all mixture distributions for each precursor ion charge
%
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
dshteynb@systemsbiology.org

*/


#include "Common/tpp_hashmap.h" // deals with different compilers
//#include "Common/ResidueMass.h" 
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef __LGPL__
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_vector.h>
#endif

#include "Validation/InterProphet/InterProphetParser/SearchHit.h"
#include "Parsers/Algorithm2XML/GradientProgram.h"
//#include <pwiz/data/msdata/MSDataFile.hpp>
//#include <pwiz/data/msdata/SpectrumInfo.hpp>
//#include <pwiz_tools/common/FullReaderList.hpp>
#include "Parsers/mzParser/mzParser.h"

using namespace std;  //MH: This is dangerous in a header. TODO: move it elsewhere


typedef TPP_STDSTRING_HASHMAP(SearchHit*) hit_hash;
typedef TPP_STDSTRING_HASHMAP(Array<SearchHit*>*) arrhit_hash;
typedef TPP_STDSTRING_HASHMAP(Array<double>*) dblarr_hash;
typedef TPP_STDSTRING_HASHMAP(dblarr_hash*) dblarr_hash_hash;
typedef TPP_STDSTRING_HASHMAP(Array<int>*) intarr_hash;
typedef TPP_STDSTRING_HASHMAP(Array<string*>*) strparr_hash;
typedef TPP_STDSTRING_HASHMAP(double) dbl_hash;
typedef TPP_STDSTRING_HASHMAP(string) str_hash;
typedef TPP_STDSTRING_HASHMAP(int) int_hash;
typedef TPP_STDSTRING_HASHMAP(int) bool_hash;
typedef TPP_STDSTRING_HASHMAP(dbl_hash*) dbl_hash_hash;
typedef TPP_STDSTRING_HASHMAP(dbl_hash_hash*) dbl_hash_hash_hash;
typedef TPP_STDSTRING_HASHMAP(int_hash*) int_hash_hash;
typedef TPP_STDSTRING_HASHMAP(bool_hash*) bool_hash_hash;

class RTInfo {
 public:
  RTInfo();
  RTInfo(int n,  double q25, double med, double q75,double mean,double stdev, double min) {
    n_ = n;
    med_ = med;
    origMed_ = med;
    q25_ = q25;
    q75_ = q75;
    siqr_  =  (q75_ - q25_)/2.;
    mean_ = mean;
    origMean_ = mean;
    stdev_ = stdev;
    min_ = min;
    origMin_ = min;
  }
  void ComputeSIQR() {
    siqr_ = (q75_ - q25_)/2.;
  }
  ~RTInfo();
  int n_;
  double origMin_;
  double origMed_;
  double origMean_;
  double min_;
  double q25_;
  double med_;
  double q75_;
  double siqr_;
  double mean_;
  double stdev_;
};

typedef TPP_STDSTRING_HASHMAP(RTInfo*) rtinfo_hash;
typedef TPP_STDSTRING_HASHMAP(rtinfo_hash*) rtinfo_hash_hash;
typedef TPP_STDSTRING_HASHMAP(string) mod_hash;



class RTCatalog {
  friend class RTCatalogParser;
  friend class RTCalculator;
 public:
  
  RTCatalog(double minProb, double minRT=0, double maxRT=10000, string* acnFile=NULL);

  RTCatalog(const char* file);
  
  ~RTCatalog();
    
  int getNumRuns();

  void progress(int tic, int step, int &tot);

  int getRTCount(string& pep);
 
  double getRTMedian(string& pep) const;

  double getRTMin(string& pep) const;
  
  double getRTSIQR(string& pep) const;
  
  double getRTMean(string& pep) const ;
  
  double getRTStdev(string& pep) const;

  RTInfo* getRTInfo(string& pep);

  RTInfo* getRTInfo(string& pep, string& run);

  void insertResult(const string& run, const string& modpep);
  void insertResult(string& run, string& spectrum,
		    double prob, Array<double>* allntt_prob, 
		    string& pepseq, string& modpep, 
		    double calcnmass, double rt, 
		    double prec_intens, double collision_eng, 
		    double matchedions_frac,  
		    string& exp_lbl, string& charge);

  void adjustResult(const int run_idx, string& modpep);
 
  bool addRun(string& run_name);
  string getPeptideString(const string& in_pep);
  void initMods();
  void sortRunNames();
  void writeRunRTStats(const char* file);

  void calcRTStatsByRun();

  void calcRTStatsCombined(ostream& out);
  void trackPeptidesReport(ostream& out, Array<string*>* peps);
  void trackPeptidesReportPV(ostream& out, Array<string*>* peps, bool iRT);

  void trackPeptidesChromatograms(ostream& out, Array<string*>* runs, dblarr_hash* peps_q1q3);
  void trackPeptidesChromatograms(Array<string*>* runs, dblarr_hash* peps_q1q3);

  bool rejectResult(double prob, double rt) ;
  void setMods(bool pv, bool cys);

  void setDisco(bool disco) { disco_ = disco; }

  void setDaltonTolerance() {
    dalTOL_ = true;
  }
  void trackPeptidesXICs(Array<string*>* runs, dblarr_hash_hash* byrun_peps_q1s, double maxTOL);
 private :

  // RTCalculator* rt_calc_;

  //  Array<dblarr_hash*>* peprts_hash_ ;  
  //Array<dblarr_hash*>* pepintens_hash_ ;  
  // array of RTs in a peptide look up table, one table for each run_idx

  mod_hash* mod_table_;

  strparr_hash* pepruns_hash_;
  // hash of runname ptrs, one for each peptide

  Array<string*>* run_names_;

  str_hash* run_files_;


  rtinfo_hash* peprtinfo_hash_;

  rtinfo_hash_hash* byrun_peprtinfo_hash_;
  dblarr_hash_hash* byrun_peprts_hash_ ; 
  dblarr_hash_hash* byrun_pepintens_hash_ ; 
  dblarr_hash_hash* byrun_pepMatchedIons_hash_ ;  
  
  
  dbl_hash_hash* byrun_pep_maxIntens_;
  dbl_hash_hash* byrun_pep_maxMatchedIons_;


  dbl_hash* byrunchrome_maxInt_hash_;
  dbl_hash* byrunchrome_intWtdStdev_hash_;
  dbl_hash* byrunchrome_intWtdMean_hash_;

  bool_hash* ismix_run_hash_;

  int num_runs_;


  double minProb_;

  double minRT_;
  double maxRT_;

  bool dalTOL_;

  bool pv_mods_;

  bool cys_cam_;

  bool disco_;

  GradientProgram* acn_gradient_;

};

#endif
