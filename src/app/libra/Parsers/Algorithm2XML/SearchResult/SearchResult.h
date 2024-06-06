#ifndef SEARCHRESULT_H
#define SEARCHRESULT_H

#include <string.h>
#include <iostream>
#include <vector>
#include <stdlib.h>

#include "Common/sysdepend.h"
#include "Common/constants.h"
#include "Parsers/Parser/Tag.h"
#include "Common/tpp_hashmap.h"

#ifdef USE_STD_MODS
#include "Common/ResidueMass/ResidueMass.h"
#include "Common/ModificationInfo/ModificationInfo.h"
#endif

//#define SIZE_BUF 8192
#define VAL_UNCERTAINTY 1


/*

Program       : SearchResult for discr_calc of PeptideProphet 
Author        : Andrew Keller <akeller@systemsbiology.org>                                                       
Date          : 11.27.02 


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
enum xlink_type_t { na, xl, loop };

class SearchResult {

 public:

  SearchResult();
  SearchResult(Array<Tag*>* tags);

  virtual ~SearchResult();

  virtual Boolean isMaldi(char* spec);
  char* strCopy(const char* orig);
  Boolean isProcessed();
  virtual char* extractDatabaseWithTags(const char* html, const char* start_tag, const char* end_tag);
  void setRunIdx(int run_idx);
  void setRunName(string* run_name);
#ifdef USE_STD_MODS
  Boolean setModificationEncoding(char* buffer, int buffer_len);
#endif
  int charge_;
  int scan_;
  char* spectrum_;
  Array<char*>* proteins_;
  char* protein_;
  char* peptide_;
  string strip_peptide_;
  char* modified_peptide_;
  Boolean maldi_; // whether or not maldi data
  int degen_; // whether or not degen
  double massdiff_;
  double adduct_;
  virtual const char* getName() = 0;
  int format_;
  double probability_; // for preexising...
  Boolean adjusted_prob_; // whether or not prob is adjusted (2/3)
  Boolean incomplete_prob_; // whether or not analysis was complete
  double neutral_mass_;
  int num_matched_ions_;
  int tot_num_ions_;
  int num_tol_term_;
  int num_missed_cl_;
  int varmod_count_;

  double unweight_spec_entropy_; //MSBooster
  double delta_RT_loess_; //MSBooster

  double CV_; // compensationn0 if not used
  double pI_; // 0 if not used
  double RT_; // 0 if not used

  xlink_type_t xlink_type_;

  //DDS: pI model
  double run_pI_diff_; // 15 if not used
  double run_RT_diff_; // 15 if not used
  int run_idx_;
  string* run_name_;
  string exp_label_;

#ifdef USE_STD_MODS
  ModificationInfo* mod_info_;
  char prev_aa_;
  char next_aa_;
#endif

  
  virtual ostream& print(ostream& os);
  virtual char* stripHTML(const char* orig);

  int getVarModCount(char aa) {
    if ((*var_mods_hash_).find(aa) != (*var_mods_hash_).end() && (*var_mods_hash_)[aa]) {
      return (*var_mods_hash_)[aa]->size();
    }
    return 0;
  }

  double getVarModMass(char aa, int i) {
    if ((*var_mods_hash_).find(aa) != (*var_mods_hash_).end() && (*var_mods_hash_)[aa]) {
      if ( i >=0 && i < (*var_mods_hash_)[aa]->size()) {
	return (*(*var_mods_hash_)[aa])[i];
      }
    }
    return 0.;
  }

 protected:

  virtual void init();
  Boolean processed_; // whether or not ok

  
  TPP_HASHMAP_T<char, vector<double>*>* var_mods_hash_;


};

#endif
