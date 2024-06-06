#ifndef PRO_DATA_COMP_H
#define PRO_DATA_COMP_H

/*

Program       : ProDataComponent                                                    
Author        : Andrew Keller <akeller@systemsbiology.org>                                                       
Date          : 11.27.02 

Primary data object holding all mixture distributions for each precursor ion charge

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

#include "SearchResult.h"
#include "SequestResult.h"
#include "MascotResult.h"
#include "CometResult.h"
#include "Common/constants.h"

struct prodatacomponent_struct {

  char xpressratio[50];

  int lightfirstscan;
  int lightlastscan;
  int heavyfirstscan;
  int heavylastscan;
  double lightmass;
  double heavymass;
  double masstol;
  int xpresslight;

  double asap_mean;
  double asap_err;
  double asap_inv_mean;
  double asap_inv_err;
  int asapratio_index;
  char score_summary[1000];

};



class ProDataComponent {

 public:
  ProDataComponent(int seq, int pk, int data, int xml_ind, int ind, Boolean heavy2light);
  ProDataComponent(int seq, int pk, int data, int xml_ind, int ind, Boolean heavy2light, int msms_run_idx);
  ~ProDataComponent();
  char* display(Array<char*>* basenames, Array<char*>* pepproph_times, Array<char*>* iproph_times, Array<char*>* dbs, Array<char*>* asap_times);
  void enter(pepDataStrct pepdata, SearchResult* result, int basename_index, int pepproph_timestamp_index, int iproph_timestamp_index,
	   int database_index, int asap_timestamp_index, prodatacomponent_struct comp_data, long scan);

//  void enter(pepDataStrct pepdata, SearchResult* result, int basename_index, int pepproph_timestamp_index,
//	     int database_index, int asap_timestamp_index, double xpress, double asap_mean, double asap_err, int asap_ind, 
//	     char* score_summary);

  Boolean match(int seq, int pk, int data);
  Boolean xml_match(int xml_ind);
#ifdef USE_STD_MODS
  void write(ostream& os, Array<char*>* inputfiles, Array<char*>* basenames, Array<char*>* pepproph_times, Array<char*>* iproph_times, Array<char*>* dbs, Array<char*>* asap_times, char* radio, Array<char*>* misc_run_conds, char* colored_aas);
#endif
#ifndef USE_STD_MODS
  void write(ostream& os, Array<char*>* inputfiles, Array<char*>* basenames, Array<char*>* pepproph_times, Array<char*>* iproph_times, Array<char*>* dbs, Array<char*>* asap_times, char* radio, Array<char*>* aa_mods, Array<char*>* term_mods, Array<char*>* misc_run_conds);
#endif

#ifdef USE_STD_MODS
  void write(ostream& os, Array<char*>* inputfiles, Array<char*>* basenames, Array<char*>* pepproph_times, Array<char*>* iproph_times, Array<char*>* dbs, Array<char*>* asap_times,  Array<Boolean>* asap_quantHighBGs, Array<Boolean>* asap_zeroBGs, Array<double>* asap_mzBounds, Array<bool>* asap_wavelets, char* radio, Array<char*>* misc_run_conds, char* colored_aas);
#endif
#ifndef USE_STD_MODS
  void write(ostream& os, Array<char*>* inputfiles, Array<char*>* basenames, Array<char*>* pepproph_times, Array<char*>* iproph_times, Array<char*>* dbs, Array<char*>* asap_times,  Array<Boolean>* asap_quantHighBGs, Array<Boolean>* asap_zeroBGs, Array<double>* asap_mzBounds,  Array<bool>* asap_wavelets, char* radio, Array<char*>* aa_mods, Array<char*>* term_mods, Array<char*>* misc_run_conds);
#endif

#ifdef USE_STD_MODS
  void writeStandardFormat(ostream& os, Array<char*>* inputfiles, Array<char*>* basenames, Array<char*>* pepproph_times, Array<char*>* iproph_times, Array<char*>* dbs, Array<char*>* asap_times, Array<Boolean>* asap_quantHighBGs, Array<Boolean>* asap_zeroBGs, Array<double>* asap_mzBounds, Array<bool>* asap_wavelets,  Boolean include_header, Array<char*>* misc_run_conds, char* colored_aas);
#endif
#ifndef USE_STD_MODS
  void writeStandardFormat(ostream& os, Array<char*>* inputfiles, Array<char*>* basenames, Array<char*>* pepproph_times, Array<char*>* iproph_times, Array<char*>* dbs, Array<char*>* asap_times, Array<Boolean>* asap_quantHighBGs, Array<Boolean>* asap_zeroBGs, Array<double>* asap_mzBounds, Array<bool>* asap_wavelets, Boolean include_header, Array<char*>* aa_mods, Array<char*>* term_mods, Array<char*>* misc_run_conds);
#endif

#ifdef USE_STD_MODS
  void writeStandardFormat(ostream& os, Array<char*>* inputfiles, Array<char*>* basenames, Array<char*>* pepproph_times, Array<char*>* iproph_times, Array<char*>* dbs, Array<char*>* asap_times, Boolean include_header, Array<char*>* misc_run_conds, char* colored_aas);
#endif
#ifndef USE_STD_MODS
  void writeStandardFormat(ostream& os, Array<char*>* inputfiles, Array<char*>* basenames, Array<char*>* pepproph_times, Array<char*>* iproph_times, Array<char*>* dbs, Array<char*>* asap_times, Boolean include_header, Array<char*>* aa_mods, Array<char*>* term_mods, Array<char*>* misc_run_conds);
#endif
  void setMassType(int val, int fragval);
  int seq_;
  int peak_;
  int data_;
  int xml_index_;
  int data_index_;
  
  int msms_run_idx_;

  pepDataStrct pepdata_;
  int basename_index_;
  long scan_;

  SearchResult* result_;


 protected:

  Boolean setModifiedPeptide(char* modpep, char* buffer, int buf_len);
  char* getColoredPeptide(const char* peptide, const char* labeled_aas, const char* starttag, const char* endtag);


  // other things nec for display....

  // prob
  //int basename_index_;
  int pepproph_timestamp_index_;
  int iproph_timestamp_index_;

  // protein
  int database_index_;
  // asap
  int asap_timestamp_index_;

  
  // xpress
  char* xpress_;

  double asap_mean_;
  double asap_error_;
  double asap_inv_mean_;
  double asap_inv_error_;
  int asap_index_;

  int LightFirstScan_;
  int LightLastScan_;
  int HeavyFirstScan_;
  int HeavyLastScan_;
  double LightMass_;
  double HeavyMass_;
  double MassTol_;
  int bXpressLight_; // 0 default, 1 forced light=1, 2 forced heavy=1
  
  Boolean heavy2light_;

  char* score_summary_;
  int masstype_; // 0 for average, 1 for monoisotopic
  int fragmasstype_; // 0 for average, 1 for monoisotopic
  
};

#endif
