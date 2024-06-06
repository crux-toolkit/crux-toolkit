#ifndef OPTION_H
#define OPTION_H

/*

Program       : Option                                                       
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

#include <string.h>
#include <iostream>
#include <stdlib.h>

#include "Common/constants.h"
#include "Parsers/mzParser/mzParser.h"
//#include "Validation/DiscriminateFunction/DiscrimValMixtureDistr.h"

enum IcatOption {
  ICAT_UNKNOWN,
  ICAT_ON,
  ICAT_OFF };

enum DeltaCnStarOption {
  DELTACN_ZERO,
  DELTACN_LEAVE,
  DELTACN_EXCLUDE };

enum MascotStarOption {
  MASCOTSTAR_LEAVE,
  MASCOTSTAR_EXCLUDE, 
  MASCOTSTAR_PENALIZE,
  MASCOT_EXPECT};

enum CometStarOption {
  COMETSTAR_LEAVE,
  COMETSTAR_EXCLUDE, 
  COMETSTAR_ZERO };

enum SearchEngine {
  SEQUEST,
  MASCOT, 
  COMET,
  PROBID};

enum SampleEnzyme {
  tryptic,
  chymotryptic, 
  gluc_bicarb,
  elastase,
  nonspecific,
  tca,
  CNBr,
  AspN,
  tryptic_CNBr };

class ModelOptions {
public:

  IcatOption icat_;
  bool silac_;
  Boolean glyc_;
  Boolean phospho_;
  bool xl_;
  bool msboostEntropy_;
  bool msboostRTloess_;
  bool matchedions_;
  bool onefval_;
  double minprob_;
  int extraitrs_;
  Boolean pI_;
  int min_pI_ntt_;
  double min_pI_prob_;
  Boolean RT_;
  Boolean CV_;
  int min_RT_ntt_;
  double min_RT_prob_;
  Boolean accMass_;
  Boolean accMassZ_;
  double massWidth_;
  Boolean ppm_;
  float bandWidthX_; //multiplicative factor on f-value bandwith in nonparam_ modeling mode
  float conservative_;
  Boolean maldi_;
  Boolean instrwarn_;
  Boolean no_neg_init_;
  Boolean no_ntt_;
  Boolean no_xlsecond_;
  Boolean no_xltopexp_;
  Boolean no_nmc_;
  Boolean no_vmc_;
  Boolean forcedistr_;
  Boolean no_mass_;
  Boolean nonparam_;
  Boolean use_expect_;
  Boolean neg_gamma_;
  Boolean low_mem_;
  int negdistr_type_;
  Boolean use_decoy_;
  int max_threads_;
  Boolean output_decoy_probs_;
  Boolean use_chg_[MAX_CHARGE];

  char enzyme_[100];
  char engine_[100];
  
  char decoy_label_[100]; 
  char* rtcat_file_;
  char* rtcoef_file_;
  char* rtann_file_;

  char* search_offsets_;

  //  SearchEngine engine_;
  //  SampleEnzyme enzyme_;
  mzParser::InstrumentStruct* massspec_;

  // HENRY: Add icat type
  int spectrast_icat_;
  // bool spectrast_delta_;
  bool optimize_fval_;
  // END HENRY

  // HENRY: Add whether to multiply by lib probability
  Boolean multiply_by_spectrast_lib_probs_;
  // END HENRY

  // Not all systems will initialize memory to zero.
  ModelOptions()
  {  memset(this, 0, sizeof(ModelOptions)); }
};

class ScoreOptions {
public:
  DeltaCnStarOption deltastar_;
  MascotStarOption mascotstar_;
  CometStarOption cometstar_;
  char inputfile_[5000];
};


struct StaticModificationCount {
  char mod;
  double mass;
  int num;
};


//void setSearchEngine(ModelOptions* opts, char* engine);
//char* getSearchEngine(ModelOptions opts);

//void setSampleEnzyme(ModelOptions* opts, char* enzyme);
//char* getSampleEnzyme(ModelOptions opts);







#endif
