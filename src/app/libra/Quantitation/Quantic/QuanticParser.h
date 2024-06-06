#ifndef _QUANTICPARSER_H_
#define _QUANTICPARSER_H_

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

#include "Quantic.h"
#include "Common/TPPVersion.h"
#include "Parsers/Parser/Parser.h"
#include "Validation/InterProphet/InterProphetParser/KDModel.h"
#include "Search/SpectraST/SpectraSTPepXMLLibImporter.hpp"
#include "cramp.hpp"
#include <sstream>
#include <ostream>

typedef vector<streampos> stream_vec;
typedef vector<cRamp*> cramp_vec;

using namespace std;

class ProbPos {
 public:
  ProbPos(double pr, int ps) {
    prob_ = pr;
    pos_ = ps;
  }
 

  ~ProbPos() {}
  double prob_;
  int pos_;
};




class QuanticParser : public Parser {

 public:
  QuanticParser(double mztol, double ppmtol, unsigned int verbose);

  ~QuanticParser() {};
  
  void parseRead(const char* c, int t = 0);
  void parseWrite(const char* c, int t = 0);
  void parseWriteDIA(const char* c, int t = 0);

  void parse(const char* c) { parseRead(c, 0); parseWrite(c,0); }


  void parseWrite(const char* c) ;
  void parseWriteUpdate(const char* c, bool update, int t = 0);
  void parseWriteUpdateDIA(const char* c, bool update, int t = 0);

  void computePTMModel();


  void parseSpecProbs(const char* c);

  void run(const char* c, const char* opts, int max_threads = 1, unsigned int top_peaks = 6);

  void runDiaMode(const char* c, const char* opts, int max_threads = 1, unsigned int top_peaks = 6);

  string pepToMolecForm(string pep);
 

  void setOutFile(string);

  void setUpdate(bool up);

  void setMinProb(double mP) {
    minProb_ = mP;
  }

  void setKeepOld(bool ko) {
    keepOld_ = ko;
  }

  void setDiaMode(bool diaMode) {
    diaMode_ = diaMode;
  }
  void setAnnotate(bool an) {
    annotate_ = an;
  } 

  void setPeptideCoupled(bool p) {
    pep_coupled_ = p;
  }

  void setModString(string& modstring) ;

  
  void setMolForms(string& molForms) {
    ptmMolForms_ = molForms;
  }


  void writeUpdatedModTags(Quantic* proph, ostream& fout, string& mod_pep_str);
  unsigned int  max_threads_;
  unsigned int  top_peaks_;
  


  unsigned int SQtot_;
    

  char* file_;
 private:

  std::string modstring_;
  ofstream fout;

  //  Peptide* pep_;
  std::string tmp_file_;
  std::string out_file_;

  std::string opts_;
  //  Quantic* ptm_proph_;
  KDModel* ptm_model_;
  KDModel* mat_model_;
  size_t total_mods_;
  size_t total_sites_;

  bool keepOld_;

  bool diaMode_;
  
  bool update_mod_tags_;

  

  vector<double> priors_;

  TPP_HASHMAP_T<char, double> stat_mods_hash_;
  TPP_HASHMAP_T<char, vector<double>*> var_mods_hash_;


  TPP_HASHMAP_T<char, double> stat_prot_termods_hash_;
  TPP_HASHMAP_T<char, vector<double>*> var_prot_termods_hash_;


  std::string ptmMolForms_;
  TPP_STDSTRING_HASHMAP(double) spec_prob_hash_;

  double mzTol_;
  
  double minProb_;

  double massOffset_;

  unsigned int verbose_;

  double ppmTol_;
  
  stream_vec SQoffsets_;

  vector<unsigned int> SQruns_;

  vector<double> massdiffs_;

  vector<string> aminoacids_;
  vector<double> massshift_;
  vector<string> ptmtypes_;
  vector<vector<double>*> neutlosses_;

   
  cramp_vec cramps_;

  bool annotate_;
  bool pep_coupled_;

  
};

static unsigned int threads_done_run_;
static unsigned int threads_start_run_;
static unsigned int current_run_;
static   unsigned int header_thread_;

#ifdef MSVC
static  HANDLE _mutex;
static  HANDLE run_mutex;
static  HANDLE hdr_mutex;
static  HANDLE fout_mutex;
static  HANDLE sq_mutex;
#else
static  pthread_mutex_t _mutex;
static  pthread_mutex_t run_mutex;
static  pthread_mutex_t hdr_mutex;
static  pthread_mutex_t fout_mutex;
static  pthread_mutex_t sq_mutex;
#endif




#endif
