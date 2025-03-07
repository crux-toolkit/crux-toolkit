#ifndef TIDE_MATCH_SET_H
#define TIDE_MATCH_SET_H

#define  BOOST_DATE_TIME_NO_LIB
#include <boost/thread.hpp>
#include <vector>
#include "raw_proteins.pb.h"
#include "tide/records.h"
#include "tide/fixed_cap_array.h"
#include "tide/peptide.h"
#include "tide/sp_scorer.h"
#include "tide/spectrum_collection.h"
#include "tide/ActivePeptideQueue.h"
#include "tide/spectrum_preprocess.h"

#include "model/Modification.h"
#include "model/PostProcessProtein.h"

using namespace std;

typedef vector<const pb::Protein*> ProteinVec;

class TideMatchSet {    
 public:
  class Scores {
   public:
    int ordinal_;
    double xcorr_score_;
    double exact_pval_;
    double refactored_xcorr_;
    int resEv_score_;
    double resEv_pval_;
    double combined_pval_;
    double tailor_; 
    int by_ion_matched_; 
    int by_ion_total_;    
    int repeat_ion_match_; 
    double sp_score_;
    double hyper_score_;
    double hyper_score_la_; 
    double delta_cn_;
    double delta_lcn_;
    bool active_;
    deque<Peptide*>::const_iterator peptide_itr_;
    Scores():ordinal_(0), xcorr_score_(0.0), exact_pval_(0.0), refactored_xcorr_(0.0), 
      resEv_score_(0.0), resEv_pval_(0.0), combined_pval_(0.0), tailor_(0.0), by_ion_matched_(0), by_ion_total_(0), 
      sp_score_(0), hyper_score_(0), hyper_score_la_(0), delta_cn_(0), delta_lcn_(0), active_(false) {}
  };
//   typedef FixedCapacityArray<Scores> PSMScores;
  typedef vector<Scores> PSMScores;
  PSMScores psm_scores_;   // This one is used to gather psms during scoring.

  int n_concat_or_target_matches_;  // concat or target
  int n_decoy_matches_;
  bool psm_scores_processed_;

  static bool cmpXcorrScore(const Scores& x, const Scores& y) {  // compare PSMs by xcorr scores. Larger comes first
    return x.xcorr_score_ < y.xcorr_score_;
  }  
  static bool cmpCombinedPvalue(const Scores& x, const Scores& y) {  // compare PSMs by P-values. smaller comes first
    // return x.refactored_xcorr_ < y.refactored_xcorr_;
    return x.combined_pval_ > y.combined_pval_;
  }  
  static bool cmpHyperScore(const Scores& x, const Scores& y) {  // compare PSMs by hyper scores. larger comes first
    return x.hyper_score_ < y.hyper_score_;
  }  

  // Define the column names and their order in the result files.
  static int XCorr_tsv_cols[];  //these are declared at the beginning of TideMatchSet.cpp
  static int Pvalues_tsv_cols[];
  static int Diameter_tsv_cols[]; 

  static int XCorr_mzTab_cols[];  //these are declared at the beginning of TideMatchSet.cpp
  static int Pvalues_mzTab_cols[];  //these are declared at the beginning of TideMatchSet.cpp

  static int XCorr_pin_cols[];  //these are declared at the beginning of TideMatchSet.cpp
  static int Pvalues_pin_cols[];

  // TideMatchSet();
  TideMatchSet(ActivePeptideQueue* active_peptide_queue, ObservedPeakSet* observed);
  ~TideMatchSet();

  static int* getColumns(TSV_OUTPUT_FORMATS_T format, size_t& numHeaders);
  static string getHeader(TSV_OUTPUT_FORMATS_T format, string tide_index_mztab_param_file); // pass filetype;
  void getReport(TSV_OUTPUT_FORMATS_T format, string spectrum_filename,
                   const SpectrumCollection::SpecCharge* sc, int spectrum_file_cnt, 
                   string &concat_or_target_report, string& decoy_report);
  void gatherTargetsDecoys();  // Additional scores are:  delta_cn, delta_lcn, tailor
  void calculateAdditionalScores(PSMScores& psm_scores, const SpectrumCollection::SpecCharge* sc);  // Additional scores are:  delta_cn, delta_lcn, tailor; 
  void printResults(TSV_OUTPUT_FORMATS_T format, string spectrum_filename, const SpectrumCollection::SpecCharge* sc, int spectrum_file_cnt, bool target, PSMScores& psm_scores, string& results,
    map<PSMScores::iterator, boost::tuple<double, double, double>>* intensity_map = NULL,
    map<PSMScores::iterator, boost::tuple<double, double, double>>* logrank_map = NULL,
    map<PSMScores::iterator, boost::tuple<double, double, double>>* coelute_map = NULL,
    map<PSMScores::iterator, boost::tuple<double, double>>* ms2pval_map = NULL,
    map<string, double>* peptide_predrt_map = NULL);

  static string GetModificationList(const pb::ModTable* mod_table, string site_prefix, string position_prefix, bool variable, int& cnt);
  /* Constants required for the tailor scoring */
  const double TAILOR_QUANTILE_TH = 0.01;
  const double TAILOR_OFFSET = 5.0 ;
  double quantile_score_;
  
  PSMScores::iterator last_psm_;
  ActivePeptideQueue* active_peptide_queue_;  
  
  // Pointer to the experimental spectrum data; This is used here to calculate 
  // the repeat_ion_match value (part of the SP scoring). Originally, the
  // repeat_ion_match value was calulated for each PSMs in TideSearch App. 
  // But it is just needed  only for top-N PSMs. So, we need to used
  // the observed spectrum vector to calculate the repeat_ion_match here.
  ObservedPeakSet* observed_;  
  
  // Global static parameters
  static SCORE_FUNCTION_T curScoreFunction_;
  static int top_matches_;
  static int decoy_num_;
  static int mass_precision_;
  static int score_precision_;
  static int mod_precision_;
  static string decoy_prefix_;
  static bool concat_;
  static int psm_id_mzTab_;
  static string fasta_file_name_;

//  private:
  PSMScores concat_or_target_psm_scores_;
  PSMScores decoy_psm_scores_;


};

#endif
