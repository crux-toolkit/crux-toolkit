#include <deque>
#include "peptides.pb.h"
#include "peptide.h"
#include "theoretical_peak_set.h"
#include "fifo_alloc.h"
#include "spectrum_collection.h"
#include "io/OutputFiles.h"
#include "IonInvertedIndex.h"

#ifndef ACTIVE_PEPTIDE_QUEUE_H
#define ACTIVE_PEPTIDE_QUEUE_H

class TheoreticalPeakCompiler;

class ActivePeptideQueue {
 public:
  ActivePeptideQueue(RecordReader* reader,
        const vector<const pb::Protein*>& proteins,
        vector<const pb::AuxLocation*>* locations=NULL, 
        bool dia_mode = false);

  ~ActivePeptideQueue();

  int SetActiveRange(double min_range, double max_range); 
  inline size_t size() const { return queue_.size(); }  

  Peptide* GetPeptide(int index) {
    return *(begin_ + index); 
  }

  SCORE_FUNCTION_T curScoreFunction_;  

  int nPeptides_;
  int nCandPeptides_;
  int CandPeptidesTarget_;
  int CandPeptidesDecoy_;

  // This vector will store the score hitograms for either the XCorr or for the HyperScore.
  // This histogram is used for score calibration, currently for Tailor.
  // Later, this histogram also can be used for the linear regression based score calibration, (E-value) like in Comet or X!Tandem
  vector<int> score_histogram_;
  // Since scores are real valued, the scores will be discretized, the bin of a scores S is 
  // score_histogram[round(score*score_scale_factor)]
  int score_scale_factor_;    
  // XCorr scores can have negative value, so we need an offset to make scores positive
  // the bin of the score_histogram[round(score+score_histogram_offset)*score_scale_factor)]
  int score_histogram_offset_;   
  int max_score_;
  int highest_bin_;
  int score_count_;
  int total_match_;

  // Score --> bin:
  // int bin = round( (score + score_histogram_offset_)*score_scale_factor_ );
  // ++score_histogram_[bin];
  // bin --> Score:
  // score = (double)end/score_scale_factor_ - score_histogram_offset_ 


  const double TAILOR_QUANTILE_TH = 0.005;
  const double TAILOR_OFFSET = 5.0 ;
  // Calculate Tailor scores. Get the 99th quantile:
  double getTailorQuantileFromHistogram();
  double getMeanMatch() {
    return (double)total_match_/(double)score_count_;
  }  
  void ResetHist();
  void AddScoreToHist(double score, int match_cnt = 0);

  // class PeptideWrapper {  // Peptide objects split into hot and cold data in order to reduce cache miss ratio. 
  // hot data:
  //   double mass;
  //   bool isDecoy;
  //   vector<unsigned int> peaks_0;   // Single charged b-y ions, in case of exact p-value, this contains only the b-ions
  //   vector<unsigned int> peaks_1;   // Double charged b-y ions
  // cold data:
    // Peptide* peptide;
  // }

  deque<Peptide*> queue_;
  deque<Peptide*>::const_iterator begin_, end_;

  int min_candidates_;
  bool dia_mode_;

  IonInvertedIndex ion_inverted_index_;

 private:

  void ComputeTheoreticalPeaksBack();    

  RecordReader* reader_;
  const vector<const pb::Protein*>& proteins_; 
  vector<const pb::AuxLocation*>* locations_;
  
  TheoreticalPeakSetBYSparse theoretical_peak_set_;
  pb::Peptide current_pb_peptide_;
};

#endif
