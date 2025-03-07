#include <deque>
#include "peptides.pb.h"
#include "peptide.h"
#include "theoretical_peak_set.h"
#include "fifo_alloc.h"
#include "spectrum_collection.h"
#include "io/OutputFiles.h"

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

  int SetActiveRange(vector<double>* min_mass, vector<double>* max_mass, 
        double min_range, double max_range); 

  Peptide* GetPeptide(int index) {
    return *(begin_ + index); 
  }

  int nPeptides_;
  int nCandPeptides_;
  int CandPeptidesTarget_;
  int CandPeptidesDecoy_;

  deque<Peptide*> queue_;
  deque<Peptide*>::const_iterator begin_, end_;  
  int min_candidates_;
  bool dia_mode_;

 private:
  bool isWithinIsotope(vector<double>* min_mass,
        vector<double>* max_mass, 
        double mass,
        int* isotope_idx);   

  void ComputeTheoreticalPeaksBack();    

  RecordReader* reader_;
  const vector<const pb::Protein*>& proteins_; 
  vector<const pb::AuxLocation*>* locations_;
  
  TheoreticalPeakSetBYSparse theoretical_peak_set_;
  pb::Peptide current_pb_peptide_;
};

#endif
