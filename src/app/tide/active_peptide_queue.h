// Benjamin Diament
//
// An ActivePeptideQueue is constructed with a file of peptides of
// non-decreasing neutral mass, and a correpsonding set of proteins. With
// successive call to SetActiveRange(min_mass, max_mass) the ActivePeptideQueue
// reads in peptides and initializes them. It also discards any peptides in 
// memory lighter than min_mass.
//
// The ActivePeptideQueue maintains in memory a window of peptides that fall
// within  a given mass range. The window to maintain (the "active peptides")
// is set by a call to SetActiveRange(). Successive values
// of min_mass and max_mass must be non-decreasing. After the call to
// SetActiveRange() the client may use the iterator interface HasNext() and
// NextPeptide() to iterate over the window. The client may also use
// GetPeptide() to get a specific peptide in the window.

#include <deque>
#include "peptides.pb.h"
#include "peptide.h"
#include "theoretical_peak_set.h"
#include "fifo_alloc.h"
#include "spectrum_collection.h"
#include "io/OutputFiles.h"

//#include "sp_scorer.h"
#ifndef ACTIVE_PEPTIDE_QUEUE_H
#define ACTIVE_PEPTIDE_QUEUE_H

class TheoreticalPeakCompiler;

class ActivePeptideQueue {
 public:
  ActivePeptideQueue(RecordReader* reader,
            const vector<const pb::Protein*>& proteins);

  ~ActivePeptideQueue();

  bool isWithinIsotope(vector<double>* min_mass, vector<double>* max_mass, double mass, int* isotope_idx);
  
  // See above for usage and .cc for implementation details.
  int SetActiveRange(vector<double>* min_mass, vector<double>* max_mass, double min_range, double max_range, vector<bool>* candidatePeptideStatus);
  int SetActiveRangeBIons(vector<double>* min_mass, vector<double>* max_mass, double min_range, double max_range, vector<bool>* candidatePeptideStatus);

  bool HasNext() const { return iter_ != end_; }
  Peptide* NextPeptide() { return *iter_; }
  const Peptide* GetPeptide(int back_index) const {
    return peptide_centric_ ? current_peptide_ : *(end_ - back_index);
  }
  void SetBinSize(double binWidth, double binOffset) {
    theoretical_b_peak_set_.binWidth_ = binWidth;
    theoretical_b_peak_set_.binOffset_ = binOffset;
  }

  deque<TheoreticalPeakSetBIons> b_ion_queue_;
  deque<TheoreticalPeakSetBIons>::const_iterator iter1_, end1_;
 
  int CountAAFrequency(double binWidth, double binOffset, double** dAAFreqN,
                       double** dAAFreqI, double** dAAFreqC, int** dAAMass);
  //Added by Andy Lin for RESIDUE_EVIDENCE_MATRIX
  int CountAAFrequencyRes(double binWidth, double binOffset, vector<double>& dAAFreqN,
                          vector<double>& dAAFreqI, vector<double>& dAAFreqC, vector<double>& dAAMass);
  
  int ActiveTargets() const { return active_targets_; }
  int ActiveDecoys() const { return active_decoys_; }

  void ReportPeptideHits(Peptide* peptide);
  void SetOutputs(OutputFiles* output_files, const vector<const pb::AuxLocation*>* locations, int top_matches,
                  bool compute_sp, ofstream* target_file, ofstream* decoy_file, double highest_mz) {
      locations_ = locations;
      output_files_ = output_files;
      top_matches_ = top_matches;
      compute_sp_ = compute_sp;
      target_file_ = target_file;
      decoy_file_ = decoy_file;
      highest_mz_ = highest_mz;
  }
  void setPeptideCentric(bool peptide_centric) {
    peptide_centric_ = peptide_centric;
  }
  
  void setElutionWindow(int elution_window) {
    elution_window_ = elution_window;
  }
  // iter_ points to the current peptide. Client access is by HasNext(),
  // GetPeptide(), and NextPeptide(). end_ points just beyond the last active
  // peptide.
  deque<Peptide*>::const_iterator iter_, end_;
  
//  const ProteinVec& proteins_;
 private:
  const vector<const pb::AuxLocation*>* locations_;
  OutputFiles* output_files_;
  int top_matches_;
  bool compute_sp_;
  ofstream* target_file_;
  ofstream* decoy_file_;
  double highest_mz_;
  Peptide* current_peptide_;
  bool exact_pval_search_;
  bool peptide_centric_;
  int elution_window_;


//  Spectrum* spectrum_;
  // IMPLEMENTATION DETAILS

  // See .cc file.
  void ComputeTheoreticalPeaksBack();
  void ComputeBTheoreticalPeaksBack();

  RecordReader* reader_;
  pb::Peptide current_pb_peptide_;

  // All amino acid sequences from which the peptides are drawn.
  const vector<const pb::Protein*>& proteins_; 

  // Workspace for computing theoretical peaks for a single peptide.
  // Gets reused for each new peptide.
  ST_TheoreticalPeakSet theoretical_peak_set_;
  TheoreticalPeakSetBIons theoretical_b_peak_set_;
  
  // The active peptides. Lighter peptides are enqueued before heavy ones.
  // queue_ maintains only the peptides that fall within the range specified
  // by the last call to SetActiveRange().
  deque<Peptide*> queue_;

  // Set by most recent call to SetActiveRange()
  double min_mass_, max_mass_;

  // While we maintain a window of active peptides, we allocate and relase them
  // on a first-in, first-out basis. We use FifoAllocators 
  // (see fifo_alloc.{h,cc}) to manage memory efficiently for this usage 
  // pattern. As we read in peptides and compute the theoretical peaks, we
  // use compiler_prog1 and compiler_prog2 to generate code on the fly for
  // taking dot products with the theoretical peak set for each peptide. 
  // FifoAllocators allow us to execute the code thus generated,
  // since they set the proper permissions. The set of theoretical peaks for 
  // "dotting" with charge 1 and charge 2 spectra, have different
  // FifoAllocators and TheoreticalPeakCompilers.
  FifoAllocator fifo_alloc_peptides_;
  FifoAllocator fifo_alloc_prog1_;
  FifoAllocator fifo_alloc_prog2_;
  TheoreticalPeakCompiler* compiler_prog1_;
  TheoreticalPeakCompiler* compiler_prog2_;

  // Number of targets and decoys in active range
  int active_targets_, active_decoys_;
};

/*
CONSIDER:
Make 3 subclasses:
  Read-as-you-go (below)
  Preread into memory (as current operation)
  Threaded reads (!!)

pb::Peptide* pb_peptide = new pb::Peptide;
CHECK(reader_->Read(pb_peptide));
!reader_->Done()
*/
#endif //ACTIVE_PEPTIDE_QUEUE_H
