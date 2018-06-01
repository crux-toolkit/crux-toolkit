// Benjamin Diament

// Peptide is a class for managing the representation of a peptide and its
// associated theoretical peaks for use at search time. It encapsulates a
// diverse set of functions, interacting closely with TheoreticalPeakSet,
// TheoreticalPeakCompiler, ActivePeptideQueue, and DotProd() (defined in
// search.cc). At indexing time this class is used with TheoreticalPeakSet to
// store peaks to a protocol buffer.
//
// TODO 255: the intractions with other classes are complicated, and the various
// functionalities here could probably be encapsulated with cleaner interfaces.
// Perhaps one way to clean this up would be to make two different Peptide
// classes: one for indexing and another for searching.
//
// At search time, ActivePeptideQueue (see active_peptide_queue.{h,cc})
// manages the task of reading Peptides from storage and keeping current
// candidate Peptides in memory.  Each Peptide knows how to take the dot
// product of its theoretical peaks with a given observed spectrum.
// Specifically, ComputeTheoreticalPeaks() generates a compiled program for
// doing so.  A different version of the generated program exists for charge 1
// and  charge 2.

#ifndef PEPTIDE_H
#define PEPTIDE_H

#include <iostream>
#include <vector>
#include "raw_proteins.pb.h"
#include "peptides.pb.h"
#include "theoretical_peak_pair.h"
#include "fifo_alloc.h"
#include "mod_coder.h"
#include "sp_scorer.h"

#include "spectrum_collection.h"
//#include "TideMatchSet.h"

using namespace std;

class FifoAllocator;
class TheoreticalPeakSet;
class TheoreticalPeakSetBIons;
class SpScorer;


// This is the TheoreticalPeakSet we are using at search time.  We may use a
// different subclass  of TheoreticalPeakSet in a final version, pending more
// tests. See Bugzilla #253.
class TheoreticalPeakSetBYSparse;
typedef TheoreticalPeakSetBYSparse ST_TheoreticalPeakSet; // ST="search time"
// This is an alternative TheoreticalPeakSet
// class TheoreticalPeakSetBYSparseOrdered;
// typedef TheoreticalPeakSetBYSparseOrdered ST_TheoreticalPeakSet;
// This alternative TheoreticalPeakSet is slow but slightly more flexible.
// class TheoreticalPeakSetMakeAll;
// typedef TheoreticalPeakSetMakeAll ST_TheoreticalPeakSet; // ST="search time"

class TheoreticalPeakCompiler;

// BIG CAUTION: At search time, you CANNOT expect even the IMPLICIT destructor
// to get called!! We actually RELY on the fact that when we use FIFO
// allocation, the destructor won't get called. We expect the destructor to 
// get called only when the object is created using the normal system memory
// allocation. This way, we can use the destructor to clean up the system 
// memory allocated for mods_ when FIFO allocation is NOT used.
class Peptide {
 public:
  // The proteins parameter is presumed to live in memory all the while the
  // Peptide exists, so that residues_ can refer to the amino acid sequence.
  Peptide(const pb::Peptide& peptide,
          const vector<const pb::Protein*>& proteins,
          FifoAllocator* fifo_alloc = NULL)
    : len_(peptide.length()), mass_(peptide.mass()), id_(peptide.id()),
    first_loc_protein_id_(peptide.first_location().protein_id()),
    first_loc_pos_(peptide.first_location().pos()), 
    has_aux_locations_index_(peptide.has_aux_locations_index()),
    aux_locations_index_(peptide.aux_locations_index()),
    mods_(NULL), num_mods_(0), decoyIdx_(peptide.has_decoy_index() ? peptide.decoy_index() : -1),
    prog1_(NULL), prog2_(NULL) {
    // Set residues_ by pointing to the first occurrence in proteins.
    residues_ = proteins[first_loc_protein_id_]->residues().data() 
                    + first_loc_pos_;
    if (peptide.modifications_size() > 0) {
      num_mods_ = peptide.modifications_size();
      if (fifo_alloc) {
        mods_ = (ModCoder::Mod*) fifo_alloc->New(sizeof(mods_[0]) * num_mods_);
      } else {
        mods_ = new ModCoder::Mod[num_mods_];
      }
      for (int i = 0; i < num_mods_; ++i)
        mods_[i] = ModCoder::Mod(peptide.modifications(i));
    }
  }
  class spectrum_matches {
   public:
      spectrum_matches(Spectrum* spectrum, double score1, double score2,
          int score3, int charge) {
          spectrum_ = spectrum;
          score1_ = score1;
          score2_ = score2;
          score3_ = score3;
          charge_ = charge;
          d_cn_ = 0.0;
          elution_score_ = 0.0;
      }
      spectrum_matches() {
          spectrum_ = NULL;
          score1_ = 0.0;
          score2_ = 0.0;
          score3_ = 0.0;
          charge_ = 0;
          d_cn_ = 0.0;
          elution_score_ = 0.0;
      }
      Spectrum* spectrum_;
      double score1_;
      double score2_;
      int score3_;
      int charge_;
      double d_cn_;
      double elution_score_;
      SpScorer::SpScoreData spData_;

      static bool compPV(const spectrum_matches &a, const spectrum_matches &b) {
          return a.score1_ < b.score1_;
      }
      static bool compSC(const spectrum_matches &a, const spectrum_matches &b) {
          return a.score1_ > b.score1_;
      }
      static bool compRT(const spectrum_matches &a, const spectrum_matches &b) {
        return a.spectrum_->RTime() < b.spectrum_->RTime();
      }
      static bool compES(const spectrum_matches &a, const spectrum_matches &b) {
        return a.elution_score_ < b.elution_score_;
      }
  };
  vector<spectrum_matches> spectrum_matches_array;
  void AddHit(Spectrum* spectrum, double score1, double score2,
          int score3, int charge) {

    spectrum_matches_array.push_back(spectrum_matches(spectrum,
                                     score1, score2, score3, charge));
  }

  // CAUTION: We do NOT expect this destructor to get called when FIFO 
  // allocation is used. It will get called only when normal system memory
  // is used for allocation. Please see the BIG CAUTION message at the top 
  // of the class definition for details.
  ~Peptide() {
    delete[] mods_;
  }

  // Allocation by FifoAllocator
  void* operator new(size_t size, FifoAllocator* fifo_alloc) {
    return fifo_alloc->New(size);
  }

  string Seq() const { return string(residues_, Len()); } // For display

  string SeqWithMods() const;

  // Compute and cache set of theoretical peaks using the provided workspace.
  // Workspace exists for efficiency: it can be reused by another Peptide
  // without reallocating memory. The pb_peptide should be the same as supplied
  // during construction, but this class doesn't own or store the underlying
  // pb::Peptide, as it is probably not FifoAllocated.
  // The first version takes a generic workspace (good for index generation and
  // testing). The second version is intended to be faster at search time: a
  // specific subclass of TheoreticalPeakSet is given, which should avoid
  // associated virtual method calls. (TODO 256: are virtual method calls indeed
  // avoided?).  The second version also produces the compiled programs for
  // taking dot products.
  void ComputeTheoreticalPeaks(TheoreticalPeakSet* workspace) const;
  void ComputeTheoreticalPeaks(ST_TheoreticalPeakSet* workspace,
                               const pb::Peptide& pb_peptide,
                               TheoreticalPeakCompiler* compiler_prog1,
                               TheoreticalPeakCompiler* compiler_prog2);
  void ComputeBTheoreticalPeaks(TheoreticalPeakSetBIons* workspace) const;

  // Return the appropriate program depending on the precursor charge.
  // TODO 257: fix the unfortunate use of max_charge.
  const void* Prog(int max_charge) const {
    return max_charge <= 2 ? prog1_ : prog2_;
  }

  void ReleaseFifo(FifoAllocator* fifo_alloc_prog1,
       FifoAllocator* fifo_alloc_prog2) {
    // TODO 258: this code should probably move to ActivePeptideQueue
    if (prog1_ != NULL) {
      fifo_alloc_prog1->Release(prog1_);
    } else {
      fifo_alloc_prog1->ReleaseAll();
    }
    if (prog2_ != NULL) {
      fifo_alloc_prog2->Release(prog2_);
    } else {
      fifo_alloc_prog2->ReleaseAll();
    }
  }

  int Len() const { return len_; }
  double Mass() const { return mass_; }
  int Id() const { return id_; }
  int FirstLocProteinId() const { return first_loc_protein_id_; }
  int FirstLocPos() const { return first_loc_pos_; }
  bool HasAuxLocationsIndex() const { return has_aux_locations_index_; }
  int AuxLocationsIndex() const { return aux_locations_index_; }
  int Mods(const ModCoder::Mod** mods) const {
    *mods = mods_;
    return num_mods_;
  }
  ModCoder::Mod* Mods() const { return mods_; }
  bool IsDecoy() const { return decoyIdx_ >= 0; }
  int DecoyIdx() const { return decoyIdx_; }
  double* getAAMasses();

 private:
  template<class W> void AddIons(W* workspace) const;
  template<class W> void AddBIonsOnly(W* workspace) const;

  void Compile(const TheoreticalPeakArr* peaks,
               const pb::Peptide& pb_peptide,
               TheoreticalPeakCompiler* compiler_prog1,
               TheoreticalPeakCompiler* compiler_prog2);
          

  void Show();

  int len_;
  double mass_;
  int id_;
  int first_loc_protein_id_;
  int first_loc_pos_;
  bool has_aux_locations_index_;
  int aux_locations_index_;
  const char* residues_;
  int num_mods_;
  ModCoder::Mod* mods_;
  int decoyIdx_;

  void* prog1_;
  void* prog2_;
};

#endif // PEPTIDE_H
