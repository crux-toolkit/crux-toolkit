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
#include "theoretical_peak_set.h"
#include "mod_coder.h"
#include "util/Params.h"

#include "spectrum_collection.h"

using namespace std;


// This is the TheoreticalPeakSet we are using at search time.  We may use a
// different subclass  of TheoreticalPeakSet in a final version, pending more
// tests. See Bugzilla #253.
class TheoreticalPeakSetBYSparse;


class TheoreticalPeakCompiler;

// BIG CAUTION: At search time, you CANNOT expect even the IMPLICIT destructor
// to get called!! We actually RELY on the fact that when we use FIFO
// allocation, the destructor won't get called. We expect the destructor to 
// get called only when the object is created using the normal system memory
// allocation. This way, we can use the destructor to clean up the system 
// memory allocated for mods_ when FIFO allocation is NOT used.
class Peptide {
 public:

  // The proteins parameter is presumed to live in memory all the time while the
  // Peptide exists, so that residues_ can refer to the amino acid sequence.
  Peptide(const pb::Peptide& peptide,
          const vector<const pb::Protein*>& proteins,
          vector<const pb::AuxLocation*>* locations = NULL);

  // CAUTION: We do NOT expect this destructor to get called when FIFO 
  // allocation is used. It will get called only when normal system memory
  // is used for allocation. Please see the BIG CAUTION message at the top 
  // of the class definition for details.
  ~Peptide() {
//    delete[] mods_;
  }

  string Seq() const { return string(residues_, Len()); } // For display
  string TargetSeq() const { return string(target_residues_, Len()); } // For display
  /**
 * Gets the protein name with the peptide position appended.
 */
  string GetLocationStr(const string& decoy_prefix);
/**
 * Gets the flanking AAs for a Tide peptide sequence
 */  
  string GetFlankingAAs();

/**
  Get the modifications from the peptide sequence
*/
  void getModifications(int mod_precision, string& mod_crux_string, string& mod_mztab_string);

/**
 * Gets the peptide sequence with modifications
 */
  string SeqWithMods(int mod_precision);

  // void DecodeMod(){  // TODO: to be removed
  //   for (int i = 0; i < num_mods_; ++i) {
  //     int index;
  //     double delta;
  //     MassConstants::DecodeMod(mods_[i], &index, &delta);
  //   }
  // }

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
  void ComputeTheoreticalPeaks(TheoreticalPeakSetBYSparse* workspace, bool dia_mode = false);
  
  int Len() const { return len_; }
  double Mass() const { return mass_; }
  int Id() const { return id_; }
  int FirstLocProteinId() const { return first_loc_protein_id_; }
  int FirstLocPos() const { return first_loc_pos_; }
  int ProteinLenth() const {return protein_length_;}
  vector<ModCoder::Mod> Mods() const {
    return mods_;
  }
//  vector<ModCoder::Mod> Mods() const { return mods_; }

  bool IsDecoy() const { return decoyIdx_ >= 0; }
  int DecoyIdx() const { return decoyIdx_; }
  vector<double> getAAMasses() const;

  static double MassToMz(double mass, int charge) { // Added for Diameter
    return mass / (charge * 1.0) + MASS_PROTON;
  }
  vector<int>& IonMzbins() { return ion_mzbins_; }
  vector<int>& BIonMzbins() { return b_ion_mzbins_; }
  vector<int>& YIonMzbins() { return y_ion_mzbins_; }
  vector<double>& IonMzs() { return ion_mzs_; } // added for debug purpose  
  
  vector<unsigned int> peaks_0;   // Single charged b-y ions, in case of exact p-value, this contains only the b-ions
  vector<unsigned int> peaks_1;   // Double charged b-y ions
  vector<unsigned int> peaks_1b;   // Single charged b ions
  vector<unsigned int> peaks_1y;   // Single charged y ions
  vector<unsigned int> peaks_2b;   // Double charged b ions
  vector<unsigned int> peaks_2y;   // Double charged y ions 

  bool active_;
  
 private:
  template<class W> void AddIons(W* workspace, bool dia_mode = false) ;

  void Compile(const TheoreticalPeakArr* peaks);
  bool find_static_mod(const pb::ModTable* mod_table, char AA, double& mod_mass, string& mod_name); // mod_mass output variable
  bool find_variable_mod(const pb::ModTable* mod_table, char AA, double mod_mass, string& mod_name); // mod_mass output variable
          
  int len_;
  double mass_;
  int id_;
  int first_loc_protein_id_;
  int first_loc_pos_;
  int protein_length_;
  vector<pb::Location> aux_locations;
  const char* residues_;
  const char* target_residues_;
  int num_mods_;
  vector<ModCoder::Mod> mods_;
  double nterm_mod_;
  double cterm_mod_;
  int decoyIdx_;
  string decoy_seq_;

  const vector<const pb::Protein*>* proteins_;

  vector<int> ion_mzbins_, b_ion_mzbins_, y_ion_mzbins_;   // Added for Diameter
  vector<double> ion_mzs_; // added for debug purpose

  string protein_id_str_;
  string flankingAAs_;
  string seq_with_mods_;
  string mod_crux_string_;
  string mod_mztab_string_;
};

#endif // PEPTIDE_H
