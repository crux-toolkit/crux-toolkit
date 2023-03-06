// Benjamin Diament

#include <iostream>
#include <limits>
#include <gflags/gflags.h>
#include "mass_constants.h"
#include "max_mz.h"
#include "fifo_alloc.h"
#include "fixed_cap_array.h"
#include "theoretical_peak_set.h"
#include "peptide.h"
#include "compiler.h"
#include "util/StringUtils.h"

#ifdef DEBUG
DEFINE_int32(debug_peptide_id, -1, "Peptide id to debug.");
#endif

#if 0
DEFINE_bool(flanks, true, "Include flanking peaks.");
DEFINE_bool(dups_ok, false, "Don't remove duplicate peaks");
#endif

string Peptide::SeqWithMods() const {
  string seq_with_mods = string(residues_, Len());
  
  if (num_mods_ > 0) {
    int mod_pos_offset = 0;
    string mod_str;
    vector<int> mod;
    
    for (int i = 0; i < num_mods_; ++i) {
      mod.push_back(mods_[i]);
    }
    
    sort(mod.begin(), mod.end());
    
    for (int i = 0; i < num_mods_; ++i) {
      int index;
      double delta;
      MassConstants::DecodeMod(mod[i], &index, &delta);
      mod_str = '[' + StringUtils::ToString(delta, mod_precision_) + ']';
      seq_with_mods.insert(index + 1 + mod_pos_offset, mod_str);
      mod_pos_offset += mod_str.length();
    }
  }
  return seq_with_mods;  
}

void Peptide::Show() {
#ifdef DEBUG
  if (Id() == FLAGS_debug_peptide_id) {
    cout << Seq() << endl;
    cout << "Charge 1 Pos" << endl;
    for (int i = 0; i < peaks_charge_1_.size(); ++i)
      cout << "Theoretical Peak[" << peaks_charge_1_[i].Bin() << "] = "
           << peaks_charge_1_[i].Type() << endl;
    cout << "Charge 1 Neg" << endl;
    for (int i = 0; i < negs_charge_1_.size(); ++i)
      cout << "Theoretical Peak[" << negs_charge_1_[i].Bin() << "] = "
           << negs_charge_1_[i].Type() << endl;
    cout << "Charge 2 Pos" << endl;
    for (int i = 0; i < peaks_charge_2_.size(); ++i)
      cout << "Theoretical Peak[" << peaks_charge_2_[i].Bin() << "] = "
             << peaks_charge_2_[i].Type() << endl;
    cout << "Charge 2 Neg" << endl;
    for (int i = 0; i < negs_charge_2_.size(); ++i)
      cout << "Theoretical Peak[" << negs_charge_2_[i].Bin() << "] = "
           << negs_charge_2_[i].Type() << endl;
    cout << endl;
  }
#endif
}

template<class W>
void Peptide::AddIons(W* workspace, bool dia_mode) {
  // Use workspace to assemble all B and Y ions. workspace will determine
  // which, if any, associated ions will be represented.
  double max_possible_peak = numeric_limits<double>::infinity();
  if (MaxBin::Global().MaxBinEnd() > 0)
    max_possible_peak = MaxBin::Global().CacheBinEnd();

  vector<double> aa_masses = getAAMasses();
  // carp(CARP_DETAILED_DEBUG, "**********aa_masses:%s", StringUtils::JoinDoubleVec(aa_masses, ',').c_str() );


  // added by Yang
  if (dia_mode) {
    ion_mzs_.clear(); ion_mzbins_.clear(); b_ion_mzbins_.clear(); y_ion_mzbins_.clear();
  }

  // Add all charge 1 B ions.
  double total = aa_masses[0];
  for (int i = 1; i < Len() && total <= max_possible_peak; ++i) {
    workspace->AddBIon(total, 1);
    if (dia_mode) {
      b_ion_mzbins_.push_back(MassConstants::mass2bin(Peptide::MassToMz(total + MassConstants::B, 1)));
      ion_mzbins_.push_back(MassConstants::mass2bin(Peptide::MassToMz(total + MassConstants::B, 1)));
      ion_mzs_.push_back(Peptide::MassToMz(total + MassConstants::B, 1));
    }
    total += aa_masses[i];
  }

  // Add all charge 1 Y ions.
  total = aa_masses[Len() - 1];
  for (int i = Len()-2; i >= 0 && total <= max_possible_peak; --i) {
    workspace->AddYIon(total, 1);
    if (dia_mode) {
      y_ion_mzbins_.push_back(MassConstants::mass2bin(Peptide::MassToMz(total + MassConstants::Y, 1)));
      ion_mzbins_.push_back(MassConstants::mass2bin(Peptide::MassToMz(total + MassConstants::Y, 1)));
      ion_mzs_.push_back(Peptide::MassToMz(total + MassConstants::Y, 1));
    }
    total += aa_masses[i];
  }

  // Add all charge 2 B ions.
  max_possible_peak = max_possible_peak*2 + 2;  //adjust for larger charge
  total = aa_masses[0];
  for (int i = 1; i < Len() && total <= max_possible_peak; ++i) {
    workspace->AddBIon(total, 2);
    total += aa_masses[i];
  }

  // Add all charge 2 Y ions.
  total = aa_masses[Len() - 1];
  for (int i = Len()-2; i >= 0 && total <= max_possible_peak; --i) {
    workspace->AddYIon(total, 2);
    total += aa_masses[i];
  }

  // added by Yang
  if (dia_mode) {
    sort(b_ion_mzbins_.begin(), b_ion_mzbins_.end());
    sort(y_ion_mzbins_.begin(), y_ion_mzbins_.end());
    sort(ion_mzbins_.begin(), ion_mzbins_.end());
    sort(ion_mzs_.begin(), ion_mzs_.end());
  }

}

template<class W>
void Peptide::AddBIonsOnly(W* workspace) const {
  // Use workspace to assemble b ions only.
  // Intended primarily to support XCorr p-value calculations.
  double max_possible_peak = numeric_limits<double>::infinity();
  if (MaxBin::Global().MaxBinEnd() > 0) {
    max_possible_peak = MaxBin::Global().CacheBinEnd();
  }
  
  vector<double> aa_masses = getAAMasses();

  // Add all charge 1 B ions.
  double total = MASS_PROTON + aa_masses[0];
  for (int i = 1; i < Len() && total <= max_possible_peak; ++i) {
    workspace -> AddBIon(total);
    total += aa_masses[i];
  }
}

#ifdef DEBUG
void DisAsm(const void* prog) {
  unsigned char* pos = (unsigned char*) prog;
  while (true) {
    switch(*pos) {
    case 3: assert(pos[1] == 130); cout << "add " << *((int*) (void*) (pos+2)) << "(%edx), %eax\n"; pos += 6; break;
    case 43: assert(pos[1] == 130); cout << "sub " << *((int*) (void*) (pos+2)) << "(%edx), %eax\n"; pos += 6; break;
    case 195: cout << "ret" << endl; goto out;
    default: assert(false);
    }
  }
 out:
  return;
}
#endif

void Peptide::Compile(const TheoreticalPeakArr* peaks,
                      const pb::Peptide& pb_peptide,
                      TheoreticalPeakCompiler* compiler_prog1,
                      TheoreticalPeakCompiler* compiler_prog2) {
#ifdef CPP_SCORING
  // Store the theoretical peak indeces for the peptide in a vector. 
  int i;
  peaks_0.clear();
  peaks_1.clear();
  peaks_0.reserve(peaks[0].size());
  peaks_1.reserve(peaks[1].size());
  
  for (i = 0; i < peaks[0].size(); ++i) {
    peaks_0.push_back(peaks[0][i]);
  }
  for (i = 0; i < peaks[1].size(); ++i) {
    peaks_1.push_back(peaks[1][i]);
  }
  prog1_ = (void*)3;  //Mark the peptide that it contains theoretical peaks to score
#else
  // Create the Assembly code for scoring. 
  int pos_size = peaks[0].size();
  prog1_ = compiler_prog1->Init(pos_size, 0);
  compiler_prog1->AddPositive(peaks[0]);
  compiler_prog1->Done();

  pos_size = peaks[0].size() + peaks[1].size();
  prog2_ = compiler_prog2->Init(pos_size, 0);
  compiler_prog2->AddPositive(peaks[0]);
  compiler_prog2->AddPositive(peaks[1]);
  compiler_prog2->Done();
#endif
}

void Peptide::ComputeTheoreticalPeaks(TheoreticalPeakSetBYSparse* workspace, bool dia_mode) {
  AddIons<TheoreticalPeakSetBYSparse>(workspace, dia_mode);   // Generic workspace
#ifdef DEBUG
  Show();
#endif
}

void Peptide::ComputeBTheoreticalPeaks(TheoreticalPeakSetBIons* workspace) const {
  AddBIonsOnly<TheoreticalPeakSetBIons>(workspace);   // workspace for b ion only peak set
#ifdef DEBUG
  Show();
#endif
}

void Peptide::ComputeTheoreticalPeaks(TheoreticalPeakSetBYSparse* workspace,
                                      const pb::Peptide& pb_peptide,
                                      TheoreticalPeakCompiler* compiler_prog1,
                                      TheoreticalPeakCompiler* compiler_prog2,
                                      bool dia_mode) {
  // Search-time fast workspace
  AddIons<TheoreticalPeakSetBYSparse>(workspace, dia_mode);

#if 0
  TheoreticalPeakArr peaks[2];
  peaks[0].Init(2000);
  peaks[1].Init(2000);
  workspace->GetPeaks(&peaks[0], NULL, &peaks[1], NULL, NULL);
  Compile(peaks, pb_peptide, compiler_prog1, compiler_prog2);
#endif

  Compile(workspace->GetPeaks(), pb_peptide, compiler_prog1, compiler_prog2);
#ifdef DEBUG
  if (Id() == FLAGS_debug_peptide_id) {
    cout << "Prog1:" << endl;
    DisAsm(prog1_);
    cout << "Prog2:" << endl;
    DisAsm(prog2_);
  }
#endif
}

// return the amino acid masses in the current peptide
vector<double> Peptide::getAAMasses() const {
  vector<double> masses_charge(Len());
  const char* residue = residues_;
  for (int i = 0; i < Len(); ++i, ++residue) {
    if (i == 0) { // nterm static pep
      if(first_loc_pos_ == 0)
        masses_charge[i] = MassConstants::nprotterm_mono_table[*residue];
      else
        masses_charge[i] = MassConstants::nterm_mono_table[*residue];
    } else if (i == Len() - 1) { // cterm static pep
      if(first_loc_pos_ + len_ == protein_length_ - 1)
        masses_charge[i] = MassConstants::cprotterm_mono_table[*residue];
      else
        masses_charge[i] = MassConstants::cterm_mono_table[*residue];
    } else { // all other mods
      masses_charge[i] = MassConstants::mono_table[*residue];
    }
  }

  for (int i = 0; i < num_mods_; ++i) {
    int index;
    double delta;
    MassConstants::DecodeMod(mods_[i], &index, &delta);
    masses_charge[index] += delta;
  }

  return masses_charge;
}
/**
 * Gets the protein name with the peptide position appended.
 */
 void  Peptide::GetLocationStr(
  const vector<const pb::Protein*>& proteins, 
  const string& decoy_prefix, 
  string& locations
) const {
  
  locations = (IsDecoy()?decoy_prefix:"") + 
    proteins[FirstLocProteinId()]->name() + 
    "(" + std::to_string(FirstLocPos()+1) + ")";
  
  for (vector<pb::Location>::const_iterator 
    loc = aux_locations.begin(); 
    loc != aux_locations.end();
    ++loc
  ) {
    const pb::Protein* protein = proteins[(*loc).protein_id()];
    int pos = (*loc).pos()+1;
    locations += "," + (IsDecoy()?decoy_prefix:"") + protein->name() + 
      "(" + std::to_string(pos) + ")";    
  }
}

/**
 * Gets the flanking AAs for a Tide peptide sequence
 */
void  Peptide::GetFlankingAAs(
  const vector<const pb::Protein*>& proteins,
  string& flankingAAs
) const {
  
  flankingAAs.clear();
  int pos = FirstLocPos();
  const string& seq = proteins[FirstLocProteinId()]->residues();

  flankingAAs = ((pos > 0) ? proteins[FirstLocProteinId()]->residues().substr(pos-1, 1) : "-") +
    ((pos+Len() <  proteins[FirstLocProteinId()]->residues().length()) ? 
    proteins[FirstLocProteinId()]->residues().substr(pos+Len(),1) : "-");
    
  for (vector<pb::Location>::const_iterator
    loc = aux_locations.begin();
    loc != aux_locations.end();
    ++loc
  ) {
    const pb::Protein* protein = proteins[(*loc).protein_id()];
    pos = (*loc).pos();
    flankingAAs += "," + ((pos > 0) ? protein->residues().substr(pos-1, 1) : "-") + 
      ((pos+Len() <  protein->residues().length()) ? 
      protein->residues().substr(pos+Len(),1) : "-");
  }
  
}

// Probably defunct, uses old calling format.
/*
int NoInlineDotProd(Peptide* peptide, const int* cache, int charge) {
  const void* prog = peptide->Prog(charge);
  int result;
#ifdef _MSC_VER
  // FIXME CEG add Windows compatible inline assembly
#else
  __asm__ __volatile__("call *%[prog]\n"
                       : "=a" (result)
                       : "d" (cache), [prog] "abcSD" (prog));
#endif
  return result;
}
*/