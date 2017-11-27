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

#ifdef DEBUG
DEFINE_int32(debug_peptide_id, -1, "Peptide id to debug.");
#endif

#if 0
DEFINE_bool(flanks, true, "Include flanking peaks.");
DEFINE_bool(dups_ok, false, "Don't remove duplicate peaks");
#endif

string Peptide::SeqWithMods() const {
  vector<char> buf(Len() + num_mods_ * 30 + 1);
  int residue_pos = 0;
  char* buf_pos = &(buf[0]);
  for (int i = 0; i < num_mods_; ++i) {
    int index;
    double delta;
    MassConstants::DecodeMod(mods_[i], &index, &delta);
    while (residue_pos <= index)
      *buf_pos++ = residues_[residue_pos++];
    buf_pos += sprintf(buf_pos, "[%s%.1f]", delta >= 0 ? "+" : "", delta);
  }
  while (residue_pos < Len())
    *buf_pos++ = residues_[residue_pos++];
  *buf_pos = '\0';
  return &(buf[0]);
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
void Peptide::AddIons(W* workspace) const {
  // Use workspace to assemble all B and Y ions. workspace will determine
  // which, if any, associated ions will be represented.
  double max_possible_peak = numeric_limits<double>::infinity();
  if (MaxBin::Global().MaxBinEnd() > 0)
    max_possible_peak = MaxBin::Global().CacheBinEnd();

  vector<double> aa_masses(Len());
  const char* residue = residues_;
  // Collect m/z values for each residue, for z = 1, 2.
  for (int i = 0; i < Len(); ++i, ++residue) {
    if (i == 0) { // nterm static pep
      aa_masses[i] = MassConstants::nterm_mono_table[*residue];
    } else if (i == Len() - 1) { // cterm static pep
      aa_masses[i] = MassConstants::cterm_mono_table[*residue];
    } else { // all other mods
      aa_masses[i] = MassConstants::mono_table[*residue];
    }
  }

  for (int i = 0; i < num_mods_; ++i) {
    int index;
    double delta;
    MassConstants::DecodeMod(mods_[i], &index, &delta);
    aa_masses[index] += delta;
  }

  // Add all charge 1 B ions.
  double total = aa_masses[0];
  for (int i = 1; i < Len() && total <= max_possible_peak; ++i) {
    workspace->AddBIon(total, 1);
    total += aa_masses[i];
  }

  // Add all charge 1 Y ions.
  total = aa_masses[Len() - 1];
  for (int i = Len()-2; i >= 0 && total <= max_possible_peak; --i) {
    workspace->AddYIon(total, 1);
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
}

template<class W>
void Peptide::AddBIonsOnly(W* workspace) const {
  // Use workspace to assemble b ions only.
  // Intended primarily to support XCorr p-value calculations.
  double max_possible_peak = numeric_limits<double>::infinity();
  if (MaxBin::Global().MaxBinEnd() > 0) {
    max_possible_peak = MaxBin::Global().CacheBinEnd();
  }
  
  vector<double> aa_masses(Len());
  const char* residue = residues_;
  // Collect m/z values for each residue, for z = 1.
  for (int i = 0; i < Len(); ++i, ++residue) {
    if (i == 0) { // nterm static pep
      aa_masses[i] = MassConstants::nterm_mono_table[*residue];
    } else if (i == Len() - 1) { // cterm static pep
      aa_masses[i] = MassConstants::cterm_mono_table[*residue];
    } else { // all other mods
      aa_masses[i] = MassConstants::mono_table[*residue];
    }
  }

  //Add modifications to amino acids
  for (int i = 0; i < num_mods_; ++i) {
    int index;
    double delta;
    MassConstants::DecodeMod(mods_[i], &index, &delta);
    aa_masses[index] += delta;
  }

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
  int pos_size = peaks[0].size();
  prog1_ = compiler_prog1->Init(pos_size, 0);
  compiler_prog1->AddPositive(peaks[0]);
//  compiler_prog1->AddPositive(pb_peptide.peak1());
//  compiler_prog1->AddNegative(pb_peptide.neg_peak1());
  compiler_prog1->Done();

  pos_size = peaks[0].size() + peaks[1].size();
  prog2_ = compiler_prog2->Init(pos_size, 0);
  compiler_prog2->AddPositive(peaks[0]);
  compiler_prog2->AddPositive(peaks[1]);
//  compiler_prog2->AddPositive(pb_peptide.peak2());
//  compiler_prog2->AddNegative(pb_peptide.neg_peak2());
  compiler_prog2->Done();
/*    cout << Seq() << endl;
    for (int i = 0; i < peaks[0].size(); ++i)
      cout << "Theoretical Peak[" << peaks[0][i].Bin() << "] = "
           << peaks[0][i].Type() << endl;
    for (int i = 0; i < peaks[1].size(); ++i)
      cout << "Theoretical Peak[" << peaks[1][i].Bin() << "] = "
           << peaks[1][i].Type() << endl;
*/
//	exit(1);  
}

void Peptide::ComputeTheoreticalPeaks(TheoreticalPeakSet* workspace) const {
  AddIons<TheoreticalPeakSet>(workspace);   // Generic workspace
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

void Peptide::ComputeTheoreticalPeaks(ST_TheoreticalPeakSet* workspace,
                                      const pb::Peptide& pb_peptide,
                                      TheoreticalPeakCompiler* compiler_prog1,
                                      TheoreticalPeakCompiler* compiler_prog2) {
  // Search-time fast workspace
  AddIons<ST_TheoreticalPeakSet>(workspace);

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
double* Peptide::getAAMasses(){
  double* masses_charge = new double[Len()];
  const char* residue = residues_;
  for (int i = 0; i < Len(); ++i, ++residue) {
    if (i == 0) { // nterm static pep
      masses_charge[i] = MassConstants::nterm_mono_table[*residue];
    } else if (i == Len() - 1) { // cterm static pep
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

// Probably defunct, uses old calling format.
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
