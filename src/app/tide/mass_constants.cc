#include <limits>
#include "io/carp.h"
#include "mass_constants.h"
#include "header.pb.h"
#include "stdio.h"

#include "util/AminoAcidUtil.h"

using namespace std;

double const MassConstants::kFixedPointScalar = 1e5; 

const double MassConstants::elts_mono[] = {
  1.007825035, // H
  12.0,        // C
  14.003074,   // N
  15.99491463, // O
  30.973762,   // P
  31.9720707,  // S
  79.9165196   // Se
};

const double MassConstants::elts_avg[] = {
  1.00794,    // H
  12.0107,    // C
  14.0067,    // N
  15.9994,    // O
  30.973761,  // P
  32.065,     // S
  78.96       // Se
};
    
double MassConstants::mono_table[256];
double MassConstants::avg_table[256];
double MassConstants::nterm_mono_table[256];
double MassConstants::cterm_mono_table[256];
double MassConstants::nterm_avg_table[256];
double MassConstants::cterm_avg_table[256];
double MassConstants::nprotterm_mono_table[256];
double MassConstants::cprotterm_mono_table[256];
double MassConstants::nprotterm_avg_table[256];
double MassConstants::cprotterm_avg_table[256];
//double* MassConstants::aa_mass_table = NULL;
//double MassConstants::aa_bin_1[256];
//double MassConstants::aa_bin_2[256];

const double MassConstants::mono_h2o = 2*MassConstants::elts_mono[0] + MassConstants::elts_mono[3];
const double MassConstants::avg_h2o = 2*MassConstants::elts_avg[0] + MassConstants::elts_avg[3];
const double MassConstants::mono_nh3 = 3*MassConstants::elts_mono[0] + MassConstants::elts_mono[2];
const double MassConstants::mono_co = MassConstants::elts_mono[1] + MassConstants::elts_mono[3];
const double MassConstants::mono_oh = MassConstants::elts_mono[0] + MassConstants::elts_mono[3];
const double MassConstants::mono_h = MassConstants::elts_mono[0] ;
const double MassConstants::A = 0 - 28.0;
//const double MassConstants::B_H2O = 0 - MassConstants::mono_h2o;
//const double MassConstants::B_NH3 = 0 - MassConstants::mono_nh3;
const double MassConstants::B = 0.0;
//const double MassConstants::Y_H2O = MassConstants::mono_h2o - 18.0;
//const double MassConstants::Y_NH3 = MassConstants::mono_h2o - 17.0;
const double MassConstants::Y = MassConstants::mono_h2o;
/*const double MassConstants::BIN_SHIFT_A_ION_CHG_1 = 28;
const double MassConstants::BIN_SHIFT_A_ION_CHG_2 = 14;
const double MassConstants::BIN_SHIFT_H2O_CHG_1 = 18;
const double MassConstants::BIN_SHIFT_H2O_CHG_2 = 9;
const double MassConstants::BIN_SHIFT_NH3_CHG_1 = 17;
const double MassConstants::BIN_SHIFT_NH3_CHG_2_CASE_A = 9;
const double MassConstants::BIN_SHIFT_NH3_CHG_2_CASE_B = 8;
*/
double MassConstants::BIN_H2O = 18;
double MassConstants::BIN_NH3 = 17;
double MassConstants::BIN_CO = 28;


//default parameter settings, will be changed during parameter parsing
double MassConstants::bin_width_ = BIN_WIDTH; // 1.0005079;
double MassConstants::bin_offset_ = BIN_OFFSET; //0.40;

FixPt MassConstants::fixp_mono_table[256];
FixPt MassConstants::fixp_avg_table[256];
FixPt MassConstants::fixp_nterm_mono_table[256];
FixPt MassConstants::fixp_cterm_mono_table[256];
FixPt MassConstants::fixp_nterm_avg_table[256];
FixPt MassConstants::fixp_cterm_avg_table[256];
FixPt MassConstants::fixp_nprotterm_mono_table[256];
FixPt MassConstants::fixp_cprotterm_mono_table[256];
FixPt MassConstants::fixp_nprotterm_avg_table[256];
FixPt MassConstants::fixp_cprotterm_avg_table[256];

const FixPt MassConstants::fixp_mono_h2o = ToFixPt(MassConstants::mono_h2o);
const FixPt MassConstants::fixp_avg_h2o  = ToFixPt(MassConstants::avg_h2o);
const FixPt MassConstants::fixp_mono_nh3 = ToFixPt(MassConstants::mono_nh3);
const FixPt MassConstants::fixp_mono_co  = ToFixPt(MassConstants::mono_co);
const FixPt MassConstants::fixp_proton   = ToFixPt(MASS_PROTON);

MassConstants::FixPtTableSet MassConstants::mono_tables{MassConstants::fixp_mono_table, 
                      MassConstants::fixp_nterm_mono_table,
                      MassConstants::fixp_cterm_mono_table,
                      MassConstants::fixp_nprotterm_mono_table,
                      MassConstants::fixp_cprotterm_mono_table  };

MassConstants::FixPtTableSet MassConstants::avg_tables{MassConstants::fixp_avg_table, 
                      MassConstants::fixp_nterm_avg_table,
                      MassConstants::fixp_cterm_avg_table,
                      MassConstants::fixp_nprotterm_avg_table,
                      MassConstants::fixp_cprotterm_avg_table  };

const pb::ModTable* MassConstants::mod_table_ = NULL; 
const pb::ModTable* MassConstants::n_mod_table_ = NULL; 
const pb::ModTable* MassConstants::c_mod_table_ = NULL; 
const pb::ModTable* MassConstants::nprot_mod_table_ = NULL; 
const pb::ModTable* MassConstants::cprot_mod_table_ = NULL;  

ModCoder MassConstants::mod_coder_;
vector<double> MassConstants::unique_deltas_;
//double* MassConstants::unique_deltas_bin_;

static bool CheckModTable(const pb::ModTable& mod_table);

void MassConstants::FillMassTable(const double* elements, double* table) {
  const double* e = elements;
  double H = e[0], C = e[1], N = e[2], O = e[3], P = e[4], S = e[5], Se = e[6];

  for (int i = 0; i < 256; ++i)
    table[i] = numeric_limits<double>::signaling_NaN();

  table['A'] = C*3  + H*5  + N   + O       ;
  table['C'] = C*3  + H*5  + N   + O   + S ;
  table['D'] = C*4  + H*5  + N   + O*3     ;
  table['E'] = C*5  + H*7  + N   + O*3     ;
  table['F'] = C*9  + H*9  + N   + O       ;
  table['G'] = C*2  + H*3  + N   + O       ;
  table['H'] = C*6  + H*7  + N*3 + O       ;
  table['I'] = C*6  + H*11 + N   + O       ;
  table['J'] = C*6  + H*11 + N   + O       ;
  table['K'] = C*6  + H*12 + N*2 + O       ;
  table['L'] = C*6  + H*11 + N   + O       ;
  table['M'] = C*5  + H*9  + N   + O   + S ;
  table['N'] = C*4  + H*6  + N*2 + O*2     ;
  table['O'] = C*12 + H*21 + N*3 + O*3     ;
  table['P'] = C*5  + H*7  + N   + O       ;
  table['Q'] = C*5  + H*8  + N*2 + O*2     ;
  table['R'] = C*6  + H*12 + N*4 + O       ;
  table['S'] = C*3  + H*5  + N   + O*2     ;
  table['T'] = C*4  + H*7  + N   + O*2     ;
  table['U'] = C*3  + H*7  + N   + O*2 + Se;
  table['V'] = C*5  + H*9  + N   + O       ;
  table['W'] = C*11 + H*10 + N*2 + O       ;
  table['Y'] = C*9  + H*9  + N   + O*2     ;

  table['w'] =        H*2        + O       ; // water
  table['a'] =        H*3  + N             ; // ammonia
  table['c'] = C                 + O       ; // carbon monoxide
}


bool MassConstants::Init(const pb::ModTable* mod_table, 
  const pb::ModTable* n_mod_table, 
  const pb::ModTable* c_mod_table, 
  const pb::ModTable* nprot_mod_table, 
  const pb::ModTable* cprot_mod_table, 
  const double bin_width, const double bin_offset) {

  if (mod_table && !CheckModTable(*mod_table))
    return false;

  if (!n_mod_table || !c_mod_table) {
    carp(CARP_FATAL, "We could not find nterm or cterm mod tables. "
    "This is a [relatively] new requirement in attempt to fix static "
    "mod discrepancies - 5/19/2015.");
  }

  carp(CARP_DEBUG, "Checking peptide N-terminal mods table.");
  if (n_mod_table && !CheckModTable(*n_mod_table))
    return false;

  carp(CARP_DEBUG, "Checking peptide C-terminal mods table.");
  if (c_mod_table && !CheckModTable(*c_mod_table))
    return false;
  
  carp(CARP_DEBUG, "Checking protein N-terminal mods table.");
  if(nprot_mod_table != nullptr && !CheckModTable(*nprot_mod_table))
    return false;

  carp(CARP_DEBUG, "Checking protein C-terminal mods table.");
  if(cprot_mod_table != nullptr && !CheckModTable(*cprot_mod_table))
    return false;

  for (int i = 0; i < 256; ++i) {
    mono_table[i] = avg_table[i] = nterm_mono_table[i] = 
    cterm_mono_table[i] = nterm_avg_table[i] = cterm_avg_table[i] = 0;
    if(nprot_mod_table != nullptr && cprot_mod_table != nullptr) {
      nprotterm_avg_table[i] = 0;
      cprotterm_avg_table[i] = 0;
      nprotterm_mono_table[i] = 0;
      cprotterm_mono_table[i] = 0;
    }
  }

  // Store the modification tables. Thay are used in reporting peptide modifications in TideMatchSet
  mod_table_ = mod_table; 
  n_mod_table_ = n_mod_table; 
  c_mod_table_ = c_mod_table; 
  nprot_mod_table_ = nprot_mod_table; 
  cprot_mod_table_ = cprot_mod_table;  


  //initialize all tables with unmodified AA MW
  FillMassTable(elts_mono, mono_table);
  FillMassTable(elts_avg, avg_table);
  FillMassTable(elts_mono, nterm_mono_table);
  FillMassTable(elts_mono, cterm_mono_table);
  FillMassTable(elts_avg, nterm_avg_table);
  FillMassTable(elts_avg, cterm_avg_table);
  FillMassTable(elts_mono, nprotterm_mono_table);
  FillMassTable(elts_mono, cprotterm_mono_table);
  FillMassTable(elts_avg, nprotterm_avg_table);
  FillMassTable(elts_avg, cprotterm_avg_table);

  if (mod_table) {
    // TODO: consider handling average and monoisotopic masses
    // differently (e.g. with different mod tables). We're deferring
    // this for now, since we haven't sorted out how to handle average
    // vs. monoisotopic masses generally.
    for (int i = 0; i < mod_table->static_mod_size(); ++i) {
      char aa = mod_table->static_mod(i).amino_acids()[0];
      double delta = mod_table->static_mod(i).delta();
      mono_table[aa] += delta;
      avg_table[aa] += delta;
      nterm_mono_table[aa] += delta;
      cterm_mono_table[aa] += delta;
      nterm_avg_table[aa] += delta;
      cterm_avg_table[aa] += delta;
      nprotterm_mono_table[aa] += delta;
      cprotterm_mono_table[aa] += delta;
      nprotterm_avg_table[aa] += delta;
      cprotterm_avg_table[aa] += delta;
    }
    carp(CARP_DEBUG, "Number of unique modification masses: %d", mod_table->unique_deltas_size());

    mod_coder_.Init(mod_table->unique_deltas_size());
    unique_deltas_.clear();
    unique_deltas_.reserve(mod_table->unique_deltas_size());
    for (int i = 0; i < mod_table->unique_deltas_size(); ++i) {
      unique_deltas_.push_back(mod_table->unique_deltas(i));
    }
  }

  //apply peptide terminal static mods on top of regular static mods
  ApplyTerminusStaticMods(n_mod_table, nterm_mono_table, nterm_avg_table);
  ApplyTerminusStaticMods(c_mod_table, cterm_mono_table, cterm_avg_table);
  //Apply peptide terminal mods to the protein terminal tables as well.
  ApplyTerminusStaticMods(n_mod_table, nprotterm_mono_table, nprotterm_avg_table);
  ApplyTerminusStaticMods(c_mod_table, cprotterm_mono_table, cprotterm_avg_table);

  //apply protein terminal static mods on top of everything else.
  if(nprot_mod_table != nullptr)
    ApplyTerminusStaticMods(nprot_mod_table, nprotterm_mono_table, nprotterm_avg_table);
  if(cprot_mod_table != nullptr)
    ApplyTerminusStaticMods(cprot_mod_table, cprotterm_mono_table, cprotterm_avg_table);
  
  SetFixPt(mono_table, avg_table, fixp_mono_table, fixp_avg_table);
  SetFixPt(nterm_mono_table, nterm_avg_table, fixp_nterm_mono_table, fixp_nterm_avg_table);
  SetFixPt(cterm_mono_table, cterm_avg_table, fixp_cterm_mono_table, fixp_cterm_avg_table);
  SetFixPt(nprotterm_mono_table, nprotterm_avg_table, fixp_nprotterm_mono_table, fixp_nprotterm_avg_table);
  SetFixPt(cprotterm_mono_table, cprotterm_avg_table, fixp_cprotterm_mono_table, fixp_cprotterm_avg_table);

  bin_width_ = bin_width;
  bin_offset_ = bin_offset;
  BIN_H2O = mass2bin(mono_h2o, 1);
  BIN_NH3 = mass2bin(mono_nh3, 1);
  BIN_CO = mass2bin(mono_co, 1);

  return true;
}

/* Updates the masses for the respective table to include the
static mod information from the mod table. */
void MassConstants::ApplyTerminusStaticMods(const pb::ModTable* mod_table, 
  double* mono_table, double* avg_table) {

  for (int i = 0; i < mod_table->static_mod_size(); ++i) {
    char aa = (mod_table->static_mod(i).amino_acids())[0];
    double delta = mod_table->static_mod(i).delta();
    if (aa == 'X') {
      const char* AA = "ACDEFGHIJKLMNOPQRSTUVWY";
      for (int currAA = 0; AA[currAA] != '\0'; currAA++) {
        mono_table[AA[currAA]] += delta;
        avg_table[AA[currAA]] += delta;
      }
    } else {
      mono_table[aa] += delta;
      avg_table[aa] += delta;
    }
  }
}

/* Sets the fixpt for the given mono and avg tables */
void MassConstants::SetFixPt(double* mono_table, double* avg_table, 
            FixPt* fixp_mono_table, FixPt* fixp_avg_table) {
  for (int i = 0; i < 256; ++i) {
    if (mono_table[i] == 0) {
      mono_table[i] = avg_table[i]/* = aa_bin_1[i] = aa_bin_2[i]*/
      = numeric_limits<double>::signaling_NaN();
      fixp_mono_table[i] = fixp_avg_table[i] = 0;
    } else {
      fixp_mono_table[i] = ToFixPt(mono_table[i]);
      fixp_avg_table[i] = ToFixPt(avg_table[i]);
    }
  }
}

static bool CheckModification(const pb::Modification& mod, bool* repeats = NULL) {
  string aa_str = mod.amino_acids();
  if (aa_str.length() != 1)
    return false;
  char aa = aa_str[0];
  const char* AA = "ACDEFGHIJKLMNOPQRSTUVWY";
  for (; (*AA != '\0') && (*AA != aa); ++AA); // find aa in AA
  if (*AA == '\0')
    return false;
  if (repeats != NULL) {
    if (repeats[aa])
      return false;
    repeats[aa] = true;
  }
  return true;
}

static bool CheckModTable(const pb::ModTable& mod_table) {
  // static mods table should not have repeated amino acids
  bool repeats[256];
  for (int i = 0; i < 256; ++i)
    repeats[i] = false;
  bool usingXMod = false;
  for (int i = 0; i < mod_table.static_mod_size(); ++i) {
    const pb::Modification& mod = mod_table.static_mod(i);
    string aa_str = mod.amino_acids();
    if (aa_str.length() != 1) {   //redundant, already checked by the parser
      return false;
    } else if (aa_str[0] == 'X') {
      usingXMod = true;
    } else if (!CheckModification(mod, repeats)) {
      carp(CARP_FATAL, "Multiple static mods on same residue detected.");
      return false;
    }
  }
  if (usingXMod && mod_table.static_mod_size() > 1) {
    carp(CARP_FATAL, "We are using X static mod, but we have detected "
    "other static mods as well. Only one static mod is allowed.");
  }

  return true;
}
