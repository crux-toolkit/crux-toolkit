// Benjamin Diament
//
// Wrapper for all mass constants. Arrays, where they appear, 
// are of length 256, and are indexed by an upper-case letter 
// symbolizing the amino acid.
//
// Unused array entries are filled with signaling_NaN, which is 
// supposed to generate an error if touched, but seems to fail silently 
// on Intel. Oh well. :(

#ifndef MASS_CONSTANTS_H
#define MASS_CONSTANTS_H

#include "mod_coder.h"

#define BIN_WIDTH 1.0005079
#define BIN_OFFSET 0.68
typedef unsigned int FixPt; // 32-bit fixed-point arithmetic
   // ALTERNATIVE unsigned long long int; // 64 bits
   // ALTERNATIVE unsigned long int; // machine register size

namespace pb { class ModTable; }

class MassConstants {
 public:
  // Elemental mass values are from www.unimod.org/masses.html.
  static const double elts_mono[];
  static const double elts_avg[];

#ifdef _MSC_VER
  // The windows compiler only allows intialization
  // of static constant integer types within a class
  static const double proton;
#else
  static const double proton = 1.00727646688;
#endif

  static bool Init(const pb::ModTable* mod_table, const double bin_width, const double bin_offset);
  
  static double mono_table[];
  static double avg_table[];
  static const double mono_h2o;
  static const double avg_h2o;
  static const double mono_nh3;
  static const double mono_co;
  static const double A;
//  static const double B_H2O;
//  static const double B_NH3;
  static const double B;
//  static const double Y_H2O;
//  static const double Y_NH3;
  static const double Y;
  static double BIN_H2O;
  static double BIN_NH3;
/*  static const double BIN_SHIFT_A_ION_CHG_1;
  static const double BIN_SHIFT_A_ION_CHG_2;
  static const double BIN_SHIFT_H2O_CHG_1;
  static const double BIN_SHIFT_H2O_CHG_2;
  static const double BIN_SHIFT_NH3_CHG_1;
  static const double BIN_SHIFT_NH3_CHG_2_CASE_A;
  static const double BIN_SHIFT_NH3_CHG_2_CASE_B;
*/  // Fixed-Point Versions
#ifdef _MSC_VER
  // The windows compiler only allows intialization
  // of static constant integer types within a class
  static const double kFixedPointScalar;
#else
  static const double kFixedPointScalar = 1e5;
#endif
  static FixPt ToFixPt(double x) {
    return FixPt(x * kFixedPointScalar + 0.5);
  }
  static double ToDouble(FixPt x) {
    return x/kFixedPointScalar;
  }

  static FixPt fixp_mono_table[];
  static FixPt fixp_avg_table[];

  static const FixPt fixp_mono_h2o;
  static const FixPt fixp_avg_h2o;

  static const FixPt fixp_mono_nh3;
  static const FixPt fixp_mono_co;

  static const FixPt fixp_proton;

  static void DecodeMod(int code, int* aa_index, double* delta) {
    int unique_delta_index;
    mod_coder_.DecodeMod(code, aa_index, &unique_delta_index);
    *delta = unique_deltas_[unique_delta_index];
  }
  static unsigned int mass2bin(double mass, double charge = 1){
    return (unsigned int)((mass + (charge - 1)*proton)/(charge*bin_width_) + 1.0 - bin_offset_);
  }

  static double bin_width_;
  static double bin_offset_;
 private:
  static void FillMassTable(const double* elements, double* table);

  static ModCoder mod_coder_;
  static double* unique_deltas_;

};

#endif
