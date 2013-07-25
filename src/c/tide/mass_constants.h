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

typedef unsigned int FixPt; // 32-bit fixed-point arithmetic
   // ALTERNATIVE unsigned long long int; // 64 bits
   // ALTERNATIVE unsigned long int; // machine register size

namespace pb { class ModTable; }

class MassConstants {
 public:
  // Elemental mass values are from www.unimod.org/masses.html.
  static const double elts_mono[];
  static const double elts_avg[];

  static const double proton = 1.007276467;
  static const double bin_width = 1.0005079;

  static bool Init(const pb::ModTable* mod_table);

  static double mono_table[];
  static double avg_table[];
  static double aa_bin_1[];
  static double aa_bin_2[];

  static const double mono_h2o;
  static const double avg_h2o;
  static const double mono_nh3;
  static const double mono_co;

  // Fixed-Point Versions
  static const double kFixedPointScalar = 1e5;
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
    *delta = unique_deltas_bin_[unique_delta_index];
  }

 private:
  static void FillMassTable(const double* elements, double* table);

  static ModCoder mod_coder_;
  static double* unique_deltas_;
  static double* unique_deltas_bin_;
};

#endif
