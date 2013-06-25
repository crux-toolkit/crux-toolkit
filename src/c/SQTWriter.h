#ifndef SQTWRITER_H
#define SQTWRITER_H

#include "carp.h"
#include "crux-file-utils.h"
#include "Index.h"
#include "parameter.h"
#include "PostProcessProtein.h"
#include "Spectrum.h"

#include <iomanip>
#include <string>

using namespace std;

class SQTWriter {

 public:
  SQTWriter();
  ~SQTWriter();
  void openFile(string filename);
  void closeFile();

  void writeHeader(
    string database,
    int num_proteins,
    bool is_decoy = false
  );

  void writeSpectrum(
    Crux::Spectrum* spectrum,
    SpectrumZState& z_state,
    int num_matches
  );

  void writePSM(
    Crux::Peptide* peptide,
    FLOAT_T xcorr_score,
    int xcorr_rank,
    FLOAT_T sp_score,
    int sp_rank,
    FLOAT_T delta_cn,
    int b_y_matched,
    int b_y_total,
    bool is_decoy
  );

 protected:
  ofstream* file_;

};

#endif // SQTWRITER_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

