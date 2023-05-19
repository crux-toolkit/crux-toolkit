#include <cstdio>
#include "app/tide/abspath.h"
#include "app/tide/records_to_vector-inl.h"

#include "io/carp.h"
#include "parameter.h"
#include "io/SpectrumRecordWriter.h"
#include "TideIndexApplication.h"
#include "TideSearchApplication.h"
#include "ParamMedicApplication.h"
#include "PSMConvertApplication.h"
#include "tide/mass_constants.h"
#include "TideMatchSet.h"
#include "util/Params.h"
#include "util/FileUtils.h"
#include "util/StringUtils.h"
#include <math.h> //Added by Andy Lin
#include <map> //Added by Andy Lin


TideSearchApplication::TideSearchApplication():
  exact_pval_search_(false), remove_index_(""), spectrum_flag_(NULL) {
}

TideSearchApplication::~TideSearchApplication() {}

int TideSearchApplication::main(int argc, char** argv) {
  return main(Params::GetStrings("tide spectra file"));
}

SpectrumCollection* TideSearchApplication::loadSpectra(const string& file) {

}

void TideSearchApplication::search(void* threadarg) {

}
vector<string> TideSearchApplication::getOptions() const {
  string arr[] = {
    "auto-mz-bin-width",
    "auto-precursor-window",
    "compute-sp",
    "concat",
    "deisotope",
    "elution-window-size",
    "exact-p-value",
    "file-column",
    "fileroot",
    "isotope-error",
    "mass-precision",
    "max-precursor-charge",
    "min-precursor-charge",
    "min-peaks",
    "mod-precision",
    "mz-bin-offset",
    "mz-bin-width",
    "mzid-output",
    "num-threads",
    "output-dir",
    "overwrite",
    "parameter-file",
    "peptide-centric-search",
    "score-function",
    "fragment-tolerance",
    "evidence-granularity",
    "pepxml-output",
    "pin-output",
    "pm-charges",
    "pm-max-frag-mz",
    "pm-max-precursor-delta-ppm",
    "pm-max-precursor-mz",
    "pm-max-scan-separation",
    "pm-min-common-frag-peaks",
    "pm-min-frag-mz",
    "pm-min-peak-pairs",
    "pm-min-precursor-mz",
    "pm-min-scan-frag-peaks",
    "pm-pair-top-n-frag-peaks",
    "pm-top-n-frag-peaks",
    "precision",
    "precursor-window",
    "precursor-window-type",
    "print-search-progress",
    "remove-precursor-peak",
    "remove-precursor-tolerance",
    "scan-number",
    "skip-preprocessing",
    "spectrum-max-mz",
    "spectrum-min-mz",
    "spectrum-parser",
    "sqt-output",
    "store-index",
    "store-spectra",
    "top-match",
    "txt-output",
    "brief-output",
    "use-flanking-peaks",
    "use-neutral-loss-peaks",
    "use-z-line",
    "use-tailor-calibration",
    "verbosity"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}
