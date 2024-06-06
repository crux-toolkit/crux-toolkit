#ifndef SEQUEST_PARAMS_H
#define SEQUEST_PARAMS_H

#include "Parsers/Algorithm2XML/SearchParams/SearchParams.h"


/*
const static residueStrct nativeAA[]


 = {
  {1, "<", 1.0078250321},
  {1, ">", 17.0027396542},
  {1, "G", 57.0214637236},
  {1, "A", 71.0371137878},
  {1, "V", 99.0684139162},
  {1, "L", 113.0840639804},
  {1, "I", 113.0840639804},
  {1, "S", 87.0320284099},
  {1, "C", 103.0091844778},
  {1, "T", 101.0476784741},
  {1, "M", 131.0404846062},
  {1, "P", 97.0527638520},
  {1, "F", 147.0684139162},
  {1, "Y", 163.0633285383},
  {1, "W", 186.0793129535},
  {1, "H", 137.0589118624},
  {1, "K", 128.0949630177},
  {1, "R", 156.1011110281},
  {1, "D", 115.0269430320},
  {1, "E", 129.0425930962},
  {1, "N", 114.0429274472},
  {1, "Q", 128.0585775114}
};

double getMonoisotopicAAMass(char* aa) {
  switch(aa) {
  case '<' : return 1.0078250321;
  case 'n' : return 1.0078250321;
  case '>' : return 17.0027396542;
  case 'c' : return 17.0027396542;
  case 'G' : return 57.0214637236;
  case 'A' : return 71.0371137878;
  case 'V' : return 99.0684139162;
  case 'L' : return 113.0840639804;
  case 'I' : return 113.0840639804;
  case 'S' : return 87.0320284099;
  case 'C' : return 103.0091844778;
  case 'T' : return 101.0476784741;
  case 'M' : return 131.0404846062;
  case 'P' : return 97.0527638520;
  case 'F' : return 147.0684139162;
  case 'Y' : return 163.0633285383;
  case 'W' : return 186.0793129535;
  case 'H' : return 137.0589118624;
  case 'K' : return 128.0949630177;
  case 'R' : return 156.1011110281;
  case 'D' : return 115.0269430320;
  case 'E' : return 129.0425930962;
  case 'N' : return 114.0429274472;
  case 'Q' : return 128.0585775114;
  case default: return 0.0;
  } // switch
}

double getAverageAAMass(char* aa) {
  switch(aa) {
  case '<' : return 1.0078250321;
  case 'n' : return 1.0078250321;
  case '>' : return 17.0027396542;
  case 'c' : return 17.0027396542;
  case 'G' : return 57.0519;
  case 'A' : return 71.0788;
  case 'V' : return 99.1326;
  case 'L' : return 113.1594;
  case 'I' : return 113.1594;
  case 'S' : return 87.0782;
  case 'C' : return 103.1388;
  case 'T' : return 101.1051;
  case 'M' : return 131.1926;
  case 'P' : return 97.1167;
  case 'F' : return 147.1766;
  case 'Y' : return 163.1760;
  case 'W' : return 186.2132;
  case 'H' : return 137.1411;
  case 'K' : return 128.1741;
  case 'R' : return 156.1875;
  case 'D' : return 115.0886;
  case 'E' : return 129.1155;
  case 'N' : return 114.1038;
  case 'Q' : return 128.1307;
  case default: return 0.0;
  } // switch
}

*/

class SequestParams : public SearchParams {
  friend class SequestOut;
  friend class Out2XML;
 public:

  SequestParams(const char* paramsfile);
  Array<Tag*>* getParameterTags();
  Array<Tag*>* getModificationTags();
  void modifyAminoacidModifications(char* peptide);
  void modifyTerminalModifications(char* peptide);
  Array<Tag*>* getSearchParamTags(const char* basename, const char* engine);
  Array<Tag*>* getSearchParamTags(const char* basename, const char* engine, const char* database);
  void writeParams(FILE* paramFile);
 protected:


  void init();
  const char* getEnzyme(int index);
  Boolean matchStart(const char* line, const char* tag);
  void setModificationMasses();

  // put sequest specific parameters here....
  double pep_mass_tol_;
  double frag_ion_tol_; // default 0.0

  Boolean use_default_params_;
  Boolean print_dups_;

  int num_output_lines_; // def 10
  int max_num_differential_AA_per_mod_; // def = 4    ; max # of modified AA per diff. mod in a peptide
  int nucleotide_reading_frame_; //def = 0           ; 0=proteinDB, 1-6, 7=forward three, 8=reverse three, 9=all six
  int remove_precursor_peak_; // default = 0              ; 0=no, 1=yes
  double ion_cutoff_percentage_; //def = 0.0            ; prelim. score cutoff % as a decimal number i.e. 0.30 for 30%
  int match_peak_count_; //def = 0                   ; number of auto-detected peaks to try matching (max 5)
  int match_peak_allowed_error_; //def = 1           ; number of allowed errors in matching auto-detected peaks
  double match_peak_tolerance_; //def = 1.0             ; mass tolerance for matching auto-detected peaks

  char protein_mass_filter_[500];
  char sequence_header_filter_[500];
  char ion_series_[150];
  char database_[4086];
};





#endif
