#ifndef _CNPRPROTEINSUMMARYHEADER_H
#define _CNPRPROTEINSUMMARYHEADER_H

#include "NeoProtXMLStructs.h"
#include "CnprProgramDetails.h"
#include <string>
#include <vector>

class CnprProteinSummaryHeader {
public:

  //Constructor
  CnprProteinSummaryHeader();

  void write(FILE* f, int tabs = -1);

  //required
  double initial_min_peptide_pro;
  double min_peptide_probability;
  double min_peptide_weight;
  int num_input_1_spectra;
  int num_input_2_spectra;
  int num_input_3_spectra;
  int num_input_4_spectra;
  int num_input_5_spectra;
  double num_predicted_correct_prots;
  std::string reference_database;
  std::string residue_substitution_list;
  std::string sample_enzyme;
  std::string source_files;
  std::string source_files_alt;
  
  //optional
  std::string organism;
  std::string source_file_xtn;
  double total_no_spectrum_ids;
  std::string win_cyg_reference_database;

  //sub elements
  CnprProgramDetails program_details;


private:

};

#endif
