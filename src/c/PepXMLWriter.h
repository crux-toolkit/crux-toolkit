/**
 * \file PepXMLWriter.h
 * \brief Writes search results in the TPP-compliant .pep.xml format.
 */
#ifndef PEPXMLWRITER_H
#define PEPXMLWRITER_H

#include <string>
#include <vector>
#include "objects.h"

class PepXMLWriter{

 public:
  PepXMLWriter();
  ~PepXMLWriter();

  bool exact_pval_search_;

  /**
   * Open a file of the given name.  Replace an existing file if
   * overwrite is true, else exit if an existing file is found.
   */
  void openFile(const char* filename, bool overwrite);

  /**
   * Close the file, if open.
   */
  void closeFile();

  /**
   * Write all of the elements that preceed the spectrum_query
   * elements.
   * Requires OpenFile has been called without CloseFile.
   */
  void writeHeader();

  /**
   * Close the msms_run_summaryand msms_pipeline_analysis tags.
   * Requires OpenFile has been called without CloseFile.
   */
  void writeFooter();

  /**
   * Define which scores will be printed.
   * Assumes that all psms will have all the same scores.  Requires that
   * scores_computed is an array with NUMBER_SCORER_TYPES elements
   * indexed by SCORER_TYPE_T.  So if scores_computed[i] == true, then
   * (SCORER_TYPE_T)i will be printed.
   * Requires OpenFile has been called without CloseFile.
   */
  void SetScoresComputed(const bool* scores_computed);

  /**
   * Write the details for a PSM to be contained in a spectrum_query
   * element.  Requires that the arrays pre_aas, post_aas,
   * protein_names, protein_descriptions are all num_proteins long.
   * Assumes scores is NUMBER_SCORER_TYPES long. 
   * Requires OpenFile has been called without CloseFile.
   */
  void writePSM(
    int spectrum_scan_number, ///< identifier for the spectrum
    const char* filename, ///< file that spectrum came from
    double spectrum_neutral_mass, ///< computed mass of the spectrum
                                  ///at this charge state
    int charge, ///< assumed charge state for the match
    
    int* PSM_rank, ///< rank of this peptide for the spectrum
    const char* unmodified_peptide_sequence, ///< sequence with no mods
    const char* modified_peptide_sequence, ///< either with symbols or masses
    double peptide_mass, ///< mass of the peptide sequence
    int num_proteins, ///< proteins matched to this peptide
    const char* flanking,  ///< "XY, AB, " X and Y are the preceeding and
                        /// following aas in the first protein 
    std::vector<std::string>& protein_names, ///<
    std::vector<std::string>& protein_descriptions, ///<
    bool* scores_computed,
    double* scores, ///< indexed by score type
    unsigned current_num_matches
  );

  /**
   * Legacy function that sets delta_cn, by_ions_matched, and by_ions_total
   */
  void writePSM(
    int spectrum_scan_number, ///< identifier for the spectrum
    const char* filename, ///< file that spectrum came from
    double spectrum_neutral_mass, ///< computed mass of the spectrum
                                  ///at this charge state
    int charge, ///< assumed charge state for the match
    
    int* PSM_rank, ///< rank of this peptide for the spectrum
    const char* unmodified_peptide_sequence, ///< sequence with no mods
    const char* modified_peptide_sequence, ///< either with symbols or masses
    double peptide_mass, ///< mass of the peptide sequence
    int num_proteins, ///< proteins matched to this peptide
    const char* flanking,  ///< "XY, AB, " X and Y are the preceeding and
                        /// following aas in the first protein 
    std::vector<std::string>& protein_names, ///<
    std::vector<std::string>& protein_descriptions, ///<
    double delta_cn, ///<
    bool* scores_computed,
    double* scores, ///< indexed by score type
    unsigned  num_matched_ions, 
    unsigned tot_matched_ions,
    unsigned current_num_matches
  );
 protected:
  void initScoreNames();
  void printSpectrumElement(int spectrum_scan_number, 
                            const char* filename,
                            double spectrum_neutral_mass, 
                            int charge);
  void closeSpectrumElement();
  std::string getSpectrumTitle(int spectrum_scan_number, 
                               const char* filename,
                               int charge);
  void printPeptideElement(int* rank,
                           const char* peptide_sequence,
                           const char* mod_peptide_sequence,
                           double peptide_mass,
                           double spectrum_mass,
                           int num_proteins,
                           const char* flanking_aas,
                           std::vector<std::string>& protein_names, 
                           std::vector<std::string>& protein_descriptions,
                           bool* scores_computed,
                           double* scores,
                           unsigned current_num_matches);

  void printScores(
    double* scores, 
    bool* scores_computed,
    int* ranks
  );

  void printAnalysis(double* scores, bool* scores_computed);

  std::string filename_;
  FILE* file_;
  bool* scores_computed_;
  std::string last_spectrum_printed_;
  int current_index_;
  int mass_precision_;
  ENZYME_T enzyme_;
  int precision_;
  std::vector<const char*> scoreNames_;
};


#endif // PEPXMLWRITER_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

