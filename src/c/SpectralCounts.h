/**
 * \file SpectralCounts.h 
 * AUTHOR: Barbara Frewen
 * CREATE DATE: 8 March 2011
 * \brief Object for running the spectral-counts command.
 */
#ifndef SPECRAL_COUNTS_H
#define SPECRAL_COUNTS_H

#include "CruxApplication.h"
#include <string>
#include <vector>
#include <map>
#include <set>
#include "utils.h"
#include "objects.h"
#include "SpectrumCollection.h"
#include "OutputFiles.h"

class SpectralCounts: public CruxApplication { 

 public:

  SpectralCounts();
  ~SpectralCounts();
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();

 private:
  // internally-used types
  /**
   * \typedef PeptideSet
   * \brief Collection of peptide objects (not a meta-peptide)
   */
  typedef std::set<PEPTIDE_T*, bool(*)(PEPTIDE_T*, PEPTIDE_T*)> PeptideSet;
  /**
   * \typedef MetaMapping
   * \brief Mapping of peptideSet to MetaProtein
   * Each entry is a set of peptides mapped to a set of proteins of which
   * all contain the set of peptides
   */
  typedef std::map<PeptideSet, MetaProtein, 
              bool(*)(PeptideSet, PeptideSet) > MetaMapping;
  /**
   * \typedef ProteinToPeptides
   * Mapping of Protein objects to a set of peptides that are part
   * of the protein sequence
   */
  typedef std::map<Protein*, PeptideSet , 
              bool(*)(Protein*, Protein*)> ProteinToPeptides;
  /**
   * \typedef MetaToScore
   * \brief Mapping of MetaProtein to the score assigned to it
   */
  typedef std::map<MetaProtein, FLOAT_T, 
              bool(*)(MetaProtein, MetaProtein)> MetaToScore;
  /**
   * \typedef ProteinToMeta
   * \brief Mapping of Protein to MetaProtein to which it belongs
   */
  typedef std::map<Protein*, MetaProtein, 
              bool(*)(Protein*, Protein*)> ProteinToMetaProtein;
  
  // private functions
  void get_parameter_values();
  void filter_matches();
  void get_peptide_scores();
  void get_protein_scores();
  void get_protein_to_peptides();
  void get_protein_to_meta_protein();
  void get_meta_mapping();
  void get_meta_ranks();
  void get_meta_scores();
  void perform_parsimony_analysis();
  void normalize_peptide_scores();
  void normalize_protein_scores();
  void make_unique_mapping();
  void getSpectra(std::map<std::pair<int,int>, Spectrum*>& spectra);
  FLOAT_T sum_match_intensity(MATCH_T* match,
                          SpectrumCollection* spectra);

  // member variables
  OutputFiles* output_;
  std::string psm_file_;
  FLOAT_T threshold_;
  std::string database_name_;
  bool unique_mapping_;
  QUANT_LEVEL_TYPE_T quantitation_;
  PARSIMONY_TYPE_T parsimony_;
  MEASURE_TYPE_T measure_;
  FLOAT_T bin_width_;
  std::set<MATCH_T*> matches_;
  PeptideToScore peptide_scores_;
  ProteinToScore protein_scores_;
  ProteinToPeptides protein_supporting_peptides_;
  ProteinToMetaProtein protein_meta_protein_;
  MetaMapping meta_mapping_;
  MetaToScore meta_protein_scores_;
  MetaToRank meta_protein_ranks_;

  // comparison function declarations
  static bool compare_peptide_sets(PeptideSet, PeptideSet);
  static bool compare_meta_proteins(MetaProtein, MetaProtein);
  static bool sets_are_equal_size(std::pair<PeptideSet, MetaProtein>,
                                  std::pair<PeptideSet, MetaProtein>);
 
}; // class


#endif
