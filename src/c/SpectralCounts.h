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
#include "MatchCollectionParser.h"
#include "objects.h"
#include "Peptide.h"
#include "SpectrumCollection.h"
#include "OutputFiles.h"

class SpectralCounts: public CruxApplication { 

 public:

  SpectralCounts();
  ~SpectralCounts();
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();
  virtual COMMAND_T getCommand();
  virtual bool needsOutputDirectory();

 private:
  // internally-used types
  /**
   * \typedef PeptideSet
   * \brief Collection of peptide objects (not a meta-peptide)
   */
  typedef std::set<Crux::Peptide*, bool(*)(Crux::Peptide*, Crux::Peptide*)> PeptideSet;
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
  typedef std::map<Crux::Protein*, PeptideSet , 
    bool(*)(Crux::Protein*, Crux::Protein*)> ProteinToPeptides;
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
  typedef std::map<Crux::Protein*, MetaProtein, 
    bool(*)(Crux::Protein*, Crux::Protein*)> ProteinToMetaProtein;
  
  // private functions
  void getParameterValues();
  void filterMatches();
  void getPeptideScores();
  
  /**
   * calculates the protein scores based upon the dNSAF metric.  First,
   * the spectral counts for the peptides unique to each protein is added,
   * then the spectral counts for the peptide shared amongst the proteins
   * are added with a distribution factor based upon the number of unique
   * counts that each protein is assigned
   */
  void getProteinScoresDNSAF();

  void getProteinScores();
  void getProteinToPeptides();
  void getProteinToMetaProtein();
  void getMetaMapping();
  void getMetaRanks();
  void getMetaScores();
  void performParsimonyAnalysis();
  void normalizePeptideScores();
  void normalizeProteinScores();
  /**
   * Checks that the normalized scores add up to one
   */
  void checkProteinNormalization();



  void computeEmpai();
  void makeUniqueMapping();
  void getSpectra(std::map<std::pair<int,int>, Crux::Spectrum*>& spectra);
  FLOAT_T sumMatchIntensity(Crux::Match* match,
                          Crux::SpectrumCollection* spectra);
  SCORER_TYPE_T get_qval_type(MatchCollection* match_collection);


  /**
   * adds all matches 
   */
  void filterMatchesNone();

  /**
   * filters matches based upon a custom threshold (not q-value)
   */
  void filterMatchesCustom();

  /**
   * filters matches based upon q-value
   */
  void filterMatchesQValue();

  /**
   * filters matches based upon a SCORER_TYPE_T
   */
  void filterMatchesScore(
    SCORER_TYPE_T scorer ///< scorer to use
  );

  /**
   * filters matches based upon a custom score that is not SCORER_TYPE_T
   */
  void filterMatchesCustomScore();

  void invalidCustomScore();

  // member variables
  OutputFiles* output_;
  std::string psm_file_;
  FLOAT_T threshold_;
  std::string database_name_;
  bool unique_mapping_;
  THRESHOLD_T threshold_type_;
  QUANT_LEVEL_TYPE_T quantitation_;
  PARSIMONY_TYPE_T parsimony_;
  MEASURE_TYPE_T measure_;
  FLOAT_T bin_width_;
  std::set<Crux::Match*> matches_;
  MatchCollection* match_collection_;
  // For custom thresholding fields
  bool threshold_min_; 
  std::string custom_threshold_name_;
  MatchCollectionParser parser_;

  PeptideToScore peptide_scores_;
  PeptideToScore peptide_scores_unique_;
  PeptideToScore peptide_scores_shared_;

  ProteinToScore protein_scores_;
  ProteinToScore protein_scores_unique_;
  ProteinToScore protein_scores_shared_;

  ProteinToPeptides protein_supporting_peptides_;
  ProteinToMetaProtein protein_meta_protein_;
  MetaMapping meta_mapping_;
  MetaToScore meta_protein_scores_;
  MetaToRank meta_protein_ranks_;

  // comparison function declarations
  static bool comparePeptideSets(PeptideSet, PeptideSet);
  static bool compareMetaProteins(MetaProtein, MetaProtein);
  static bool setsAreEqualSize(std::pair<PeptideSet, MetaProtein>,
                               std::pair<PeptideSet, MetaProtein>);
 
}; // class


#endif
