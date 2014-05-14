/**
 * \file ProteinMatchCollection.h
 * $Revision: 1.00 $
 * \brief Object for holding a collection of ProteinMatches, PeptideMatches, and SpectrumMatches.
 ************************************************************************************************/
#ifndef PROTEINMATCHCOLLECTION_H_
#define PROTEINMATCHCOLLECTION_H_
#include "PeptideMatch.h"
#include "ProteinMatch.h"
#include "objects.h"
#include "match_objects.h"

#include <deque>
#include <set>
#include <string>

class ProteinMatchCollection {

 protected:

  struct cmpSeq {
    bool operator() (const MODIFIED_AA_T* lhs, const MODIFIED_AA_T* rhs) {
      size_t idx = -1;
      MODIFIED_AA_T lhs_aa, rhs_aa;
      do {
        ++idx;
        lhs_aa = lhs[idx];
        rhs_aa = rhs[idx];
        if (lhs_aa != rhs_aa) {
          return lhs_aa < rhs_aa;
        }
      } while (lhs_aa != MOD_SEQ_NULL && rhs_aa != MOD_SEQ_NULL);
      return rhs_aa != MOD_SEQ_NULL;
    }
  };

  std::map<std::string, ProteinMatch*> protein_match_map_;
  std::map<MODIFIED_AA_T*, PeptideMatch*, cmpSeq> peptide_match_map_;
  std::deque<ProteinMatch*> protein_matches_; ///< All protein matches
  std::deque<PeptideMatch*> peptide_matches_; ///< All peptide matches
  std::deque<SpectrumMatch*> spectrum_matches_; ///< All spectrum matches

  std::map<std::pair<int, int>, int> spectrum_counts_; ///< matches/spectrum
  bool distinct_matches_; ///< are matches distinct?
 public:

  /**
   * \returns a blank ProteinMatchCollection
   */
  ProteinMatchCollection();

  /**
   * \returns a ProteinMatchCollection using a MatchCollection
   * TODO - remove this later
   */
  ProteinMatchCollection(
    MatchCollection* match_collection ///< matches to add
  );
  
  /**
   * Default destructor
   */
  virtual~ ProteinMatchCollection();

  /**
   * \returns the begin iterator for the ProteinMatch objects
   */
  ProteinMatchIterator proteinMatchBegin();

  /**
   * \returns the end iterator for the ProteinMatch objects
   */
  ProteinMatchIterator proteinMatchEnd();

  /**
   * \returns the begin iterator for the PeptideMatch objects
   */
  PeptideMatchIterator peptideMatchBegin();
  
  /**
   * \returns the end iterator for all of the PeptideMatch objects
   */
  PeptideMatchIterator peptideMatchEnd();

  /**
   * \returns the begin iterator for all of the SpectrumMatch objects
   */
  SpectrumMatchIterator spectrumMatchBegin();

  /**
   * \returns the end iterator for all of the SpectrumMatch objects
   */
  SpectrumMatchIterator spectrumMatchEnd();

  /**
   * \returns the ProteinMatch for a Protein, creates a new
   * one if it is not found and create is true
   */  
  ProteinMatch* getProteinMatch(
    Crux::Protein* protein,  ///< Protein for which to find the protein match
    bool create = true ///< Create the ProteinMatch if it doesn't exist
  );

  /**
   * \returns the ProteinMatch for a Protein, null if it doesn't exist
   */
  ProteinMatch* getProteinMatch(
    const std::string& id ///< id of the protein
  );

  /**
   * \returns a PeptideMatch for a Peptide object, creating it
   * if it doesn't exist and create is true
   */
  PeptideMatch* getPeptideMatch(
    Crux::Peptide* peptide,  ///< peptide to find  
    bool create=true ///< create if it doesn't exist?
  );

  /**
   * \returns the PeptideMatch for the sequence, null if it doesn't exist
   */
  PeptideMatch* getPeptideMatch(
    MODIFIED_AA_T* mod_seq ///< modified sequence to find
    );

  /**
   * Helper method that adds a Crux match to the ProteinCollection,
   * creating the SpectrumMatch, PeptideMatch, and ProteinMatch objects.
   */
  void addMatch(
    MatchCollection* match_collection, ///<Collection from where the match came from
    Crux::Match* match ///< Match to add
  );

  /**
   * Helper method that adds a Crux match to the ProteinMatchCollection,
   * Adding all SpectrumMatch, PeptideMatch, and ProteinMatch objects.
   */
  void addMatches(MatchCollection* match_collection);

  /**
   * Get the matches/spectrum as a map, where the key is <scan, charge>
   */
  const std::map<std::pair<int, int>, int>& getMatchesSpectrum();
  
  /**
   * \returns whether matches are distinct are not
   */
  bool hasDistinctMatches();

};

#endif // PROTEINMATCHCOLLECTION_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

