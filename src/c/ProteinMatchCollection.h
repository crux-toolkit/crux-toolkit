/**
 * \file ProteinMatchCollection.h
 * $Revision: 1.00 $
 * \brief Object for holding a collection of ProteinMatches, PeptideMatches, and SpectrumMatches.
 ************************************************************************************************/
#ifndef PROTEINMATCHCOLLECTION_H_
#define PROTEINMATCHCOLLECTION_H_
#include "objects.h"
#include "match_objects.h"

#include <string>



class ProteinMatchCollection {

 protected:
  std::vector<ProteinMatch*> protein_matches_; ///< All protein matches
  std::vector<PeptideMatch*> peptide_matches_; ///< All peptide matches
  std::vector<SpectrumMatch*> spectrum_matches_; ///< All spectrum matches

 public:

  /**
   * \returns a blank ProteinMatchCollection
   */
  ProteinMatchCollection();

  /**
   * \returns a ProteinMatchCollection using a MatchCollection
   * TOOD - remove this later
   */
  ProteinMatchCollection(MatchCollection* match_collection);
  
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
    const std::string& sequence ///< sequence to find
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

};

#endif
