/**
 * \file PeptideMatch.h
 * $Revision: 1.00 $
 * \brief Object for holding peptide scores
 ******************************************************/
#ifndef PEPTIDEMATCH_H_
#define PEPTIDEMATCH_H_

#include "AbstractMatch.h"
#include "match_objects.h"
#include "Peptide.h"

class PeptideMatch : public AbstractMatch {

 protected:  
  std::deque<ProteinMatch*> protein_matches_; ///< protein matches associated with this peptide match
  std::map<ProteinMatch*, PeptideSrc*> protein_match_to_peptide_src_; ///< peptide srcs for the protein match
  std::deque<SpectrumMatch*> spectrum_matches_; ///< spectrum matches for this peptide
  Crux::Peptide* peptide_; ///< peptide itself

 public:
  /**
   * \returns an empty PeptideMatch
   */
  PeptideMatch();

  /**
   * \returns a PeptideMatch with the peptide set
   */
  PeptideMatch(
    Crux::Peptide* peptide ///< Peptide to set
  );

  /**
   * Default destructor
   */
  virtual ~PeptideMatch();

  /**
   * Sets the peptide for this match
   */
  void setPeptide(
    Crux::Peptide* peptide ///< peptide to set
  );

  /**
   * \returns the peptide for the match
   */
  Crux::Peptide* getPeptide();
   
  /**
   * adds a ProteinMatch to this PeptideMatch
   */
  void addProteinMatch(
    ProteinMatch* protein_match, ///< ProteinMatch to set
    PeptideSrc* src ///< Location in protein
  );

  /**
   * \returns the associated PeptideSrc from the ProteinMatch
   */
  PeptideSrc* getSrc(
    ProteinMatch* protein_match ///< gets the associated src 
  );

  /**
   * Adds a spectrumMatch to this PeptideMatch
   */
  void addSpectrumMatch(
    SpectrumMatch* spectrum_match ///<SpectrumMatch to add
  );
  
  /**
   * \returns the begin iterator of the spectrummatches
   */
  SpectrumMatchIterator spectrumMatchBegin();

  /**
   * \returns the end iterator of the spectrum matches
   */
  SpectrumMatchIterator spectrumMatchEnd();

  /**
   * \returns the begin iterator of the spectrum matches
   */
  ProteinMatchIterator proteinMatchBegin();

  /**
   * \returns the end iterator of the spectrum matches
   */
  ProteinMatchIterator proteinMatchEnd();
 
};

#endif //PEPTIDEMATCH_H_

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
