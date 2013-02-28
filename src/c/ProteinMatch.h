/**
 * \file ProteinMatch.h
 * $Revision: 1.00 $
 * \brief Object for holding protein scores
 ******************************************************/
#ifndef PROTEINMATCH_H_
#define PROTEINMATCH_H_

#include "match_objects.h"
#include "AbstractMatch.h"

#include <map>
#include <string>
#include <vector>

class ProteinMatch : public AbstractMatch {

 protected:
  std::vector<PeptideMatch*> peptide_matches_; ///< PeptideMatches associated with this ProteinMatch
  Crux::Protein* protein_; ///< Protein 
  
 public:

  /**
   * \returns an empty ProteinMatch object
   */
  ProteinMatch();
  
  /**
   * \returns a ProteinMatch object with the protein set
   */
  ProteinMatch(
    Crux::Protein* protein ///< Protein to set
  );
  
  /**
   * Destructor
   */
  virtual ~ProteinMatch();

  /** 
   * \returns the peptide for the match
   */
  Crux::Protein* getProtein();

  /**
   * \returns the begining of the peptidematches vector
   */
  PeptideMatchIterator peptideMatchBegin();

  /**
   * \returns the end of the peptide matches vector
   */
  PeptideMatchIterator peptideMatchEnd();

  /**
   * adds a peptide match to the protein match
   */
  void addPeptideMatch(
    PeptideMatch* peptide_match ///< peptide_match to add
  );

  /**
   * \returns the id for this protein
   */
  std::string getId();

};

#endif //PROTEINMATCH_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
