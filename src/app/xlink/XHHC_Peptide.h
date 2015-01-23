/**
 * \file XHHC_Peptide.h 
 * AUTHOR: Sean McIlwain and Paul Draghicescu
 * CREATE DATE: 19 January 2012
 * \brief Object for keeping track of the links on a
 * cross-linkable peptide.
 *****************************************************************************/
#ifndef XHHC_PEPTIDE_H
#define XHHC_PEPTIDE_H

#include "objects.h"


#include <string>
#include <vector>

class XHHC_Peptide {
 protected:
  bool mass_calculated_[NUMBER_MASS_TYPES]; ///< MONO or Average
  int num_links_; ///< the number of links this peptide has
  std::vector<bool> links_; ///< vector bools for the aa position of the link
  std::string sequence_; ///< sequence of the peptide
  int length_; ///< sequence length
  FLOAT_T mass_[NUMBER_MASS_TYPES]; ///< mass of the peptide

 public:
  /**
    * \returns a blank XHHC_Peptide object
    */
  XHHC_Peptide();

  /**
   * \returns XHHC_Peptide object with the sequence initialized
   */
  XHHC_Peptide(
    std::string sequence ///< the peptide sequence
  );


  /**
   * Destructor
   */
  virtual ~XHHC_Peptide();

  /**
   * cleaves the peptide at an index, adding b and y ions
   * skips if a cleave site on self loop between the link
   */
  void splitAt(
    int index, ///< index to cleave 
    std::vector<std::pair<LinkedPeptide, LinkedPeptide> >& pairs, ///< the b-y ion pairs
    int charge, ///< charge of the peptide
    XHHC_Peptide& other, ///< the other peptide that we are linked to
    bool is_loop ///< Is this a self loop?
    );

  /**
   * \returns whether a link exists at the index
   */
  bool hasLinkAt(
    int index ///<amino acid index
  );

  /**
   * \returns the sequence length of the peptide
   */
  int getLength();

  /**
   * returns whether a link exists
   */
  bool hasLink();

  /**
   * sets the peptide sequence
   */
  void setSequence(
    std::string sequence
  );

  /**
   * \returns the peptide sequence
   */
  std::string getSequence();

  /**
   * \returns whether we have initialized the peptide
   */
  bool isEmpty();
          
  /**
   * \returns the number of links
   */
  int getNumLinks();

  /**
   * add a link at a amino acid index
   */
  void addLink(
    int index ///< the amino acid index
  );

  /**
   * removes a link from the peptide
   */
  void removeLink(
    int index ///< the amino acid index
  );
  
  /**
   * \returns the first index of a link, -1 if none exist
   */
  int linkSite();

  /**
   * \returns the mass of the peptide
   */
  FLOAT_T getMass(
    MASS_TYPE_T mass_type ///< MONO or AVERAGE
  );

  /**
   * \returns a shuffled peptide, preserving any links
   */
  XHHC_Peptide shuffle();
};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
