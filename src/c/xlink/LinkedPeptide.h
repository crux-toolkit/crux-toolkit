/**
 * \file LinkedPeptide.h 
 * AUTHOR: Sean McIlwain and Paul Draghicescu
 * CREATE DATE: 18 January 2012
 * \brief Object for keeping track of a (non-)crosslinked peptide.
 *****************************************************************************/
#ifndef LinkedPeptide_H
#define LinkedPeptide_H

#include "objects.h"
#include "xhhc.h"


class LinkedPeptide {
 protected:
  bool mass_calculated_[NUMBER_MASS_TYPES]; ///<MONO or AVERAGE.
  int charge_; ///< charge of the LinkedPeptide
  bool decoy_; ///< indicates that this LinkedPeptide is a decoy
  ION_TYPE_T type_; ///<B_ION or Y_ION
  FLOAT_T mass_[NUMBER_MASS_TYPES]; ///<mass of the LinkedPeptide
  FLOAT_T mz_[NUMBER_MASS_TYPES]; ///<mz of the LinkedPeptide
  std::vector<XHHC_Peptide> peptides_; ///<Peptides involved in the LinkedPeptide
  static FLOAT_T linker_mass_; ///<Mass of the linker

 public:
  
  /**
   * Initializes the object
   */
  void init();

  /**
   * \returns a blank LinkedPeptide object
   */
  LinkedPeptide(); 

  /**
   * \returns a linked peptide object.
   * constructor for a linked peptide. If A or B null, then
   * a self link will be created. If an index is -1, a link to nothing
   * will be created.
   */
  LinkedPeptide(
    char* peptide_A, ///< First peptide
    char* peptide_B, ///< Second peptide (can be NULL)
    int posA, ///< Index of link on first peptide
    int posB, ///< Index of link on second peptide
    int charge ///< charge of product
  );

  /**
   * \returns a blank LinkedPeptide object with 
   * the charge set
   */
  LinkedPeptide(
    int charge ///< charge of product
  );

  /**
   * Destructor
   */
  virtual ~LinkedPeptide();

  /**
   * adds a XHHC_Peptide to the list of peptide
   * in this LinkedPeptide object
   */
  void addPeptide(
    XHHC_Peptide& peptide ///< the XHHC_Peptide
  );
  
  /**
   * \returns a reference to the internal
   * peptides vector
   */
  std::vector<XHHC_Peptide>& getPeptides();

  /**
   * \returns the charge of the LinkedPeptide
   */
  int getCharge();

  /**
   * sets the charge of the LinkedPeptide
   */
  void setCharge(
    int charge ///< the charge of the LinkedPeptide
    );

  /**
   * sets the IonType for the LinkedPeptide
   */
  void setIonType(
    ION_TYPE_T type ///<The ion type
    );

  /**
   * \returns the ion type for this LinkedPeptide
   */
  ION_TYPE_T getIonType();

  /**
   * \returns the number of peptides in this linked
   * Peptide (should be 1 or 2)
   */
  int size();
  
  /**
   * indicates that this LinkedPeptide is a
   * decoy
   */
  void setDecoy();

  /**
   * \returns whether the LinkedPeptide is a
   * decoy
   */
  bool isDecoy();

  /**
   * \returns whether the LinkedPeptide is a
   * crosslinked peptide
   */
  bool isCrossLinked();

  /**
   * \returns whether the Peptide is a deadend
   */
  bool isDeadEnd();

  /**
   * \returns whether the Peptide is a self-loop
   */
  bool isSelfLoop();

  /**
   * \returns whether the Peptide is a normal/linear peptide
   */
  bool isLinear();

  /**
   * \returns the m/z of the LinkedPeptide
   */
  FLOAT_T getMZ(
    MASS_TYPE_T mass_type ///< MONO or AVERAGE
    );

  /**
   * calculates the mass of the LinkedPeptide
   */
  void calculateMass(
    MASS_TYPE_T mass_type ///< MONO or AVERAGE
    );
  
  /**
   * /returns the mass of the LinkedPeptide
   */
  FLOAT_T getMass(
    MASS_TYPE_T mass_type //< Average or Monoisotopic
  );
  
  /**
   * splits the first LinkedPeptide into b-y ions
   */
  void splitA(
    std::vector<std::pair<LinkedPeptide, LinkedPeptide> >& ion_pairs ///<the b-y ion pairs
  );

  /**
   * splits the second LinkedPeptide into b-y ions
   */
  void splitB(
    std::vector<std::pair<LinkedPeptide, LinkedPeptide> >& ion_pairs ///< the b-y ion pairs
  );

  /**
   * splits the LinkedPeptide into b-y ions
   */
  void split(
    std::vector<std::pair<LinkedPeptide, LinkedPeptide> >& ion_pairs ///< the b-y ion pairs
  );

  /**
   * prints the LinkedPeptide to a stream
   */
  friend std::ostream &operator<< (
    std::ostream& os, 
    LinkedPeptide& lp); 

  /**
   * used for sorting LinkedPeptides
   */
  friend bool operator< (
    const LinkedPeptide &lp1, 
    const LinkedPeptide &lp2
  );
 
  friend bool compareMassMono(
    const LinkedPeptide& lp1, 
    const LinkedPeptide& lp2
  );

  friend bool compareMassAverage(
    const LinkedPeptide& lp1, 
    const LinkedPeptide& lp2
  );

  friend bool compareMZMono(
    const LinkedPeptide& lp1,
    const LinkedPeptide& lp2
  );

  friend bool compareMZAverage(
    const LinkedPeptide& lp1,
    const LinkedPeptide& lp2
  );

  /**
   * Sorts a vector of LinkedPeptide by mass
   */
  static void sortByMass(
    std::vector<LinkedPeptide>& linked_peptides, ///< The LinkedPeptides to sort
    MASS_TYPE_T mass_type=MONO ///< MONO or AVERAGE
  );

  /**
   * Sorts a vector of LinkedPeptide by m/z
   */
  static void sortByMZ(
    std::vector<LinkedPeptide>& linked_peptides, ///< The LinkedPeptides to sort
    MASS_TYPE_T mass_type=MONO ///< MONO or AVERAGE
  );

  /**
   * Sets the linker_mass_ static variable
   */
  static void setLinkerMass(
    FLOAT_T linker_mass ///<the linker mass
    );

  /**
   * \returns the linker_mass_ static variable
   */
  static FLOAT_T getLinkerMass();

};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
