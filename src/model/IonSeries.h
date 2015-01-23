/**
 * \file IonSeries.h 
 * AUTHOR: Chris Park
 * CREATE DATE: 28 June 2006
 * \brief Object for a series of ions.
 *****************************************************************************/
#ifndef IONSERIES_H
#define IONSERIES_H

#include <stdio.h>
#include <vector>
#include "objects.h"
#include "Peptide.h"
#include "Ion.h"
#include "IonConstraint.h"

static const int MAX_NUM_ION_TYPE = 8; // number of different ion_types

/**
 * \class IonSeries
 * \brief An object to represent a series of ions, and organize them.
 * For which additional data structures will be created as needed 
 * loss_limit can be equal to NULL, thus if need to use should always
 * check that it is not NULL.
 */
class IonSeries {
  friend class XLinkPeptide;
 protected:
  // TODO change name to unmodified_char_seq
  char* peptide_; ///< The peptide sequence for this ion series
  MODIFIED_AA_T* modified_aa_seq_; ///< sequence of the peptide
  FLOAT_T peptide_mass_; ///< The peptide neutral mass. For efficiency. 
  int charge_; ///< /<The charge state of the peptide for this ion series
  IonConstraint* constraint_; ///< The constraints which these ions obey
  std::vector<Ion*> ions_; ///< The ions in this series
  bool is_predicted_; ///< has this ion_series been predicted already?
  std::vector<Ion*> specific_ions_[MAX_NUM_ION_TYPE]; 
    ///< specific ions in the series, reference to master array of ions
  LOSS_LIMIT_T* loss_limit_; 
    ///< nh3, h2o loss limit for a given cleavage index, 
    ///< before using this array should always sheck if not NULL
  int peptide_length_;   ///< the length of the peptide

  // ??? what is the difference between peptide_length and num_ions

  /**
   * \brief Creates an array in which element i is the sum of the masses
   * of amino acids 0 to (i-1).  At i=0 is stored the length of the
   * peptide.  
   * \returns an array of ion masses for all sub sequences
   */
  FLOAT_T* createIonMassMatrix(
    //char* peptide, ///< The peptide for this ion series. -in
    MODIFIED_AA_T* modified_seq, ///< the sequence
    MASS_TYPE_T mass_type, ///< the mass_type to use MONO|AVERAGE
    int peptide_length ///< the length of the peptide
    );

  /**
   * user must ensure that there is enough space for this ion
   * adds ion to ion_series' master ion_array and if B|Y ion to the specific ion_array
   */
  void addIon(
    Ion* ion ///< ion to add -in
  );

  /**
   * helper function: add_ions
   * add all the ions to ion_series up to the max charge
   *\returns TRUE if successfully adds all ions, else FALSE
   */
  bool addIonsByCharge(
    FLOAT_T mass, ///< the base mass of the ion to add
    int cleavage_idx, ///< the absolute cleavage index (A,B,C from left X,Y,Z from right)
    ION_TYPE_T ion_type ///< the ion type of the ions to be added
    );

  /**
   * The modification depends on the loss/add && if the ion contains RKQN or STED
   * The number of losses possible cannot exceed the number of RKQN or STED in the ion
   * The loss_limit array in the ion_series must be populated prior to this method call
   *\returns TRUE if the ion can lose the mod_type modification, else FALSE
   */
  bool canIonGenerateModification(
    Ion* ion, ///< the ion to check if can lose nh3 -in
    ION_MODIFICATION_T mod_type, ///< generate ions of this modification_type -in/out
    int increment  ///< the add/loss of the modification
    );

  /**
   * Initializes an (empty) ion_series object, only called by constructors.
   */
  void init();

 public:



  /**
   * \returns An (empty) ion_series object.
   */
  IonSeries();

  /**
   * copies in the peptide sequence
   * Use this method to create ion_series only when few are needed,
   * because the memory allocation process is expensive.
   * If need a repeated new ion-series for different peptides, 
   * use "new_ion_series_generic" & "update_ion_series" combination, thus only allocate one 
   * ion_seires object.
   *\returns Instantiates a new ion_series object from the given peptide sequence and charge
   */
  IonSeries(
    const char* peptide, ///< The peptide for this ion series. -in
    int charge, ///< The charge for this ion series -in
    IonConstraint* constraint ///< The constraints which the ions in this series obey.
    );

  /**
   * Creates a heap allocated generic ion_series object that must be updated by "update_ion_series" method
   * to transform the object into a ion-series for a specific instance of a peptide sequence and charge.
   *\returns Instantiates a new generic ion_series object that must be updated for each peptide instance
   */
  IonSeries(
    IonConstraint* constraint, ///< The constraints which the ions in this series obey.
    int charge ///< The charge for this ion series -in
    );

  
/**
 * Updates an ion_series to a specific instance of a peptide sequence.
 * If the ion_series has been already generated its ions, will free ions up.
 * Copies in the peptide sequence.
 * and re-initialize for the new peptide sequence.
 */
 void update(
   char* peptide, ///< The peptide sequence for this ion series. -in
  MODIFIED_AA_T* mod_seq ///< modified version of seq -in
  );


  /**
   * Frees an allocated ion_series object.
   */
  virtual ~IonSeries();

  /**
   *Iterator access
   */
  IonIterator begin();
  IonIterator end();


  /**
   * Prints a ion_series object to file.
   */
  void print(
    FILE* file ///< file for output -out
    );

  /**
   * Prints a ion_series object to file, in GMTK single-ion format.
   */
  void printSingleGmtk(
    IonConstraint* ion_constraint, ///< ion_constraint to obey -in 
    FILE* file, ///< file -out
    int sentence_idx
    );

  /**
   * Prints a ion_series object to file, in GMTK paired-ion format.
   */
  void printPairedGmtk(
    IonConstraint* first_ion_constraint, ///< ion_constraint to obey -in 
    IonConstraint* second_ion_constraint, ///< ion_constraint to obey -in 
    FILE* file, ///< file output
    int sentence_idx
    );

  /**
   * \brief Find instances of amino acid which can incur neutral
   * losses: H2O (S|T|E|D), NH3(R|K|Q|N).  
   * Set the count of those observed so far for each cleavage index.
   * If no instance of amino acid, the count is assigned to 0
   * The information is used to determine if how many nh3 or h2o neutral
   * losses are possible. 
   */
  void scanForAAForNeutralLoss();

  /**
   * Creates all the ions with no modifications up to the max charge
   * Adds each ion to ion_series
   *\returns TRUE if successfully generates all the ions, else FALSE
   */
  bool generateIonsNoModification(
    FLOAT_T* mass_matrix ///< the mass matrix that stores the mass
    );

  /**
   * creates all the ions with specific modifications up to the max charge
   * copies all the existing ions that can be modified,
   * then applies the different modifications then adds the new modified ions to ion_series
   *\returns TRUE if successfully generates all the ions with modifications, else FALSE
   */
  bool generateIons(
    ION_MODIFICATION_T mod_type ///< generate ions of this modification_type -in/out
    );

  /**
   * creates all the flanking ions up to the max charge
   * can only create flanking ions that are B|Y ions and don't have modification
   * assumes the ions with no modification all are at the begining of the ion[] in ion_series
   * copies all the existing ions that can be modified,
   * then applies the different modifications then adds the new modified ions to ion_series
   *\returns TRUE if successfully generates all the ions with modifications, else FALSE
   */
  bool generateIonsFlank();

  /**
   * Predict ion series
   */
  void predictIons();

  /**
   * Copies ion_series object from src to dest.
   *  must pass in a memory allocated ION_SERIES_T dest
   */
  static void copy(
    IonSeries* src,///< ion to copy from -in
    IonSeries* dest///< ion to copy to -out
    );

  /**
   * remove an ion from IonSeries, does not free ion.
   */
  void removeIon(
    Ion* ion ///<ion to remove
  );


  /*************************************
   * ION_SERIES_T: get and set methods
   ************************************/

  /**
   * \returns the ion that meets the constraint or NULL
   */
  Ion* getIon(
    IonConstraint* ion_constraint,
    int cleavage_idx
    );

  /**
   * User should not free the peptide sequence seperate from the ion_series
   *\returns a pointer to the original parent peptide sequence of the ion_series object
   */
  char* getPeptide();

  /**
   *\returns the peptide length of which the ions are made
   */
  int getPeptideLength();

  /**
   * copies in the peptide sequence to heap allocated sequence.
   * set the parent peptide sequence of the ion_series object
   */
  void setPeptide(
    char* peptide///< the peptide sequence to set -in
  );

  /**
   *\returns the charge of the ion_series object
   */
  int getCharge();

  /**
   * set the charge of the ion_series object
   */
  void setCharge(
    int charge///< the charge of the ion -in
    );

  /**
   * get the is_predicted field of the ion_series object
   */
  bool getIsPredicted();

  /**
   * set the is_predicted field of the ion_series object
   */
  void setIsPredicted(
    bool is_predicted///< the is_predicted field -in
    );

  /**
   * get the modified_aa_seq of the ion_series object
   */
  MODIFIED_AA_T* getModifiedAASeq();

  /**
   *\returns the constraint of the ion_series object
   */
  IonConstraint* getIonConstraint();

  /**
   * set the of the ion_series object
   */
  void setIonConstraint(
    IonConstraint* constraint///<  -in
  );

  /**
   *\returns the total number of ions in the ion_series object
   */
  int getNumIons();

  /**
   *\returns the total number of ion_type in the ion_series object
   */
  int getNumIonsOneType(
    ION_TYPE_T ion_type ///< the type of ions -in
  );

  /**
   *\returns the ions corresponding to a specific ion type.
   */
  std::vector<Ion*>& getSpecificIons(
    ION_TYPE_T ion_type /// < the type of ions -in
  );
};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
