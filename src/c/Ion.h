/**
 * \file Ion.h
 * $Revision: 1.15 $
 * \brief Object for representing one ion in an ion_series.
 *
 */
#ifndef ION_H
#define ION_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "objects.h"
#include "utils.h"
#include "mass.h"
#include "objects.h"
#include "Peak.h"

static const int MAX_MODIFICATIONS = 4; ///< maximum modifications allowed per ion

class Ion {

 protected:
  ION_TYPE_T type_;  ///< type of the ion 
  int cleavage_idx_; ///< index of peptide amide that fragments to form this ion, starting from the N-term end 
  // N.b. this is different than the b1,y1 index, in that it always starts
  // from the N-term
  int charge_; ///< the ion charge
  char* peptide_sequence_; ///< the peptide sequence that fragments to form this ion
  FLOAT_T peptide_mass_; ///< the mass of the peptide. For efficiency
  int modification_counts_[MAX_MODIFICATIONS]; ///< an array of the number of different ion modifications
  FLOAT_T ion_mass_z_;   ///< The mass/z of the ion. 
  Peak * peak_;  ///< The assigned peak. NULL if no peak // TODO add ptr count
  int pointer_count_; ///< count of number of references to this ion.

  static const int MZ_INT_MAX = 10;
  static const int MZ_INT_MIN = 0;
  static const int PAIRED_ION_INTS = 12;
  static const int PAIRED_ION_FLOATS = 6;
  static const int SINGLE_ION_FLOATS = 3;
  static const int SINGLE_ION_INTS = 9;

  /**
   * helper function
   * only copies the pointer to the peptide sequence
   * creates an heap allocated ion, ion mass is not calculated.
   * \returns an ION_T object
   */
  void initBasicIon(
    ION_TYPE_T type,   ///< intensity for the new ion -in 
    int cleavage_idx, ///< index into the peptide amide bonds of this ion
    int charge, ///< charge of the ion
    char* peptide ///< location for the new ion -in
    ); 

  
  /**
   * initializes the mass array
   */
  static void initializeModificationMasses(
    MASS_TYPE_T mass_type ///< mass type (average, mono) -in
  );

 /**
  *\return the modified mass_z accodring to the modification type
  */
  FLOAT_T modifyMassZ(
    FLOAT_T mass_z, ///< the pre-modified ion mass_z -in
    int modification_count, ///< number of times more the modification occurs -in
    ION_MODIFICATION_T ion_modification, ///< the type of modification -in
    int charge, ///< charge of the ion
    MASS_TYPE_T mass_type ///< mass type (average, mono) -in
    );

  /**
   *\return the modified mass accodring to the modification type
   */
  FLOAT_T modifyMass(
    FLOAT_T mass, ///< the pre-modified ion mass -in
    int modification_count, ///< number of times more the modification occurs -in
    ION_MODIFICATION_T ion_modification, ///< the type of modification -in
    MASS_TYPE_T mass_type ///< mass type (average, mono) -in
    );

  /**
   *\returns the ion's AA mass added all up
   */
  FLOAT_T getMass(
    MASS_TYPE_T mass_type ///< mass type (average, mono) -in
    );    

 public:

  /**
   * initializes an Ion object.
   */
  void init();


  /**
   * \returns An (empty) ion object.
   */
  Ion();
  
 /**
  * The peptide sequence is copied as a pointer.
  * only copies the pointer to the peptide sequence
  * cleavage index starts from 0...n (ex 0:A:1:K:2:V:3:..n:L)
  * \returns an ION_T object
  */
  Ion(
    ION_TYPE_T type,   ///< intensity for the new ion -in 
    int cleavage_idx, ///< index into the peptide amide bonds of this ion
    int charge, ///< charge of the ion
    char* peptide, ///< location for the new ion -in
    MASS_TYPE_T mass_type ///< mass type (average, mono) -in
    ); 

 /**
  * only copies the pointer to the peptide sequence
  * inputs a array of all the modification counts
  * \returns an ION_T object
  */
  Ion(
    ION_TYPE_T type,   ///< intensity for the new ion -in 
    int cleavage_idx, ///< index into the peptide amide bonds of this ion
    int charge, ///< charge of the ion
    char* peptide, ///< location for the new ion -in
    MASS_TYPE_T mass_type, ///< mass type (average, mono) -in
    int* modification_counts ///< an array of modification counts for each modification -in
    );

  /**
   * only copies the pointer to the peptide sequence
   * inputs a array of all the modification counts
   * inputs the pre modified mass, of just all AA mass summed up.
   * \returns an ION_T object
   */
  Ion(
    ION_TYPE_T type,   ///< intensity for the new ion -in 
    int cleavage_idx, ///< index into the peptide amide bonds of this ion
    int charge, ///< charge of the ion
    char* peptide, ///< location for the new ion -in
    MASS_TYPE_T mass_type, ///< mass type (average, mono) -in
    FLOAT_T base_mass, ///< the base mass of the ion -in
    int* modification_counts ///< an array of modification counts for each modification -in
    );

  /**
   * only copies the pointer to the peptide sequence
   * inputs the pre modified mass, of just all AA mass summed up.
   * \returns an ION_T object
   */
  Ion(
    ION_TYPE_T type,   ///< intensity for the new ion -in 
    int cleavage_idx, ///< index into the peptide amide bonds of this ion
    int charge, ///< charge of the ion
    char* peptide, ///< location for the new ion -in
    MASS_TYPE_T mass_type, ///< mass type (average, mono) -in
    FLOAT_T base_mass ///< the base mass of the ion -in
    );

  /**
   * frees A ION_T object
   */
  ~Ion();

  /**
   * decrements the pointer and free the ion
   */
  static void freeIon(
    Ion* ion ///< ion to free
  );

  /**
   * Increments the pointer count
   */
  void incrementPointerCount();

  /**
   * prints the location and fields of ION_T object to the file, in the
   * following format:
   *
   * m/z \\t mass \\t charge \\t ion-series \\t  ...
   *  peptide-bond-index \\t modifications \n
   *
   * Where:
   *
   * m/z - is the ion's mass-to-charge
   * mass - is the ion's (charged) mass
   * charge - is the ion's charge e.g. 1,2,3
   * ion-series - is one of (b,y,p)
   * bond-index - is in [1...n), where n is peptide length
   * modifications - is one of (none|nh3|h2o)
   *
   * if the ion has more than one modification, each will be printed on a
   * separate line, with the necessary number of tabs to right justify
   */
  void print(
    FILE* file ///< to this file -in
  );

  /**
   * prints the location and fields of ION_T object to the file, in the
   * following format for GMTK single-ion models:
   *
   * m/z \\t mass \\t charge \\t ion-series \\t  ...
   *  peptide-bond-index \\t modifications \n
   *
   * Where:
   *
   * m/z - is the ion's mass-to-charge
   * mass - is the ion's (charged) mass
   * charge - is the ion's charge e.g. 1,2,3
   * ion-series - is one of (b,y,p)
   * bond-index - is in [1...n), where n is peptide length
   * modifications - is one of (none|nh3|h2o)
   *
   * if the ion has more than one modification, each will be printed on a
   * separate line, with the necessary number of tabs to right justify
   */
  void printGmtkSingle(
    FILE* file ///< to this file -in
  );

  /**
   * As print_ion_gmtk_single above, but in binary.
   */
  void printGmtkSingleBinary(
    FILE* file, ///< to this file -in
    int sentence_idx,
    int frame_idx
  );

  /**
   * Some hack routines to print out a null ion if there are none in the series.
   * For using neutral losses with GMTK.
   * Come in both binary and ascii versions.
   */
  void printNullGmtkSingleBinary(
    FILE* file,
    int sentence_idx,
    int frame_idx
    );

  void printNullGmtkPairedBinary(
    FILE* file,
    int sentence_idx,
    int frame_idx
    );

  /**
   * A hack routine to print out a null ion if there are none in the series.
   * For using neutral losses with GMTK.
   * Come in both binary and ascii versions.
   */
  void printNullGmtkSingle(
    FILE* file
    );

  /**
   * prints the location and fields of the two ION_T objects to the file
   * following format for GMTK paired-ion models:
   *
   * ints 
   *
   * 1.  m/z ratio int (from N-term)
   * 2.  m/z ratio int (from C-term)
   * 3.  peptide idx (from N-term)
   * 4.  peptide idx (from C-term)
   * 5.  aa Id (N-term)
   * 6.  aa Id (C-term)
   * 7.  first possible
   * 8.  first detectable
   * 9.  first detected
   * 10. second possible
   * 11. second observable
   * 12. second detected
   *
   * floats
   *
   * 1. m/z ratio FLOAT_T (from N-term)
   * 2. m/z ratio FLOAT_T (from C-term)
   * 3. first raw
   * 4. second raw
   * 5. first rank
   * 6. second rank
   *
   */
  static void printGmtkPairedBinary(
    Ion* first_ion, ///< print this ion -in
    Ion* second_ion, ///< print this ion -in
    FILE* file, ///< to this file -in
    int sentence_idx,
    int frame_idx
    );

  /**
   * Adds the given ION_MODIFICATION to this ion
   */
  void addModification(
    ION_MODIFICATION_T modification, ///< add this modification to the ion -in
    int modification_count,  ///< the number of modifications
    MASS_TYPE_T mass_type ///< mass type (average, mono) -in
    );

  /**
   * Adds the given ION_MODIFICATION to this ion
   */
  void setPeak(
    Peak * peak ///< peak to add to this ion -in
    );

  /**
   * Copies ion object from src to dest.
   * must pass in a seperate pointer peptide sequence from its own ion_series object
   * must pass in a memory allocated ION_T* dest
   */
  static void copy(
    Ion* src,///< ion to copy from -in
    Ion* dest,///< ion to copy to -out
    char* peptide_sequence ///< the peptide sequence that the dest should refer to -in
    );

  /**
   *\returns TRUE if forward ion_type(A,B,C), else reverse ion_type(X,Y,Z) FALSE
   */
  bool isForwardType();

  /**
   *\returns TRUE if the ion has modifications, else FALSE
   */
  bool isModified();

  /*********************************
   * get, set methods for ion fields
   *********************************/

  /**
   * \returns the location of ION_T object
   */
  FLOAT_T getMassZ();

  /**
   * sets the mass/z of the ION_T object
   * while this can be calculated from the char*, cleavage_idx and
   * modifications, it allows some optimizations if we allow it to be set
   * instead
   */
  void setMassZ(
    FLOAT_T mass_z ///< the m/z location -in
  );

  FLOAT_T getMassFromMassZ();

  void setMassZFromMass(
    FLOAT_T mass
  );


  /**
   * return the cleavage_idx of the ion object
   */
  int getCleavageIdx();

  /**
   * set the cleavage_idx of the ion object
   */
  void setCleavageIdx(
    int cleavage_idx ///< the cleavage index in the peptide -in
  );

  /**
   * return the charge of the ion object
   */
  int getCharge();

  /**
   * set the charge of the ion object
   */
  void setCharge(
    int charge ///< the charge of this ion -in
  );

  /**
   * return the ION_TYPE_T of the ion object
   */
  ION_TYPE_T getType();

  /**
   * set the ION_TYPE_T of the ion object
   */
  void setType(
    ION_TYPE_T ion_type ///< the ion type of this ion -in 
    );

  /**
   * return the parent peptide sequence of the ion object
   * returns a pointer to the sequence, should not free
   */
  char* getPeptideSequence();

  /**
   * return a pointer to the modification_count array of the ion object
   */
  int* getModificationCounts();

  /**
   * return the count of in the modification_count array of the ion object
   */
  int getSingleModificationCount(
    ION_MODIFICATION_T mod_type ///< the modification count wanted -in
    );

  /**
   * return the total number of modifications this ion has.
   */
  int getTotalModificationCount();

  /**
   * set the parent peptide_sequence of the ion object
   */
  void setPeptideSequence(
    char* peptide_sequence ///< the parent peptide's sequence of this ion -in 
    );

  /**
   * is_modified, indiciates if there are any modification to the ion
   * speeds up the proccess if FLASE.
   *\returns TRUE if successfully computes the mass/z of the ion, else FALSE
   */
  bool calcMassZ(
    MASS_TYPE_T mass_type, ///< mass type (average, mono) -in
    bool is_modified ///< are there any modifications for this ion? -in
    );

  /**
   * is_modified, indiciates if there are any modification to the ion
   * speeds up the proccess if FLASE.
   *\returns TRUE if successfully computes the mass/z of the ion, else FALSE
   */
  bool calcMassZWithMass(
    MASS_TYPE_T mass_type, ///< mass type (average, mono) -in
    FLOAT_T mass, ///< the basic mass of the ion -in
    bool is_modified ///< are there any modifications for this ion? -in
    );

};

/**
 *
 */
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
