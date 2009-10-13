/**
 * \file ion_series.h 
 * AUTHOR: Chris Park
 * CREATE DATE: 28 June 2006
 * $Revision: 1.23 $
 * \brief Object for a series of ions.
 *****************************************************************************/
#ifndef ION_SERIES_H
#define ION_SERIES_H

#include <stdio.h>
#include "objects.h"
#include "peptide.h"
#include "ion.h"
#include "ion_series.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \returns An (empty) ion_series object.
 */
ION_SERIES_T* allocate_ion_series(void);

/**
 * copies in the peptide sequence
 * Use this method to create ion_series only when few are needed,
 * because the memory allocation process is expensive.
 * If need a repeated new ion-series for different peptides, 
 * use "new_ion_series_generic" & "update_ion_series" combination, thus only allocate one 
 * ion_seires object.
 *\returns Instantiates a new ion_series object from the given peptide sequence and charge
 */
ION_SERIES_T* new_ion_series(
  const char* peptide, ///< The peptide for this ion series. -in
  int charge, ///< The charge for this ion series -in
  ION_CONSTRAINT_T* constraint ///< The constraints which the ions in this series obey.
  );


/**
 * Creates a heap allocated generic ion_series object that must be updated by "update_ion_series" method
 * to transform the object into a ion-series for a specific instance of a peptide sequence and charge.
 *\returns Instantiates a new generic ion_series object that must be updated for each peptide instance
 */
ION_SERIES_T* new_ion_series_generic(
  ION_CONSTRAINT_T* constraint, ///< The constraints which the ions in this series obey.
  int charge ///< The charge for this ion series -in
  );


/**
 * Updates an ion_series to a specific instance of a peptide sequence.
 * If the ion_series has been already generated its ions, will free ions up.
 * Copies in the peptide sequence.
 * and re-initialize for the new peptide sequence.
 */
void update_ion_series(
  ION_SERIES_T* ion_series, ///< the working ion_series -in
  char* peptide, ///< The peptide sequence for this ion series. -in
  MODIFIED_AA_T* mod_seq ///< modified version of seq -in
  );

/**
 * Frees an allocated ion_series object.
 */
void free_ion_series(
  ION_SERIES_T* ion_series ///< the ion collection to free - in
);

/**
 * Prints a ion_series object to file.
 */
void print_ion_series(
  ION_SERIES_T* ion_series, ///< ion_series to print -in 
  FILE* file ///< file for output -out
  );

/**
 * Prints a ion_series object to file, in GMTK single-ion format.
 */
void print_ion_series_single_gmtk(
  ION_SERIES_T* ion_series, ///< ion_series to print -in 
  ION_CONSTRAINT_T* ion_constraint, ///< ion_constraint to obey -in 
  FILE* file, ///< file -out
  int sentence_idx
  );

/**
 * Prints a ion_series object to file, in GMTK paired-ion format.
 */
void print_ion_series_paired_gmtk(
  ION_SERIES_T* ion_series, ///< ion_series to print -in 
  ION_CONSTRAINT_T* first_ion_constraint, ///< ion_constraint to obey -in 
  ION_CONSTRAINT_T* second_ion_constraint, ///< ion_constraint to obey -in 
  FILE* file, ///< file output
  int sentence_idx
  );

/**
 * Predict ion series
 */
void predict_ions(
  ION_SERIES_T* ion_series ///< the ion series to predict ions for
);

/**
 * Assign peaks to the nearest ions, within a tolerance (set in param file)
 */
void ion_series_assign_nearest_peaks(
    ION_SERIES_T* ion_series, 
    SPECTRUM_T* spectrum);

/**
 * Copies ion_series object from src to dest.
 *  must pass in a memory allocated ION_SERIES_T dest
 */
void copy_ion_series(
  ION_SERIES_T* src,///< ion to copy from -in
  ION_SERIES_T* dest///< ion to copy to -out
  );

/** 
 * Access routines of the form get_<object>_<field> and set_<object>_<field>. 
 */

/*************************************
 * ION_SERIES_T: get and set methods
 ************************************/

/**
 * \returns the ion that meets the constraint or NULL
 */
ION_T* get_ion_series_ion(
  ION_SERIES_T* ion_series, ///< the ion_series -in                          
  ION_CONSTRAINT_T* ion_constraint,
  int cleavage_idx
  );

/**
 * User should not free the peptide sequence seperate from the ion_series
 *\returns a pointer to the original parent peptide sequence of the ion_series object
 */
char* get_ion_series_peptide(
  ION_SERIES_T* ion_series ///< the working ion_series -in                          
  );

/**
 *\returns the peptide length of which the ions are made
 */
int get_ion_series_peptide_length(
  ION_SERIES_T* ion_series ///< the working ion_series -in                          
  );

/**
 * copies in the peptide sequence to heap allocated sequence.
 * set the parent peptide sequence of the ion_series object
 */
void set_ion_series_peptide(
  ION_SERIES_T* ion_series, ///< the working ion_series -in
  char* peptide///< the peptide sequence to set -in
  );

/**
 *\returns the charge of the ion_series object
 */
int get_ion_series_charge(
  ION_SERIES_T* ion_series ///< the working ion_series -in                          
  );

/**
 * set the charge of the ion_series object
 */
void set_ion_series_charge(
  ION_SERIES_T* ion_series, ///< the working ion_series -in
  int charge///< the charge of the ion -in
  );

/**
 *\returns the constraint of the ion_series object
 */
ION_CONSTRAINT_T* get_ion_series_ion_constraint(
  ION_SERIES_T* ion_series ///< the working ion_series -in                          
  );

/**
 * set the of the ion_series object
 */
void set_ion_series_ion_constraint(
  ION_SERIES_T* ion_series, ///< the working ion_series -in
  ION_CONSTRAINT_T* constraint///<  -in
  );

/**
 *\returns the total number of ions in the ion_series object
 */
int get_ion_series_num_ions(
  ION_SERIES_T* ion_series ///< the working ion_series -in                          
  );

/**
 *\returns the total number of ion_type in the ion_series object
 */
int get_ion_series_num_ions_one_type(
  ION_SERIES_T* ion_series, ///< the working ion_series -in                          
  ION_TYPE_T ion_type ///< the type of ions -in
  );


/******************************/

/*************************
 * ION_CONSTRAINT methods
 *************************/

/**
 *\returns an empty heap allocated ion_constraint
 */
ION_CONSTRAINT_T* allocate_ion_constraint(void);

/**
 * modification, all modifications 0
 * add more modifications as needed using the set_ion_constraint_modification
 *\returns a new heap allocated ion_constraint
 */
ION_CONSTRAINT_T* new_ion_constraint(
  MASS_TYPE_T mass_type, ///< the mass_type to use MONO|AVERAGE
  int max_charge, ///< the maximum charge of the ions, cannot exceed the parent peptide's charge
  ION_TYPE_T ion_type, ///< the ion types the peptide series should include
  BOOLEAN_T precursor_ion  ///< should include precursor ion?
  );

/**
 * \brief Create a new ion constraint based on the score type and the
 * charge of the peptide to be modeled.  Uses other
 * new_ion_constraint_ methods for some types.
 *
 * \returns A newly allocated ion constraint.
 */
ION_CONSTRAINT_T* new_ion_constraint_smart(
  SCORER_TYPE_T score_type,
  int charge
);
/**
 * modification, sets all fields for sequest settings
 *\returns a new heap allocated ion_constraint
 */
ION_CONSTRAINT_T* new_ion_constraint_sequest(
  MASS_TYPE_T mass_type, ///< the mass_type to use MONO|AVERAGE
  int max_charge, ///< the maximum charge of the ions, cannot exceed the parent peptide's charge
  ION_TYPE_T ion_type, ///< the ion types the peptide series should include
  BOOLEAN_T precursor_ion ///< should include precursor ion?
  );

/**
 * modification, sets all fields for GMTK settings
 *\returns a new heap allocated ion_constraint
 */
ION_CONSTRAINT_T* new_ion_constraint_gmtk(
  int charge ///< the charge of the peptide for which to predict ions
  );


/**
 * modification, sets all fields for sequest Sp scoring settings
 *\returns a new heap allocated ion_constraint
 */
ION_CONSTRAINT_T* new_ion_constraint_sequest_sp(
  int max_charge ///< the maximum charge of the ions, cannot exceed the parent peptide's charge
  );

/**
 * modification, sets all fields for Sequest Xcorr scoring settings
 * make B, Y, A type ions
 *\returns a new heap allocated ion_constraint
 */
ION_CONSTRAINT_T* new_ion_constraint_sequest_xcorr(
  int charge ///< the maximum charge of the ions, cannot exceed the parent peptide's charge
  );

/**
 * Frees an allocated ion_constraint object.
 */
void free_ion_constraint(
  ION_CONSTRAINT_T* ion_constraint///< the ion constraints to enforce -in
  );

/**
 * copies the ion_constraint pointer
 */
ION_CONSTRAINT_T* copy_ion_constraint_ptr(
  ION_CONSTRAINT_T* ion_constraint
  );

/**
 * copies ion_constraint object from src to dest
 * must pass in a memory allocated ION_CONSTRAINT_T dest
 */
void copy_ion_constraint(
  ION_CONSTRAINT_T* src,///< ion_constraint to copy from -in
  ION_CONSTRAINT_T* dest///< ion_constraint to copy to -out
);

/** 
 * Determines if a ion satisfies a ion_constraint.
 * \returns TRUE if the constraint is satisified. FALSE if not.
 */
BOOLEAN_T ion_constraint_is_satisfied(
   ION_CONSTRAINT_T* ion_constraint,///< the ion constraints to enforce -in
   ION_T* ion ///< query ion -in
   );

/**
 * Sets the modification count
 * can only add isotopes
 */
void set_ion_constraint_modification(
  ION_CONSTRAINT_T* ion_constraint,///< the ion constraints to enforce -in
  ION_MODIFICATION_T mod_type, ///< ion modification type -in
  int count  ///< the count of the modification -in  
  );

/**
 * Sets the exact modification boolean to exactness criteria
 * and if exactness is true sets min_charge = max_charge.
 * In other words, the constraint is now exact, in that it refers to a
 * particular ion series, charge states, and modification state, as opposed
 * to e.g. b-ions of charge state +1 or +2, or with or without NH3 loss
 */
void set_ion_constraint_exactness(
  ION_CONSTRAINT_T* ion_constraint,///< the ion constraints to enforce -in
  BOOLEAN_T exactness ///< whether to be exact or not -in
  );

/**
 * gets the modification count for specific mod_type
 */
int get_ion_constraint_modification(
  ION_CONSTRAINT_T* ion_constraint,///< the ion constraints to enforce -in
  ION_MODIFICATION_T mod_type ///< ion modification type -in
  );

/**
 * gets the mass type of the ion_constraint
 */
MASS_TYPE_T get_ion_constraint_mass_type(
  ION_CONSTRAINT_T* ion_constraint///< the ion constraints to enforce -in
  );

/**************************
 *  ION_ITERATOR_T object
 **************************/

/**
 * Instantiates a new ion_iterator object from ion_series.
 * \returns a ION_ITERATOR_T object.
 */
ION_ITERATOR_T* new_ion_iterator(
  ION_SERIES_T* ion_series ///< ion_series to iterate -in
  );        

/**
 * does not free ion
 * Frees an allocated ion_iterator object.
 */
void free_ion_iterator(
  ION_ITERATOR_T* ion_iterator///< free ion_iterator -in
);

/**
 * The basic iterator function has_next.
 */
BOOLEAN_T ion_iterator_has_next(
  ION_ITERATOR_T* ion_iterator///< is there a next ion? -in
);

/**
 * The basic iterator function next.
 */
ION_T* ion_iterator_next(
  ION_ITERATOR_T* ion_iterator///< return the next ion -in
);

/**********************************
 * ION_FILTERED_ITERATOR_T object
 **********************************/

/**
 * Only copies in the constraint as pointer
 * Instantiates a new ion_filtered_iterator object from ion_series.
 * \returns a ION_FILTERED_ITERATOR_T object.
 */
ION_FILTERED_ITERATOR_T* new_ion_filtered_iterator(
  ION_SERIES_T* ion_series, ///< ion_series to iterate -in
  ION_CONSTRAINT_T* constraint  ///< ion_constraint which returned ions satisfy
  );        

/**
 * The constraint is NOT freed from the iterator.
 * Frees an allocated ion_filtered_iterator object.
 */
void free_ion_filtered_iterator(
  ION_FILTERED_ITERATOR_T* ion_iterator///< free ion_iterator -in
);

/**
 * The basic iterator function has_next.
 */
BOOLEAN_T ion_filtered_iterator_has_next(
  ION_FILTERED_ITERATOR_T* ion_iterator///< is there a next ion? -in
);

/**
 * The basic iterator function next.
 */
ION_T* ion_filtered_iterator_next(
  ION_FILTERED_ITERATOR_T* ion_iterator///< return the next ion -in
);

#ifdef __cplusplus
}
#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
