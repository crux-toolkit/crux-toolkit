#ifndef OBJECTS_H
#define OBJECTS_H

/**
 * \typedef SPECTRUM_T 
 */
typedef struct spectrum SPECTRUM_T;

/**
 * \typedef PEAK_ITERATOR_T An object to iterate over the peaks in a
 * spectrum object.
 */
typedef struct peak_iterator PEAK_ITERATOR_T;

/**
 * \typedef SPECTRUM_COLLECTION_T 
 */
typedef struct spectrum_collection SPECTRUM_COLLECTION_T;

/**
 * \typedef SPECTRUM_ITERATOR_T An object to iterate over the spectra in a
 * spectrum_collection object.
 */
typedef struct spectrum_iterator SPECTRUM_ITERATOR_T;

/**
 * \typedef PROTEIN_T
 */
typedef struct protein PROTEIN_T;

/**
 * \typedef PEPTIDE_ITERATOR_T
 */
typedef struct peptide_iterator PEPTIDE_ITERATOR_T;

/**
 * \typedef PEPTIDE_T
 */
typedef struct peptide PEPTIDE_T;

/**
 * \typedef RESIDUE_ITERATOR_T An object to iterate over the residues
 * in a peptide.
 */
typedef struct residue_iterator RESIDUE_ITERATOR_T;

/**
 * \typedef PEAK_T 
 */
typedef struct peak PEAK_T;


#endif
