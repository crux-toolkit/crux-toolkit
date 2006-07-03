#ifndef OBJECTS_H
#define OBJECTS_H

/**
 * \typedef SPECTRUM_T 
 * A spectrum
 */
typedef struct spectrum SPECTRUM_T;

/**
 * \typedef PEAK_ITERATOR_T 
 * An object to iterate over the peaks in a spectrum
 */
typedef struct peak_iterator PEAK_ITERATOR_T;

/**
 * \typedef SPECTRUM_COLLECTION_T 
 * A collection of spectrum objects
 */
typedef struct spectrum_collection SPECTRUM_COLLECTION_T;

/**
 * \typedef SPECTRUM_ITERATOR_T 
 * An object to iterate over the spectrum objects in a spectrum_collection
 */
typedef struct spectrum_iterator SPECTRUM_ITERATOR_T;

/**
 * \typedef PROTEIN_T
 * A protein sequence
 */
typedef struct protein PROTEIN_T;

/**
 * \typedef PEPTIDE_ITERATOR_T
 * An objects to iterate over the peptides in a protein sequence
 */
typedef struct peptide_iterator PEPTIDE_ITERATOR_T;

/**
 * \typedef PEPTIDE_T
 * A peptide subsequence of a protein
 */
typedef struct peptide PEPTIDE_T;

/**
 * \typedef RESIDUE_ITERATOR_T 
 * An object to iterate over the residues in a peptide
 */
typedef struct residue_iterator RESIDUE_ITERATOR_T;

/**
 * \typedef PEAK_T 
 * A peak in a spectrum
 */
typedef struct peak PEAK_T;


#endif
