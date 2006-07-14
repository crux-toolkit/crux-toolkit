/**
 * \file objects.h 
 * $Revision: 1.8 $
 * \brief The defined objects
 *****************************************************************************/
#ifndef OBJECTS_H 
#define OBJECTS_H

#include <stdio.h>

/**
 * \typedef PEAK_T 
 * A peak in a spectrum
 */
typedef struct peak PEAK_T;

/**
 * \typedef SPECTRUM_T 
 * \brief A spectrum
 */
typedef struct spectrum SPECTRUM_T;

/**
 * The enum for spectrum type (MS1, MS2, MS3)
 */
enum _spectrum_type { MS1, MS2, MS3 };

/**
 * \typedef SPECTRUM_TYPE_T 
 * \brief The typedef for spectrum type (MS1, MS2, MS3)
 */
typedef enum _spectrum_type SPECTRUM_TYPE_T;

/**
 * \typedef PEAK_ITERATOR_T 
 * \brief An object to iterate over the peaks in a spectrum
 */
typedef struct peak_iterator PEAK_ITERATOR_T;

/**
 * \typedef SPECTRUM_COLLECTION_T 
 * \brief A collection of spectra
 */
typedef struct spectrum_collection SPECTRUM_COLLECTION_T;

/**
 * \typedef SPECTRUM_ITERATOR_T 
 * \brief An object to iterate over the spectra in a spectrum_collection
 */
typedef struct spectrum_iterator SPECTRUM_ITERATOR_T;

/**
 * \typedef PEPTIDE_T
 * \brief A peptide subsequence of a protein
 */
typedef struct peptide PEPTIDE_T;

/**
 * \typedef PEPTIDE_CONSTRAINT_T
 * \brief An object representing constraints which a peptide may or may not
 * satisfy.
 */
typedef struct peptide_constraint PEPTIDE_CONSTRAINT_T;

/**
 * \typedef RESIDUE_ITERATOR_T 
 * \brief An object to iterate over the residues in a peptide
 */
typedef struct residue_iterator RESIDUE_ITERATOR_T;

/**
 * \typedef PROTEIN_PEPTIDE_ASSOCIATION_ITERATOR_T 
 * \brief An object to iterate over the protein peptide associations in a peptide
 */
typedef struct protein_peptide_association_iterator PROTEIN_PEPTIDE_ASSOCIATION_ITERATOR_T;


/**
 * \brief The enum for peptide type, with regard to trypticity.
 */
enum _peptide_type { TRYPTIC, PARTIALLY_TRYPTIC, NOT_TRYPTIC, ANY_TRYPTIC}; 

/**
 * \typedef PEPTIDE_TYPE_T 
 * \brief The typedef for peptide type, with regard to trypticity.
 */
typedef enum _peptide_type PEPTIDE_TYPE_T;

/**
 * \typedef PROTEIN_PEPTIDE_ASSOCIATION_T
 * \brief object for mapping a peptide to it's parent protein.
 */
typedef struct protein_peptide_association PROTEIN_PEPTIDE_ASSOCIATION_T;


/**
 * \typedef PROTEIN_T
 * \brief A protein sequence
 */
typedef struct protein PROTEIN_T;

/**
 * \typedef PROTEIN_PEPTIDE_ITERATOR_T
 * \brief An object to iterate over the peptides in a protein sequence
 */
typedef struct protein_peptide_iterator PROTEIN_PEPTIDE_ITERATOR_T;

/**
 * \typedef DATABASE_T
 * \brief A database of protein sequences.
 */
typedef struct database DATABASE_T;

/**
 * \typedef DATABASE_PROTEIN_ITERATOR_T
 * \brief An object to iterate over the proteins in a database 
 */
typedef struct database_protein_iterator DATABASE_PROTEIN_ITERATOR_T;

/**
 * \typedef DATABASE_PEPTIDE_ITERATOR_T
 * \brief An object to iterate over the peptides in a database 
 */
typedef struct database_peptide_iterator DATABASE_PEPTIDE_ITERATOR_T;

#endif
