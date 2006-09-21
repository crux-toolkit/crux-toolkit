/**
 * \file objects.h 
 * $Revision: 1.17 $
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
 * \typedef PEPTIDE_SRC_ITERATOR_T 
 * \brief An object to iterate over the protein peptide associations in a peptide
 */
typedef struct peptide_src_iterator PEPTIDE_SRC_ITERATOR_T;


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
 * The enum for isotopic mass type (average, mono)
 */
enum _mass_type {AVERAGE, MONO };

/**
 * \typedef MASS_TYPE_T
 * \brief The typedef for mass type (average, mono);
 */
typedef enum _mass_type MASS_TYPE_T;

/**
 * \typedef PEPTIDE_SRC_T
 * \brief object for mapping a peptide to it's parent protein.
 */
typedef struct peptide_src PEPTIDE_SRC_T;


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


/**
 * The enum for sort type (mass, length, lexical, none)
 */
enum _sort_type {MASS, LENGTH, LEXICAL, NONE};

/**
 * \typedef SORT_TYPE_T
 * \brief The typedef for sort type (mass, length)
 */
typedef enum _sort_type SORT_TYPE_T;

/**
 * \typedef DATABASE_SORTED_PEPTIDE_ITERATOR_T
 * \brief An object to iterate over the peptides in a database in sorted order 
 */
typedef struct database_sorted_peptide_iterator DATABASE_SORTED_PEPTIDE_ITERATOR_T;

/**
 * \typedef PEPTIDE_WRAPPER_T
 * \brief An object to wrap a peptide, allowing a linked list of peptides.
 */
typedef struct peptide_wrapper PEPTIDE_WRAPPER_T;

/**
 * \typedef INDEX_T
 * \brief An index of a database 
 */
typedef struct index INDEX_T;

/**
 * \typedef INDEX_PEPTIDE_ITERATOR_T
 * \brief An object to iterate over the peptides in an index
 */
typedef struct index_peptide_iterator INDEX_PEPTIDE_ITERATOR_T;


/**
 * \struct index_filtered_peptide_iterator
 * \brief An iterator to filter out the peptides wanted from the index_peptide_iterator
 */
typedef struct index_filtered_peptide_iterator INDEX_FILTERED_PEPTIDE_ITERATOR_T;

/**
 * \typedef SORTED_PEPTIDE_ITERATOR_T
 * \brief An object to iterate over the peptides in sorted order 
 */
typedef struct sorted_peptide_iterator SORTED_PEPTIDE_ITERATOR_T;


/**
 * The enum for index type
 */
enum _index_type {DB_INDEX, BIN_INDEX};

/**
 * \typedef INDEX_TYPE_T
 * \brief The typedef for index type (db_index, bin_index)
 */
typedef enum _index_type INDEX_TYPE_T;

/**
 * \struct bin_peptide_iterator
 * \brief An iterator to iterate over the peptides in a bin( one file handler)
 */
typedef struct bin_peptide_iterator BIN_PEPTIDE_ITERATOR_T;

/**
 * \struct bin_sorted_peptide_iterator
 * \brief Object to iterate over the peptides within a bin, in an
 * sort in mass
 */
typedef struct bin_sorted_peptide_iterator BIN_SORTED_PEPTIDE_ITERATOR_T;


/**
 * \struct protein_index
 * \brief Object to store the protein relation to the fasta file
 */
typedef struct protein_index PROTEIN_INDEX_T;

/**
 * \struct protein_index_iterator
 * \brief Object to iterate over the protein index in the protein index file
 */
typedef struct protein_index_iterator PROTEIN_INDEX_ITERATOR_T;

#endif
