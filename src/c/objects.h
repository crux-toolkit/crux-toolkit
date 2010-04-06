/**
 * \file objects.h 
 * $Revision: 1.62 $
 * \brief The defined objects
 *****************************************************************************/
#ifndef OBJECTS_H 
#define OBJECTS_H

#include <stdio.h>

#define QSORT_COMPARE_METHOD int(*)(const void*, const void*)


#ifdef __cplusplus
class DelimitedFile;
#endif

/**
 * \typedef PEAK_T 
 * A peak in a spectrum
 */
typedef struct peak PEAK_T;

/**
 * The enum for peak sort type(_PEAK_LOCATION, _PEAK_INTENSITY)
 */
enum _peak_sort_type {_PEAK_LOCATION, _PEAK_INTENSITY};

/**
 * \typedef  PEAK_SORT_TYPE_T 
 * \brief The typedef for peak sort type(_PEAK_LOCATION, _PEAK_INTENSITY)
 */
typedef enum _peak_sort_type PEAK_SORT_TYPE_T;

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
 * \typedef FILTERED_SPECTRUM_CHARGE_ITERATOR_T 
 * \brief An object to iterate over the spectra in a spectrum_collection
 */
typedef struct filtered_spectrum_charge_iterator FILTERED_SPECTRUM_CHARGE_ITERATOR_T;

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

// REPLACE PEPTIDE_TYPE_T with DIGEST_T and ENZYME_T
/**
 * \enum _digest_type
 * The rule governing how a peptide was cleaved from its source
 * protein sequence. 
 */
enum _digest_type {
  INVALID_DIGEST,      ///< required invalid value for the enum
  FULL_DIGEST,         ///< c- AND n-term specific to ENZYME_T
  PARTIAL_DIGEST,      ///< c- OR n-term specific to ENZYME_T
  NON_SPECIFIC_DIGEST, ///< not specific to any enzyme cleavage rules
};
#define NUMBER_DIGEST_TYPES 4
/**
 * \typedef DIGEST_T
 * \brief The rule governing how a peptide was digested.  Used in
 * conjunction with ENZYME_T to define how peptides are generated.
 */
typedef enum _digest_type DIGEST_T;

/**
 * \enum _enzyme_type
 */
enum _enzyme_type {
  INVALID_ENZYME,        ///< required invalid value for the enum
  NO_ENZYME,             ///< cleave anywhere
  TRYPSIN,               ///< cleave after K or R, not before P
  CHYMOTRYPSIN,          ///< cleave after FWY, not before P
  ELASTASE,              ///< cleave after ALIV, not before P
  CLOSTRIPAIN,           ///< cleave after R
  CYANOGEN_BROMIDE,      ///< cleave after M
  IODOSOBENZOATE,        ///< cleave after W
  PROLINE_ENDOPEPTIDASE, ///< cleave after P
  STAPH_PROTEASE,        ///< cleave after E
  ASPN,                  ///< cleave before D
  MODIFIED_CHYMOTRYPSIN, ///< cleave after FWYL, not before P
  ELASTASE_TRYPSIN_CHYMOTRYPSIN, ///< cleave after ALIVKRWFY, not before P
  CUSTOM_ENZYME    ///< cleave after/before user-defined residues
};
#define NUMBER_ENZYME_TYPES 14
/**
 * \typedef ENZYME_T
 * \brief The enzyme with which a peptide was digested.  Used in
 * conjunction with DIGEST_T to define how peptides are generated.
 */
typedef enum _enzyme_type ENZYME_T;

/**
 * Being PARTIALLY_TRYPTIC is N or C terminus tryptic
 * \brief The enum for peptide type, with regard to trypticity.
 */
/*
#define NUMBER_PEPTIDE_TYPES 6
enum _peptide_type { TRYPTIC, PARTIALLY_TRYPTIC, N_TRYPTIC, C_TRYPTIC, NOT_TRYPTIC, ANY_TRYPTIC}; 
*/
/**
 * \typedef PEPTIDE_TYPE_T 
 * \brief The typedef for peptide type, with regard to trypticity.
 */
//typedef enum _peptide_type PEPTIDE_TYPE_T;

/**
 * The enum for isotopic mass type (average, mono)
 */
enum _mass_type {AVERAGE, MONO };
#define NUMBER_MASS_TYPES 2

/**
 * \typedef MASS_TYPE_T
 * \brief The typedef for mass type (average, mono);
 */
typedef enum _mass_type MASS_TYPE_T;


/**
 * The enum for window type for selecting peptides or assigning ions.
 */
enum _window_type {
  WINDOW_INVALID,
  WINDOW_MASS, 
  WINDOW_MZ, 
  WINDOW_PPM 
};
#define NUMBER_WINDOW_TYPES 4
/**
 * \typedef WINDOW_TYPE_T
 * \brief The typedef for window type (mass, mz, ppm);
 */
typedef enum _window_type WINDOW_TYPE_T;

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
enum _sort_type {SORT_NONE, SORT_MASS, SORT_LENGTH, SORT_LEXICAL};
#define NUMBER_SORT_TYPES 4

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
 * \typedef INDEX_FILTERED_PEPTIDE_ITERATOR_T
 * \brief An iterator to filter out the peptides wanted from the index_peptide_iterator
 */
typedef struct index_filtered_peptide_iterator INDEX_FILTERED_PEPTIDE_ITERATOR_T;

/**
 * \typedef SORTED_PEPTIDE_ITERATOR_T
 * \brief An object to iterate over the peptides in sorted order 
 */
typedef struct sorted_peptide_iterator SORTED_PEPTIDE_ITERATOR_T;

/**
 * \typedef ION_T 
 * \brief An object to represent a (fragment) ion of a peptide
 */
typedef struct ion ION_T;

/**
 * \typedef ION_SERIES_T 
 * \brief An object to represent a series of ions
 */
typedef struct ion_series ION_SERIES_T;

/**
 * \typedef ION_CONSTRAINT_T
 * \brief An object to represent a constraint to be applied to ions
 */
typedef struct ion_constraint ION_CONSTRAINT_T;

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
 * The enum for an ion type (P_ion is the precursor ion)
 * BY_ION(B & Y ion), BYA_ION(B & Y & A ion)
 */
enum _ion_type {A_ION, B_ION, C_ION, X_ION, Y_ION, Z_ION, 
  P_ION, BY_ION, BYA_ION, ALL_ION};
#define NUMBER_ION_TYPES 10

/**
 * \typedef ION_TYPE_T
 * \brief The typedef for ion type (a,b,c,x,y,z ions)
 */
typedef enum _ion_type ION_TYPE_T;

/**
 * The enum for an ion modification
 */
enum _ion_modification {NH3, H2O, ISOTOPE, FLANK, ALL_MODIFICATION}; 

/**
 * \typedef ION_MODIFICATION_T
 * \brief The typedef for ion modification type (NH3, H2O etc.)
 */
typedef enum _ion_modification ION_MODIFICATION_T;

/**
 * \typedef BIN_PEPTIDE_ITERATOR_T
 * \brief An iterator to iterate over the peptides in a bin( one file handler)
 */
typedef struct bin_peptide_iterator BIN_PEPTIDE_ITERATOR_T;

/**
 * \typedef BIN_SORTED_PEPTIDE_ITERATOR_T
 * \brief Object to iterate over the peptides within a bin, in an
 * sort in mass
 */
typedef struct bin_sorted_peptide_iterator BIN_SORTED_PEPTIDE_ITERATOR_T;

/**
 * \typedef  PROTEIN_INDEX_T
 * \brief Object to store the protein relation to the fasta file
 */
typedef struct protein_index PROTEIN_INDEX_T;

/**
 * \typedef PROTEIN_INDEX_ITERATOR_T
 * \brief Object to iterate over the protein index in the protein index file
 */
typedef struct protein_index_iterator PROTEIN_INDEX_ITERATOR_T;

/**
 * \typedef ION_ITERATOR_T
 * \brief An object to iterate over all ion objects in the ion_series
 */
typedef struct ion_iterator ION_ITERATOR_T;

/**
 * \typedef ION_FILTERED_ITERATOR_T
 * \brief An object to iterate over ion objects that meet constraint in the ion_series
 */
typedef struct ion_filtered_iterator ION_FILTERED_ITERATOR_T;

/**
 *\typedef LOSS_LIMIT_T
 *\brief An object that specifies the max amount of neutral loss possible at a given cleavage index
 * all numbers are for forward ions(A,B,C) subtract from total to get reverse limit
 */
typedef struct loss_limit LOSS_LIMIT_T;

/**
 * \typedef SCORER_T
 * \brief An object to score a spectrum v. ion_series or spectrum v. spectrum
 */
typedef struct scorer SCORER_T;

/**
 * The enum for scorer type
 * The QRANKER scores were added after Crux had been released for a while.
 * We don't want to change the CSM file layout at this point, so when
 * reading and writing scores from the CSM we omit the last two scores
 * in the score type enum.
 */
#define NUMBER_SCORER_TYPES 16 //BF added for consistant naming
//enum _scorer_type { SP, XCORR, DOTP, LOGP_EXP_SP, LOGP_BONF_EXP_SP, LOGP_EVD_XCORR, LOGP_BONF_EVD_XCORR, LOGP_WEIBULL_SP, LOGP_BONF_WEIBULL_SP, LOGP_WEIBULL_XCORR, LOGP_BONF_WEIBULL_XCORR, Q_VALUE, PERCOLATOR_SCORE, LOGP_QVALUE_WEIBULL_XCORR};
enum _scorer_type { 
  SP,                  ///< SEQUEST preliminary score
  XCORR,               ///< SEQUEST primary score
  DOTP,                ///< not yet implemented
  LOGP_EXP_SP,                     // this spot hijacked for zscore
  //ZSCORE,            ///< z-score (mean-max)/stdev
  //LOGP_BONF_EXP_SP,              // this spot hijacked for decoy-x-qval
  DECOY_XCORR_QVALUE,  ///< Benjamini-Hochberg q-value from xcorrs
  //LOGP_EVD_XCORR,               // this spot hijacked for decoy-p-qval
  DECOY_PVALUE_QVALUE, ///< Benjamini-Hochberg q-value from Weibull p-vals
  LOGP_BONF_EVD_XCORR,
  LOGP_WEIBULL_SP,
  LOGP_BONF_WEIBULL_SP,
  LOGP_WEIBULL_XCORR,
  LOGP_BONF_WEIBULL_XCORR,
  Q_VALUE,
  PERCOLATOR_SCORE,
  LOGP_QVALUE_WEIBULL_XCORR,
  QRANKER_SCORE,
  QRANKER_Q_VALUE
};

/*

enum _scorer_type { SP, XCORR, DOTP, 
LOGP_BONF_WEIBULL_SP, sp-logp

LOGP_BONF_WEIBULL_XCORR, xcorr-logp
};
*/
#define _SCORE_TYPE_NUM 16 ///< the number of different score types

/**
 * \typedef SCORER_TYPE_T
 * \brief The typedef for scorer type (SP, XCORR, DOTP)
 */
typedef enum _scorer_type SCORER_TYPE_T;

/**
 *\typedef GENERATE_PEPTIDES_ITERATOR_T
 *\brief An object that navigates the options and selects the correct peptide iterator to use
 */
typedef struct generate_peptides_iterator_t GENERATE_PEPTIDES_ITERATOR_T;

/**
 *\typedef HIT_T
 *\brief An object that contains the a protein and its score. 
 */
typedef struct hit HIT_T;

/**
 *\typedef HIT_COLLECTION_T
 *\brief An object that contains multiple hit objects
 */
typedef struct hit_collection HIT_COLLECTION_T;

/**
 *\typedef HIT_ITERATOR_T
 *\brief An object that navigates the hits in a hit collection
 */
typedef struct hit_iterator HIT_ITERATOR_T;

/**
 * The enum for protein scorer type
 */
#define NUMBER_PROTEIN_SCORER_TYPES 2
enum _protein_scorer_type { PROTEIN_SCORER_PVALUE, PROTEIN_SCORER_OLIVER };

/**
 * \typedef PROTEIN_SCORER_TYPE_T
 * \brief The typedef for protein scorer type
 */
typedef enum _protein_scorer_type PROTEIN_SCORER_TYPE_T;

/**
 *\typedef MATCH_T
 *\brief An object that contains the information of a peptide and the scoring of multiple types
 */
typedef struct match MATCH_T;

/**
 *\typedef MATCH_COLLECTION_T
 *\brief An object that contains mutiple match objects
 */
typedef struct match_collection MATCH_COLLECTION_T;

/**
 *\typedef MATCH_ITERATOR_T
 *\brief An object that navigates the matches
 */
typedef struct match_iterator MATCH_ITERATOR_T;

/**
 *\typedef MATCH_COLLECTION_ITERATOR_T
 *\brief An object that navigates the match_collection objects
 */
typedef struct match_collection_iterator MATCH_COLLECTION_ITERATOR_T;

#define NUMBER_ALGORITHM_TYPES 6 //BF added for consistant naming
/**
 * The enum for algorithm type (PERCOLATOR, CZAR, ALL)
 */
enum _algorithm {PERCOLATOR_ALGORITHM, RCZAR_ALGORITHM, QVALUE_ALGORITHM, NO_ALGORITHM, ALL_ALGORITHM, QRANKER_ALGORITHM};

/**
 * \typedef ALGORITHM_TYPE_T
 * \brief The typedef for _algorithm (PERCOLATOR, CZAR, ALL)
 */
typedef enum _algorithm ALGORITHM_TYPE_T;

/**
 * One value for each command that can be passed to crux
 * (e.g. search-for-matches, sequest-search, percolator).
 */
enum _command {
  INVALID_COMMAND,      ///< required by coding standards
  INDEX_COMMAND,        ///< create-index
  SEARCH_COMMAND,       ///< search-for-matches
  SEQUEST_COMMAND,      ///< sequest-search
  QVALUE_COMMAND,       ///< compute-q-values
  PERCOLATOR_COMMAND,   ///< percolator
  QRANKER_COMMAND,      ///< q-ranker
  PROCESS_SPEC_COMMAND, ///< print-processed-spectra
  XLINK_SEARCH_COMMAND, ///< search-for-xlinks	

  NUMBER_COMMAND_TYPES  ///< always keep this last so the value
                        /// changes as cmds are added
};

typedef enum _command COMMAND_T;
/**
 * \typedef RECORD_T
 * \brief RECORD_T for each value/key pair
 */
typedef struct record RECORD_T;

/**
 * \typedef HASH_T
 * \brief HASH_T hash table, contains the records
 */
typedef struct hash HASH_T;

/**
 * \typedef HASH_ITERATOR_T
 * \brief HASH_ITERATOR_T iterator for keys in a hash
 */
typedef struct hash_iterator HASH_ITERATOR_T;



/**
 * Identifying which set the PSM belongs to
 */
enum  _set_type {SET_TARGET=0,SET_DECOY1,SET_DECOY2,SET_DECOY3};

/**
 * \typedef SET_TYPE_T
 * \brief the typedef for set types for match type TARGET, DECOY1,
 * DECOY2, DECOY3 
 */
typedef enum _set_type SET_TYPE_T;

/**
 * \typedef MODIFIED_AA_T
 * \brief The alternate type for encoding a peptide sequence (instead
 * of char).  Allows modifications to be added to each AA.  See
 * modifications.h for more details.
 */
// why doesn't this work when I put it in modifications.h????
typedef unsigned short MODIFIED_AA_T; ///< letters in the expanded peptide

/**
 * \typedef AA_MOD_T
 * \brief The struct _aa_mod is typdefed as AA_MOD_T
 */
typedef struct _aa_mod AA_MOD_T;

/**
 * \typedef PEPTIDE_MOD_T
 * \brief The struct _peptide_mod is typdefed as PEPTIDE_MOD_T
 */
typedef struct _peptide_mod PEPTIDE_MOD_T;

/**
 * \enum _mod_position (typedefed as MOD_POSITION_T)
 * \brief An indication of where an AA_MOD may occur within a peptide.
 * Default is ANY_POSITION.
 */
enum _mod_position{ 
  ANY_POSITION, ///< at any position in any peptide
  C_TERM, ///< only c-terminus of peptide, seq[0]
  N_TERM  ///< only n-terminus of peptide, seq[len]
};

/**
 * \typedef MOD_POSITION_T
 * \brief The typedef of the indicator for where an amino acid
 * modification can occur within a peptide and/or protein.
 */
typedef enum _mod_position MOD_POSITION_T;

/**
 * \typedef _linked_list_head is typdefed as LINKED_LIST_T*
 * All list actions can be performed with an object of this type.
 */
typedef struct _linked_list_head LINKED_LIST_T;

/**
 * \typedef _linked_list_node is typdefed as LINKED_LIST_T*
 * This is an element of a list.  Can be used for adding to the end of
 * list or walking through a list.  Cannot allocate a new one.
 */
typedef struct _linked_list_node LIST_POINTER_T;

/**
 * \typedef modified_peptides_iterator_t is typedefed as
 * MODIFIED_PEPTIDES_ITERATOR_T 
 */
typedef struct modified_peptides_iterator_t MODIFIED_PEPTIDES_ITERATOR_T;

#endif

