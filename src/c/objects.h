
/**
 * \file objects.h 
 * $Revision: 1.62 $
 * \brief The defined objects
 *****************************************************************************/
#ifndef OBJECTS_H 
#define OBJECTS_H

#include <stdio.h>
#include <set>
#include <map>
#include "utils.h"

#define QSORT_COMPARE_METHOD int(*)(const void*, const void*)


class DelimitedFile;
class DelimitedFileReader;
class IonConstraint;
class IonFilteredIterator;
class IonSeries;
class MatchFileReader;
class Spectrum;
class SpectrumZState;


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
 * \class SpectrumCollection
 * \brief A collection of spectra
 */
class SpectrumCollection;

/**
 * \typedef SPECTRUM_ITERATOR_T 
 * \brief An object to iterate over the spectra in a spectrum_collection
 */
typedef struct spectrum_iterator SPECTRUM_ITERATOR_T;

/**
 * \class FilteredSpectrumChargeIterator 
 * \brief An object to iterate over the spectra in a spectrum_collection
 */
class FilteredSpectrumChargeIterator;

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
  NUMBER_DIGEST_TYPES  ///< keep last, number of types
};

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
  CUSTOM_ENZYME,         ///< cleave after/before user-defined residues
  NUMBER_ENZYME_TYPES    ///< leave last, number of types
};

/**
 * \typedef ENZYME_T
 * \brief The enzyme with which a peptide was digested.  Used in
 * conjunction with DIGEST_T to define how peptides are generated.
 */
typedef enum _enzyme_type ENZYME_T;

/**
 * The enum for isotopic mass type (average, mono)
 */
enum _mass_type {AVERAGE, MONO, NUMBER_MASS_TYPES };

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
  WINDOW_PPM,
  NUMBER_WINDOW_TYPES  
};


/**
 * The enum for measure type for spectral counts
 */
enum _measure_type {
  MEASURE_INVALID,
  MEASURE_SIN,
  MEASURE_NSAF,
  MEASURE_EMPAI,
  NUMBER_MEASURE_TYPES
};

/**
 * \typedef MEASURE_TYPE_T
 * \brief The typedef for measure type (sin, nsaf)
 */
typedef enum _measure_type MEASURE_TYPE_T;

/**
 * The quantification level type for spectral counts
 */
enum _quant_level_type {
  QUANT_LEVEL_INVALID,
  PEPTIDE_QUANT_LEVEL,
  PROTEIN_QUANT_LEVEL,
  NUMBER_QUANT_LEVEL_TYPES
};

/**
 * \typedef QUANT_LEVEL_TYPE_T
 * \brief The typdef for quantificaiton level (peptide, protein)
 */
typedef enum _quant_level_type QUANT_LEVEL_TYPE_T;


/**
 * The enum for parsimony type for spectral counts
 */
enum _parsimony_type {
  PARSIMONY_INVALID,
  PARSIMONY_SIMPLE,
  PARSIMONY_GREEDY,
  PARSIMONY_NONE,
  NUMBER_PARSIMONY_TYPES
};

/*
 * \typedef PARSIMONY_TYPE_T
 * \brief The typedef for parsimony type (simple, greedy, none)
 */
typedef enum _parsimony_type PARSIMONY_TYPE_T;


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
 * \class Protein
 * \brief A protein sequence
 */
class Protein;


/**
 * \class ProteinPeptideIterator
 * \brief An object to iterate over the peptides in a protein sequence
 */
class ProteinPeptideIterator;

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
enum _sort_type {SORT_NONE, 
                 SORT_MASS, 
                 SORT_LENGTH, 
                 SORT_LEXICAL, 
                 NUMBER_SORT_TYPES };

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
                P_ION, BY_ION, BYA_ION, ALL_ION, NUMBER_ION_TYPES };

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
 */
enum _scorer_type { 
  SP,                  ///< SEQUEST preliminary score
  XCORR,               ///< SEQUEST primary score

  DECOY_XCORR_QVALUE,  ///< q-value derived from empirical null (decoys)
  DECOY_XCORR_PEPTIDE_QVALUE,

  LOGP_WEIBULL_XCORR,
  LOGP_BONF_WEIBULL_XCORR,
  LOGP_QVALUE_WEIBULL_XCORR,
  LOGP_PEPTIDE_QVALUE_WEIBULL,

  PERCOLATOR_SCORE,
  PERCOLATOR_QVALUE,
  PERCOLATOR_PEPTIDE_QVALUE,

  QRANKER_SCORE,
  QRANKER_QVALUE,
  QRANKER_PEPTIDE_QVALUE,

  NUMBER_SCORER_TYPES,
  INVALID_SCORER_TYPE
};

/**
 * \typedef SCORER_TYPE_T
 * \brief The typedef for scorer type.
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
enum _protein_scorer_type { PROTEIN_SCORER_PVALUE, 
                            PROTEIN_SCORER_OLIVER, 
                            NUMBER_PROTEIN_SCORER_TYPES };

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

/**
 * The enum for algorithm type (PERCOLATOR, CZAR, ALL)
 */
enum _algorithm {PERCOLATOR_ALGORITHM, 
                 RCZAR_ALGORITHM, 
                 QVALUE_ALGORITHM, 
                 NO_ALGORITHM, 
                 ALL_ALGORITHM, 
                 QRANKER_ALGORITHM, 
                 NUMBER_ALGORITHM_TYPES };

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
  SPECTRAL_COUNTS_COMMAND, ///< spectral counts
  QRANKER_COMMAND,      ///< q-ranker
  PROCESS_SPEC_COMMAND, ///< print-processed-spectra
  XLINK_SEARCH_COMMAND, ///< search-for-xlinks
  VERSION_COMMAND,      ///< just print the version number

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

/**
 * \enum _xlink_site_t (typedefed as XLINK_SITE_T)
 * \brief An indication of the type of the crosslinking site that
 * may occur in a peptide.
 * Default is UNKNOWN.
 */
enum XLINK_SITE_T{
  XLINKSITE_UNKNOWN,
  XLINKSITE_NTERM,
  XLINKSITE_ALL,
  XLINKSITE_AA,
  NUMBER_XLINKSITES
};

/**
 * \enum COMPARISON_T
 */
enum COMPARISON_T{
  COMPARISON_INVALID,
  COMPARISON_LT,
  COMPARISON_LTE,
  COMPARISON_EQ,
  COMPARISON_GTE,
  COMPARISON_GT,
  COMPARISON_NEQ,
  NUMBER_COMPARISONS,
};

/**
 * \enum COLTYPE_T
 */
enum COLTYPE_T{
  COLTYPE_INVALID,
  COLTYPE_INT,
  COLTYPE_REAL,
  COLTYPE_STRING,
  NUMBER_COLTYPES
};


/**
 * \typedef peptideToScore
 * \brief Mapping of peptide object to scores
 */
typedef std::map<PEPTIDE_T*, FLOAT_T, bool(*)(PEPTIDE_T*, PEPTIDE_T*) > PeptideToScore;

/**
 * \typedef ProteinToScore
 * \brief Mapping of protein object to scores
 */
typedef std::map<Protein*, FLOAT_T, bool(*)(Protein*, Protein*) > ProteinToScore;

/**
 * \typedef MetaProtein
 * \brief Collection of protein objects which contain exactly the same
 * set of peptides.
 */
typedef std::set<Protein*, bool(*)(Protein*, Protein*) > MetaProtein;

/**
 * \typedef ProteinToMeta
 * \brief Mapping of Protein to MetaProtein to which it belongs
 */
typedef std::map<Protein*, MetaProtein, bool(*)(Protein*, Protein*) > ProteinToMetaProtein;

/**
 * \typedef MetaToRank
 * \brief Mapping of MetaProtein to ranks to the rank asigned to it
 */
typedef std::map<MetaProtein, int, bool(*)(MetaProtein, MetaProtein) > MetaToRank;

#endif

