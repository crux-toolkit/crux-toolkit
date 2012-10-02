/**
 * \file index.h 
 * $Revision: 1.25 $
 * \brief Object for representing an index of a index
 *****************************************************************************/
#ifndef INDEX_H 
#define INDEX_H

#include <stdio.h>

#include "utils.h"
#include "objects.h"
#include "Peptide.h"
#include "Protein.h"
#include "carp.h"
#include "PeptideConstraint.h"
#include "Database.h"
#include "Alphabet.h"
#include "IndexMap.h"


static const int MAX_INDEX_FILES = 1024;


/**
 * \class index
 * \brief A index of a database
 *
 * An index consists of three parts: a file with protein sequences, 
 * a set of files that point into the sequence file, and a map file
 * listing which files contain pointers for which peptide masses.  All
 * files are all contained in a single directory.
 *
 * The directory is also the name of the index given by the user.
 * There are no assumptions made about the name of the index
 * (directory) and the name of the fasta file from which it was
 * created.  The binary version of the fasta is named [fasta
 * filename]-binary-fasta.  The file with a list of index files and
 * the masses they contain is crux_index_map.
 *
 * The index has two peptide constraints: one to define ALL the
 * peptides that are indexed and one to define the peptides being
 * searched for.  When an index is being created, the search
 * constraint is NULL.  When an index is being loaded from disk, the
 * disk_constraint is set according to the values in the header of
 * crux_index_map.  Each time a peptide generator is created for the
 * index, the search_constraint is set (any existing ones being
 * deleted).  
 *
 * Bounds checking for a search is done at two points.  When the index
 * is loaded from disk, the global constraints for the search are
 * compared to the constraints for the index and the process is
 * terminated if they are not compatible. (see
 * check_index_constraints)  When a peptide iterator is created and a 
 * new search_constraint is set, only the min and max mass is checked
 * (as no other parts of the constraint change in the course of a
 * search).  If the iterator requests a mass range partly outside the
 * bounds of the index, a warning is printed (one for over max and one
 * for under min).
 */
class Index {

 protected:
  int num_pointers_; ///< The number of pointers to this index.
  Database* database_; ///< The database that has been indexed.
  Database* decoy_database_; ///< The decoy verision of the database.
  IndexMap* index_map_;
  IndexMap* decoy_index_map_;
  char* directory_; ///< The directory containing the indexed files
  PeptideConstraint* disk_constraint_;///< Defines peptides on disk
  PeptideConstraint* search_constraint_;///< Defines peptides being searched
  bool on_disk_; ///< Does this index exist on disk yet?
  FLOAT_T mass_range_;  ///< the range of masses in each index file -in
  bool is_unique_; ///< only unique peptides? -in
  DECOY_TYPE_T decoys_; ///< the type of decoys stored
  FLOAT_T* static_mods_; ///< static mod for each amino acid

  /* Private Functions */

  /**
   * initializes an Index object.
   */
  void init();

  /**
   * \brief Initializes a new Index object by setting its member
   * variables to the values in the index on disk.
   *
   * Assumes that the name of the index (directory) has been set by
   * caller.  Checks that the directory exists.  Reads the mapfile
   * therein and gets the constraint values from the header.  Returns
   * true if all of the above are true, else prints an error and returns
   * false. 
   *
   * \returns true if fields were successfully set, false if index on
   * disk could not be used to initialize the object.
   */
  bool setFieldsFromDisk();

  /**
   * \brief Private function to take a header line from index_map and
   * set the appropriate value in the index.
   */
  void setFieldFromMap(char* line);

  /**
   * \brief Private function to confirm that the attributes of the index
   * will work for this run.
   * Requires that the range from min to max length and mass include the
   * range requested in parameter.c.  The index must have at least as
   * many missed cleavages as requested.  The mass type must be the
   * same. Requires index and request be both unique or both not 
   * unique.  The peptide cleavage type of the index must be no more
   * restrictive than requested (where TRYPTIC is most restrictive and
   * ANY_TRYPTIC is least).
   * \returns true if all constraints will work with those in parameter.c.
   */
  bool checkConstraints();

  /**
   * \brief The initialization step in creating an index after
   * allocation and before parsing database.  Private function called
   * only by new_index().  No assumptions about fasta
   * or index filenames or locations relative to eachother.  Does not
   * initialize the database to be associated with this index.
   *
   * \returns void
   */
  void setFields(
    const char* output_dir,      ///< The name of the new index
    PeptideConstraint* constraint,  
    ///< Constraint which these peptides satisfy -in
    FLOAT_T mass_range,  
    ///< the range of mass that each index file should be partitioned into -in
    bool is_unique, ///< only unique peptides? -in
    DECOY_TYPE_T decoys ///< the kind of decoys to store
    );
  
  /**
   * write to the file stream various information of the
   * index files created
   */
  bool writeHeader(
    FILE* file ///< out put stream for crux_index_map -in
    );

  /**
   * write to the file stream various information of the
   * index files created in human readable format
   *\returns true if successfully creates README file, else false
   */
  bool writeReadmeFile(
    FILE* file ///< out put stream for README file -in
    );

  /**
   * Calculate the total number of bins( file handlers) that will be needed
   * \returns the total number of bins needed
   */
  long getNumBinsNeeded(
    int* mass_limits  ///< an array that holds the min/max mass limit -in
    );

  /**
   * \brief For one file bin, reparses the peptides, sorts them, and
   * reprints them to the crux_index file 
   *
   * \returns the sorted bin
   */
  FILE* sortBin(
    FILE* file, ///< the working file handle to the bin -in
    long bin_idx, ///< bin index in the file array -in
    unsigned int peptide_count, ///< the total peptide count in the bin -in
    FILE* text_file, ///< optional file to write sequences to -in
    const char* file_prefix, ///< beginning of index file names -in
    Database* database ///< target or decoy database to use -in
    );

  /**
   * The steps of creating an index that are repeated for the target and
   * decoy databases.  For the decoy, the info_file and text_file can be NULL.
   */
  void index_database(
    Database* database, ///< the database to index
    const char* file_prefix, ///< name for the index files
    FILE* info_file, ///< index map
    FILE* text_file); ///< optional peptides file

 public:
  // Member constants
  static const char* index_file_prefix;
  static const char* decoy_index_file_prefix;

  /**
   * \returns An (empty) index object.
   */
  Index();

  /**
   * Constructor for creating a new, empty index to be
   * populated from a fasta file. 
   */
  Index(
    const char* fasta_filename,  ///< The fasta file
    const char* output_dir,      ///< The name of the new index
    PeptideConstraint* constraint,///< Constraint which these peptides satisfy
    FLOAT_T mass_range,  ///< the range of masses contained in each index file
    DECOY_TYPE_T decoys ///< the type of decoys to store
    );         

  /**
   * Constructor for opening an existing index for searching.
   */
  Index(
    const char* fasta_filename  ///< The fasta file
    );

  /**
   * Frees an allocated index object.
   */
  static void free(Index* index);
  virtual ~Index();

  /**
   * Increments the index pointer count.
   */
  Index* copyPtr();

  /**
   * Scans the index on disk to populate fields in index.
   * \returns True if success, false if failure.
   */
  bool parse();

  /**
   * The main index method. Does all the heavy lifting, creating files
   * serializing peptides, etc.
   *
   * \returns TRUE if success. FALSE if failure.
   */
  bool create(
    bool create_text_file ///< Should an ASCII text file be create? -in
    );

  /**
   * \brief Looks in given directory for a file ending in
   * "-binary-fasta" and returns a heap-allocated string of the full
   * name including the index directory.
   *
   * Exits with error if index_name does not exist, no file
   * *-binary-fasta exists, or more than one *binary-fasta file exists.
   * \returns A string with the name of the existing binary fasta file
   * for this index.
   */
  static char* getBinaryFastaName(
    const char* index_name ///< directory containing index
    );

  /**
   * \brief Looks in given directory for a file ending in
   * "-decoy-binary-fasta" and returns a heap-allocated string of the
   * full name including the index directory.
   *
   * Exits with error if index_name does not exist, no file
   * *-decoy-binary-fasta exists, or more than one *binary-fasta file exists.
   * \returns A string with the name of the existing binary fasta file
   * for this index.
   */
  static char* getDecoyBinaryFastaName(
    const char* index_name ///< directory containing index
    );


  /*********************************************
   * set and get methods for the object fields
   *********************************************/

  /**
   * \returns A const pointer to the directory name.
   */
  const char* getDirectory();

  /**
   * Sets the directory of the index.
   * index->directory must been initiailized.
   */
  void setDirectory(const char* directory); ///< the directory to add -in
  
  /**
   * \returns A pointer to the database.
   */
  Database* getDatabase(bool is_decoy); ///< return target or decoy 

  /**
   *\returns a pointer to the peptides constraint
   */
  PeptideConstraint* getSearchConstraint();

  /**
   * \brief Sets the peptide search constraint to be used by the
   * generate_peptides_iterator.  Makes a copy of the constraint pointer.
   * Deletes any existing search constraint. 
   */
  void setSearchConstraint(
    PeptideConstraint* constraint ///< Constraint for the next iterator
    );

  /**
   * \returns the disk constraint pointer from index
   */
  PeptideConstraint* getDiskConstraint();

  /**
   * \returns the IndexMap for the target or decoy index
   */
  IndexMap* getIndexMap(
    bool decoy ///< return decoy index?
  );

  /**
   *\returns TRUE if only allow unique peptides else FALSE
   */
  bool getIsUnique();
                           

  /**
   * \returns The number of proteins in the index.
   */
  int getNumProteins();

  /**
   * \returns The type of decoys stored in the database.
   */
  DECOY_TYPE_T getDecoyType();
};



/***********************************************
 * Iterators index
 ***********************************************/

/**********************************************************************
 * wrapper, for generate_peptides_iterator, cast back to original type
 ***********************************************************************/
/**
 * Frees an allocated index_peptide_iterator object.
 */
void void_free_index_peptide_iterator(
    void* index_peptide_iterator  ///< the iterator to free -in
    );

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over,
 * FALSE if not.
 */
bool void_index_peptide_iterator_has_next(
    void* index_peptide_iterator ///< the iterator of interest -in
    );

/**
 * \returns The next peptide in the index.
 */
Crux::Peptide* void_index_peptide_iterator_next(
    void* index_peptide_iterator ///< the iterator of interest -in
    );



/***********************************************
 * Iterators BIN
 ***********************************************/


/**
 * Instantiates a new bin_peptide_iterator from a gvien bin file handler.
 * \returns a new heap allocated bin_peptide_iterator object
 */
BIN_PEPTIDE_ITERATOR_T* new_bin_peptide_iterator(
  Index* index, ///< The index object which we are iterating over -in
  FILE* file ///< the bin to parse peptides
  );

/**
 *  The basic iterator functions.
 * \returns The next peptide in the index.
 */
Crux::Peptide* bin_peptide_iterator_next(
  BIN_PEPTIDE_ITERATOR_T* bin_peptide_iterator 
   ///< the bin_peptide_iterator to get peptide -in
  );

/**
 * The basic iterator functions.
 * check to see if the bin_peptide_iterator has more peptides to return
 * \returns TRUE if there are additional peptides to iterate over,
 * FALSE if not.
 */
bool bin_peptide_iterator_has_next(
  BIN_PEPTIDE_ITERATOR_T* bin_peptide_iterator 
  ///< the bin_peptide_iterator to initialize -in
  );

/**
 * Frees an allocated bin_peptide_iterator object.
 */
void free_bin_peptide_iterator(
  BIN_PEPTIDE_ITERATOR_T* bin_peptide_iterator ///< the iterator to free -in
  );


/***********************************************
 * Iterators sorted BIN
 ***********************************************/

/**
 * Instantiates a new sorted_bin_peptide_iterator from a gvien bin
 * file handle.
 * \returns A new heap allocated sorted_bin_peptide_iterator object.
 */
BIN_SORTED_PEPTIDE_ITERATOR_T* new_bin_sorted_peptide_iterator(
  Index* index, ///< The index object which we are iterating over -in
  FILE* file,///< the working file handler to the bin -in
  unsigned int peptide_count, ///< the total peptide count in the bin -in
  Database* database ///< target or decoy database
  );

/**
 *  The basic iterator functions.
 * \returns The next peptide in the index.
 */
Crux::Peptide* bin_sorted_peptide_iterator_next(
  BIN_SORTED_PEPTIDE_ITERATOR_T* bin_sorted_peptide_iterator 
  ///< the bin_peptide_iterator to get peptide -in
  );

/**
 * The basic iterator functions.
 * Check to see if the bin_sorted_peptide_iterator has more peptides
 * to return.
 *\returns TRUE if there are additional peptides to iterate over,
 * FALSE if not.
 */
bool bin_sorted_peptide_iterator_has_next(
  BIN_SORTED_PEPTIDE_ITERATOR_T* bin_sorted_peptide_iterator 
  ///< the bin_peptide_iterator to initialize -in
  );

/**
 * Frees an allocated bin_peptide_iterator object.
 */
void free_bin_sorted_peptide_iterator(
  BIN_SORTED_PEPTIDE_ITERATOR_T* bin_sorted_peptide_iterator 
  ///< the iterator to free -in
  );


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
