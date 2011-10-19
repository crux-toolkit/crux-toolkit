/**
 * \file MatchCollectionIterator.cpp 
 * $Revision: 1.38 $
 * \brief An object that iterates over the MatchCollection objects in
 * the specified directory of serialized MatchCollections 
 */

#include "MatchCollectionIterator.h"

/******************************
 * match_collection_iterator
 ******************************/
     
/**
 * \brief Finds the next match_collection in directory and prepares
 * the iterator to hand it off when 'next' called.
 *
 * When no more match_collections (i.e. psm files) are available, set
 * match_collection_iterator->has_another_collection to false
 * \returns void
 */
void MatchCollectionIterator::setup()
{
  // are there any more match_collections to return?
  if(collection_idx_ < number_collections_){
    // then go parse the match_collection
    match_collection_ = 
      new MatchCollection(this, (SET_TYPE_T)collection_idx_);

    // we have another match_collection to return
    has_another_collection_ = true;
    
    // let's move on to the next one next time
    ++collection_idx_;

    // reset directory
    rewinddir(working_directory_);
  }
  else{
    // we're done, no more match_collections to return
    has_another_collection_ = false;
  }
}

/**
 * Create a match_collection iterator from a directory of serialized files.
 * Only handles up to one target and three decoy sets per folder.
 *\returns match_collection iterator instantiated from a result folder
 */
MatchCollectionIterator::MatchCollectionIterator(
  const char* output_file_directory, 
    ///< the directory path where the PSM output files are located -in
  const char* fasta_file, 
    ///< The name of the fasta file for peptides for match_collections. -in
  int* decoy_count
) :
  working_directory_(NULL), directory_name_(NULL), database_(NULL),
  decoy_database_(NULL), number_collections_(0), collection_idx_(-1),
  match_collection_(NULL), has_another_collection_(false)
{
  carp(CARP_DEBUG, 
       "Creating match collection iterator for dir %s and protein database %s",
       output_file_directory, fasta_file);

  struct dirent* directory_entry = NULL;
  bool use_index = is_directory(fasta_file);

  /*
    BF: I think that this step is to count how many decoys there are
    per target file.  This is prone to errors as all it really does is
    check for the presence of a file with *decoy_1*, and one with
    *decoy_2* and *decoy_3*.  In fact, the three files could be from
    different targets.  Nothing was being done with the check for a
    target file.  There must be a better way to do this.
   */


  // do we have these files in the directory
  bool boolean_result = false;
  bool decoy_1 = false;
  bool decoy_2 = false;
  bool decoy_3 = false;

  // open PSM file directory
  working_directory_ = opendir(output_file_directory);
  
  if (working_directory_ == NULL) {
    carp(CARP_FATAL, "Failed to open PSM file directory: %s", 
        output_file_directory);
  }
  
  // determine how many decoy sets we have
  while((directory_entry = readdir(working_directory_))){
    
    if(suffix_compare(directory_entry->d_name, "decoy-1.txt")) {
      carp(CARP_DEBUG, "Found decoy file %s", directory_entry->d_name);
      decoy_1 = true;
    }
    else if(suffix_compare(directory_entry->d_name, "decoy.txt")) {
      decoy_1 = true;
    }
    else if(suffix_compare(directory_entry->d_name, "decoy-2.txt")) {
      decoy_2 = true;
    }
    else if(suffix_compare(directory_entry->d_name, "decoy-3.txt")) {
      decoy_3 = true;
    }    
    else if(suffix_compare(directory_entry->d_name, ".txt")){
      carp(CARP_DEBUG, "Found target file %s", directory_entry->d_name);
      boolean_result = true;
    }
    if (boolean_result && decoy_1 && decoy_2 && decoy_3) {
      break; // We've found all the files we can use.
    }
  }
  
  // set total_sets count
  int total_sets = 0;

  if(decoy_3){
    total_sets = 4; // 3 decoys + 1 target
    *decoy_count = 3;
  }
  else if(decoy_2){
    total_sets = 3; // 2 decoys + 1 target
    *decoy_count = 2;
  }
  else if(decoy_1){
    total_sets = 2; // 1 decoys + 1 target
    *decoy_count = 1;
  }
  else{
    total_sets = 1;
    *decoy_count = 0;
    carp(CARP_INFO, "No decoy sets exist in directory: %s", 
        output_file_directory);
  }
  if(!boolean_result){
    carp(CARP_FATAL, "No PSM files found in directory '%s'", 
         output_file_directory);
  }

  // get binary fasta file name with path to crux directory 
  if (use_index == true){ 
    char* binary_fasta = Index::getBinaryFastaName(fasta_file);
    database_ = new Database(binary_fasta, true);// is memmapped
    free(binary_fasta);
    binary_fasta = Index::getDecoyBinaryFastaName(fasta_file);
    if( binary_fasta != NULL ){
      decoy_database_ = new Database(binary_fasta, true);// is memmapped
      decoy_database_->parse();
      free(binary_fasta);
    }
  } else {
    database_ = new Database(fasta_file, false);// not memmapped
    database_->transformTextToMemmap(".", true);// is temp
    decoy_database_ = NULL;
  }
  database_->parse();

  // reset directory
  rewinddir(working_directory_);
 

  // set match_collection_iterator fields
  collection_idx_ = 0;
  number_collections_ = total_sets;
  directory_name_ = 
    my_copy_string(output_file_directory);
  has_another_collection_ = false;

  carp(CARP_DETAILED_DEBUG,"num collections:%d",number_collections_);

  cols_in_file_ = new vector<bool>();

  // setup the match collection iterator for iteration
  // here it will go parse files to construct match collections
  setup();
  carp(CARP_DEBUG, "Done creating match collection iterator");
}

/**
 *\returns true, if there's another match_collection to return, else return false
 */
bool MatchCollectionIterator::hasNext()
{
  // Do we have another match_collection to return
  return has_another_collection_;
}

/**
 * free match_collection_iterator
 */
MatchCollectionIterator::~MatchCollectionIterator()
{
  // free unclaimed match_collection
  if(match_collection_ != NULL){
    delete match_collection_;
  }
  
  // free up all match_collection_iteratory.
  free(directory_name_);
  Database::freeDatabase(database_);
  Database::freeDatabase(decoy_database_);
  closedir(working_directory_); 
  delete cols_in_file_;
}

/**
 * \brief Fetches the next match collection object and prepares for
 * the next iteration 
 *\returns The next match collection object
 */
MatchCollection* MatchCollectionIterator::next()
{
  MatchCollection* match_collection = NULL;
  
  if(has_another_collection_){
    match_collection = match_collection_;
    match_collection_ = NULL;
    setup();
    return match_collection;
  }
  else{
    carp(CARP_ERROR, "No match_collection to return");
    return NULL;
  }
}

/**
 *\returns the database
 */
Database* MatchCollectionIterator::getDatabase() {
  return database_;
}
    
/**
 *\returns the decoy database
 */
Database* MatchCollectionIterator::getDecoyDatabase() {
  return decoy_database_;
}
    

/**
 *\returns the total number of match_collections to return
 */
int MatchCollectionIterator::getNumberCollections()
{
  return number_collections_;
}

/**
 * \brief Get the name of the directory the match_collection_iterator
 * is working in.
 * \returns A const pointer to the directory name.
 */
const char* MatchCollectionIterator::getDirectoryName()
{
  return directory_name_;
}

/**
 * \brief Get the working directory
 */
DIR* MatchCollectionIterator::getWorkingDirectory() {

  return working_directory_;
}

vector<bool>& MatchCollectionIterator::getColsInFile(){

  return *cols_in_file_;
}
