/**
 * \file MatchSearch.cpp
 * BASED ON: original_match_search.c
 * DATE: Aug 19, 2008
 * AUTHOR: Barbara Frewen
 * \brief Main file for crux-search-for-matches.
 *
 * Given an ms2 file and a fasta file or index, compare all spectra to
 * peptides in the fasta file/index and return high scoring matches.
 * Peptides are determined by parameters for length, mass, mass
 * tolerance, cleavages, modifications. Score each spectrum with
 * respect to all candidates, and rank by score. Output in binary csm
 * file format or text sqt file format.
 */
/*
 * Here is the outline for how the new search should work

   for each spectrum
     for each charge state
      for each peptide modification
        create a peptide iterator
        for each peptide
         score peptide/spectra
      if passes criteria, print results and move on
      else next peptide modification  
 */

#include "MatchSearch.h"
#include "FilteredSpectrumChargeIterator.h"
#include "SearchProgress.h"
#include "SpectrumCollectionFactory.h"
#include "ModifiedPeptidesIterator.h"

using namespace std;
using namespace Crux;
/**
 * \returns a blank MatchSearch object
 */
MatchSearch::MatchSearch() {
}

/**
 * Destructor
 */
MatchSearch::~MatchSearch() {
}

/**
 * \brief Look at matches and search parameters to determine if a
 * sufficient number PSMs have been found.  Returns true if the
 * maximum number of modifications per peptide have been considered.
 * In the future, implement and option and test for a minimum score.
 * \returns true if no more PSMs need be searched.
 */
bool MatchSearch::isSearchComplete(
  MatchCollection* matches, ///< matches to consider
  int mods_per_peptide){ ///< modifications per peptide searched


  if( matches == NULL ){
    return false;
  }

  // keep searching if no limits on how many mods per peptide
  if( get_int_parameter("max-mods") == MAX_PEPTIDE_LENGTH ){
    return false;
  }
  // stop searching if at max mods per peptide
  if( mods_per_peptide == get_int_parameter("max-mods") ){ 
    return true;
  }

  // test for minimun score found

  return false;
  
}


/**
 * \brief Search the database OR index with up to num_peptide_mods from
 * the list for matches to the spectrum. 
 * Scored PSMs are added to the match_collection, possibly truncating
 * the collection and deleting existing matches in the collection.
 * After searching with each peptide mod, assess if there exists a
 * "good enough" match and end the search if there is, returning the
 * number of peptide mods that were searched.
 * \return The number of peptide mods searched.
 */
int MatchSearch::searchPepMods(
  MatchCollection* match_collection, ///< store PSMs here
  bool is_decoy,   ///< generate decoy peptides from index/db
  Index* index,       ///< index to use for generating peptides
  Database* database, ///< db to use for generating peptides
  Spectrum* spectrum, ///< spectrum to search
  SpectrumZState& zstate, ///< seach spectrum at this z-state
  PEPTIDE_MOD_T** peptide_mods, ///< list of peptide mods to apply
  int num_peptide_mods, ///< how many p_mods to use from the list
  bool compute_sp,  ///< compute sp scores
  bool store_scores///< save all scores for p-value estimation
  ){

  // set match_collection charge
  match_collection->setZState(zstate);

  // get spectrum precursor mz
  double mz = spectrum->getPrecursorMz();

  int mod_idx = 0;

  // assess scores after all pmods with x amods have been searched
  int cur_aa_mods = 0;

  // for each peptide mod
  for(mod_idx=0; mod_idx<num_peptide_mods; mod_idx++){
    // get peptide mod
    PEPTIDE_MOD_T* peptide_mod = peptide_mods[mod_idx];

    // is it time to assess matches?
    int this_aa_mods = peptide_mod_get_num_aa_mods(peptide_mod);
    
    if( this_aa_mods > cur_aa_mods ){
      carp(CARP_DEBUG, "Finished searching %i mods", cur_aa_mods);
      bool passes = isSearchComplete(match_collection, cur_aa_mods);
      if( passes ){
        carp(CARP_DETAILED_DEBUG, 
             "Ending search with %i modifications per peptide", cur_aa_mods);
        break;
      }// else, search with more mods
      cur_aa_mods = this_aa_mods;
    }
      
    //TODO SJM:  Figure out why this code gives different results for the sequest 
    //smoke test (this was changed in Rev. 2006).
    //      20014c20014
    //< S     21134   21134   3       0.00    server  2140.03 0.00    0.00    213
    //---
    //> S     21134   21134   3       0.00    server  2140.03 0.00    0.00    207
    
    
    // get peptide iterator
    ModifiedPeptidesIterator* peptide_iterator = new
      ModifiedPeptidesIterator(mz, zstate, peptide_mod, is_decoy,
                               index, database); 
    
    // score peptides
    int added = match_collection->addMatches(spectrum, 
                            zstate, 
                            peptide_iterator,
                            is_decoy,
                            store_scores,
                            compute_sp,
                            false // don't filtery by Sp
                            );
    
    carp(CARP_DEBUG, "Added %i matches", added);

    delete peptide_iterator;
    
  }//next peptide mod

  return mod_idx;
}

/**
 * Print the target and decoy match collections to their respective
 * target and decoy files.
 *
 * Three possibilities: 1. combine the target and all decoy
 * collections and print to target file.  2. print targets to target
 * file and combine all decoys and print to one decoy file.  3. print
 * each collection to a separate file.
 * Possible side effectos: Collections may be merged and re-ranked.
 */
void MatchSearch::printSpectrumMatches(
  OutputFiles& output_files, ///< files to print to
  MatchCollection* target_psms, ///< target psms to print
  vector<MatchCollection*>& decoy_psms, ///< decoy psms to print
  Spectrum* spectrum, ///< spectrum for all psms
  bool combine_target_decoy, ///< print targets and decoys to one file
  int num_decoy_files ///< number of decoy files to print
  ){

  // now print matches to one, two or several files
  if( combine_target_decoy ){
    // merge all collections
    MatchCollection* all_psms = target_psms;
    for(size_t decoy_idx = 0; decoy_idx < decoy_psms.size(); decoy_idx++){
      MatchCollection::merge(decoy_psms[decoy_idx], all_psms);
    }
    
    // sort and rank
    if( all_psms->getScoredType(SP) ){
      all_psms->populateMatchRank(SP);
      all_psms->saveTopSpMatch();
    }
    all_psms->populateMatchRank(XCORR);

    vector<MatchCollection*> empty_list;
    output_files.writeMatches(all_psms, // target matches
                              empty_list, // decoy matches
                              XCORR, spectrum); 
    
  }else{ // targets and decoys in separate files
    
    // if decoys in one file
    if( num_decoy_files == 1 ){
      // merge decoys
      MatchCollection* merged_decoy_psms = decoy_psms[0];
      for(size_t decoy_idx = 1; decoy_idx < decoy_psms.size(); decoy_idx++){
        MatchCollection::merge(decoy_psms[decoy_idx],
                                merged_decoy_psms);
      }
      
      // sort and rank
      // NOTE (BF 09-14-10): since the multiple decoy collections have already
      // been truncated, the merged ranks aren't accurate for the total space
      // of decoys searched
      merged_decoy_psms->populateMatchRank(XCORR);
      
      vector<MatchCollection*> decoy_list(1, merged_decoy_psms);
      output_files.writeMatches(target_psms, 
                                decoy_list, 
                                XCORR, spectrum);
      
    }else{
      // already sorted and ranked
      output_files.writeMatches(target_psms, decoy_psms, 
                                XCORR, spectrum);
    }
  }
}

// TODO this should be in match_collection
/**
 * Search the given database or index using shuffled peptides and the
 * spectrum/charge in the target psm match collection.  Add those
 * scores to the target psm match collection for use in weibull
 * parameter estimation but do not save the matches.  Repeat the
 * search with all peptide mods in the list.
 */
void MatchSearch::addDecoyScores(
  MatchCollection* target_psms, ///< add scores to these matches
  Spectrum* spectrum, ///< spectrum to score
  SpectrumZState& zstate, ///< charge/mass to use for spectrum
  Index* index, ///< search this index if not null
  Database* database, ///< search this database if not null
  PEPTIDE_MOD_T** peptide_mods, ///< list of peptide mods to search
  int num_peptide_mods ///< number of mods in the above array
){

  int mod_idx = 0;
  // for each peptide mod in the list
  for(mod_idx = 0; mod_idx < num_peptide_mods; mod_idx++){
    ModifiedPeptidesIterator* peptide_iterator = 
      new ModifiedPeptidesIterator(spectrum->getPrecursorMz(),
                                   zstate,
                                   peptide_mods[mod_idx],
                                   true, // is decoy
                                   index,
                                   database);
    target_psms->addDecoyScores(spectrum, 
                                zstate,
                                peptide_iterator);  
    delete peptide_iterator;
  }


}

int MatchSearch::main(int argc, char** argv){

  /* Define optional command line arguments */
  const char* option_list[] = {
    "fileroot",
    "output-dir",
    "overwrite",
    "sqt-output",
    "mzid-output",
    "pepxml-output",
    "txt-output",
    "num-decoys-per-target",
    "decoys",
    "decoy-location",
    "compute-sp",
    "compute-p-values",
    "spectrum-min-mz",
    "spectrum-max-mz",
    "spectrum-charge",
    "max-ion-charge",
    "scan-number",
    "mz-bin-width",
    "mz-bin-offset",
    "spectrum-parser",
    "parameter-file",
    "verbosity"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */
  const char* argument_list[] = {"ms2 file", "protein-database"};
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  initialize(argument_list, num_arguments,
                 option_list, num_options, argc, argv);

  // Get input: ms2 file 
  const char* ms2_file = get_string_parameter_pointer("ms2 file");

  // open ms2 file
  Crux::SpectrumCollection* spectra = SpectrumCollectionFactory::create(ms2_file);

  // parse the ms2 file for spectra
  carp(CARP_INFO, "Reading in ms2 file %s", ms2_file);
  if(!spectra->parse()){
    carp(CARP_FATAL, "Failed to parse ms2 file: %s", ms2_file);
  }
  
  carp(CARP_DEBUG, "There were %i spectra found in the ms2 file",
       spectra->getNumSpectra());

  /* Get input: protein file */
  char* input_file = get_string_parameter("protein-database");

  /* Prepare input, fasta or index */
  Index* index = NULL;
  Database* database = NULL;
  int num_proteins = prepare_protein_input(input_file, &index, &database); 
  free(input_file);

  carp(CARP_DEBUG, "Found %i proteins", num_proteins);
  if( num_proteins == 0 ){
    carp(CARP_FATAL, "No proteins were found in the protein source.");
  }
  
  /* Prepare output files */

  bool combine_target_decoy = get_boolean_parameter("tdc");
  OutputFiles output_files(this); 
  output_files.writeHeaders(num_proteins, combine_target_decoy);
  //output_files.writeHeader();
  // TODO (BF oct-21-09): consider adding pvalue file to OutputFiles
  FILE* decoy_pvalue_file = NULL;
  if( get_boolean_parameter("decoy-p-values") ){
    carp(CARP_DEBUG, "Opening decoy p-value file.");
    char* decoy_pvalue_filename 
      = get_string_parameter("search-decoy-pvalue-file");
    prefix_fileroot_to_name(&decoy_pvalue_filename);
    char* output_directory = get_string_parameter("output-dir");
    decoy_pvalue_file = create_file_in_path(decoy_pvalue_filename, 
                                            output_directory, 
                                            get_boolean_parameter("overwrite"));
    free(decoy_pvalue_filename);
    free(output_directory);
  }

  /* Perform search: loop over spectra*/

  // create spectrum iterator
  FilteredSpectrumChargeIterator* spectrum_iterator =
    new FilteredSpectrumChargeIterator(spectra);

  // get search parameters for match_collection
  bool compute_sp = get_boolean_parameter("compute-sp");
  if (get_boolean_parameter("sqt-output") && !compute_sp){
    compute_sp = true;
    carp(CARP_INFO, "Enabling parameter compute-sp since SQT output is enabled "
                    " (this will increase runtime).");
  }
  bool compute_pvalues = get_boolean_parameter("compute-p-values");
  int num_decoy_files = get_int_parameter("num-decoy-files");

  // For remembering and reporting number of searches
  SearchProgress progress(spectra->getNumChargedSpectra());

  // get list of mods
  PEPTIDE_MOD_T** peptide_mods = NULL;
  int num_peptide_mods = generate_peptide_mod_list( &peptide_mods );

  // do we need decoys
  bool have_index = (index != NULL);
  int num_decoy_collections = get_num_decoys(have_index);

  // for each spectrum
  while(spectrum_iterator->hasNext()) {
    SpectrumZState zstate;
    Spectrum* spectrum = spectrum_iterator->next(zstate);
    int charge = zstate.getCharge();
    
    bool is_decoy = false;

    progress.report(spectrum->getFirstScan(), charge);

    // with the target database decide how many peptide mods to use
    MatchCollection* target_psms = new MatchCollection(is_decoy); 
    int max_pep_mods = searchPepMods( target_psms, 
                                        is_decoy,   
                                        index,       
                                        database, 
                                        spectrum, 
                                        zstate,
                                        peptide_mods, 
                                        num_peptide_mods,
                                        compute_sp,
                                        compute_pvalues); 
 
    // are there any matches?
    if( target_psms->getMatchTotal() == 0 ){
      // don't print and don't search decoys
      carp(CARP_WARNING, "No matches found for spectrum %i, charge %i",
           spectrum->getFirstScan(), charge);
      delete target_psms;
      progress.increment(false);
      continue; // next spectrum
    }
    
    // now search decoys with the same number of mods
    is_decoy = true;
    // create separate decoy match_collections 
    vector<MatchCollection*> decoy_collection_list;

    int decoy_idx = 0;
    for(decoy_idx = 0; decoy_idx < num_decoy_collections; decoy_idx++){

      MatchCollection* decoy_psms = new MatchCollection(is_decoy);
      decoy_psms->setTargetExperimentSize(target_psms->getExperimentSize());
      decoy_collection_list.push_back(decoy_psms);

      searchPepMods(decoy_psms, 
                      is_decoy, 
                      index, 
                      database, 
                      spectrum, 
                      zstate, 
                      peptide_mods, 
                      max_pep_mods,
                      compute_sp,
                      compute_pvalues);
    }

    // calculate p-values for each collection of PSMs separately
    // use targets to get Weibull parameters, use same params for decoys
    if( compute_pvalues == true ){

      carp(CARP_DEBUG, "Estimating Weibull parameters.");
      while( ! target_psms->hasEnoughWeibullPoints() ){
        // generate more scores from new decoys if there are not enough
        addDecoyScores(target_psms, spectrum, zstate, index, 
                         database, peptide_mods, max_pep_mods);
        
      }
      target_psms->estimateWeibullParametersFromXcorrs(spectrum,
                                              charge);
      target_psms->computePValues(NULL);

      // use same params for each decoy set
      int decoy_idx = 0;
      for(decoy_idx = 0; decoy_idx < num_decoy_collections; decoy_idx++){
        MatchCollection* cur_collection = decoy_collection_list[decoy_idx];

        MatchCollection::transferWeibull(target_psms, cur_collection);

        carp(CARP_DEBUG, "Calculating p-values.");
        cur_collection->computePValues(decoy_pvalue_file);
      
      }// next collection
    }

    printSpectrumMatches(output_files, 
                           target_psms, 
                           decoy_collection_list,
                           spectrum, 
                           combine_target_decoy, 
                           num_decoy_files);

    progress.increment(true);

    // clean up
    delete target_psms;
    for(decoy_idx = 0; decoy_idx < num_decoy_collections; decoy_idx++){
      delete decoy_collection_list[decoy_idx];
    }

  }// next spectrum
  output_files.writeFooters();

  // clean up
  delete spectrum_iterator; 
  delete spectra;
  for(int mod_idx = 0; mod_idx < num_peptide_mods; mod_idx++){
    free_peptide_mod(peptide_mods[mod_idx]);
  }
  free(peptide_mods);
  Index::free(index);
  Database::freeDatabase(database);

  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6);
  carp(CARP_INFO, "Finished crux-search-for-matches");

  return 0;
}// end main


/**
 * \returns the command name for MatchSearch
 */
string MatchSearch::getName() {
  return "search-for-matches";
}

/**
 * \returns the description for MatchSearch
 */
string MatchSearch::getDescription() {
  return 
  "Search a collection of spectra against a sequence "
  "database, returning a collection of peptide-spectrum "
  "matches (PSMs) scored by XCorr.";

}

/**
 * \returns the file stem of the application, default getName.
 */
string MatchSearch::getFileStem() {
  return "search";
}

/**
 * \returns the enum of the application, default MISC_COMMAND
 */
COMMAND_T MatchSearch::getCommand() {
  return SEARCH_COMMAND;
}

/**
 * \returns whether the application needs the output directory or not. (default false).
 */
bool MatchSearch::needsOutputDirectory() {
  return true;
}

bool MatchSearch::hidden() {
  return true;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
