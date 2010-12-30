/**
 * \file match_search.cpp
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
#include "match_collection.h"
#include "carp.h"
#include "crux-utils.h"
#include "parameter.h"
#include "SpectrumCollection.h"
#include "FilteredSpectrumChargeIterator.h"
#include <errno.h>
#include "OutputFiles.h"
#include "SearchProgress.h"

/* Private functions */
int search_pep_mods(
  MATCH_COLLECTION_T* match_collection, ///< store PSMs here
  BOOLEAN_T is_decoy,   ///< generate decoy peptides from index/db
  INDEX_T* index,       ///< index to use for generating peptides
  DATABASE_T* database, ///< db to use for generating peptides
  Spectrum* spectrum, ///< spectrum to search
  int charge,           ///< seach spectrum at this charge state
  PEPTIDE_MOD_T** pep_mod_list, ///< list of peptide mods to apply
  int num_peptide_mods, ///< how many p_mods to use from the list
  BOOLEAN_T store_scores///< keep all scores for p-value estimation
);
void add_decoy_scores(
  MATCH_COLLECTION_T* target_psms, ///< add scores to these matches
  Spectrum* spectrum, ///<
  int charge, ///< 
  INDEX_T* index, ///< search this index if not null
  DATABASE_T* database, ///< search this database if not null
  PEPTIDE_MOD_T** peptitde_mods, ///< list of peptide mods to search
  int num_peptide_mods ///< number of mods in the above array
);
BOOLEAN_T is_search_complete(MATCH_COLLECTION_T* matches, 
                             int mods_per_peptide);
void print_spectrum_matches(
  OutputFiles& output_files,       
  MATCH_COLLECTION_T* target_psms, 
  vector<MATCH_COLLECTION_T*>& decoy_psms,
  Spectrum* spectrum,             
  BOOLEAN_T combine_target_decoy,
  int num_decoy_files
                   );

#ifdef SEARCH_ENABLED // Discard this code in open source release
int search_main(int argc, char** argv){

  /* Define optional command line arguments */
  const char* option_list[] = {
    "fileroot",
    "output-dir",
    "overwrite",
    "num-decoys-per-target",
    "decoy-location",
    "compute-sp",
    "compute-p-values",
    "spectrum-min-mass",
    "spectrum-max-mass",
    "spectrum-charge",
    "scan-number",
    "xcorr-var-bin",
    "mz-bin-width",
    "mz-bin-offset",
    "parameter-file",
    "verbosity"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */
  const char* argument_list[] = {"ms2 file", "protein database"};
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  initialize_run(SEARCH_COMMAND, argument_list, num_arguments,
                 option_list, num_options, argc, argv);

  // Get input: ms2 file 
  const char* ms2_file = get_string_parameter_pointer("ms2 file");

  // open ms2 file
  SpectrumCollection* spectra = new SpectrumCollection(ms2_file);

  // parse the ms2 file for spectra
  carp(CARP_INFO, "Reading in ms2 file %s", ms2_file);
  if(!spectra->parse()){
    carp(CARP_FATAL, "Failed to parse ms2 file: %s", ms2_file);
  }
  
  carp(CARP_DEBUG, "There were %i spectra found in the ms2 file",
       spectra->getNumSpectra());

  /* Get input: protein file */
  char* input_file = get_string_parameter("protein database");

  /* Prepare input, fasta or index */
  INDEX_T* index = NULL;
  DATABASE_T* database = NULL;
  int num_proteins = prepare_protein_input(input_file, &index, &database); 
  free(input_file);

  carp(CARP_DEBUG, "Found %i proteins", num_proteins);
  if( num_proteins == 0 ){
    carp(CARP_FATAL, "No proteins were found in the protein source.");
  }
  
  /* Prepare output files */
  OutputFiles output_files(SEARCH_COMMAND); 
  output_files.writeHeaders(num_proteins);
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
  BOOLEAN_T compute_pvalues = get_boolean_parameter("compute-p-values");
  BOOLEAN_T combine_target_decoy = get_boolean_parameter("tdc");
  int num_decoy_files = get_int_parameter("num-decoy-files");

  // For remembering and reporting number of searches
  SearchProgress progress;

  // get list of mods
  PEPTIDE_MOD_T** peptide_mods = NULL;
  int num_peptide_mods = generate_peptide_mod_list( &peptide_mods );

  // for each spectrum
  while(spectrum_iterator->hasNext()) {
    int charge = 0;
    Spectrum* spectrum = spectrum_iterator->next(&charge);
    BOOLEAN_T is_decoy = FALSE;

    progress.report(spectrum->getFirstScan(), charge);

    // with the target database decide how many peptide mods to use
    MATCH_COLLECTION_T* target_psms = new_empty_match_collection(is_decoy); 
    int max_pep_mods = search_pep_mods( target_psms, 
                                        is_decoy,   
                                        index,       
                                        database, 
                                        spectrum, 
                                        charge,
                                        peptide_mods, 
                                        num_peptide_mods,
                                        compute_pvalues); 
 
    // are there any matches?
    if( get_match_collection_match_total(target_psms) == 0 ){
      // don't print and don't search decoys
      carp(CARP_WARNING, "No matches found for spectrum %i, charge %i",
           spectrum->getFirstScan(), charge);
      free_match_collection(target_psms);
      progress.increment(FALSE);
      continue; // next spectrum
    }
    
    // now search decoys with the same number of mods
    is_decoy = TRUE;
    // create separate decoy match_collections 
    int num_decoy_collections = get_int_parameter("num-decoys-per-target"); 
    vector<MATCH_COLLECTION_T*> decoy_collection_list;

    int decoy_idx = 0;
    for(decoy_idx = 0; decoy_idx < num_decoy_collections; decoy_idx++){

      MATCH_COLLECTION_T* decoy_psms = new_empty_match_collection(is_decoy);
      decoy_collection_list.push_back(decoy_psms);

      search_pep_mods(decoy_psms, 
                      is_decoy, 
                      index, 
                      database, 
                      spectrum, 
                      charge, 
                      peptide_mods, 
                      max_pep_mods,
                      compute_pvalues);
    }

    // calculate p-values for each collection of PSMs separately
    // use targets to get Weibull parameters, use same params for decoys
    if( compute_pvalues == TRUE ){

      carp(CARP_DEBUG, "Estimating Weibull parameters.");
      while( ! has_enough_weibull_points(target_psms) ){
        // generate more scores from new decoys if there are not enough
        add_decoy_scores(target_psms, spectrum, charge, index, 
                         database, peptide_mods, max_pep_mods);
        
      }
      estimate_weibull_parameters_from_xcorrs(target_psms,
                                              spectrum,
                                              charge);
      compute_p_values(target_psms, NULL);

      // use same params for each decoy set
      int decoy_idx = 0;
      for(decoy_idx = 0; decoy_idx < num_decoy_collections; decoy_idx++){
        MATCH_COLLECTION_T* cur_collection = decoy_collection_list[decoy_idx];

        transfer_match_collection_weibull(target_psms, cur_collection);

        carp(CARP_DEBUG, "Calculating p-values.");
        compute_p_values(cur_collection, decoy_pvalue_file);
      
      }// next collection
    }

    print_spectrum_matches(output_files, 
                           target_psms, 
                           decoy_collection_list,
                           spectrum, 
                           combine_target_decoy, 
                           num_decoy_files);

    progress.increment(TRUE);

    // clean up
    free_match_collection(target_psms);
    for(decoy_idx = 0; decoy_idx < num_decoy_collections; decoy_idx++){
      free_match_collection(decoy_collection_list[decoy_idx]);
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
  free_index(index);
  free_database(database);

  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6);
  carp(CARP_INFO, "Finished crux-search-for-matches");

  return 0;
}// end main

#else // SEARCH_ENABLED not defined
int search_main(int argc, char **argv){
  (void) argc;
  (void) argv;
  fputs(
    "You are using the open source version of Crux. Due to intellectual\n"
    "property issues, we are unable to provide database search functionality\n"
    "in this version. To obtain a licence for the full functional version of\n"
    "Crux that includes the database search tools, please visit the following URL:\n"
    "\nhttp://depts.washington.edu/ventures/UW_Technology/Express_Licenses/crux.php\n",
    stderr
  );
  return 1;
}
#endif // SEARCH_ENABLED



/* Private function definitions */

/**
 * \brief Look at matches and search parameters to determine if a
 * sufficient number PSMs have been found.  Returns TRUE if the
 * maximum number of modifications per peptide have been considered.
 * In the future, implement and option and test for a minimum score.
 * \returns TRUE if no more PSMs need be searched.
 */
BOOLEAN_T is_search_complete(MATCH_COLLECTION_T* matches, 
                             int mods_per_peptide){


  if( matches == NULL ){
    return FALSE;
  }

  // keep searching if no limits on how many mods per peptide
  if( get_int_parameter("max-mods") == MAX_PEPTIDE_LENGTH ){
    return FALSE;
  }
  // stop searching if at max mods per peptide
  if( mods_per_peptide == get_int_parameter("max-mods") ){ 
    return TRUE;
  }

  // test for minimun score found

  return FALSE;
  
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
int search_pep_mods(
  MATCH_COLLECTION_T* match_collection, ///< store PSMs here
  BOOLEAN_T is_decoy,   ///< generate decoy peptides from index/db
  INDEX_T* index,       ///< index to use for generating peptides
  DATABASE_T* database, ///< db to use for generating peptides
  Spectrum* spectrum, ///< spectrum to search
  int charge,           ///< seach spectrum at this charge state
  PEPTIDE_MOD_T** peptide_mods, ///< list of peptide mods to apply
  int num_peptide_mods, ///< how many p_mods to use from the list
  BOOLEAN_T store_scores///< save all scores for p-value estimation
){

  // set match_collection charge
  set_match_collection_charge(match_collection, charge);

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
      BOOLEAN_T passes = is_search_complete(match_collection, cur_aa_mods);
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
    MODIFIED_PEPTIDES_ITERATOR_T* peptide_iterator =
      new_modified_peptides_iterator_from_mz(spectrum->getPrecursorMz(),
                                             charge,
                                             peptide_mod, 
                                             is_decoy,
                                             index,
                                             database);
    
    
    // score peptides
    int added = add_matches(match_collection, 
                            spectrum, 
                            charge, 
                            peptide_iterator,
                            is_decoy,
                            store_scores,
                            get_boolean_parameter("compute-sp"),
                            FALSE // don't filtery by Sp
                            );
    
    carp(CARP_DEBUG, "Added %i matches", added);

    free_modified_peptides_iterator(peptide_iterator);
    
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
void print_spectrum_matches(
  OutputFiles& output_files,       
  MATCH_COLLECTION_T* target_psms, 
  vector<MATCH_COLLECTION_T*>& decoy_psms,
  Spectrum* spectrum,             
  BOOLEAN_T combine_target_decoy,
  int num_decoy_files
                   ){

  // now print matches to one, two or several files
  if( combine_target_decoy == TRUE ){
    // merge all collections
    MATCH_COLLECTION_T* all_psms = target_psms;
    for(size_t decoy_idx = 0; decoy_idx < decoy_psms.size(); decoy_idx++){
      merge_match_collections(decoy_psms[decoy_idx], all_psms);
    }
    
    // sort and rank
    if( get_match_collection_scored_type(all_psms, SP) == TRUE ){
      populate_match_rank_match_collection(all_psms, SP);
    }
    populate_match_rank_match_collection(all_psms, XCORR);

    vector<MATCH_COLLECTION_T*> empty_list;
    output_files.writeMatches(all_psms, // target matches
                              empty_list, // decoy matches
                              XCORR, spectrum); 
    
  }else{ // targets and decoys in separate files
    
    // if decoys in one file
    if( num_decoy_files == 1 ){
      // merge decoys
      MATCH_COLLECTION_T* merged_decoy_psms = decoy_psms[0];
      for(size_t decoy_idx = 1; decoy_idx < decoy_psms.size(); decoy_idx++){
        merge_match_collections(decoy_psms[decoy_idx],
                                merged_decoy_psms);
      }
      
      // sort and rank
      // NOTE (BF 09-14-10): since the multiple decoy collections have already
      // been truncated, the merged ranks aren't accurate for the total space
      // of decoys searched
      populate_match_rank_match_collection(merged_decoy_psms, XCORR);
      
      vector<MATCH_COLLECTION_T*> decoy_list(1, merged_decoy_psms);
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
void add_decoy_scores(
  MATCH_COLLECTION_T* target_psms, ///< add scores to these matches
  Spectrum* spectrum, ///<
  int charge, ///< 
  INDEX_T* index, ///< search this index if not null
  DATABASE_T* database, ///< search this database if not null
  PEPTIDE_MOD_T** peptide_mods, ///< list of peptide mods to search
  int num_peptide_mods ///< number of mods in the above array
){

  int mod_idx = 0;
  // for each peptide mod in the list
  for(mod_idx = 0; mod_idx < num_peptide_mods; mod_idx++){
    MODIFIED_PEPTIDES_ITERATOR_T* peptide_iterator = 
      new_modified_peptides_iterator_from_mz(
                                          spectrum->getPrecursorMz(),
                                          charge,
                                          peptide_mods[mod_idx],
                                          TRUE, // is decoy
                                          index,
                                          database);
    add_decoy_scores_match_collection(target_psms, 
                                      spectrum, 
                                      charge, 
                                      peptide_iterator);  

    free_modified_peptides_iterator(peptide_iterator);
  }


}

