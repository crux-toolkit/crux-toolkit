/**
 * \file sequest-search.cpp
 * AUTHOR: Barbara Frewen
 * CREATE DATE: Oct 2, 2009
 * PROJECT: crux
 * \brief The crux search routine that emulates SEQUEST.
 *
 * Scores all candidate peptides with Sp, deletes all but the 500
 * top-scoring candidates, scores remaining 500 with xcorr, sorts
 * results by xcorr and returns the top 5 plus the match with the best
 * Sp score.  Writes results to .sqt, .txt, and .csm files.  Does not
 * compute p-values.
 */

#include "sequest-search.h"

// Private functions, commented below at definition
void print_matches(
  OutputFiles& output_files,       
  MATCH_COLLECTION_T* target_psms, 
  vector<MATCH_COLLECTION_T*>& decoy_psms,
  Spectrum* spectrum,             
  BOOLEAN_T combine_target_decoy,
  int num_decoy_files
                   );

/**
 * \brief The starting point for the crux sequest-search command.
 *
 * After parsing command line and opening input and output files,
 * iterates over every spectrum at each charge state.  Generates
 * matches for each and prints.  Matches to decoy spectra can also be
 * generated.  If num-decoys-per-target is n, there are n decoy match
 * collections which may be merged for printing.  
 * If the user gave the command line
 *
 * crux sequest-search [options] ms2file proteins
 * 
 * then this is passed everything after the 'crux' token
 */
#ifdef SEARCH_ENABLED // Discard this code in open source release
int sequest_search_main(int argc,   ///< number of cmd line tokens
                        char** argv)///< array of cmd line tokens
{
  const char* option_list[] = {
    "verbosity",
    "parameter-file",
    "overwrite",
    "spectrum-min-mass",
    "spectrum-max-mass",
    "spectrum-charge",
    "output-dir",
    "scan-number",
    "fileroot",
    "num-decoys-per-target",
    "decoy-location"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  // Define required command line arguments 
  const char* argument_list[] = {"ms2 file", "protein database"};
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  initialize_run(SEQUEST_COMMAND, argument_list, num_arguments,
                 option_list, num_options, argc, argv);

  // Get input: protein file
  char* input_file = get_string_parameter("protein database");

  // Prepare input, fasta or index 
  INDEX_T* index = NULL;
  DATABASE_T* database = NULL;
  int num_proteins = prepare_protein_input(input_file, &index, &database); 
  free(input_file);

  carp(CARP_DEBUG, "Found %i proteins", num_proteins);
  if( num_proteins == 0 ){
    carp(CARP_FATAL, "No proteins were found in the protein source.");
  }
  
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

  // Prepare output files 
  OutputFiles output_files(SEQUEST_COMMAND); 
  output_files.writeHeaders(num_proteins);

  // get search parameters for match_collection
  BOOLEAN_T combine_target_decoy = get_boolean_parameter("tdc");
  int num_decoy_files = get_int_parameter("num-decoy-files");
  int num_decoys_per_target = get_int_parameter("num-decoys-per-target");

  SearchProgress progress;

  // get list of mods
  PEPTIDE_MOD_T** peptide_mods = NULL;
  int num_peptide_mods = generate_peptide_mod_list( &peptide_mods );

  // create spectrum iterator
  FilteredSpectrumChargeIterator* spectrum_iterator = 
    new FilteredSpectrumChargeIterator(spectra);

  // Perform search on each spectrum
  while(spectrum_iterator->hasNext()){
    int charge = 0;
    Spectrum* spectrum = 
      spectrum_iterator->next(&charge);

    progress.report(spectrum->getFirstScan(), charge);

    // create empty match collections to store results in
    MATCH_COLLECTION_T* target_psms = new_empty_match_collection(FALSE); 
    set_match_collection_charge(target_psms, charge);

    vector<MATCH_COLLECTION_T*> decoy_psm_collections;
    for(int decoy_idx=0; decoy_idx < num_decoys_per_target; decoy_idx++){
      MATCH_COLLECTION_T* psms = new_empty_match_collection(TRUE);
      set_match_collection_charge(psms, charge);
      decoy_psm_collections.push_back(psms);
    }

    // search with one peptide iterator for each peptide mod
    for(int mod_idx = 0; mod_idx < num_peptide_mods; mod_idx++){

      // get peptide mod
      PEPTIDE_MOD_T* peptide_mod = peptide_mods[mod_idx];

      // get peptide iterator
      FLOAT_T precursor_mz = spectrum->getPrecursorMz();
      MODIFIED_PEPTIDES_ITERATOR_T* peptide_iterator =
        new_modified_peptides_iterator_from_mz(precursor_mz,
                                               charge,
                                               peptide_mod, 
                                               FALSE, // not decoy
                                               index,
                                               database);

      // add matches to targets
      int added = add_matches(target_psms,
                              spectrum,
                              charge,
                              peptide_iterator,
                              FALSE, // not decoy
                              FALSE, // don't save scores for p-values
                              TRUE,  // do preliminary Sp scoring
                              TRUE   // filter by Sp
                              ); 

      // add matches to each decoy
      for(int decoy_idx = 0; decoy_idx < num_decoys_per_target; decoy_idx++){

        // get new peptide iterator
        free_modified_peptides_iterator(peptide_iterator);
        peptide_iterator =
          new_modified_peptides_iterator_from_mz(precursor_mz,
                                                 charge,
                                                 peptide_mod, 
                                                 TRUE,  // is decoy
                                                 index,
                                                 database);
        // add matches
        MATCH_COLLECTION_T* cur_decoys = decoy_psm_collections.at(decoy_idx);
        add_matches(cur_decoys,
                    spectrum,
                    charge,
                    peptide_iterator,
                    TRUE,  // is decoy
                    FALSE, // don't save scores for p-values
                    TRUE,  // do preliminary Sp scoring
                    TRUE   // filter by Sp
                    ); 
      }

      carp(CARP_DEBUG, "Found %d peptides.", added);

      // clean up for next peptide mod
      free_modified_peptides_iterator(peptide_iterator);

    } // next peptide mod

    // print matches
    int total_matches = get_match_collection_match_total(target_psms);

    if( total_matches == 0 ){
      carp(CARP_WARNING, "No matches found for spectrum %i, charge %i.",
           spectrum->getFirstScan(), charge);
      progress.increment(FALSE);

    }else{  
      print_matches(output_files, target_psms, decoy_psm_collections,
                    spectrum, combine_target_decoy, num_decoy_files);
      progress.increment(TRUE);
    }

    // clean up
    free_match_collection(target_psms);
    vector<MATCH_COLLECTION_T*>::iterator it = decoy_psm_collections.begin();
    for(; it < decoy_psm_collections.end(); ++it){
      MATCH_COLLECTION_T* psms = *it;
      free_match_collection(psms);
    }

  } // next spectrum

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
  carp(CARP_INFO, "Finished crux sequest-search");
  
  return 0;
}// end main

#else // SEARCH_ENABLED not defined
int sequest_search_main(int argc, char **argv){
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
 * \brief Pring the target and decoy match collections to their
 * respective target and decoy files.
 *
 * Three possibilities: 1. combine the target and all decoy
 * collections and print to target file.  2. print targets to target
 * file and combine all decoys and print to one decoy file.  3. print
 * each collection to a separate file.
 * Possible side effectos: Collections may be merged and re-ranked.
 */
void print_matches(
  OutputFiles& output_files,       ///< files to print to
  MATCH_COLLECTION_T* target_psms, ///< target psms to print
  vector<MATCH_COLLECTION_T*>& decoy_psms,///< decoy psms to print
  Spectrum* spectrum,            ///< all matches are to this spec
  BOOLEAN_T combine_target_decoy,  ///< merge targets and decoys?
  int num_decoy_files              ///< merge decoys?
){ 

  if( combine_target_decoy ){

    // merge all collections
    MATCH_COLLECTION_T* all_psms = target_psms;
    for(size_t decoy_idx = 0; decoy_idx < decoy_psms.size(); decoy_idx++){
      merge_match_collections(decoy_psms.at(decoy_idx), all_psms);
    }

    // sort and rank
    populate_match_rank_match_collection(all_psms, SP);
    save_top_sp_match(all_psms);
    populate_match_rank_match_collection(all_psms, XCORR);
    vector<MATCH_COLLECTION_T*> empty_list;
    output_files.writeMatches(all_psms, // target matches
                              empty_list,     // decoy matches
                              XCORR,    // use XCORR rank for cutoff
                              spectrum); 


  }else{ // targets and decoys in separate files

    if( num_decoy_files == 1 ){ // combine decoys

      // merge decoys
      MATCH_COLLECTION_T* merged_decoy_psms = decoy_psms.at(0);
      for(size_t decoy_idx = 1; decoy_idx < decoy_psms.size(); decoy_idx++){
        merge_match_collections(decoy_psms.at(decoy_idx),
                                merged_decoy_psms);
      }
      
      // sort and rank
      populate_match_rank_match_collection(merged_decoy_psms, SP);
      save_top_sp_match(merged_decoy_psms);
      populate_match_rank_match_collection(merged_decoy_psms, XCORR);

      vector<MATCH_COLLECTION_T*> decoy_list(1, merged_decoy_psms);
      output_files.writeMatches(target_psms, decoy_list, 
                                XCORR, spectrum);
      
    }else{   // write targets and decoys to separate files
      // TODO write a version of OutputFiles::writeMatches that takes
      // a vector of decoy match collections
      int num_decoys = decoy_psms.size();
      vector<MATCH_COLLECTION_T*> decoy_psm_array;
      for(int decoy_idx = 0; decoy_idx < num_decoys; decoy_idx++){
        decoy_psm_array.push_back(decoy_psms.at(decoy_idx));
      }

      output_files.writeMatches(target_psms, decoy_psm_array, 
                                XCORR, spectrum);
    }
  }


}








/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
