/**
 * \file match_search.c
 * BASED ON: original_match_search.c
 * DATE: Aug 19, 2008
 * AUTHOR: Barbara Frewen
 * DESCRIPTION: Main file for crux-search-for-matches.  Given an ms2
 * file and a fasta file or index, compare all spectra to peptides in
 * the fasta file/index and return high scoring matches.  Peptides are
 * determined by parameters for length, mass, mass tolerance, cleavages,
 * modifications. Score first by a preliminary method, keep only the
 * top ranking matches, score those with a second method and re-rank
 * by second score.  Output in binary csm file format or text sqt file
 * format. 
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
#include "carp.h"
#include "crux-utils.h"
#include "parameter.h"
#include "spectrum_collection.h"
#include "match_collection.h"
#include <errno.h>

#define NUM_SEARCH_OPTIONS 12
#define NUM_SEARCH_ARGS 2
#define PARAM_ESTIMATION_SAMPLE_COUNT 500

/* Private functions */
int prepare_protein_input(char* input_file, 
                          INDEX_T** index, 
                          DATABASE_T** database);
void open_output_files(char *output_directory,
                       BOOLEAN_T overwrite,
                       FILE*** binary_filehandle_array, 
                       FILE** sqt_filehandle,
                       FILE** decoy_sqt_filehandle,
                       FILE** tab_file,
                       FILE** decoy_tab_file);
BOOLEAN_T is_search_complete(MATCH_COLLECTION_T* matches, 
                             int mods_per_peptide);

//int main(int argc, char** argv){
int search_main(int argc, char** argv){

  /* Verbosity level for set-up/command line reading */
  set_verbosity_level(CARP_ERROR);

  /* Define optional command line arguments */
  int num_options = NUM_SEARCH_OPTIONS;
  char* option_list[NUM_SEARCH_OPTIONS] = {
    "verbosity",
    "version",
    "parameter-file",
    "overwrite",
    "compute-p-values",
    "spectrum-min-mass",
    "spectrum-max-mass",
    "spectrum-charge",
    "output-dir",
    "fileroot",
    "num-decoys-per-target",
    "decoy-location"
  };

  /* Define required command line arguments */
  int num_arguments = NUM_SEARCH_ARGS;
  char* argument_list[NUM_SEARCH_ARGS] = {"ms2 file", "protein input"};

  /* Initialize parameter.c and set default values*/
  initialize_parameters();

  /* Define optional and required arguments */
  select_cmd_line_options(option_list, num_options);
  select_cmd_line_arguments(argument_list, num_arguments);

  /* Parse the command line, including optional params file
     Includes syntax, type, and bounds checking, dies on error */
  parse_cmd_line_into_params_hash(argc, argv, "crux search-for-matches");

  /* Set seed for random number generation */
  if(strcmp(get_string_parameter_pointer("seed"), "time")== 0){
    time_t seconds; // use current time to seed
    time(&seconds); // Get value from sys clock and set seconds variable.
    srand((unsigned int) seconds); // Convert seconds to a unsigned int
  }
  else{
    srand((unsigned int)atoi(get_string_parameter_pointer("seed")));
  }
  
  /* Create output directory */ 
  char* output_folder = get_string_parameter("output-dir");
  BOOLEAN_T overwrite = get_boolean_parameter("overwrite");
  int result = create_output_directory(
    output_folder, 
    TRUE // Allow existing directory
  );
  if( result == -1 ){
    carp(CARP_FATAL, "Unable to create output directory %s.", output_folder);
  }

  /* Open the log file to record carp messages */
  char* log_file_name = get_string_parameter("search-log-file");
  open_log_file(&log_file_name);
  free(log_file_name);
  log_command_line(argc, argv);

  carp(CARP_INFO, "Beginning crux search-for-matches");

  // Write the parameter file
  char* param_file_name = get_string_parameter("search-param-file");
  print_parameter_file(&param_file_name);
  free(param_file_name);

  /* Get input: ms2 file */
  char* ms2_file = get_string_parameter_pointer("ms2 file");

  // open ms2 file
  SPECTRUM_COLLECTION_T* spectra = new_spectrum_collection(ms2_file);

  // parse the ms2 file for spectra
  carp(CARP_INFO, "Reading in ms2 file %s", ms2_file);
  if(!parse_spectrum_collection(spectra)){
    carp(CARP_FATAL, "Failed to parse ms2 file: %s", ms2_file);
  }
  
  carp(CARP_DEBUG, "There were %i spectra found in the ms2 file",
       get_spectrum_collection_num_spectra(spectra));

  /* Get input: protein file */
  char* input_file = get_string_parameter("protein input");

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

  FILE** psm_file_array = NULL; //file handle array
  FILE* sqt_file = NULL;
  FILE* decoy_sqt_file  = NULL;
  FILE* tab_file = NULL;
  FILE* decoy_tab_file  = NULL;

  open_output_files(
    output_folder,
    overwrite,
    &psm_file_array, 
    &sqt_file, 
    &decoy_sqt_file, 
    &tab_file,
    &decoy_tab_file
  );
  free(output_folder);

  //print headers
  serialize_headers(psm_file_array);
  print_sqt_header(sqt_file, "target", num_proteins, FALSE);// !analyze-matches
  print_sqt_header(decoy_sqt_file, "decoy", num_proteins, FALSE);
  print_tab_header(tab_file);
  print_tab_header(decoy_tab_file);
  /* Perform search: loop over spectra*/

  // create spectrum iterator
  FILTERED_SPECTRUM_CHARGE_ITERATOR_T* spectrum_iterator = 
    new_filtered_spectrum_charge_iterator(spectra);

  // get search parameters for match_collection
  BOOLEAN_T compute_pvalues = get_boolean_parameter("compute-p-values");
  BOOLEAN_T combine_target_decoy = get_boolean_parameter("tdc");

  // The total number of searches attempted.
  // This is the value that gets reported to the user on stderr.
  int spectrum_searches_counter = 0; 

  // The number of searches that found at least one candidate.
  // This is the value that goes into the .csm header.
  int num_successful_searches = 0;

  // flags and counters for loop
  int mod_idx = 0;
  int num_decoys = get_int_parameter("number-decoy-set");
  int progress_increment = get_int_parameter("print-search-progress");
  if( progress_increment == 0 ){
    progress_increment = BILLION;
  }

  // get list of mods
  PEPTIDE_MOD_T** peptide_mods = NULL;
  int num_peptide_mods = generate_peptide_mod_list( &peptide_mods );

  // for each spectrum
  while(filtered_spectrum_charge_iterator_has_next(spectrum_iterator)){
    int charge = 0;
    SPECTRUM_T* spectrum = 
      filtered_spectrum_charge_iterator_next(spectrum_iterator, &charge);
    double mass = get_spectrum_neutral_mass(spectrum, charge);

    if( ((spectrum_searches_counter+1) % progress_increment) == 0 ){
      carp(CARP_INFO, 
           "Searching spectrum number %i, charge %i, search number %i",
           get_spectrum_first_scan(spectrum), charge,
           spectrum_searches_counter+1 );
    }

    // with just the target database decide how many peptide mods to use
    // create an empty match collection 
    MATCH_COLLECTION_T* match_collection = 
      new_empty_match_collection( FALSE ); // is decoy = false

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

      // get peptide iterator
      MODIFIED_PEPTIDES_ITERATOR_T* peptide_iterator = 
        new_modified_peptides_iterator_from_mass(mass,
                                                 peptide_mod,
                                                 index,
                                                 database);
      // score peptides
      int added = add_matches(match_collection, 
                              spectrum, 
                              charge, 
                              peptide_iterator,
                              FALSE // is decoy
                              );

      carp(CARP_DEBUG, "Added %i matches", added);

      free_modified_peptides_iterator(peptide_iterator);

    }//next peptide mod

    // in case we searched all mods, do we need to assess again?

    // are there any matches?
    if( get_match_collection_match_total(match_collection) == 0 ){
      // don't print and don't search decoys
      carp(CARP_WARNING, "No matches found for spectrum %i, charge %i",
           get_spectrum_first_scan(spectrum), charge);
      free_match_collection(match_collection);
      spectrum_searches_counter++;
      continue; // next spectrum
    }
    
    // calculate p-values
    if( compute_pvalues == TRUE ){
      carp(CARP_DEBUG, "Estimating Weibull parameters for target.");
      if( estimate_weibull_parameters_from_xcorrs(match_collection,
                                                          spectrum,
                                                          charge) ){
        carp(CARP_DEBUG, "Calculating p-values.");
        compute_p_values(match_collection);
      }else{
        set_p_values_as_unscored(match_collection);
      }
    }

    if( combine_target_decoy == FALSE ){
      // print matches
      carp(CARP_DEBUG, "About to print target matches");
      print_matches(
                    match_collection, 
                    spectrum, 
                    FALSE,// is decoy
                    psm_file_array[0], 
                    sqt_file, 
                    decoy_sqt_file, 
                    tab_file,
                    decoy_tab_file
                    );
      
      //does this free all the matches, all the spectra and all the peptides?
      free_match_collection(match_collection);
      match_collection = NULL;
    }

    // now score same number of mods for decoys
    int max_mods = mod_idx;

    // for num_decoys  and num_decoy_repeats
    int decoy_idx = 0;
    int repeat_idx = 0;
    int num_decoy_repeats = get_int_parameter("num-decoys-per-target");
    // for tdc, the number of decoy files is 0 but we want to do at least one decoy search
    num_decoys = (combine_target_decoy) ? 1 : num_decoys;

    carp(CARP_DEBUG, "num_decoys (decoy set) %i, num_decoy_repeats (decoy per target) %i, max mods %i", num_decoys, num_decoy_repeats, max_mods);

    for(decoy_idx = 0; decoy_idx < num_decoys; decoy_idx++ ){
      carp(CARP_DETAILED_DEBUG, "Searching decoy %i", decoy_idx+1);

      if( combine_target_decoy == FALSE ){// create an empty match collection 
        match_collection = new_empty_match_collection( TRUE ); // is decoy
      }// else add to match_collection from above

      for(mod_idx = 0; mod_idx < max_mods; mod_idx++){
          
        // get peptide mod
        PEPTIDE_MOD_T* peptide_mod = peptide_mods[mod_idx];
        
        // for multiple decoy searches written to one file, repeat here
        for(repeat_idx=0; repeat_idx < num_decoy_repeats; repeat_idx++){
          // get peptide iterator
          MODIFIED_PEPTIDES_ITERATOR_T* peptide_iterator = 
            new_modified_peptides_iterator_from_mass(mass,
                                                     peptide_mod,
                                                     index,
                                                     database);
          // score peptides
          int added = add_matches(match_collection, 
                                  spectrum, 
                                  charge, 
                                  peptide_iterator,
                                  //0, // no sampling for param estimation
                                  TRUE);// is decoy
          carp(CARP_DEBUG, "Added %i matches", added);
          
          free_modified_peptides_iterator(peptide_iterator);
        }// next repeat
        
      }// last mod

      // calculate p-values
      if( compute_pvalues == TRUE ){
        carp(CARP_DEBUG, "Estimating Weibull parameters for decoys.");
        if( estimate_weibull_parameters_from_xcorrs(match_collection, 
                                                    spectrum, 
                                                    charge) ){
          carp(CARP_DEBUG, "Calculating p-values for decoys.");
          compute_p_values(match_collection);

        }else{ // there were too few scores to do estimation
          set_p_values_as_unscored(match_collection);
        }
      }

      // print matches
      // only print first decoy to sqt
      FILE* tmp_decoy_sqt_file = decoy_sqt_file;
      FILE* tmp_decoy_tab_file = decoy_tab_file;
      if( decoy_idx > 0 ){ 
        tmp_decoy_sqt_file = NULL; 
        tmp_decoy_tab_file = NULL; 
      }
      
      if( combine_target_decoy == FALSE ){ // print to decoy file
        carp(CARP_DEBUG, "About to print decoy matches");
        print_matches(
                      match_collection, 
                      spectrum, 
                      TRUE,// is decoy
                      psm_file_array[1+decoy_idx], 
                      sqt_file, 
                      tmp_decoy_sqt_file,
                      tab_file, 
                      tmp_decoy_tab_file
                      );
        
        free_match_collection(match_collection);
      }else{ // print all to target file
        carp(CARP_DEBUG, "About to print target and decoy matches");
        print_matches(
                      match_collection, 
                      spectrum, 
                      FALSE,// is decoy
                      psm_file_array[0], 
                      sqt_file, 
                      decoy_sqt_file, 
                      tab_file,
                      decoy_tab_file
                      );
        // free_match_collection?
      }

    }// last decoy

    spectrum_searches_counter++;
    num_successful_searches++;

    // clean up
    //free_match_collection(match_collection); ?? why commented out??

  }// next spectrum

  // finished searching!

  // fix headers in csm files
  int file_idx;
  for(file_idx=0; file_idx < num_decoys + 1; file_idx++){
    carp(CARP_DEBUG, "Changing csm header to have %i spectrum searches",
         num_successful_searches);
    serialize_total_number_of_spectra(num_successful_searches,
                                      psm_file_array[file_idx]);
  }
  // clean up memory

  carp(CARP_INFO, "Finished crux-search-for-matches");
  exit(0);
}// end main




/* Private function definitions */
/**
 * \brief Open either the index or fasta file and prepare it for
 * searching.  Die if the input file cannot be found or read.
 * \returns the number of proteins in the file/index
 */
int prepare_protein_input(char* input_file, 
                          INDEX_T** index, 
                          DATABASE_T** database){

  int num_proteins = 0;
  BOOLEAN_T use_index = is_directory(input_file);

  if (use_index == TRUE){
    carp(CARP_INFO, "Preparing protein index %s", input_file);
    *index = new_index_from_disk(input_file);

    if (index == NULL){
      carp(CARP_FATAL, "Could not create index from disk for %s", input_file);
    }
    num_proteins = get_index_num_proteins(*index);

  } else {
    carp(CARP_INFO, "Preparing protein fasta file %s", input_file);
    *database = new_database(input_file, FALSE);         
    if( database == NULL ){
      carp(CARP_FATAL, "Could not create protein database");
    } 

    if(!parse_database(*database)){
      carp(CARP_FATAL, "Error with protein input");
    } 
    num_proteins = get_database_num_proteins(*database);
  }
  return num_proteins;
}

/**
 * \brief A private function for crux-search-for-matches to prepare
 * binary psm, tab-delimited text, and sqt files.
 *
 * Opens psm file(s) if requested, setting a given
 * pointer to the array of filehandles.  Opens sqt file(s) if
 * requested, setting the given pointers to each file handle.  If
 * binary files not requested, creates an array of NULL pointers.  If
 * sqt files not requested, sets given pointers to NULL. 
 *
 * \returns void.  Sets given arguments to newly created filehandles.
 */
void open_output_files(
  char *output_directory, ///< name of output directory -in
  BOOLEAN_T overwrite,     ///< overwrite existing files -in
  FILE*** psm_file_array, ///< put binary psm filehandles here -out
  FILE** sqt_file,        ///< put text sqt filehandle here -out
  FILE** decoy_sqt_file,  ///< put decoy sqt filehandle here -out
  FILE** tab_file,        ///< put text sqt filehandle here -out
  FILE** decoy_tab_file)  ///< put decoy sqt filehandle here -out
{
  // create binary psm files (allocate memory, even if not used)
  *psm_file_array = create_psm_files();

  //create sqt file handles
  carp(CARP_DEBUG, "Opening sqt files");
  char* sqt_filename = get_string_parameter("search-sqt-output-file");
  prefix_fileroot_to_name(&sqt_filename);
  *sqt_file = create_file_in_path(sqt_filename, 
                                  output_directory, 
                                  overwrite);
  free(sqt_filename);
  if( get_int_parameter("number-decoy-set") > 0 ){
    char* decoy_sqt_filename = get_string_parameter("decoy-sqt-output-file");
    prefix_fileroot_to_name(&decoy_sqt_filename);
    *decoy_sqt_file = create_file_in_path(decoy_sqt_filename,
                                          output_directory,
                                          overwrite);
    free(decoy_sqt_filename);
  }

  //create tab-delimited file handles
  carp(CARP_DEBUG, "Opening tab delimited files");
  char* tab_filename = get_string_parameter("search-tab-output-file");
  prefix_fileroot_to_name(&tab_filename);
  *tab_file = create_file_in_path(tab_filename, 
                                  output_directory, 
                                  overwrite);
  free(tab_filename);
  char* decoy_tab_filename = get_string_parameter("decoy-tab-output-file");
  prefix_fileroot_to_name(&decoy_tab_filename);
  if( get_int_parameter("number-decoy-set") > 0 ){
    *decoy_tab_file = create_file_in_path(decoy_tab_filename,
                                          output_directory,
                                          overwrite);
    free(decoy_tab_filename);
  }

  carp(CARP_DEBUG, "Finished opening output files");
}

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
