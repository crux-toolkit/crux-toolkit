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
#include "parameter.h"
#include "spectrum_collection.h"
#include "match_collection.h"

#define NUM_SEARCH_OPTIONS 16
#define NUM_SEARCH_ARGS 2
#define PARAM_ESTIMATION_SAMPLE_COUNT 500

/* Private functions */
int prepare_protein_input(char* input_file, 
                          INDEX_T** index, 
                          DATABASE_T** database);
void open_output_files(FILE*** binary_filehandle_array, 
                       FILE** sqt_filehandle,
                       FILE** decoy_sqt_filehandle,
                       FILE** tab_file,
                       FILE** decoy_tab_file);
BOOLEAN_T is_search_complete(MATCH_COLLECTION_T* matches, 
                             int mods_per_peptide);

int search_main(int argc, char** argv){

  /* Verbosity level for set-up/command line reading */
  set_verbosity_level(CARP_ERROR);

  /* Define optional command line arguments */
  int num_options = NUM_SEARCH_OPTIONS;
  char* option_list[NUM_SEARCH_OPTIONS] = {
    "verbosity",
    "version",
    "parameter-file",
    "write-parameter-file",
    "overwrite",
    "use-index",
    "compute-p-values",
    "spectrum-min-mass",
    "spectrum-max-mass",
    "spectrum-charge",
    "match-output-folder",
    "output-mode",
    "sqt-output-file",
    "tab-output-file",
    "decoy-sqt-output-file",
    "number-decoy-set"
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
  
  carp(CARP_INFO, "Beginning crux search-for-matches");

  /* Get input: ms2 file */
  char* ms2_file = get_string_parameter_pointer("ms2 file");

  // open ms2 file
  SPECTRUM_COLLECTION_T* spectra = new_spectrum_collection(ms2_file);

  // parse the ms2 file for spectra
  carp(CARP_INFO, "Reading in ms2 file %s", ms2_file);
  if(!parse_spectrum_collection(spectra)){
    carp(CARP_FATAL, "Failed to parse ms2 file: %s", ms2_file);
    free_spectrum_collection(spectra);
    exit(1);
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
    exit(1);
  }
  
  /* Prepare output files */

  FILE** psm_file_array = NULL; //file handle array
  FILE* sqt_file = NULL;
  FILE* decoy_sqt_file  = NULL;
  FILE* tab_file = NULL;
  FILE* decoy_tab_file  = NULL;

  open_output_files(
    &psm_file_array, 
    &sqt_file, 
    &decoy_sqt_file, 
    &tab_file,
    &decoy_tab_file
  );

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
  //  int sample_count = (compute_pvalues) ? PARAM_ESTIMATION_SAMPLE_COUNT : 0;
  BOOLEAN_T combine_target_decoy = get_boolean_parameter("tdc");

  // flags and counters for loop
  int spectrum_searches_counter = 0; //for psm file header, sum(spec*charges)
  int mod_idx = 0;
  int num_decoys = get_int_parameter("number-decoy-set");
  int progress_increment = get_int_parameter("print-search-progress");
  if( progress_increment == 0 ){
    progress_increment = BILLION;
  }

  // get list of mods
  PEPTIDE_MOD_T** peptide_mods = NULL;
  int num_peptide_mods = generate_peptide_mod_list( &peptide_mods );

  // DELETE ME
  // for estimating params for p-values, randomly select a total of 
  //    sample_count matches, a constant fraction from each peptide mod
  //int sample_per_pep_mod =  sample_count / num_peptide_mods;
  //carp(CARP_DEBUG, "Got %d peptide mods, sample %i per", 
  //     num_peptide_mods, sample_per_pep_mod);

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
      continue; // next spectrum
    }
    
    // calculate p-values
    if( compute_pvalues == TRUE ){
      carp(CARP_DEBUG, "Estimating Weibull parameters.");
      //if( estimate_weibull_parameters_from_sample_matches(match_collection,
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
    // decoy files is 0 but we want to do at least one decoy searc for tdc
    num_decoys = (combine_target_decoy) ? 1 : num_decoys;

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
      }



    }// last decoy

    spectrum_searches_counter++;

    // clean up
    //    free_match_collection(match_collection);
  }// next spectrum

  // fix headers in csm files
  int file_idx;
  for(file_idx=0; file_idx < num_decoys + 1; file_idx++){
    carp(CARP_DEBUG, "Changing csm header to have %i spectrum searches",
         spectrum_searches_counter);
    serialize_total_number_of_spectra(spectrum_searches_counter,
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
  BOOLEAN_T use_index = get_boolean_parameter("use-index");

  if (use_index == TRUE){
    carp(CARP_INFO, "Preparing protein index %s", input_file);
    *index = new_index_from_disk(input_file);

    if (index == NULL){
      carp(CARP_FATAL, "Could not create index from disk for %s", input_file);
      exit(1);
    }
    num_proteins = get_index_num_proteins(*index);

  } else {
    carp(CARP_INFO, "Preparing protein fasta file %s", input_file);
    *database = new_database(input_file, FALSE);         
    if( database == NULL ){
      carp(CARP_FATAL, "Could not create protein database");
      exit(1);
    } 

    if(!parse_database(*database)){
      carp(CARP_FATAL, "Error with protein input");
      exit(1);
    } 
    num_proteins = get_database_num_proteins(*database);
  }
  return num_proteins;
}

/**
 * \brief A private function for crux-search-for-matches to prepare
 * binary psm and text sqt files.
 *
 * Reads the --overwrite and --output-mode values from
 * parameter.c. Opens psm file(s) if requested, setting a given
 * pointer to the array of filehandles.  Opens sqt file(s) if
 * requested, setting the given pointers to each file handle.  If
 * binary files not requested, creates an array of NULL pointers.  If
 * sqt files not requested, sets given pointers to NULL. 
 *
 * \returns void.  Sets given arguments to newly created filehandles.
 */
void open_output_files(
  FILE*** psm_file_array, ///< put binary psm filehandles here -out
  FILE** sqt_file,        ///< put text sqt filehandle here -out
  FILE** decoy_sqt_file,  ///< put decoy sqt filehandle here -out
  FILE** tab_file,        ///< put text sqt filehandle here -out
  FILE** decoy_tab_file)  ///< put decoy sqt filehandle here -out
{
  char* match_output_folder = get_string_parameter("match-output-folder");
  MATCH_SEARCH_OUTPUT_MODE_T output_type = get_output_type_parameter(
                                                    "output-mode");
  BOOLEAN_T overwrite = get_boolean_parameter("overwrite");
  carp(CARP_DEBUG, "The output type is %d (binary, sqt, tab, all)" \
       " and overwrite is '%d'", (int)output_type, (int)overwrite);


  // create binary psm files (allocate memory, even if not used)
  *psm_file_array = create_psm_files();

  if(output_type == SQT_OUTPUT || output_type == ALL_OUTPUT){

    //create sqt file handles
    carp(CARP_DEBUG, "Opening sqt files");
    char* sqt_filename = get_string_parameter_pointer("sqt-output-file");
    *sqt_file = create_file_in_path(sqt_filename, 
                                    match_output_folder, 
                                    overwrite);
    char* decoy_sqt_filename = get_string_parameter_pointer(
                                                    "decoy-sqt-output-file");
    if( get_int_parameter("number-decoy-set") > 0 ){
      *decoy_sqt_file = create_file_in_path(decoy_sqt_filename,
                                            match_output_folder,
                                            overwrite);
    }

    if(sqt_file == NULL || decoy_sqt_file == NULL){
      carp(CARP_DEBUG, "sqt file or decoy is null");
    }
  }

  if(output_type == TAB_OUTPUT || output_type == ALL_OUTPUT){

    //create sqt file handles
    carp(CARP_DEBUG, "Opening tab delimited files");
    char* tab_filename = get_string_parameter_pointer("tab-output-file");
    *tab_file = create_file_in_path(tab_filename, 
                                    match_output_folder, 
                                    overwrite);
    char* decoy_tab_filename = get_string_parameter_pointer(
                                                    "decoy-tab-output-file");
    if( get_int_parameter("number-decoy-set") > 0 ){
      *decoy_tab_file = create_file_in_path(decoy_tab_filename,
                                            match_output_folder,
                                            overwrite);
    }

    if(tab_file == NULL || decoy_tab_file == NULL){
      carp(CARP_DEBUG, "tab file or decoy tab file is null");
    }

  }

  free(match_output_folder);
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
