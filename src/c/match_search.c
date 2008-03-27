/*************************************************************************//**
 * \file match_search.c
 * AUTHOR: Chris Park
 * CREATE DATE: 6/18/2007
 * DESCRIPTION: Given as input an ms2 file, a sequence database, and
 * an optional parameter file, search all the spectrum against the
 * peptides in the sequence database, and return high scoring
 * peptides. ouput as binary ouput and optional sqt file format
 * REVISION: 
 ****************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "carp.h"
#include "peptide.h"
#include "protein.h"
#include "parse_arguments.h"
#include "parameter.h"
#include "spectrum.h"
#include "spectrum_collection.h"
#include "generate_peptides_iterator.h"
#include "crux-utils.h"
#include "scorer.h"
#include "objects.h"
#include "match.h"
#include "match_collection.h"

#define NUM_SEARCH_OPTIONS 14
#define NUM_SEARCH_ARGS 2

/* Private functions */
int prepare_protein_input(char* input_file, 
                          INDEX_T** index, 
                          DATABASE_T** database);
void open_output_files(FILE*** binary_filehandle_array, 
                       FILE** sqt_filehandle,
                       FILE** decoy_sqt_filehandle);

int main(int argc, char** argv){

  /* Define optional command line arguments */
  int num_options = NUM_SEARCH_OPTIONS;
  char* option_list[NUM_SEARCH_OPTIONS] = {
    "verbosity",
    "parameter-file",
    "overwrite",
    "use-index",
    "prelim-score-type",
    "score-type",
    "spectrum-min-mass",
    "spectrum-max-mass",
    "spectrum-charge",
    "match-output-folder",
    "output-mode",
    "sqt-output-file",
    "decoy-sqt-output-file",
    "number-decoy-set"
  };

  /* Define required command line arguments */
  int num_arguments = NUM_SEARCH_ARGS;
  char* argument_list[NUM_SEARCH_ARGS] = {"ms2 file", "protein input"};

  /* for debugging of parameter processing */
  // TODO change to a make flag
  //set_verbosity_level(CARP_DETAILED_DEBUG);
  set_verbosity_level(CARP_ERROR);

  /* Initialize parameter.c and set default values*/
  initialize_parameters();

  /* Define optional and required arguments */
  select_cmd_line_options(option_list, num_options);
  select_cmd_line_arguments(argument_list, num_arguments);


  /* Parse the command line, including optional params file
     Includes syntax, type, and bounds checking, dies on error */
  parse_cmd_line_into_params_hash(argc, argv, "crux-search-for-matches");

  /* Set verbosity */
  verbosity = get_int_parameter("verbosity");
  set_verbosity_level(verbosity);

  /* Get input: ms2 file */
  char* ms2_file = get_string_parameter_pointer("ms2 file");

  // open ms2 file
  SPECTRUM_COLLECTION_T* collection = new_spectrum_collection(ms2_file);
  
  // parse the ms2 file for spectra
  carp(CARP_INFO, "Reading in ms2 file %s", ms2_file);
  if(!parse_spectrum_collection(collection)){
    carp(CARP_FATAL, "Failed to parse ms2 file: %s", ms2_file);
    free_spectrum_collection(collection);
    exit(1);
  }
  
  carp(CARP_DEBUG, "There were %i spectra found in the ms2 file",
       get_spectrum_collection_num_spectra(collection));

  /* Get input: protein file */
  char* input_file = get_string_parameter_pointer("protein input");

  /* Prepare input, fasta or index */
  INDEX_T* index = NULL;
  DATABASE_T* database = NULL;
  int num_proteins = prepare_protein_input(input_file, &index, &database); 

  carp(CARP_DEBUG, "Found %i proteins", num_proteins);

  /* Prepare output files */

  FILE** psm_file_array = NULL; //file handle array
  FILE* sqt_file = NULL;
  FILE* decoy_sqt_file  = NULL;

  open_output_files(&psm_file_array, &sqt_file, &decoy_sqt_file);

  //print headers
  serialize_headers(psm_file_array);
  print_sqt_header(sqt_file, "target", num_proteins);
  print_sqt_header(decoy_sqt_file, "decoy", num_proteins);

  /* Perform search: loop over spectra*/

  carp(CARP_INFO, "Searching spectra");

  // create spectrum iterator
  SPECTRUM_ITERATOR_T* spectrum_iterator = new_spectrum_iterator(collection);

  // get search parameters (most could be done within new_match_collection)
  long int max_rank_preliminary = get_int_parameter("max-rank-preliminary");
  SCORER_TYPE_T prelim_score = get_scorer_type_parameter("prelim-score-type");
  SCORER_TYPE_T main_score = get_scorer_type_parameter("score-type");
  double spectrum_min_mass = get_double_parameter("spectrum-min-mass");
  double spectrum_max_mass =  get_double_parameter("spectrum-max-mass");

  // get list of mods
  PEPTIDE_MOD_T* peptide_mods = NULL;
  // uses aa_mods in parameter.c
  int num_peptide_mods = generate_peptide_mod_list( &peptide_mods );
  // so it will compile
  carp(CARP_DEBUG, "Got %d peptide mods", num_peptide_mods);

  // flags and counters for loop
  int spectrum_counter = 0;
  int spectrum_searches_counter = 0; //for psm file header, spec*charges
  int z_i = 0;
  int file_i = 0;
  int total_files = get_int_parameter("number-decoy-set") + 1;

  // find matches for each spectrum
  while(spectrum_iterator_has_next(spectrum_iterator)){

    SPECTRUM_T* spectrum = spectrum_iterator_next(spectrum_iterator);

    // select spectra that are within m/z target interval
    if(get_spectrum_precursor_mz(spectrum) <  spectrum_min_mass ||
       get_spectrum_precursor_mz(spectrum) >= spectrum_max_mass){
      continue; //skip this spectrum, search next
      }
    
    //for each charge, search spectrum
    int* charge_array = NULL;
    int num_charges = get_charges_to_search(spectrum, &charge_array);

    for(z_i=0; z_i < num_charges; z_i++){
      int charge = charge_array[z_i];
      carp(CARP_DETAILED_DEBUG, 
           "Searching spectrum number %i, charge %i, search number %i",
           get_spectrum_first_scan(spectrum), charge,
           spectrum_searches_counter+1 );

      // for each database (real/rand), search spectrum
      BOOLEAN_T is_decoy = FALSE; //first target, then decoys

      for(file_i=0; file_i < total_files; file_i++){

        carp(CARP_DETAILED_DEBUG, "csm alloc");
        MATCH_COLLECTION_T* match_collection = 
          new_match_collection_from_spectrum( spectrum,
                                              charge,
                                              max_rank_preliminary,
                                              prelim_score,
                                              main_score,
                                              0,//mass_offset,
                                              is_decoy,
                                              index,
                                              database);
        if(match_collection == NULL){
          file_i = total_files; //don't search decoys
          spectrum_searches_counter--;  //don't count this search
          continue;
        }
        carp(CARP_DETAILED_DEBUG, "about to print matches");
        
        print_matches(match_collection, spectrum, is_decoy,
                      psm_file_array[file_i], sqt_file, decoy_sqt_file);


        free_match_collection(match_collection);
        is_decoy = TRUE;
        //exit(1); // to get gmon.out for a single spectrum uncomment this line
      }// next set (target, decoy, decoy...)

      spectrum_searches_counter++;
    }// next charge state, same spectrum

    free(charge_array);
    
    spectrum_counter++;
    if( spectrum_counter %1000 == 0 ){
      carp(CARP_INFO, "Searched %d spectra", spectrum_counter);
    }
  }// next spectrum    
  carp(CARP_DEBUG, "finished searching");

  /* Post-search processing */

  // update header with number of searches
  int file_idx;
  for(file_idx=0; file_idx < total_files; ++file_idx){
    carp(CARP_DEBUG, "Updating header with %d searches", 
         spectrum_searches_counter);
    serialize_total_number_of_spectra(spectrum_searches_counter, 
                                      psm_file_array[file_idx]);
  }
  carp(CARP_DEBUG, "Fixed headers");

  // clean up files
  for(file_idx = 0; file_idx < total_files; ++file_idx){
    if( psm_file_array[file_idx] != NULL){
      fclose(psm_file_array[file_idx]);
    }
  }
  if( sqt_file != NULL ){
    fclose(sqt_file);
  }
  if( decoy_sqt_file != NULL ){
    fclose(decoy_sqt_file);
  }

  // clean up memory
  free(psm_file_array);
  free_spectrum_iterator(spectrum_iterator);
  free_spectrum_collection(collection);
  free_parameters();

  carp(CARP_INFO, "crux-search-for-matches finished");
  exit(0);


  //where is this used?  Should it be put back?
  
  // seed for random rnumber generation
  //hide this from user?
  if(strcmp(get_string_parameter_pointer("seed"), "time")== 0){
    time_t seconds; // use current time to seed
    time(&seconds); // Get value from sys clock and set seconds variable.
    srand((unsigned int) seconds); // Convert seconds to a unsigned int
  }
  else{
    srand((unsigned int)atoi(get_string_parameter_pointer("seed")));
  }
  
 
}


/* Private function definitions */

/*
  an alternative is to create a type.  It could include things
  like 2or3
 */
/*int get_selected_charge_states(){
  int charge_state = 0;

  char* charge_str = get_string_parameter_pointer("spectrum-charge");

  if( strcmp( charge_str, "all") == 0){
    return charge_state;
  }

  charge_state = atoi(charge_str);

  if( (charge_state < 1) || (charge_state > 3) ){
    carp(CARP_FATAL, "spectrum-charge option must be 1,2,3, or 'all'.  " \
    "%s is not valid", charge_str);
    exit(1);
  }
  return charge_state;
}
*/

int prepare_protein_input(char* input_file, 
                          INDEX_T** index, 
                          DATABASE_T** database){

  int num_proteins = 0;
  BOOLEAN_T use_index = get_boolean_parameter("use-index");
  BOOLEAN_T is_unique = get_boolean_parameter("unique-peptides");

  if (use_index == TRUE){
    carp(CARP_INFO, "Preparing protein index %s", input_file);
    *index = new_index_from_disk(input_file, is_unique);

    if (index == NULL){
      carp(CARP_FATAL, "Could not create index from disk for %s", input_file);
      exit(1);
    }
    num_proteins = get_index_num_proteins(*index);

  } else {
    carp(CARP_INFO, "Preparing protein fasta file %s", input_file);
    *database = new_database(input_file, FALSE);         
    if( database == NULL ){
      carp(CARP_FATAL, "Could not read fasta file %s", input_file);
      exit(1);
    } 
    //BF added this, might not be correct
    parse_database(*database);
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
  FILE*** psm_file_array, ///< will put binary psm filehandles here -out
  FILE** sqt_file,        ///< will put text sqt filehandle here -out
  FILE** decoy_sqt_file)  ///< will put decoy sqt filehandle here -out
{
  char* match_output_folder = get_string_parameter_pointer(
                                                    "match-output-folder");
  MATCH_SEARCH_OUTPUT_MODE_T output_type = get_output_type_parameter(
                                                    "output-mode");
  BOOLEAN_T overwrite = get_boolean_parameter("overwrite");
  carp(CARP_DEBUG, "The output type is %d (binary, sqt, all)" \
       " and overwrite is '%d'", (int)output_type, (int)overwrite);


  // create binary psm files (allocate memory, even if not used)
  *psm_file_array = create_psm_files();

  if( output_type != BINARY_OUTPUT ){ //ie sqt or all

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

  carp(CARP_DEBUG, "Finished opening output files");
}


/* NOTES */
// generate sqt ouput file if not set by user
//TODO move the generation of file name to where name is used
//       I don't think it's working anyway
/*    if(strcmp(
      get_string_parameter_pointer("sqt-output-file"), "target.sqt") ==0){
      sqt_output_file =generate_name(ms2_file, "-target.sqt", ".ms2", NULL);
      decoy_sqt_output_file = 
      generate_name(ms2_file, "-decoy.sqt", ".ms2", NULL);
      set_string_parameter("sqt-output-file", sqt_output_file);
      set_string_parameter("decoy-sqt-output-file", decoy_sqt_output_file);
      }
*/


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

