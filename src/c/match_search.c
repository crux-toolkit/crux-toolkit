/*****************************************************************************
 * \file match_search.c
 * AUTHOR: Chris Park
 * CREATE DATE: 6/18/2007
 * DESCRIPTION: Given as input an ms2 file, a sequence database, and an optional parameter file, 
 * search all the spectrum against the peptides in the sequence database, and return high scoring peptides. 
 * ouput as binary ouput and optional sqt file format
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

#define NUM_SEARCH_OPTIONS 13
#define NUM_SEARCH_ARGS 2

/* Private functions */
int get_selected_charge_states();

int main(int argc, char** argv){

  /* Declarations */
  // command line options
  int verbosity;
  BOOLEAN_T use_index;
  double spectrum_min_mass; 
  double spectrum_max_mass; 
  char* match_output_folder = NULL; 
  //char* sqt_output_file = NULL;
  char* sqt_filename = NULL;
  //char* decoy_sqt_output_file = NULL;
  char* decoy_sqt_filename = NULL;
  int number_decoy_set;
  SCORER_TYPE_T main_score;
  SCORER_TYPE_T prelim_score;
  MATCH_SEARCH_OUTPUT_MODE_T output_type;

  // required arguments
  char* ms2_file = NULL;
  char* input_file = NULL;
  
  /* Define optional command line arguments */
  int num_options = NUM_SEARCH_OPTIONS;
  char* option_list[NUM_SEARCH_OPTIONS] = {
    "verbosity",
    "parameter-file",
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

  // parameter file options
  long int max_rank_preliminary = 500;
  long int max_rank_result = 500;
  int top_match = 1;
  BOOLEAN_T run_all_charges = TRUE;
  int spectrum_charge_to_run = 0;
  float mass_offset = 0; //delete me?


  /* for debugging of parameter processing */
  // TODO change to a make flag
  set_verbosity_level(CARP_DETAILED_DEBUG);
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

  
  /* Get parameter values */
  // inputs
  ms2_file = get_string_parameter_pointer("ms2 file");
  input_file = get_string_parameter_pointer("protein input");
  use_index = get_boolean_parameter("use-index");

  //outputs
  output_type = get_output_type_parameter("output-mode");
  match_output_folder = get_string_parameter_pointer("match-output-folder");
  sqt_filename = get_string_parameter_pointer("sqt-output-file");
  decoy_sqt_filename = get_string_parameter_pointer(
						"decoy-sqt-output-file");
  //TODO generate ms2-target.sqt file names if default is set, but do
  // it when file names are used,  see notes at end of file

  //searching  
  number_decoy_set = get_int_parameter("number-decoy-set");
  spectrum_charge_to_run = get_selected_charge_states();
  if( spectrum_charge_to_run > 0){
    run_all_charges = FALSE;
  }
  spectrum_min_mass = get_double_parameter("spectrum-min-mass");
  spectrum_max_mass =  get_double_parameter("spectrum-max-mass");

  //scoring
  main_score = get_scorer_type_parameter("score-type");
  prelim_score = get_scorer_type_parameter("prelim-score-type");
  max_rank_preliminary = get_int_parameter("max-rank-preliminary");

  //results
  max_rank_result = get_int_parameter("max-rank-result");//print to sqt
  // set max number of matches to be serialized per spectrum 
  top_match = get_int_parameter("top-match");
  
  // get mass offset from precursor mass to search for candidate peptides
  //delete me?
  mass_offset = get_double_parameter("mass-offset");    
  
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
  
  
  /************** Finished parameter setting **************/
  
  SPECTRUM_COLLECTION_T* collection = NULL; ///<spectrum collection
  SPECTRUM_ITERATOR_T* spectrum_iterator = NULL;
  MATCH_COLLECTION_T* match_collection = NULL;
  int possible_charge = 0;
  int* possible_charge_array = NULL;

  char** psm_result_filenames = NULL;
  //FILE** psm_result_file = NULL; //file handle array
  FILE** psm_file_array = NULL; //file handle array
  //FILE* psm_result_file_sqt = NULL;
  FILE* sqt_file = NULL;
  //  FILE* decoy_result_file_sqt  = NULL;
  FILE* decoy_sqt_file  = NULL;
  int total_files = number_decoy_set + 1; // plus one for target file
  int file_idx = 0;
  
  // read ms2 file
  collection = new_spectrum_collection(ms2_file);
  
  // parse the ms2 file for spectra
  if(!parse_spectrum_collection(collection)){
    carp(CARP_ERROR, "Failed to parse ms2 file: %s", ms2_file);
    free_spectrum_collection(collection);
    exit(1);
  }
  
  // create spectrum iterator
  spectrum_iterator = new_spectrum_iterator(collection);
  
  carp(CARP_DEBUG, "There were %i spectra found in the ms2 file",
       get_spectrum_collection_num_spectra(collection));

  // get psm result file handle array
  // this includes ones for the target and for the decoys
  psm_file_array = 
    get_spectrum_collection_psm_result_filenames( collection,
						  match_output_folder,
						  &psm_result_filenames,
						  number_decoy_set,
						  ".ms2"
                                                  );
  
  /* for debugging */
  int tempi = 0;
  for(tempi=0; tempi<total_files; tempi++){
    carp(CARP_DETAILED_DEBUG, "Result file name is %s", 
	 psm_result_filenames[tempi]);
  }
  
  // get psm_result sqt file handle if needed
  carp(CARP_DETAILED_DEBUG, "sqt output file is %s", sqt_filename);

  if(output_type != BINARY_OUTPUT ){ //ie binary only
    sqt_file= create_file_in_path(sqt_filename, match_output_folder);
    decoy_sqt_file = 
      create_file_in_path(decoy_sqt_filename, match_output_folder);
  }

  //do we really need sqts to stdout???
  //TODO turn this into a function
  /*  if(output_type == SQT_OUTPUT || output_type == ALL_OUTPUT){
    if (strcmp(sqt_filename, "STDOUT") == 0){
      sqt_file= stdout;
    } else {
      sqt_file= 
	create_file_in_path(sqt_filename, match_output_folder);
    }
    if (strcmp(decoy_sqt_filename, "STDOUT") == 0){
      decoy_sqt_file = stdout;
    } else {
      decoy_sqt_file = 
	create_file_in_path(decoy_sqt_filename, match_output_folder);
    }
  }
  */

  // check for at least one file handle for results
  if(psm_file_array[0] == NULL ||
     //((output_type == SQT_OUTPUT || output_type == ALL_OUTPUT) &&
     ((output_type != BINARY_OUTPUT) && sqt_file== NULL)){
    carp(CARP_FATAL, "Failed to create file handles for results");
    exit(1);
  }
  
  /**
   * General order of serialization is, 
   * - serialize_header
   * - serialize_psm_features
   * - serialize_total_number_of_spectra
   */
  
  /* Write headers to files */
  for(file_idx=0; file_idx < total_files; ++file_idx){
    serialize_header(collection, input_file, psm_file_array[file_idx]);
/********  TODO write this function
    write_sqt_header(sqt_file, decoy_file);

//Here's an eg header
H       SQTGenerator SEQUEST
H       SQTGeneratorVersion     2.7
H       Comment SEQUEST was written by J Eng and JR Yates, III
H       Comment SEQUEST ref. J. Am. Soc. Mass Spectrom., 1994, v. 4, p. 976
H       Comment SEQUEST ref. Eng,J.K.; McCormack A.L.; Yates J.R.
H       Comment SEQUEST is licensed to Finnigan Corp.
H       Comment Paralellization Program is run_ms2
H       Comment run_ms2 was written by Rovshan Sadygov
H       StartTime 05/18/2007, 05:02 AM
H       EndTime 05/18/2007, 10:47 AM
H       Database        /home/maccoss/dbase/human-050906-contam.fasta
H       DBSeqLength     16478009
H       DBLocusCount    34224
H       PrecursorMasses AVG
H       FragmentMasses  MONO
H       Alg-PreMassTol  3.000
H       Alg-FragMassTol 0.0
H       Alg-XCorrMode   0
H       StaticMod       C=160.139
H       Alg-DisplayTop  5
H       Alg-IonSeries   0 1 1 0.0 1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0
H       EnzymeSpec      No_Enzyme

*/
  }
  
  carp(CARP_DETAILED_DEBUG, "Headers written to output files");
  
  /* Prepare input, fasta or index */
  
  INDEX_T* index = NULL;
  DATABASE_T* database = NULL;
  BOOLEAN_T is_unique = get_boolean_parameter("unique-peptides");
  //todo make this a helper function
  if (use_index == TRUE){
    carp(CARP_DETAILED_DEBUG, "Using existing index");
    index = new_index_from_disk(input_file, is_unique);

    if (index == NULL){
      carp(CARP_FATAL, "Could not create index from disk for %s", input_file);
      exit(1);
    }
  } else {
    carp(CARP_DETAILED_DEBUG, "Using non-indexed fasta file");
    database = new_database(input_file, FALSE);         
    if( database == NULL ){
      carp(CARP_FATAL, "Could not read fasta file %s", input_file);
      exit(1);
    } 
    //BF added this, might not be correct
    parse_database(database);
  }
  
  /* Perform search: iterate over all spectra in ms2 file and score */
  int spectra_idx = 0; //superfluous
  int spectrum_searches_counter = 0; //spec (z-specific) w/ psms
  int bf_loop_i = 0;   //superfluous
  SPECTRUM_T* spectrum = NULL;
  int charge_index = 0;
  BOOLEAN_T is_decoy = FALSE;

  while(spectrum_iterator_has_next(spectrum_iterator)){
    bf_loop_i++;
    // get next spectrum
    spectrum = spectrum_iterator_next(spectrum_iterator);
    
    carp(CARP_DETAILED_DEBUG, 
	 "Searching spectrum number %i, loop number %i, search number %i", 
	 get_spectrum_first_scan(spectrum), bf_loop_i, spectra_idx+1);
    
    // select spectra that are within m/z target interval
    if(get_spectrum_precursor_mz(spectrum) <  spectrum_min_mass ||
       get_spectrum_precursor_mz(spectrum) >= spectrum_max_mass)
      {
	continue; //skip this spectrum, search next
      }
    
    // get possible charge state
    //TODO pass charge_to_run this function and eliminate first if()
    possible_charge = get_spectrum_num_possible_z(spectrum);
    possible_charge_array = get_spectrum_possible_z_pointer(spectrum);
    
    // iterate over all possible charge states for each spectrum
    for(charge_index = 0; charge_index < possible_charge; ++charge_index){
      
      // skip spectra that are not in the charge state to be run
      if(!run_all_charges && 
	      spectrum_charge_to_run != possible_charge_array[charge_index]){
	      continue;  //get next charge state
      }
      
      ++spectra_idx;
      spectrum_searches_counter++;
      
      // iterate over first for target next and for all decoy sets
      for(file_idx = 0; file_idx < total_files; ++file_idx){
	      if(file_idx == 0){ 	// is it target ?
	        is_decoy = FALSE;
	      }
	      else{
	        is_decoy = TRUE;
	      }
	
	// get match collection with scored, ranked match collection
	match_collection = 
	  new_match_collection_from_spectrum(
					     spectrum, 
					 possible_charge_array[charge_index], 
					     max_rank_preliminary, 
					     prelim_score, 
					     main_score, 
					     mass_offset, 
					     is_decoy,
					     index,
					     database
					     );

	if ((match_collection == NULL) && !is_decoy){
	  //if no targets found, don't search decoys
	  file_idx = total_files;
	  spectrum_searches_counter--;
	  continue;
	}

	// serialize the psm features from rank 1 to 'top_match'
	serialize_psm_features(match_collection, psm_file_array[file_idx], 
			       top_match, prelim_score, main_score);
	
	// write to SQT files
	// FIXME ONLY one header
	if(output_type == SQT_OUTPUT || output_type == ALL_OUTPUT){

	  if (file_idx == 0){
	    print_match_collection_sqt(sqt_file, max_rank_result,
				       match_collection, spectrum, 
				       prelim_score, main_score);
	  } else if (file_idx == 1){
	    print_match_collection_sqt(decoy_sqt_file, max_rank_result,
				       match_collection, spectrum, 
				       prelim_score, main_score);
	  } 
	} //next target/decoy        
	
	// free up match_collection
	free_match_collection(match_collection);          
      }
    }
  } //end while iterator has spectra
  
  // Modify the header serialized information for all files(target & decoy)
  // Set the total number of spectra serialized in the PSM result files
  for(file_idx=0; file_idx < total_files; ++file_idx){
    serialize_total_number_of_spectra(spectra_idx, psm_file_array[file_idx]);
    serialize_total_number_of_spectra(spectrum_searches_counter, psm_file_array[file_idx]);
  }
  
  // DEBUG
  carp(CARP_DEBUG, "Total spectrum searches: %d, total with psms: %i", 
       spectra_idx, spectrum_searches_counter);
  
  if(output_type == SQT_OUTPUT || output_type == ALL_OUTPUT){
    fclose(sqt_file);
    fclose(decoy_sqt_file);
  }
  
  // ok, now close all psm_result_files and free filenames
  for(file_idx = 0; file_idx < total_files; ++file_idx){
    fclose(psm_file_array[file_idx]);
    free(psm_result_filenames[file_idx]);
  }
  
  if (use_index == TRUE){
    free_index(index);
  } else {
    free_database(database);
  }
  free(psm_result_filenames);
  free(psm_file_array);
  free_spectrum_iterator(spectrum_iterator);
  free_spectrum_collection(collection);
  free_parameters();

  carp(CARP_INFO, "crux-search-for-matches finished");
  exit(0);
}


/* Private function definitions */

/*
  an alternative is to create a type.  It could include things
  like 2or3
 */
int get_selected_charge_states(){
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

