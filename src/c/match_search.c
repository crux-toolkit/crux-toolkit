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
  //char* spectrum_charge_str = NULL;
  //  double number_runs;
  char* match_output_folder = NULL; 
  char* sqt_output_file = NULL;
  char* decoy_sqt_output_file = NULL;
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
    //"number-runs",       //delete this
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
  float mass_offset = 0;


  /* for debugging of parameter processing */
  // TODO change to a make flag
  set_verbosity_level(CARP_DETAILED_DEBUG);

  /* Initialize parameter.c and set default values*/
  initialize_parameters();

  /* Define optional and required arguments */
  select_cmd_line_options(option_list, num_options);
  select_cmd_line_arguments(argument_list, num_arguments);


  /* Parse the command line, including optional params file
     Includes syntax, type, and bounds checking, dies on error */
  parse_cmd_line_into_params_hash(argc, argv);

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
  sqt_output_file = get_string_parameter_pointer("sqt-output-file");
  decoy_sqt_output_file = get_string_parameter_pointer(
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
  // set max number of matches to be serialized per spectrum ??
  top_match = get_int_parameter("top-match");
  
  // get mass offset from precursor mass to search for candidate peptides
  //get rid of this?
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
  FILE** psm_result_file = NULL; //file handle array
  FILE* psm_result_file_sqt = NULL;
  FILE* decoy_result_file_sqt  = NULL;
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
  
  carp(CARP_DETAILED_DEBUG, "There were %i spectra found in the ms2 file",
       get_spectrum_collection_num_spectra(collection));
  // get psm_result file handler array
  // this includes one for the target and for the decoys
  psm_result_file = 
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
  //do we really need sqts to stdout???
  carp(CARP_DETAILED_DEBUG, "sqt output file is %s", sqt_output_file);
  if(output_type == SQT_OUTPUT || output_type == ALL_OUTPUT){
    if (strcmp(sqt_output_file, "STDOUT") == 0){
      psm_result_file_sqt = stdout;
    } else {
      psm_result_file_sqt = 
	create_file_in_path(sqt_output_file, match_output_folder);
    }
    if (strcmp(decoy_sqt_output_file, "STDOUT") == 0){
      decoy_result_file_sqt = stdout;
    } else {
      decoy_result_file_sqt = 
	create_file_in_path(decoy_sqt_output_file, match_output_folder);
    }
  }
  
  // check for at least one file handle for results
  if(psm_result_file[0] == NULL ||
     ((output_type == SQT_OUTPUT || output_type == ALL_OUTPUT) &&
      psm_result_file_sqt == NULL)){
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
    serialize_header(collection, input_file, psm_result_file[file_idx]);
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
  int spectra_idx = 0;
  int bf_spectrum_i = 0;
  SPECTRUM_T* spectrum = NULL;
  int charge_index = 0;
  BOOLEAN_T is_decoy = FALSE;

  while(spectrum_iterator_has_next(spectrum_iterator)){
    bf_spectrum_i++;
    carp(CARP_DETAILED_DEBUG, 
	 "Searching spectrum number %i, search number %i", 
	 bf_spectrum_i, spectra_idx+1);
    
    // check if total runs exceed limit user defined
    //TODO disable this
    /*    if(number_runs <= spectra_idx){
      break;
    }
    */

    // get next spectrum
    spectrum = spectrum_iterator_next(spectrum_iterator);
    
    // select spectra that are within m/z target interval
    if(get_spectrum_precursor_mz(spectrum) <  spectrum_min_mass ||
       get_spectrum_precursor_mz(spectrum) >= spectrum_max_mass)
      {
	continue; //get next spectrum
      }
    
    // get possible charge state
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
      
      // iterate over first for target next and for all decoy sets
      for(file_idx = 0; file_idx < total_files; ++file_idx){
	// is it target ?
	if(file_idx == 0){
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
	if (match_collection == NULL){
	  continue;
	}

	// serialize the psm features from rank 1 to 'top_match'
	serialize_psm_features(match_collection, psm_result_file[file_idx], 
			       top_match, prelim_score, main_score);
	
	// write to SQT files
	// Output only for the target set (bf why??)
	// FIXME ONLY one header
	if(output_type == SQT_OUTPUT || output_type == ALL_OUTPUT){
	  // only output the first and second decoy sets
	  if (file_idx == 0){
	    print_match_collection_sqt(psm_result_file_sqt, max_rank_result,
				       match_collection, spectrum, 
				       prelim_score, main_score);
	  } else if (file_idx == 1){
	    print_match_collection_sqt(decoy_result_file_sqt, max_rank_result,
				       match_collection, spectrum, 
				       prelim_score, main_score);
	  } 
	}        
	
	// free up match_collection
	free_match_collection(match_collection);          
      }
    }
  } //end while iterator has spectra
  
  // Modify the header serialized information for all files(target & decoy)
  // Set the total number of spectra serialized in the PSM result files
  for(file_idx=0; file_idx < total_files; ++file_idx){
    serialize_total_number_of_spectra(spectra_idx, psm_result_file[file_idx]);
  }
  
  // DEBUG
  carp(CARP_DEBUG, "Total spectrum searches: %d", spectra_idx);
  
  if(output_type == SQT_OUTPUT || output_type == ALL_OUTPUT){
    fclose(psm_result_file_sqt);
    fclose(decoy_result_file_sqt);
  }
  
  // ok, now close all psm_result_files and free filenames
  for(file_idx = 0; file_idx < total_files; ++file_idx){
    fclose(psm_result_file[file_idx]);
    free(psm_result_filenames[file_idx]);
  }
  
  if (use_index == TRUE){
    free_index(index);
  } else {
    free_database(database);
  }
  free(psm_result_filenames);
  free(psm_result_file);
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

