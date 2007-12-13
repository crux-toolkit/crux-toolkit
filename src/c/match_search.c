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

#define NUM_SEARCH_OPTIONS 14
#define NUM_SEARCH_ARGS 2
/**
 * When wrong command is seen carp, and exit
 */
void wrong_command(char* arg, char* comment){
  char* usage = parse_arguments_get_usage("search_spectra");
  carp(CARP_FATAL, "incorrect argument: %s", arg);

  // print comment if given
  if(comment != NULL){
    carp(CARP_FATAL, "%s", comment);
  }

  fprintf(stderr, "%s", usage);
  free(usage);
  exit(1);
}

int main(int argc, char** argv){

  /* Declarations */
  int verbosity;
  double spectrum_min_mass; 
  double spectrum_max_mass; 
  char* spectrum_charge_str = NULL;
  double number_runs;
  char* match_output_folder = NULL; 
  char* sqt_output_file = NULL;
  char* decoy_sqt_output_file = NULL;
  int number_decoy_set;
  
  // required
  char* ms2_file = NULL;
  char* fasta_file = NULL;
  
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
    "number-runs",       //delete this
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
  // change to a make flag
  set_verbosity_level(CARP_DETAILED_DEBUG);

  /* Initialize parameter.c and set default values*/
  initialize_parameters();

  /* Define optional and required arguments */
  select_cmd_line_options(option_list, num_options);
  select_cmd_line_arguments(argument_list, num_arguments);


  /* Parse the command line, including optional params file
     Includes syntax, type, and bounds checking, dies on error */
  parse_cmd_line_into_params_hash(argc, argv);

  SPECTRUM_T* spectrum = NULL;
  SPECTRUM_COLLECTION_T* collection = NULL; ///<spectrum collection
  SPECTRUM_ITERATOR_T* spectrum_iterator = NULL;
  MATCH_COLLECTION_T* match_collection = NULL;
  int possible_charge = 0;
  int* possible_charge_array = NULL;
  int charge_index = 0;
  long int max_rank_preliminary = 500;
  long int max_rank_result = 500;
  int top_match = 1;
  float mass_offset = 0;
  BOOLEAN_T run_all_charges = TRUE;
  int spectrum_charge_to_run = 0;
  
  /* Set verbosity */
  verbosity = get_int_parameter("verbosity");
  set_verbosity_level(verbosity);

  
  //TODO move the generation of file name to where name is used
  //       I don't think it's working anyway
  // generate sqt ouput file if not set by user
  /*    if(strcmp(
	get_string_parameter_pointer("sqt-output-file"), "target.sqt") ==0){
	sqt_output_file =generate_name(ms2_file, "-target.sqt", ".ms2", NULL);
	decoy_sqt_output_file = 
        generate_name(ms2_file, "-decoy.sqt", ".ms2", NULL);
	set_string_parameter("sqt-output-file", sqt_output_file);
	set_string_parameter("decoy-sqt-output-file", decoy_sqt_output_file);
	}
  */
  
  /* Get parameter values */
  ms2_file = get_string_parameter_pointer("ms2 file");
  fasta_file = get_string_parameter_pointer("protein input");
  match_output_folder = get_string_parameter_pointer("match-output-folder");
  sqt_output_file = get_string_parameter_pointer("sqt-output-file");
  decoy_sqt_output_file = get_string_parameter_pointer(
						"decoy-sqt-output-file");
    
  // how many runs of search to perform
  //what is this??? a max number of spec to search.  Disable this
  number_runs = get_double_parameter("number-runs");
  
  // what charge state of spectra to search
  spectrum_charge_str = get_string_parameter_pointer("spectrum-charge");
  spectrum_charge_to_run = atoi(spectrum_charge_str);

  if( spectrum_charge_to_run < 1 ){
    //assert that it was 'all' and not other string
    run_all_charges = TRUE;
  }else if( spectrum_charge_to_run > 3 ){
    carp(CARP_FATAL, "spectrum-charge option must be 1,2,3, or 'all'.  " \
	 "%s is not valid", spectrum_charge_str);
    exit(1);
  }

  // number_decoy_set
  number_decoy_set = get_int_parameter("number-decoy-set");
  
  // main score type
  //FIXME
  SCORER_TYPE_T main_score;//get_scorer_type_parameter("score-type");
  char* score_type = get_string_parameter_pointer("score-type");
  string_to_scorer_type(score_type, &main_score);
  SCORER_TYPE_T prelim_score = get_scorer_type_parameter("prelim-score-type");

  // get output-mode
  MATCH_SEARCH_OUTPUT_MODE_T output_type = get_output_type_parameter(
                                             "output-mode");

  spectrum_min_mass = get_double_parameter("spectrum-min-mass");
  spectrum_max_mass =  get_double_parameter("spectrum-max-mass");
  // set max number of preliminary scored peptides to use for final scoring
  max_rank_preliminary = get_int_parameter("max-rank-preliminary");
  
  // set max number of final scoring matches to print as output in sqt
  max_rank_result = get_int_parameter("max-rank-result");
  
  // set max number of matches to be serialized per spectrum
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
  
  char** psm_result_filenames = NULL;
  FILE** psm_result_file = NULL; //file handle array
  FILE* psm_result_file_sqt = NULL;
  FILE* decoy_result_file_sqt  = NULL;
  int total_files = number_decoy_set + 1; // plus one for target file
  int file_idx = 0;
  BOOLEAN_T is_decoy = FALSE;
  
  // read ms2 file
  collection = new_spectrum_collection(ms2_file);
  
  // parse the ms2 file for spectra
  if(!parse_spectrum_collection(collection)){
    carp(CARP_ERROR, "failed to parse ms2 file: %s", ms2_file);
    // free, exit
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
  
  // did we get the file handles?
  // check that there's at least one for result
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
  
  // serialize the header information for all files(target & decoy)
  for(file_idx=0; file_idx < total_files; ++file_idx){
    serialize_header(collection, fasta_file, psm_result_file[file_idx]);
  }
  
  carp(CARP_DETAILED_DEBUG, "Headers written to output files");
  
  BOOLEAN_T use_index = get_boolean_parameter("use-index");
  char* in_file = get_string_parameter_pointer("protein input");
  
  INDEX_T* index = NULL;
  DATABASE_T* database = NULL;
  BOOLEAN_T is_unique = get_boolean_parameter("unique-peptides");
  if (use_index == TRUE){
    carp(CARP_DETAILED_DEBUG, "Using existing index");
    if ((index = new_index_from_disk(in_file, is_unique)) == NULL){
      carp(CARP_FATAL, "Could not create index from disk for %s", in_file);
      exit(1);
    }
  } else {
    carp(CARP_DETAILED_DEBUG, "Using non-indexed fasta file");
    database = new_database(in_file, FALSE);         
    //BF added this, might not be correct
    parse_database(database);
  }
  
  int spectra_idx = 0;
  int bf_spectrum_i = 0;
  // iterate over all spectrum in ms2 file and score
  while(spectrum_iterator_has_next(spectrum_iterator)){
    bf_spectrum_i++;
    carp(CARP_DETAILED_DEBUG, 
	 "Searching spectrum number %i, search number %i", 
	 bf_spectrum_i, spectra_idx+1);
    
    // check if total runs exceed limit user defined
    //TODO disable this
    if(number_runs <= spectra_idx){
      break;
    }
    
    // get next spectrum
    spectrum = spectrum_iterator_next(spectrum_iterator);
    
    // select spectra that are within m/z target interval
    if(get_spectrum_precursor_mz(spectrum) <  spectrum_min_mass ||
       get_spectrum_precursor_mz(spectrum) >= spectrum_max_mass)
      {
	continue;
      }
    
    // get possible charge state
    possible_charge = get_spectrum_num_possible_z(spectrum);
    possible_charge_array = get_spectrum_possible_z_pointer(spectrum);
    
    // iterate over all possible charge states for each spectrum
    for(charge_index = 0; charge_index < possible_charge; ++charge_index){
      
      // skip spectra that are not in the charge state to be run
      if(!run_all_charges && 
	 spectrum_charge_to_run != possible_charge_array[charge_index]){
	continue;
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
	// serialize the psm features to ouput file upto 'top_match' number of 
	// top peptides among the match_collection
	// carp(CARP_WARNING, "Outputting to %s\n", psm_result_file[file_idx]);
	
	serialize_psm_features(match_collection, psm_result_file[file_idx], 
			       top_match, prelim_score, main_score);
	
	// should I output the match_collection result as a SQT file? // Output only for the target set
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
  }
  
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
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
