/**
 * \file percolator.c
 */
/*
 * AUTHOR: Barbara Frewen
 * CREATE DATE: November 25, 2008
 * DESCRIPTION: Copied from match_analysis.c with only the percolator
 *         functionality kept.
 *         Given as input a directory containing binary psm files and
 *         a protein database, run percolator and return an sqt file
 *         with results.
 *
 *         Handles at most 4 files (target and decoy).  Looks for .csm
 *         files in the input directory and for corresponding
 *         -decoy[123].csm files.  Multiple target files in the given
 *         directory are concatinated together and presumed to be
 *         non-overlaping parts of the same ms2 file. 
 * 
 * $Revision: 1.2 $
 ****************************************************************************/
#include "percolator.h"

#ifdef PERCOLATOR
#define NUM_PERCOLATOR_OPTIONS 8
#define NUM_PERCOLATOR_ARGUMENTS 2
/* 
 * Private function declarations.  Details below
 */
MATCH_COLLECTION_T* run_percolator(
  char* psm_result_folder, 
  char* fasta_file, 
  char* feature_file); 


void print_sqt_file_perc( 
  MATCH_COLLECTION_T* match_collection,
  SCORER_TYPE_T scorer_type,
  SCORER_TYPE_T second_scorer_type
  );

#endif
  
/**
 * \brief crux-analyze-matches: takes in a directory containing binary
 * psm files and a protein index and analyzes the psms.
 */
int percolator_main(int argc, char** argv){

#ifdef PERCOLATOR
  /* Define command line arguments */
  int num_options = NUM_PERCOLATOR_OPTIONS;
  char* option_list[NUM_PERCOLATOR_OPTIONS] = {
    "version",
    "verbosity",
    "parameter-file",
    "write-parameter-file",
    "feature-file",
    "use-index",
    "overwrite",
    "sqt-output-file"
  };

  int num_arguments = NUM_PERCOLATOR_ARGUMENTS;
  char* argument_list[NUM_PERCOLATOR_ARGUMENTS] = {
    "psm-folder",
    "protein input",
  };

  /* for debugging handling of parameters*/
  set_verbosity_level(CARP_ERROR);

  /* Set up parameters and set defaults in parameter.c */
  initialize_parameters();

  /* Define optional and required arguments in parameter.c */
  select_cmd_line_options(option_list, num_options );
  select_cmd_line_arguments(argument_list, num_arguments);

  /* Parse the command line and optional paramter file
     does sytnax, type, and bounds checking and dies on error */
  parse_cmd_line_into_params_hash(argc, argv, "crux-analyze-matches");

  /* Get arguments */
  char* psm_file = get_string_parameter("psm-folder");
  char* protein_input_name = get_string_parameter("protein input");
  char* feature_file = get_string_parameter("feature-file");

  /* Get options */
  SCORER_TYPE_T scorer_type = PERCOLATOR_SCORE;
  SCORER_TYPE_T second_scorer_type = Q_VALUE;
  MATCH_COLLECTION_T* match_collection = NULL;

  /* Perform the analysis */
  carp(CARP_INFO, "Running percolator");
  match_collection = run_percolator(psm_file,
                                    protein_input_name,
                                    feature_file);
  scorer_type = PERCOLATOR_SCORE;
  second_scorer_type = Q_VALUE;
    
  carp(CARP_INFO, "Outputting matches.");
  print_sqt_file_perc(match_collection, scorer_type, second_scorer_type);

  // MEMLEAK below causes seg fault (or used to)
  // free_match_collection(match_collection);

  // clean up
  free(psm_file);
  free(protein_input_name);
  free(feature_file);


  carp(CARP_INFO, "crux percolator finished.");
  exit(0);

#else
  // Silence warnings about unused variables
  // We could use __attribute_unused__, but that is gcc specific
  int dummy = argc;
  dummy = strlen(argv[0]);
  dummy = 0;
  fprintf(stderr, "Percolator functionality not included.\n");
  exit(1);
#endif
}

#ifdef PERCOLATOR
/*  ****************** Subroutines ****************/


/*
 */
void print_sqt_file_perc(
  MATCH_COLLECTION_T* match_collection,
  SCORER_TYPE_T scorer,
  SCORER_TYPE_T second_scorer
  ){

  // get filename and open file
  char* sqt_filename = get_string_parameter("sqt-output-file");
  BOOLEAN_T overwrite = get_boolean_parameter("overwrite");
  FILE* sqt_file = create_file_in_path( sqt_filename, NULL, overwrite );

  // print header
  int num_proteins = get_match_collection_num_proteins(match_collection);
  print_sqt_header( sqt_file, "target", num_proteins, TRUE);

  fprintf(sqt_file, "H\tComment\tmatches analyzed by percolator\n");

  // get match iterator sorted by spectrum
  MATCH_ITERATOR_T* match_iterator = 
    new_match_iterator_spectrum_sorted(match_collection, scorer);

  // print each spectrum only once, keep track of which was last printed
  int cur_spectrum_num = -1;
  int cur_charge = 0;
  int match_counter = 0;
  int max_matches = get_int_parameter("max-sqt-result");

  // for all matches
  while( match_iterator_has_next(match_iterator) ){

    // get match and spectrum
    MATCH_T* match = match_iterator_next(match_iterator);
    SPECTRUM_T* spectrum = get_match_spectrum(match);
    int this_spectrum_num = get_spectrum_first_scan(spectrum);
    int charge = get_match_charge(match);

    carp(CARP_DETAILED_DEBUG, 
         "SQT printing scan %i (current %i), charge %i (current %i)", 
         this_spectrum_num, cur_spectrum_num, charge, cur_charge);

    // if this spectrum has not been printed...
    if( cur_spectrum_num != this_spectrum_num
        || cur_charge != charge){

      carp(CARP_DETAILED_DEBUG, "Printing new S line");
      // print S line to sqt file
      cur_spectrum_num = this_spectrum_num;
      cur_charge = charge;
      int num_peptides = get_match_ln_experiment_size(match);
      num_peptides = expf(num_peptides);

      print_spectrum_sqt(spectrum, sqt_file, num_peptides, charge);

      // print match to sqt file
      print_match_sqt(match, sqt_file, scorer, second_scorer);
      match_counter = 1;
    }
    // if this spectrum has been printed
    else{  
      if( match_counter < max_matches ){
        print_match_sqt(match, sqt_file, scorer, second_scorer);
        match_counter++;
      }
    }

  }// next match
  free_match_iterator(match_iterator);
  free(sqt_filename);

}

/**
 * \brief Analyze matches using the percolator algorithm
 * 
 * Runs Lukas Kall's percolator on the PSMs in the psm_result_folder
 * for a search against the sequence database fasta_file. Optionally 
 * puts the percolator PSM feature vectors into feature_file, if it is 
 * not NULL.
 * \returns a pointer to a MATCH_COLLECTION_T object
 * \callgraph
 */
MATCH_COLLECTION_T* run_percolator(
  char* psm_result_folder, 
  char* fasta_file, 
  char* feature_file){ 

  ALGORITHM_TYPE_T algorithm = PERCOLATOR_ALGORITHM;
  unsigned int number_features = 20;
  double* features = NULL;    
  double* results_q = NULL;
  double* results_score = NULL;
  double pi0 = get_double_parameter("pi0");
  char** feature_names = generate_feature_name_array(algorithm);
  MATCH_ITERATOR_T* match_iterator = NULL;
  MATCH_COLLECTION_T* match_collection = NULL;
  MATCH_COLLECTION_T* target_match_collection = NULL;
  MATCH_T* match = NULL;
  FILE* feature_fh = NULL;
  int set_idx = 0;
  
  // optional feature_file
  if(feature_file != NULL){  
    if((feature_fh = fopen(feature_file, "w")) == NULL){
      carp(CARP_FATAL, "Problem opening output file %s", feature_file);
      return NULL;
    }
  }

  carp(CARP_DETAILED_DEBUG, "Created feature file");

  // create MATCH_COLLECTION_ITERATOR_T object
  // which will read in the serialized output PSM results and return
  // first the match_collection of TARGET followed by 
  // the DECOY* match_collections.
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator =
    new_match_collection_iterator(psm_result_folder, fasta_file);

  if( match_collection_iterator == NULL ){
    carp(CARP_FATAL, "Failed to create a match collection iterator");
    exit(1);
  }
  carp(CARP_DETAILED_DEBUG, "Created the match collection iterator");

  // iterate over each, TARGET, DECOY 1..3 match_collection sets
  int iterations = 0;
  while(match_collection_iterator_has_next(match_collection_iterator)){
    
    carp(CARP_DEBUG, "Match collection iteration: %i" , iterations++);

    // get the next match_collection
    match_collection = 
      match_collection_iterator_next(match_collection_iterator);
    
    // intialize percolator, using information from first match_collection
    if(set_idx == 0){
      // the first match_collection is the target_match_collection
      target_match_collection = match_collection;

      // result array that stores the algorithm scores
      results_q = (double*)mycalloc(
          get_match_collection_match_total(match_collection), sizeof(double));
      results_score = (double*)mycalloc(
          get_match_collection_match_total(match_collection), sizeof(double));
      
      // Call that initiates percolator
      pcInitiate(
          (NSet)get_match_collection_iterator_number_collections(
                  match_collection_iterator), 
          number_features, 
          get_match_collection_match_total(match_collection), 
          feature_names, 
          pi0);
      
      // Call that sets verbosity level
      // 0 is quiet, 2 is default, 5 is more than you want
      if(verbosity < CARP_ERROR){
        pcSetVerbosity(0);
      }    
      else if(verbosity < CARP_INFO){
        pcSetVerbosity(1);
      }
      else{
        pcSetVerbosity(5);
      }
    }

    // create iterator, to register each PSM feature to Percolator
    match_iterator = new_match_iterator(match_collection, XCORR, FALSE);
    
    while(match_iterator_has_next(match_iterator)){
      match = match_iterator_next(match_iterator);
      // Register PSM with features to Percolator    
      features = get_match_percolator_features(match, match_collection);

      if (feature_fh != NULL){
        
        fprintf(feature_fh, "%i\t",
            get_spectrum_first_scan(get_match_spectrum(match))
            );
        if (get_match_null_peptide(match) == FALSE){
          fprintf(feature_fh, "1\t");
        } else { 
          fprintf(feature_fh, "-1\t");
        };

        unsigned int feature_idx;
        for (feature_idx = 0; feature_idx < number_features; feature_idx++){
          if (feature_idx < number_features - 1){
            fprintf(feature_fh, "%.4f\t", features[feature_idx]);
          } else {
            fprintf(feature_fh, "%.4f\n", features[feature_idx]);
          }
        }
      }
      
      pcRegisterPSM((SetType)set_idx, 
                    NULL, // no sequence used
                    features);
      
      free(features);
    }

    // ok free & update for next set
    // MEMLEAK 
    free_match_iterator(match_iterator);

    // don't free the target_match_collection
    if(set_idx != 0){
      free_match_collection(match_collection);
    }

    ++set_idx;
  } // end iteratation over each, TARGET, DECOY 1..3 match_collection sets

  if (feature_fh != NULL){
    fclose(feature_fh);
  }
  
  /***** PERCOLATOR run *********/

    carp(CARP_DETAILED_DEBUG, "got to here");
  // Start processing
  pcExecute(); 
  
  /* Retrieving target scores and qvalues after 
   * processing, the array should be numSpectra long and will be filled in 
   * the same order as the features were inserted */
  pcGetScores(results_score, results_q); 
       
  // fill results for Q_VALUE
  fill_result_to_match_collection(
      target_match_collection, results_q, Q_VALUE, TRUE);
  
  // fill results for PERCOLATOR_SCORE
  fill_result_to_match_collection(
      target_match_collection, results_score, PERCOLATOR_SCORE, FALSE);
   
  // Function that should be called after processing finished
  pcCleanUp();
  
  // TODO put free back in. took out because claimed it was double free
  // free names
  unsigned int name_idx;
  for(name_idx=0; name_idx < number_features; ++name_idx){
    free(feature_names[name_idx]);
  }
  free(feature_names);
  

  free(results_q);
  free(results_score);
  free_match_collection_iterator(match_collection_iterator);

  // TODO put free back in. took out because glibc claimed it was corrupted
  // double linked list
  // free_parameters();
  return target_match_collection;
}

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
