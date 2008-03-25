#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "check-peak.h"
#include "check-match.h"
#include "../spectrum.h"
#include "../spectrum_collection.h"
#include "../peak.h"
#include "../crux-utils.h"
#include "../scorer.h"
#include "../objects.h"
#include "../parameter.h"
#include "../ion_series.h"
#include "../generate_peptides_iterator.h"
#include "../match.h"
#include "../match_collection.h"

#define scan_num 16
#define ms2_file "test.ms2"
#define parameter_file "test_parameter_file"

//THE parameter system will not work if set CK_FORK=no
START_TEST (test_create){
  SPECTRUM_T* spectrum = NULL;
  SPECTRUM_COLLECTION_T* collection = NULL; ///<spectrum collection
  MATCH_COLLECTION_T* match_collection = NULL;
  MATCH_ITERATOR_T* match_iterator = NULL;
  MATCH_T* match = NULL;
    
  /********** comment this parameter section out, when using CK_FORK=no, valgrind ******/
  // Parameters have been already confirmed in check_scorer.
  
  //set verbbosity level
  int  verbosity = CARP_INFO;
  set_verbosity_level(verbosity);

  //parse paramter file
  //  parse_update_parameters(parameter_file);
  initialize_parameters();

  //add fasta file parameter_file fasta-file
  //add_parameter("fasta-file", "fasta_file");
  
  //parameters has been confirmed
  //parameters_confirmed();
  /***************************************************/

  //read ms2 file
  collection = new_spectrum_collection(ms2_file);
  spectrum = allocate_spectrum();

  DATABASE_T* database = new_database("fasta-file", FALSE);

  
  //search for spectrum with correct scan number
  fail_unless(get_spectrum_collection_spectrum(collection, scan_num, spectrum), "failed to find scan_num in ms3 file");
  
  //get match collection with perliminary score of SP, and main score of XCORR
  match_collection = new_match_collection_from_spectrum(spectrum, 1, 500, SP, XCORR, 0, FALSE, NULL, database);
  
  fail_unless(get_match_collection_scored_type(match_collection, SP), "failed to set match_collection scored type, SP");
  fail_unless(get_match_collection_scored_type(match_collection, XCORR), "failed to set match_collection scored type, SP");

  //LOGP_EXP_SP should not be scored yet
  fail_unless(!get_match_collection_scored_type(match_collection, LOGP_EXP_SP), "failed to set match_collection scored type, xcorr");
  fail_unless(!get_match_collection_iterator_lock(match_collection), "match_collection lock is not set correctly"); 
  
  //create match iterator
  match_iterator = new_match_iterator(match_collection, SP, TRUE);
  
  //match_collection should be locked now..
  fail_unless(get_match_collection_iterator_lock(match_collection), "match_collection lock is not set correctly"); 
  
  //iterate over all matches
  while(match_iterator_has_next(match_iterator)){
    match = match_iterator_next(match_iterator);
    print_match(match, stdout, TRUE, SP);
  }

  //free match iterator
  free_match_iterator(match_iterator);
  
  //should be unlocked
  fail_unless(!get_match_collection_iterator_lock(match_collection), "match_collection lock is not set correctly"); 

  free_match_collection(match_collection);
  free_spectrum_collection(collection);
  free_spectrum(spectrum);
  free_parameters();
}
END_TEST


Suite *match_suite(void){
  Suite *s = suite_create("match");
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  return s;
}
