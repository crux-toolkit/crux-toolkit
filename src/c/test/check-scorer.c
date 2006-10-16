#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "check-peak.h"
#include "../spectrum.h"
#include "../spectrum_collection.h"
#include "../peak.h"
#include "../crux-utils.h"
#include "../scorer.h"
#include "../objects.h"
#include "../parameter.h"
#include "../ion_series.h"

#define scan_num 16
#define ms2_file "test.ms2"
#define parameter_file "test_parameter_file"

START_TEST (test_create){
  SPECTRUM_T* spectrum = NULL;
  SPECTRUM_COLLECTION_T * collection = NULL; ///<spectrum collection
  ION_SERIES_T* ion_series = NULL;
  SCORER_T* scorer = NULL;
  float score = 0;

  //parse paramter file
  parse_parameter_file(parameter_file);

  //set ion constraint to sequest settings
  ION_CONSTRAINT_T* ion_constraint = new_ion_constraint_sequest_sp();  
  
  //create new ion series
  ion_series = new_ion_series("AKLVKNMT", 2, ion_constraint);

  //now predict ions
  predict_ions(ion_series);

  //read ms2 file
  collection = new_spectrum_collection(ms2_file);
  spectrum = allocate_spectrum();
  
  //search for spectrum with correct scan number
  fail_unless(get_spectrum_collection_spectrum(collection, scan_num, spectrum), "failed to find scan_num in ms3 file");

  //create new scorer
  scorer = new_scorer(SP);  
  
  //calculates the Sp score
  score = score_spectrum_v_ion_series(scorer, spectrum, ion_series);

  //print the Sp score
  printf("Sp score is: %.2f\n", score);

  //free heap
  free_scorer(scorer);
  free_ion_constraint(ion_constraint);
  free_ion_series(ion_series);
  free_spectrum_collection(collection);
  free_spectrum(spectrum);
}
END_TEST

Suite *scorer_suite(void){
  Suite *s = suite_create("Scorer");
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  return s;
}
