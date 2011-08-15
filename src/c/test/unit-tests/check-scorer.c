#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "check-peak.h"
#include "Spectrum.h"
#include "SpectrumCollectionFactory.h"
#include "Peak.h"
#include "crux-utils.h"
#include "Scorer.h"
#include "objects.h"
#include "parameter.h"
#include "IonSeries.h"
#include "IonConstraint.h"

#define scan_num 16
#define ms2_file "test.ms2"
#define parameter_file "test_parameter_file"

START_TEST (test_create){
  Spectrum* spectrum = NULL;
  SpectrumCollection* sp_collection = NULL; ///<spectrum collection
  IonSeries* ion_series = NULL;
  Scorer* scorer = NULL;
  float score = 0;

  //parse paramter file
  //  parse_update_parameters(parameter_file);
  
  //set parameter for fasta_file, although not used here...
  //set_string_parameter("fasta-file", "fasta_file");
  
  //parameters has been confirmed
  //parameters_confirmed();
  
  initialize_parameters();
  int peptide_charge = get_int_parameter("charge");

  //set ion constraint to sequest settings
  IonConstraint* ion_constraint = IonConstraint::newIonConstraintSequestSp(peptide_charge);  
  
  //create new ion series
  ion_series = new IonSeries("AKLVKNMT", 2, ion_constraint);

  //now predict ions
  ion_series->predictIons();

  //read ms2 file
  sp_collection = SpectrumCollectionFactory::create(ms2_file);
  spectrum = new Spectrum();
  
  //search for spectrum with correct scan number
  fail_unless(sp_collection->getSpectrum(scan_num) == NULL, "failed to find scan_num in ms3 file");

  //create new scorer
  scorer = new Scorer(SP);  

  //check if scorer has been set right
  fail_unless(scorer->getType() == SP, "failed to set scorer type");
  fail_unless(compare_float(scorer->getSpBeta(), 0.075) == 0, "failed to set beta");
  fail_unless(compare_float(scorer->getSpMaxMz(), 4000) == 0, "failed to set max mz");

  //calculates the Sp score
  score = scorer->scoreSpectrumVIonSeries(spectrum, ion_series);

  //print the Sp score
  printf("Sp score is: %.2f\n", score);
  
  fail_unless(compare_float(score, 5.35885334014892578125) == 0, "sp score does not match the expected value");
  
  //free heap
  delete scorer;
  delete ion_constraint;
  delete ion_series;
  delete sp_collection;
  delete spectrum;
}
END_TEST

Suite *scorer_suite(void){
  Suite *s = suite_create("Scorer");
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  return s;
}
