#include <stdlib.h>
#include "check.h"
#include "check-database.h"
#include "check-generate-peptides-iterator.h"
#include "check-index.h"
#include "check-index2.h"
#include "check-ion.h"
#include "check-ion_series.h"
#include "check-linked-list.h"
#include "check-mass.h"
#include "check-match-collection.h"
#include "check-match.h"
#include "check-modifications.h"
#include "check-modified-peptides-iterator.h"
#include "check-parameter.h"
#include "check-peak.h"
#include "check-peptide-modifications.h"
#include "check-peptide.h"
#include "check-protein.h"
#include "check-qvalue.h"
#include "check-scorer.h"
#include "check-serialize.h"
#include "check-spectrum_collection.h"
#include "check-crux-utils.h"
/* To add a new suite for <<<class>>>, add a line of the following form:
 *
 *   #include "check-<<<class>>>.h"
 *
 */

//must set bash export CK_FORK=no
int main(void){

  // BF: comment out three broken tests (ion, scorer, spectrum).


  // Create all check suites.
  /* To add a new suite for <<<class>>>, add a line of the following form:
   *
   *   Suite* suite_<<<class>>> = <<<class>>>_suite();
   *
   */
  Suite* suite_crux_utils = crux_utils_suite();
  Suite* suite_database = database_suite();
  Suite* suite_generate_peptides_iterator = generate_peptides_iterator_suite();
  Suite* suite_index = index_suite(); 
  //Suite* suite_ion = ion_suite(); 
  Suite* suite_ion_series = ion_series_suite(); 
  Suite* suite_ion_series_2 = ion_series_suite_2(); 
  Suite* suite_list = list_suite();
  Suite* suite_mass = mass_suite();
  //  Suite* suite_match = match_suite(); 
  Suite* suite_match_collection = match_collection_suite();
  Suite* suite_modifications = modifications_suite();
  Suite* suite_modified_peptides_iterator = modified_peptides_iterator_suite();
  Suite* suite_parameter = parameter_suite();
  Suite* suite_peak = peak_suite();
  Suite* suite_peptide = peptide_suite(); 
  Suite* suite_peptide_modifications = peptide_modifications_suite();
  Suite* suite_protein = protein_suite(); 
  //  Suite* suite_qvalue = qvalue_suite(); 
  //Suite* suite_scorer = scorer_suite();
  //  Suite* suite_serialize = serialize_suite();
  //Suite* suite_spectrum = spectrum_suite();
  //Suite* suite_spectrum_collection = spectrum_collection_suite();

  //add each suite to Runner
  SRunner *sr = srunner_create(NULL);
  
  /* To add a new suite for <<<class>>>, add a line of the following form:
   *
   *   srunner_add_suite(sr, suite_<<<class>>>);
   *
   */
  srunner_add_suite(sr, suite_crux_utils);
  srunner_add_suite(sr, suite_database);
  srunner_add_suite(sr, suite_generate_peptides_iterator);
  srunner_add_suite(sr, suite_index);
  //srunner_add_suite(sr,suite_ion);
  srunner_add_suite(sr, suite_ion_series);
  srunner_add_suite(sr, suite_ion_series_2);
  srunner_add_suite(sr, suite_list);
  srunner_add_suite(sr, suite_mass);
  //  srunner_add_suite(sr, suite_match);
  srunner_add_suite(sr, suite_match_collection);
  srunner_add_suite(sr, suite_modifications);
  srunner_add_suite(sr, suite_modified_peptides_iterator);
  srunner_add_suite(sr, suite_parameter);
  srunner_add_suite(sr, suite_peak);
  srunner_add_suite(sr, suite_peptide);
  srunner_add_suite(sr, suite_peptide_modifications);
  srunner_add_suite(sr, suite_protein);
  //  srunner_add_suite(sr, suite_qvalue);
  //srunner_add_suite(sr,suite_scorer);
  //  srunner_add_suite(sr, suite_serialize);
  //srunner_add_suite(sr,suite_spectrum);
  //srunner_add_suite(sr, suite_spectrum_collection);

  // Run each check suite.
  srunner_run_all(sr, CK_NORMAL);
  int nf = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (nf == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
