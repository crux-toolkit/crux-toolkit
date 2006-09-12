#include <stdlib.h>
#include <check.h>
#include "check-peak.h"
#include "check-spectrum.h"
#include "check-spectrum_collection.h"
#include "check-peptide.h"
#include "check-protein.h"
#include "check-index.h"
#include "check-database.h"

//must set bash export CK_FORK=no
int main(void){
  int nf;
  //create all check suite
  Suite* suite_peak = peak_suite();
  Suite* suite_spectrum = spectrum_suite();
  Suite* suite_spectrum_collection = spectrum_collection_suite();
  Suite* suite_peptide = peptide_suite(); 
  Suite* suite_protein = protein_suite(); 
  Suite* suite_database = database_suite();
  Suite* suite_index = index_suite(); 
  
  
  //add each suite to Runner
  SRunner *sr = srunner_create(NULL);
  srunner_add_suite(sr,suite_peak);
  srunner_add_suite(sr,suite_spectrum);
  srunner_add_suite(sr,suite_spectrum_collection);
  srunner_add_suite(sr,suite_peptide);
  srunner_add_suite(sr,suite_protein);
  srunner_add_suite(sr,suite_database);
  srunner_add_suite(sr,suite_index);


  //run each check suite
  srunner_run_all(sr, CK_NORMAL);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (nf == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
