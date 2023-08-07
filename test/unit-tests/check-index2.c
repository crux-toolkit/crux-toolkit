#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "check-peak.h"
#include "mass.h"
#include "objects.h"
#include "Spectrum.h"
#include "Peak.h"
#include "Peptide.h"
#include "PeptideSrc.h"
#include "Protein.h"
#include "Database.h"
#include "Index.h"

// declare things to set up

void index_setup(){

}

void index_teardown(){
}

START_TEST (test_next_file){

}
END_TEST


START_TEST (test_create){
}
END_TEST

Suite* index_suite(void){
  Suite *s = suite_create("Index");
  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_next_file);
  tcase_add_test(tc_core, test_create);
  tcase_add_checked_fixture(tc_core, index_setup, index_teardown);

  suite_add_tcase(s, tc_core);

  // boundry conditions
  /*
  TCase *tc_limits = tcase_create("Limits");
  tcase_add_test(tc_limits, test_too_many_mods);
  tcase_add_checked_fixture(tc_limits, mod_setup, mod_teardown);
  suite_add_tcase(s, tc_limits);
  */

  return s;
}
