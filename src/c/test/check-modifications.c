#include "check-modifications.h"
#include "../modifications.h"
//#include "../<other-needed-headers>.h"

// declare things to set up
AA_MOD_T *amod1, *amod2, *amod3;
AA_MOD_T* amod_list[3];
PEPTIDE_MOD_T* pmod1, *pmod2;

void mod_setup(){
  // create stuff
  // item1 = new_<NAME>();
}

void mod_teardown(){
  // clean up
  // free(item1);
}

/*
START_TEST(test_aa_init){
  // initialize an aa_mod try more than one index
  // use getters to test fields
}
END_TEST

START_TEST(test_aa_setters){
  // change mass, max, position, distance
  // test with getters
}
END_TEST

START_TEST(test_aa_list){
  // get pointer to aa list
  // change some values
  // test with ??
}
END_TEST
*/

START_TEST(test_create){
  // set stuff up (if it wasn't done in setup)

  // do stuff and test it
  fail_unless( 1 == 1, "this should have worked");
  //  fail_unless( item1->something() == something, "failure message");

}
END_TEST

START_TEST(test_somelimit){
  // think of a boundry condition to test
}
END_TEST

Suite* modifications_suite(){
  Suite* s = suite_create("Modifications\n");
  // Test basic features
  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_create);
  tcase_add_checked_fixture(tc_core, mod_setup, mod_teardown);
  suite_add_tcase(s, tc_core);

  // Test boundry conditions
  TCase *tc_limits = tcase_create("Limits");
  tcase_add_test(tc_limits, test_somelimit);
  tcase_add_checked_fixture(tc_core, mod_setup, mod_teardown);
  suite_add_tcase(s, tc_limits);

  return s;
}













