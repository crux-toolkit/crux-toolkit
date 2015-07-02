#include "check-<<class>>.h"
#include "<<class>>.h"

// declare things to set up
int myint1, myint2, *myintptr;

void <<class>>_setup(){
}

void <<class>>_teardown(){
}

START_TEST(test_create){
}
END_TEST


/* Boundry conditions test suite */
START_TEST(test_null){
}
END_TEST

Suite* <<class>>_suite(){
  Suite* s = suite_create("<<class>>");
  // Test basic features
  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_create);
  //  tcase_add_test(tc_core, test_...);

  tcase_add_checked_fixture(tc_core, <<class>>_setup, <<class>>_teardown);
  suite_add_tcase(s, tc_core);

  // Test boundry conditions
  TCase *tc_limits = tcase_create("Limits");
  tcase_add_test(tc_limits, test_null);
  tcase_add_checked_fixture(tc_limits, <<class>>_setup, <<class>>_teardown);
  suite_add_tcase(s, tc_limits);

  return s;
}













