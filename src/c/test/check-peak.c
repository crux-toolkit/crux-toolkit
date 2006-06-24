#include <math.h>
#include <stdlib.h>
#include "check-peak.h"
#include "../peak.h"

PEAK_T * test_peak;
void setup(void){
  test_peak = new_peak(20.0,3.0); //intensity and m/z
}
void teardown(void){
  free_peak(test_peak);
}

START_TEST (test_create){
  fail_unless( peak_intensity(test_peak) == 20.1, "Peak height not setup correctly to 20.0");
  //fail_unless( peak_mz(p) == 3, "Peak m/z not setup correctly to 3");
}
END_TEST

/*
START_TEST (test_create2){
  fail_unless( peak_height(p) == 20, "Peak height not setup correctly to 20");
  fail_unless( peak_location(test_peak) == 3.5, "Peak m/z not setup correctly to 3.0");
}
END_TEST
*/

Suite *peak_suite(void){
  Suite *s = suite_create("Peak");
  TCase *tc_core = tcase_create("Core");
  //TCase *tc_core2 = tcase_create("Core2");

  suite_add_tcase(s, tc_core);
  //suite_add_tcase(s, tc_core2);

  tcase_add_test(tc_core, test_create);
  //tcase_add_test(tc_core2, test_create2);
  tcase_add_checked_fixture(tc_core, setup, teardown);
  //tcase_add_checked_fixture(tc_core2, setup, teardown);
  return s;

}
