#include <math.h>
#include <stdlib.h>
#include "check-peak.h"
#include "../peak.h"
#include "../crux-utils.h"

PEAK_T * test_peak;

void setup(void){
  test_peak = new_peak(20.7,3.7); //intensity and m/z
}
void teardown(void){
  free_peak(test_peak);
}

START_TEST (test_create){
  fail_unless( compare_float(get_peak_intensity(test_peak), 20.7) == 0, "Peak height not setup correctly to 20.7");
  fail_unless( compare_float(get_peak_location(test_peak), 3.7) ==0 , "Peak m/z not setup correctly to 3.7");
}
END_TEST

Suite *peak_suite(void){
  Suite *s = suite_create("Peak");
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  tcase_add_checked_fixture(tc_core, setup, teardown);
  return s;
}
