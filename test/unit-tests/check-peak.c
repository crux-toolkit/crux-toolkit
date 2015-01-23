#include <math.h>
#include <stdlib.h>
#include "check-peak.h"
#include "Peak.h"
#include "crux-utils.h"

Peak* test_peak;

void setup(void){
  test_peak = new Peak(20.7,3.7); //intensity and m/z
}
void teardown(void){
  delete (test_peak);
}

START_TEST (test_create){
  fail_unless( compare_float(test_peak->getIntensity(), 20.7) == 0, "Peak height not setup correctly to 20.7");
  fail_unless( compare_float(test_peak->getLocation(), 3.7) ==0 , "Peak m/z not setup correctly to 3.7");
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
