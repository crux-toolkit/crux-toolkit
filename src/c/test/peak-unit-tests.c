#include <stdlib.h>
#include <check.h>
#include "check-peak.h"

int main(void){
  int nf;
  Suite *s = peak_suite();
  SRunner *sr = srunner_create(s);
  srunner_run_all(sr, CK_NORMAL);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (nf == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
