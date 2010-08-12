#include "check-qvalue.h"
#include "crux-utils.h"
#include "q-value.h"

/**
 * Test if two values are close to each other.  I needed this because
 * "compare_float" was not forgiving enough.
 */ 
static BOOLEAN_T approximately_equal
(FLOAT_T a,
 FLOAT_T b,
 FLOAT_T epsilon)
{
  if (fabs(a - b) < epsilon) {
    return(TRUE);
  }
  return(FALSE);
}

START_TEST (test_create){

  // A simple test of the decoy scoring routine.
  {
    FLOAT_T scores1[] = {10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
    FLOAT_T scores2[] = {5.5, 4.5, 3.5, 2.5, 1.5, 0, -1, -2, -3, -4};
    FLOAT_T* qvalues = compute_decoy_qvalues(scores1, 10, scores2, 10, 1.0);
    fail_unless(compare_float(qvalues[0], 0.0)==0);
    fail_unless(compare_float(qvalues[1], 0.0)==0);
    fail_unless(compare_float(qvalues[2], 0.0)==0);
    fail_unless(compare_float(qvalues[3], 0.0)==0);
    fail_unless(compare_float(qvalues[4], 0.0)==0);
    fail_unless(compare_float(qvalues[5], (1.0 / 6.0))==0);
    fail_unless(compare_float(qvalues[6], (2.0 / 7.0))==0);
    fail_unless(compare_float(qvalues[7], (3.0 / 8.0))==0);
    fail_unless(compare_float(qvalues[8], (4.0 / 9.0))==0);
    fail_unless(compare_float(qvalues[9], (5.0 / 10.0))==0);
    free(qvalues);
  }

  // A simple test of the p-value to q-value conversion.
  // (Cribbed from meme/trunk/src/unit-tests/check-qvalue.c).
  {
    FLOAT_T pvalues[] = {0.01, 0.10, 0.50};
    pvalues[0] = -log(pvalues[0]);
    pvalues[1] = -log(pvalues[1]);
    pvalues[2] = -log(pvalues[2]);
    FLOAT_T* qvalues = compute_qvalues_from_pvalues(pvalues, 3, 1.0);
    // 0.01 * 3 / 1 = 0.03
    fail_unless(approximately_equal(qvalues[0], 0.03, 1e-8));
    // 0.10 * 3 / 2 = 0.15
    fail_unless(approximately_equal(qvalues[1], 0.15, 1e-7));
    // 0.50 * 3 / 3 = 0.50
    fail_unless(approximately_equal(qvalues[2], 0.50, 1e-8));
    free(qvalues);
  }

  // Similar to above, but FDRs are non-monotonic.
  {
    FLOAT_T pvalues[] = {0.02, 0.03, 0.04};
    pvalues[0] = -log(pvalues[0]);
    pvalues[1] = -log(pvalues[1]);
    pvalues[2] = -log(pvalues[2]);
    FLOAT_T* qvalues = compute_qvalues_from_pvalues(pvalues, 3, 1.0);
    // 0.02 * 3 / 1 = 0.06  -> 0.04
    fail_unless(approximately_equal(qvalues[0], 0.04, 1e-8));
    // 0.03 * 3 / 2 = 0.045 -> 0.04
    fail_unless(approximately_equal(qvalues[1], 0.04, 1e-8));
    // 0.04 * 3 / 3 = 0.04  -> 0.04
    fail_unless(approximately_equal(qvalues[2], 0.04, 1e-8));
    free(qvalues);
  }

}
END_TEST

Suite *qvalue_suite(void){
  Suite *s = suite_create("qvalue");
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  return s;
}
