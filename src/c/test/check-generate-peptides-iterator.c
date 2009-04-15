#include "check-generate-peptides-iterator.h"
#include "../parameter.h"
#include "../generate_peptides_iterator.h"

// private testing functions in generate_peptides_iterator.c


// declare things to set up
GENERATE_PEPTIDES_ITERATOR_T* gpiter1;
DATABASE_T* dbase;

void gpiter_setup(){
  initialize_parameters();
  dbase = new_database("input-data/test.fasta", FALSE);
  gpiter1 = new_generate_peptides_iterator_from_mass(1268, NULL, dbase);
  //gpiter1 = new_generate_peptides_iterator_from_mass(1566, index, NULL);
}

void gpiter_teardown(){
  free_database(dbase);
  free_generate_peptides_iterator(gpiter1);
}

START_TEST(test_create){
  fail_unless( gpiter1 != NULL,
               "Got NULL instead of peptide iterator.");


  fail_unless( generate_peptides_iterator_has_next(gpiter1) == TRUE,
               "Iterator should have peptide as set up");

  PEPTIDE_T* next_p = generate_peptides_iterator_next(gpiter1);
  fail_unless( next_p != NULL, "Next returned a null peptide");
  char* seq = get_peptide_sequence(next_p);
  fail_unless( strcmp(seq, "ITNHLVAMIEK") == 0,
               "First peptide should be ITNHLVAMIEK but is %s", seq);
  fail_unless( get_double_parameter("mass-window") == 3,
               "Default mass window should be 3.");
  fail_unless( generate_peptides_iterator_has_next(gpiter1) == TRUE,
               "Iterator should have second peptide");

  next_p = generate_peptides_iterator_next(gpiter1);
  free(seq);
  seq = get_peptide_sequence(next_p);
  fail_unless( strcmp(seq, "QGQVATVLSAPAK") == 0,
               "Second peptide should be QGQVATVLSAPAK but is %s", seq);
  free(seq);
}
END_TEST

// tests for a range with no peptides, for a larger range, 
// with different constraints

START_TEST(test_somethingelse){
}
END_TEST

Suite* generate_peptides_iterator_suite(){
  Suite* s = suite_create("Generate-peptides-iterator");
  // Test basic features
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  tcase_add_checked_fixture(tc_core, gpiter_setup, gpiter_teardown);

  // Test boundry conditions
  TCase *tc_limits = tcase_create("Limits");
  suite_add_tcase(s, tc_limits);
  tcase_add_test(tc_limits, test_somethingelse);
  tcase_add_checked_fixture(tc_limits, gpiter_setup, gpiter_teardown);

  return s;
}













