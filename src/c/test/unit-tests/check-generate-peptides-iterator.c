#include "check-generate-peptides-iterator.h"
#include "parameter.h"
#include "GeneratePeptidesIterator.h"

using namespace std;
// private testing functions in GeneratePeptidesIterator.c


// declare things to set up
GeneratePeptidesIterator* gpiter1;
Database* dbase;

void gpiter_setup(){
  initialize_parameters();
  dbase = new Database("input-data/test.fasta", false);
  gpiter1 = new GeneratePeptidesIterator(
                       pair<FLOAT_T,FLOAT_T>(1268 - 3, 1268 + 3),
                       false, dbase, NULL);
}

void gpiter_teardown(){
  Database::freeDatabase(dbase);
  delete gpiter1;
}

START_TEST(test_create){
  fail_unless( gpiter1 != NULL,
               "Got NULL instead of peptide iterator.");


  fail_unless( gpiter1->hasNext() == true,
               "Iterator should have peptide as set up");

  Peptide* next_p = gpiter1->next();
  fail_unless( next_p != NULL, "Next returned a null peptide");
  char* seq = next_p->getSequence();
  fail_unless( strcmp(seq, "ITNHLVAMIEK") == 0,
               "First peptide should be ITNHLVAMIEK but is %s", seq);
  fail_unless( get_double_parameter("precursor-window") == 3,
               "Default mass window should be 3.");
  fail_unless( gpiter1->hasNext() == true,
               "Iterator should have second peptide");

  next_p = gpiter1->next();
  free(seq);
  seq = next_p->getSequence();
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













