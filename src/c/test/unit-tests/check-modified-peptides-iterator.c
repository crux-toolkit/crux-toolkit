#include "check-modified-peptides-iterator.h"
#include "modified_peptides_iterator.h"
#include "parameter.h"
#include "peptide_modifications.h"
//#include "modifications.h"
// from parameter.c
void force_set_aa_mod_list(AA_MOD_T** amod_list, int num_mods);


// declare things to set up
static MODIFIED_PEPTIDES_ITERATOR_T *iter1, *iter3;
static PEPTIDE_MOD_T* pmod1;
static AA_MOD_T *amod1, *amod2, *amod3;
static AA_MOD_T* amod_list[3];
static DATABASE_T* dbase;

void mpi_setup(){
  initialize_parameters();

  // set up inputs for iterator, use default parameter settings
  amod1 = new_aa_mod(0);
  amod2 = new_aa_mod(1);
  amod3 = new_aa_mod(2);

  amod_list[0] = amod1;
  amod_list[1] = amod2;
  amod_list[2] = amod3;

  force_set_aa_mod_list(amod_list, 3);

  pmod1 = new_peptide_mod();
  dbase = new_database("input-data/test.fasta", FALSE);

  //iter1 = new_modified_peptides_iterator_from_mass(1566, pmod1, NULL, dbase);
  //there seems to be a bug so that the two peptides in this window aren't being returned
  iter1 = new_modified_peptides_iterator_from_mass(1268, pmod1, 
                                                   false, // not decoy
                                                   NULL,  // no index
                                                   dbase);

  //iter2 = new_modified_peptides_iterator_from_mass(1566, pmod1, index);
}

void mpi_teardown(){
  free_modified_peptides_iterator(iter1);
  free_database(dbase);
  free_peptide_mod(pmod1);
}

// unmodified peptide mod
START_TEST(test_has_next_unmod){
  // there should be at least one peptide for iter1 as set up
  fail_unless( iter1 != NULL, "Iterator was set up as NULL." );
  fail_unless( modified_peptides_iterator_has_next(iter1) == TRUE, 
               "Iterator as set up should have a peptide");

  // get the first peptide
  PEPTIDE_T* next_p = modified_peptides_iterator_next(iter1);
  fail_unless( next_p != NULL, "Next returned a null peptide");
  char* seq = get_peptide_sequence(next_p);
  fail_unless( strcmp(get_peptide_sequence(next_p), "QGQVATVLSAPAK") == 0,
               "First peptide should be QGQVATVLSAPAK");

  fail_unless( modified_peptides_iterator_has_next(iter1) == TRUE, 
               "Iterator should have a second peptide");

  free_peptide(next_p);
  free(seq);

  next_p = modified_peptides_iterator_next(iter1);
  fail_unless( next_p != NULL, "Next returned a null second peptide");
  fail_unless( strcmp(get_peptide_sequence(next_p), "ITNHLVAMIEK") == 0,
               "Second peptide should be ITNHLVAMIEK");
  free_peptide(next_p);
  
}
END_TEST

// peptide mod with max 1 mod of one aa
START_TEST(test_has_next_one_mod){
  // create the peptide mod
  aa_mod_set_mass_change(amod1, 10);
  aa_mod_set_max_per_peptide(amod1, 1);
  BOOLEAN_T* aas = aa_mod_get_aa_list(amod1);
  aas['Q' - 'A'] = TRUE;
  //aamod1 should have max 1 +10 on Q
  peptide_mod_add_aa_mod(pmod1, 0, 1); // aamod is index 0, 1 copy
  iter3 = new_modified_peptides_iterator_from_mass(1268+10, pmod1,
                                                   false, // not decoy
                                                   NULL, // no index
                                                   dbase);

  // test if the iterator has two modified peptides
  fail_unless( modified_peptides_iterator_has_next(iter3) == TRUE,
               "Iterator with one Q mod should have first peptide");
  PEPTIDE_T* pep = modified_peptides_iterator_next(iter3);
  fail_unless( pep != NULL, "Iterator returned a NULL peptide");
  // check unmodified sequence
  fail_unless( strcmp(get_peptide_sequence(pep), "QGQVATVLSAPAK") == 0,
               "First peptide should be QGQVATVLSAPAK");
  // check modified sequence
  char* mod_seq = get_peptide_modified_sequence_with_symbols(pep);
  fail_unless( strcmp(mod_seq, "Q*GQVATVLSAPAK") == 0,
               "First peptide should be Q*GQVATVLSAPAK but is %s", mod_seq);

  // test for second peptide
  fail_unless( modified_peptides_iterator_has_next(iter3) == TRUE,
               "Iterator with one Q mod should have second peptide");
  free_peptide(pep);
  pep = NULL;
  pep = modified_peptides_iterator_next(iter3);
  fail_unless( pep != NULL, "Iterator returned a NULL peptide");
  // check unmodified sequence
  fail_unless( strcmp(get_peptide_sequence(pep), "QGQVATVLSAPAK") == 0,
               "Second peptide should be QGQVATVLSAPAK");
  // check modified sequence
  mod_seq = get_peptide_modified_sequence_with_symbols(pep);
  fail_unless( strcmp(mod_seq, "QGQ*VATVLSAPAK") == 0,
               "Second peptide should be QGQ*VATVLSAPAK but is %s", mod_seq);

  // check that there are no more modified peptides
  free_peptide(pep);
  fail_unless( ! modified_peptides_iterator_has_next(iter3),
               "Iterator should have no more peptides");
  fail_unless( NULL == modified_peptides_iterator_next(iter3),
               "Empty iterator should return NULL.");
}
END_TEST

START_TEST(test_all_pep){
  aa_mod_set_mass_change(amod1, -89);
  aa_mod_set_max_per_peptide(amod1, 1);
  BOOLEAN_T* aas = aa_mod_get_aa_list(amod1);
  aas['Q' - 'A'] = TRUE;
  peptide_mod_add_aa_mod(pmod1, 0, 1); // aamod is index 0, 1 copy

  iter3 = new_modified_peptides_iterator(pmod1, NULL, dbase);

  /*
  while( modified_peptides_iterator_has_next(iter3) ){
    PEPTIDE_T* pep = modified_peptides_iterator_next(iter3);
    char* mod_seq = modified_aa_string_to_string(get_peptide_modified_sequence(pep));
    printf("%s\n", mod_seq);
  }
  */
}
END_TEST


// test has next with one mod on first
// with one mod on second
// with multi mods on first, none on second
// with multi mods on second, none on firs
// with multi mods on both

START_TEST(test_looksee){
  // try a number of combinations just to see what is being produced

  // basic steps: select AAs changed by AA_MOD, add AA_MOD(s) to P_MOD
  //              create iterator, iterate over all peptides

  /*
  // one change per peptide on one of three aas
  BOOLEAN_T* aas = aa_mod_get_aa_list(amod2);
  aas['Q' - 'A'] = TRUE;
  aas['V' - 'A'] = TRUE;
  aas['L' - 'A'] = TRUE;
  peptide_mod_add_aa_mod(pmod1, 1, 1); // aamod is index 1, 1 copy
  printf("One copy of @ on Q,V,L\n");
  iter3 = new_modified_peptides_iterator_from_mass(1268, pmod1, NULL, dbase);

  while( modified_peptides_iterator_has_next(iter3)){
    PEPTIDE_T* pep = modified_peptides_iterator_next(iter3);
    char* mod_seq = modified_aa_string_to_string(get_peptide_modified_sequence(pep));
    printf("%s\n", mod_seq);
  }

  printf("One copy of @ on Q,V,L and two copies of * on T\n");
  aas = aa_mod_get_aa_list(amod1);
  aas['T' - 'A'] = TRUE;
  peptide_mod_add_aa_mod(pmod1, 0, 2); // aamod is index 1, 1 copy
  iter3 = new_modified_peptides_iterator_from_mass(1268, pmod1, NULL, dbase);
  while( modified_peptides_iterator_has_next(iter3)){
    PEPTIDE_T* pep = modified_peptides_iterator_next(iter3);
    char* mod_seq = modified_aa_string_to_string(get_peptide_modified_sequence(pep));
    printf("%s\n", mod_seq);
  }
  */
}
END_TEST

START_TEST(test_null){
  MODIFIED_PEPTIDES_ITERATOR_T* null_iter = NULL;
  fail_unless( modified_peptides_iterator_has_next(null_iter) == FALSE,
               "NULL iterator should not have next");
  fail_unless( modified_peptides_iterator_next(null_iter) == NULL,
               "NULL iterator should return NULL for next");
}
END_TEST

Suite* modified_peptides_iterator_suite(){
  //Suite* s = suite_create("Modified-peptides-iterator\n");
  Suite* s = suite_create("Modified-peptides-iterator");
  // Test basic features
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_has_next_unmod);
  tcase_add_test(tc_core, test_has_next_one_mod);
  tcase_add_test(tc_core, test_all_pep);
  tcase_add_test(tc_core, test_looksee);
  tcase_add_checked_fixture(tc_core, mpi_setup, mpi_teardown);

  // Test boundry conditions
  TCase *tc_limits = tcase_create("Limits");
  suite_add_tcase(s, tc_limits);
  tcase_add_test(tc_limits, test_null);
  tcase_add_checked_fixture(tc_limits, mpi_setup, mpi_teardown);

  return s;
}













