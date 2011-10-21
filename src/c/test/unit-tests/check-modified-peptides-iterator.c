#include "check-modified-peptides-iterator.h"
#include "ModifiedPeptidesIterator.h"
#include "parameter.h"
#include "peptide_modifications.h"
// from parameter.c
void force_set_aa_mod_list(AA_MOD_T** amod_list, int num_mods);


// declare things to set up
static ModifiedPeptidesIterator *iter1, *iter3;
static PEPTIDE_MOD_T* pmod1;
static SpectrumZState zstate;
static AA_MOD_T *amod1, *amod2, *amod3;
static AA_MOD_T* amod_list[3];
static Database* dbase;

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
  dbase = new Database("input-data/test.fasta", false);
  zstate.setNeutralMass(1268, 1);

  //iter1 = new_modified_peptides_iterator_from_mass(1566, pmod1, NULL, dbase);
  //there seems to be a bug so that the two peptides in this window aren't being returned
  iter1 = new ModifiedPeptidesIterator(1268, zstate, pmod1, 
                                       false, // not decoy
                                       NULL,  // no index
                                       dbase);

  //iter2 = new_modified_peptides_iterator_from_mass(1566, pmod1, index);
}

void mpi_teardown(){
  delete iter1;
  Database::freeDatabase(dbase);
  free_peptide_mod(pmod1);
}

// unmodified peptide mod
START_TEST(test_has_next_unmod){
  // there should be at least one peptide for iter1 as set up
  fail_unless( iter1 != NULL, "Iterator was set up as NULL." );
  fail_unless( iter1->hasNext() == true, 
               "Iterator as set up should have a peptide");

  // get the first peptide
  Peptide* next_p = iter1->next();
  fail_unless( next_p != NULL, "Next returned a null peptide");
  char* seq = next_p->getSequence();
  fail_unless( strcmp(next_p->getSequence(), "ITNHLVAMIEK") == 0,
               "First peptide should be ITNHLVAMIEK");

  fail_unless( iter1->hasNext() == true, 
               "Iterator should have a second peptide");

  delete next_p;
  free(seq);

  next_p = iter1->next();
  fail_unless( next_p != NULL, "Next returned a null second peptide");
  fail_unless( strcmp(next_p->getSequence(), "QGQVATVLSAPAK") == 0,
               "Second peptide should be QGQVATVLSAPAK");
  delete next_p;
  
}
END_TEST

// peptide mod with max 1 mod of one aa
START_TEST(test_has_next_one_mod){
  // create the peptide mod
  aa_mod_set_mass_change(amod1, 10);
  aa_mod_set_max_per_peptide(amod1, 1);
  bool* aas = aa_mod_get_aa_list(amod1);
  aas['Q' - 'A'] = true;
  //aamod1 should have max 1 +10 on Q
  peptide_mod_add_aa_mod(pmod1, 0, 1); // aamod is index 0, 1 copy
  zstate.setNeutralMass(1268+10, 1);
  iter3 = new ModifiedPeptidesIterator(1268+10, zstate, pmod1,
                                                   false, // not decoy
                                                   NULL, // no index
                                                   dbase);

  // test if the iterator has two modified peptides
  fail_unless( iter3->hasNext() == true,
               "Iterator with one Q mod should have first peptide");
  Peptide* pep = iter3->next();
  fail_unless( pep != NULL, "Iterator returned a NULL peptide");
  // check unmodified sequence
  fail_unless( strcmp(pep->getSequence(), "QGQVATVLSAPAK") == 0,
               "First peptide should be QGQVATVLSAPAK");
  // check modified sequence
  char* mod_seq = pep->getModifiedSequenceWithSymbols();
  fail_unless( strcmp(mod_seq, "Q*GQVATVLSAPAK") == 0,
               "First peptide should be Q*GQVATVLSAPAK but is %s", mod_seq);

  // test for second peptide
  fail_unless( iter3->hasNext() == true,
               "Iterator with one Q mod should have second peptide");
  delete pep;
  pep = NULL;
  pep = iter3->next();
  fail_unless( pep != NULL, "Iterator returned a NULL peptide");
  // check unmodified sequence
  fail_unless( strcmp(pep->getSequence(), "QGQVATVLSAPAK") == 0,
               "Second peptide should be QGQVATVLSAPAK");
  // check modified sequence
  mod_seq = pep->getModifiedSequenceWithSymbols();
  fail_unless( strcmp(mod_seq, "QGQ*VATVLSAPAK") == 0,
               "Second peptide should be QGQ*VATVLSAPAK but is %s", mod_seq);

  // check that there are no more modified peptides
 delete pep;
  fail_unless( ! iter3->hasNext(),
               "Iterator should have no more peptides");
  fail_unless( NULL == iter3->next(),
               "Empty iterator should return NULL.");
}
END_TEST

START_TEST(test_all_pep){
  aa_mod_set_mass_change(amod1, -89);
  aa_mod_set_max_per_peptide(amod1, 1);
  bool* aas = aa_mod_get_aa_list(amod1);
  aas['Q' - 'A'] = true;
  peptide_mod_add_aa_mod(pmod1, 0, 1); // aamod is index 0, 1 copy

  iter3 = new ModifiedPeptidesIterator(pmod1, NULL, dbase);

  /*
  while( modified_peptides_iterator_has_next(iter3) ){
    Peptide* pep = modified_peptides_iterator_next(iter3);
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
  bool* aas = aa_mod_get_aa_list(amod2);
  aas['Q' - 'A'] = true;
  aas['V' - 'A'] = true;
  aas['L' - 'A'] = true;
  peptide_mod_add_aa_mod(pmod1, 1, 1); // aamod is index 1, 1 copy
  printf("One copy of @ on Q,V,L\n");
  iter3 = new_modified_peptides_iterator_from_mass(1268, pmod1, NULL, dbase);

  while( modified_peptides_iterator_has_next(iter3)){
    Peptide* pep = modified_peptides_iterator_next(iter3);
    char* mod_seq = modified_aa_string_to_string(get_peptide_modified_sequence(pep));
    printf("%s\n", mod_seq);
  }

  printf("One copy of @ on Q,V,L and two copies of * on T\n");
  aas = aa_mod_get_aa_list(amod1);
  aas['T' - 'A'] = true;
  peptide_mod_add_aa_mod(pmod1, 0, 2); // aamod is index 1, 1 copy
  iter3 = new_modified_peptides_iterator_from_mass(1268, pmod1, NULL, dbase);
  while( modified_peptides_iterator_has_next(iter3)){
    Peptide* pep = modified_peptides_iterator_next(iter3);
    char* mod_seq = modified_aa_string_to_string(pep->getModifiedSequence());
    printf("%s\n", mod_seq);
  }
  */
}
END_TEST

Suite* modified_peptides_iterator_suite(){
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
  tcase_add_checked_fixture(tc_limits, mpi_setup, mpi_teardown);

  return s;
}













