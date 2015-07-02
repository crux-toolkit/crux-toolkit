#include <stdlib.h>
#include "check-serialize.h"
#include "modifications.h"
#include "parameter.h"
#include "database.h"
#include "protein.h"
#include "peptide.h"
#include "match_collection.h"

// also from parameter.c
void force_set_aa_mod_list(AA_MOD_T** amod_list, int num_mods);

// declare things to set up
static char pep_filename[128] = "test-peptides.bin";
static PEPTIDE_T *pep1, *pep2, *pep3, *pep4;
static PROTEIN_T *protein1;
static DATABASE_T* db;
static AA_MOD_T *amod1, *amod2;
static AA_MOD_T* amod_list[2];
static PEPTIDE_MOD_T *pmod1, *pmod2;

void ser_setup(){
  db = new_database("input-data/protseq1.fasta", FALSE);
  parse_database(db);

  protein1 = get_database_protein_at_idx(db, 0);
  pep1 = new_peptide( 10, 1087.20, protein1, 20);//VADILESNAR
  pep2 = new_peptide( 16, 1736.02, protein1, 1);//MRVLKFGGTSVANAER
  pep3 = new_peptide( 4, 547.57, protein1, 487);//QQSW

  amod1 = new_aa_mod(0);
  amod2 = new_aa_mod(1);

  amod_list[0] = amod1;
  amod_list[1] = amod2;

  aa_mod_set_mass_change(amod1, 100);
  aa_mod_set_mass_change(amod2, 200);

  // initialize mods in parameter.c
  initialize_parameters();
  force_set_aa_mod_list(amod_list, 2);

  pmod1 = new_peptide_mod();
  pmod2 = new_peptide_mod();
  peptide_mod_add_aa_mod(pmod1, 0, 2); // add 1cpy of amod1
  peptide_mod_add_aa_mod(pmod2, 0, 1); // add 1cpy of amod1
  peptide_mod_add_aa_mod(pmod2, 1, 1); // add 1cpy of amod2

}

void ser_teardown(){
  free_database(db);
  free_peptide(pep1);
  free_peptide(pep2);
  free_peptide(pep3);
}

START_TEST(test_setup){
  fail_unless( TRUE == get_database_is_parsed(db),
               "Failed to parse database");
  char* seq = get_peptide_sequence(pep1);
  fail_unless(strcmp("VADILESNAR", seq)==0,
              "peptide1 should be initialized to VADILESNAR but is %s", seq ); 
  free(seq);
  seq = get_peptide_sequence(pep2);
  fail_unless(strcmp("MRVLKFGGTSVANAER", seq)==0,
              "peptide2 should be initialized to MRVLKFGGTSVANAER but is %s",
              seq);
  free(seq);
  seq = get_peptide_sequence(pep3);
  fail_unless(strcmp("QQSW", seq)==0,
              "peptide3 should be initialized to QQSW but is %s", seq); 
}
END_TEST

START_TEST(test_serialize){
  // open file for writing
  FILE* file = fopen(pep_filename, "w");
  // serialize three peptides
  fail_unless( serialize_peptide(pep1, file, NULL) == TRUE,
               "Failed to serialize peptide1");
  fail_unless( serialize_peptide(pep2, file, NULL) == TRUE,
               "Failed to serialize peptide2");
  fail_unless( serialize_peptide(pep3, file, NULL) == TRUE,
               "Failed to serialize peptide3");
  // close file
  fclose(file);

  // open file for reading
  file = fopen(pep_filename, "r");

  // read in three peptides
  pep4 = parse_peptide(file, db, FALSE); // don't use src array
  // test length, mass, seq
  int len = get_peptide_length(pep4);
  double mass = get_peptide_peptide_mass(pep4);
  char* seq = get_peptide_sequence(pep4);
  fail_unless( len == 10, 
               "First peptide length is %d but should be 10", len );
  fail_unless( (int)mass == 1087, //1087.2, comes back as .199999
               "First peptide mass is %f but should be 1087.20", mass );
  fail_unless( strcmp(seq, "VADILESNAR") == 0, 
               "First peptide length is %s but should be VADILESNAR", seq );
  free(seq);

  pep4 = parse_peptide(file, db, FALSE); // don't use src array
  len = get_peptide_length(pep4);
  mass = get_peptide_peptide_mass(pep4);
  seq = get_peptide_sequence(pep4);
  fail_unless( len == 16, 
               "Second peptide length is %d but should be 16", len );
  fail_unless( (int)mass == 1736, // comes back 1736.020020
               "Second peptide mass is %f but should be 1736.02", mass );
  fail_unless( strcmp(seq, "MRVLKFGGTSVANAER") == 0, 
               "Second peptide length is %s but should be MRVLKFGGTSVANAER", 
               seq );
  free(seq);

  pep4 = parse_peptide(file, db, FALSE); // don't use src array
  len = get_peptide_length(pep4);
  mass = get_peptide_peptide_mass(pep4);
  seq = get_peptide_sequence(pep4);
  fail_unless( len == 4, 
               "Third peptide length is %d but should be 4", len );
  fail_unless( (int)mass == 547, // comes back 547.570007
               "Third peptide mass is %f but should be 547.57", mass );
  fail_unless( strcmp(seq, "QQSW") == 0, 
               "Third peptide length is %s but should be QQSW", seq );
  // close file

}
END_TEST

// !!! add a peptide that starts with A
START_TEST(test_serialize_mod){
  // modify peptides 1 and 3

  // create mod seq
  char* pep_seq = get_peptide_sequence(pep1);
  MODIFIED_AA_T* mod_seq = NULL;
  convert_to_mod_aa_seq(pep_seq, &mod_seq);
  modify_aa(&mod_seq[2], amod1);
  modify_aa(&mod_seq[9], amod1);
  // add mod to peptide
  set_peptide_mod(pep1, mod_seq, pmod1);
  // check mod seq
  char* pep1_mod_str = get_peptide_modified_sequence_with_symbols(pep1);
  fail_unless( strcmp(pep1_mod_str, "VAD*ILESNAR*") == 0,
               "Modified sequence of pep1 should be VAD*ILESNAR* but is %s",
               pep1_mod_str);

  // repeat for pep3
  free(pep_seq);
  pep_seq = get_peptide_sequence(pep3);
  convert_to_mod_aa_seq(pep_seq, &mod_seq);
  modify_aa(&mod_seq[0], amod1);
  modify_aa(&mod_seq[1], amod2);
  // add mod to peptide
  set_peptide_mod(pep3, mod_seq, pmod2);
  // check mod seq
  char* pep3_mod_str = get_peptide_modified_sequence_with_symbols(pep3);
  fail_unless( strcmp(pep3_mod_str, "Q*Q#SW") == 0,
               "Modified sequence of pep3 should be Q*Q#SW but is %s",
               pep3_mod_str);

  // serialize all three peptides
  FILE* file = fopen(pep_filename, "w");
  // serialize three peptides
  fail_unless( serialize_peptide(pep1, file, NULL) == TRUE,
               "Failed to serialize peptide1");
  fail_unless( serialize_peptide(pep2, file, NULL) == TRUE,
               "Failed to serialize peptide2");
  fail_unless( serialize_peptide(pep3, file, NULL) == TRUE,
               "Failed to serialize peptide3");
  // close file
  fclose(file);

  // read in three peptides
  file = fopen(pep_filename, "r");

  pep4 = parse_peptide(file, db, FALSE); // don't use src array
  char* pep4_seq = get_peptide_modified_sequence_with_symbols(pep4);
  fail_unless( strcmp(pep4_seq, pep1_mod_str) == 0,
               "First peptide sequence is %s but should be %s",
               pep4_seq, pep1_mod_str);

  pep4 = parse_peptide(file, db, FALSE); // don't use src array
  pep4_seq = get_peptide_modified_sequence_with_symbols(pep4);
  char* pep2_seq = get_peptide_modified_sequence_with_symbols(pep2);
  fail_unless( strcmp(pep4_seq, pep2_seq) == 0,
               "Second peptide sequence is %s but should be %s",
               pep4_seq, pep2_seq);

  pep4 = parse_peptide(file, db, FALSE); // don't use src array
  pep4_seq = get_peptide_modified_sequence_with_symbols(pep4);
  fail_unless( strcmp(pep4_seq, pep3_mod_str) == 0,
               "Third peptide sequence is %s but should be %s",
               pep4_seq, pep3_mod_str);

  fclose(file);
}
END_TEST


Suite* serialize_suite(){
  Suite* s = suite_create("Serialization");
  // Test basic features
  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_setup);
  tcase_add_test(tc_core, test_serialize);
  tcase_add_test(tc_core, test_serialize_mod);

  tcase_add_checked_fixture(tc_core, ser_setup, ser_teardown);
  suite_add_tcase(s, tc_core);

  // Test boundry conditions
  /*
  TCase *tc_limits = tcase_create("Limits");
  tcase_add_test(tc_limits, test_too_many_mods);
  tcase_add_checked_fixture(tc_limits, mod_setup, mod_teardown);
  suite_add_tcase(s, tc_limits);
  */
  return s;
}













