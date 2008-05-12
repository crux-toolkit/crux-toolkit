#include <stdlib.h>
#include "check-modifications.h"
#include "../modifications.h"

// declare things to set up
AA_MOD_T *amod1, *amod2, *amod3;
AA_MOD_T* amod_list[3];

void mod_setup(){
  amod1 = new_aa_mod(0);
  amod2 = new_aa_mod(1);
  amod3 = new_aa_mod(2);

  amod_list[0] = amod1;
  amod_list[1] = amod2;
  amod_list[2] = amod3;
}


void mod_teardown(){
  free_aa_mod(amod1); 
  free_aa_mod(amod2); 
  free_aa_mod(amod3); 
}

START_TEST(test_create){
  // check default values

  fail_unless( aa_mod_get_mass_change(amod1) == 0, 
               "amod1 should have had mass change 0 but had %.2f",
               aa_mod_get_mass_change(amod1));

  fail_unless( aa_mod_get_max_per_peptide(amod1) == 0, 
               "amod1 should have had max per peptide 0 but had %d",
               aa_mod_get_max_per_peptide(amod1));

  fail_unless( aa_mod_get_max_distance(amod1) == 40000,//MAX_PROTEIN_SEQ_LENGTH
               "amod1 should have had mass change 40000 but had %d",
               aa_mod_get_max_distance(amod1));

  fail_unless( aa_mod_get_position(amod1) == ANY_POSITION, 
               "amod1 should have had position %d but had %.2f",
               ANY_POSITION, aa_mod_get_position(amod1));

  fail_unless( aa_mod_get_symbol(amod1) == '*', 
               "amod1 should have had symbol * but had %c",
               aa_mod_get_symbol(amod1));
  // FIXME: when these actually get used
  fail_unless( aa_mod_get_identifier(amod1) == 1, 
               "amod1 should have had identifier 1 but had %d",
               aa_mod_get_identifier(amod1));

  fail_unless( aa_mod_get_symbol(amod2) == '@', 
               "amod2 should have had symbol @ but had %c",
               aa_mod_get_symbol(amod2));
  // FIXME: when these actually get used
  fail_unless( aa_mod_get_identifier(amod2) == 2, 
               "amod2 should have had identifier 2 but had %d",
               aa_mod_get_identifier(amod2));
  fail_unless( aa_mod_get_symbol(amod3) == '#', 
               "amod3 should have had symbol # but had %c",
               aa_mod_get_symbol(amod3));
  // FIXME: when these actually get used
  fail_unless( aa_mod_get_identifier(amod3) == 3, 
               "amod3 should have had identifier 3 but had %d",
               aa_mod_get_identifier(amod3));
}
END_TEST

START_TEST(test_set){
  // check default values
  aa_mod_set_mass_change(amod1, 45.6);
  fail_unless( aa_mod_get_mass_change(amod1) == 45.6, 
               "amod1 should have had mass change 45.6 but had %.2f",
               aa_mod_get_mass_change(amod1));

  aa_mod_set_max_per_peptide(amod1, 3);
  fail_unless( aa_mod_get_max_per_peptide(amod1) == 3, 
               "amod1 should have had max per peptide 3 but had %d",
               aa_mod_get_max_per_peptide(amod1));

  aa_mod_set_max_distance(amod1, 1);
  fail_unless( aa_mod_get_max_distance(amod1) == 1,
               "amod1 should have had mass change 1 but had %d",
               aa_mod_get_max_distance(amod1));

  aa_mod_set_max_distance(amod1, -1);
  fail_unless( aa_mod_get_max_distance(amod1) == 40000,//MAX_PROTEIN_SEQ_LENGTH
               "amod1 should have had mass change 40000 but had %d",
               aa_mod_get_max_distance(amod1));

  aa_mod_set_position(amod1, C_TERM);
  fail_unless( aa_mod_get_position(amod1) == C_TERM, 
               "amod1 should have had position %d but had %.2f",
               C_TERM, aa_mod_get_position(amod1));
}
END_TEST

START_TEST(test_char_to_mod){
  char* seq = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  MODIFIED_AA_T* converted = convert_to_mod_aa_seq(seq);

  int i = 0;
  for(i = 0; i < strlen(seq); i++){ 
    fail_unless( converted[i] == i, 
                 "AA %c should have been converted to %d, instead is %d",
                 seq[i], i, converted[i]);
  }
  free(converted);
}
END_TEST

START_TEST(test_copy_mod_seq){
  char* seq = "GBBKATRM"; 
  MODIFIED_AA_T* mod_seq = convert_to_mod_aa_seq(seq);

  MODIFIED_AA_T* copied_seq = copy_mod_aa_seq( mod_seq );
  int i = 0;
  for(i = 0; i<strlen(seq); i++){
    fail_unless( copied_seq[i] == mod_seq[i],
                 "Copied seq %s letter %d should be %d and is %d", 
                 seq, i, mod_seq[i], copied_seq[i]);
  }
  free(mod_seq);
  free(copied_seq);
}
END_TEST

/* Boundry conditions test suite */
START_TEST(test_too_many_mods){
  /* This fails as it should, but doesn't get caught nicely
  int i =0;
  for(i=0; i< MAX_AA_MODS +2; i++){
    AA_MOD_T* m = new_aa_mod(i);
    fail_unless( m != NULL, "Could not create the %dth aa mod", i);
    free(m);
  }
  */
}
END_TEST

  // test_aa_null

Suite* modifications_suite(){
  Suite* s = suite_create("Modifications");
  // Test basic features
  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_create);
  tcase_add_test(tc_core, test_set);
  tcase_add_test(tc_core, test_char_to_mod);
  tcase_add_test(tc_core, test_copy_mod_seq);
  tcase_add_checked_fixture(tc_core, mod_setup, mod_teardown);
  suite_add_tcase(s, tc_core);

  // Test boundry conditions
  TCase *tc_limits = tcase_create("Limits");
  tcase_add_test(tc_limits, test_too_many_mods);
  tcase_add_checked_fixture(tc_limits, mod_setup, mod_teardown);
  suite_add_tcase(s, tc_limits);

  return s;
}













