#include <stdlib.h>
#include "check-modifications.h"
#include "../modifications.h"
#include "../parameter.h"
// also in parameter.c
void force_set_aa_mod_list(AA_MOD_T** amod_list, int num_mods);

// declare things to set up
static AA_MOD_T *amod1, *amod2, *amod3;
static AA_MOD_T* amod_list[3];
static MODIFIED_AA_T mod_aa_D, id_3, id_max;

void mod_setup(){
  // assigns identifiers and symbols to each aamod
  amod1 = new_aa_mod(0);
  amod2 = new_aa_mod(1);
  amod3 = new_aa_mod(2);

  amod_list[0] = amod1;
  amod_list[1] = amod2;
  amod_list[2] = amod3;

  mod_aa_D = 3; // D
  id_3   = 0x0080;  // bitmask for 3rd
  id_max = 0x8000;  // bitmask for 11th mod
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

  fail_unless( aa_mod_get_identifier(amod1) == 0x0020, 
               "amod1 should have had identifier 0x0020 but had %d",
               aa_mod_get_identifier(amod1));

  fail_unless( aa_mod_get_symbol(amod2) == '#', 
               "amod2 should have had symbol # but had %c",
               aa_mod_get_symbol(amod2));

  fail_unless( aa_mod_get_identifier(amod2) == 0x0040, 
               "amod2 should have had identifier 0x0040 but had 0x%x",
               aa_mod_get_identifier(amod2));
  fail_unless( aa_mod_get_symbol(amod3) == '@', 
               "amod3 should have had symbol @ but had %c",
               aa_mod_get_symbol(amod3));

  fail_unless( aa_mod_get_identifier(amod3) == 0x0080, 
               "amod3 should have had identifier 0x0080 but had 0x%x",
               aa_mod_get_identifier(amod3));
}
END_TEST

START_TEST(test_set){
  // set values and check 
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

START_TEST(test_is_modified){
  mod_aa_D = mod_aa_D | id_3;

  fail_unless( is_aa_modified(mod_aa_D, amod3),
               "AA should be modified by mod3.");

  AA_MOD_T* amod11 = new_aa_mod(MAX_AA_MODS-1);
  mod_aa_D = mod_aa_D | id_max;
  fail_unless( is_aa_modified(mod_aa_D, amod11),
               "AA should be modified by mod11.");

  fail_unless( ! is_aa_modified(mod_aa_D, amod1),
               "AA should NOT be modified by mod1.");

  free_aa_mod(amod11);
}
END_TEST


START_TEST(test_char_to_mod){
  MODIFIED_AA_T mod_aa = char_aa_to_modified('C');
  fail_unless( mod_aa == 2,  
               "C should have been converted to 2, but was %d", mod_aa);

  const char* seq = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  MODIFIED_AA_T* converted = convert_to_mod_aa_seq(seq);

  int i = 0;
  for(i = 0; i < (int)strlen(seq); i++){ 
    fail_unless( converted[i] == i, 
                 "AA %c should have been converted to %d, instead is %d",
                 seq[i], i, converted[i]);
  }

  free(converted);
}
END_TEST

START_TEST(test_mod_to_char){
  fail_unless( modified_aa_to_char(mod_aa_D) == 'D',
               "mod aa D should convert to D, but it is %c",
               modified_aa_to_char(mod_aa_D));

  const char* seq = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  MODIFIED_AA_T* converted = convert_to_mod_aa_seq(seq);

  // now modify some of the aas and check that it still is converted correctly
  // test modify function later
  converted[3] = converted[3] | 0x8000; // first aa mod in param list
  converted[4] = converted[4] | 0x4000; // second aa mod in param list
  converted[5] = converted[5] | 0x2000; // third aa mod in param list
  converted[6] = converted[6] | 0x0400; // fifth aa mod in param list
  converted[7] = converted[7] | 0x0020; // last aa mod in param list

  int i = 0;
  for(i = 0; i < (int)strlen(seq); i++){ 
    fail_unless( modified_aa_to_char(converted[i]) == seq[i], 
                 "AA %c returned char %c",
                 seq[i], i, modified_aa_to_char(converted[i]));
  }
}
END_TEST

START_TEST(test_mod_to_string){
  // init params including mods
  initialize_parameters();
  force_set_aa_mod_list(amod_list, 3);

  // add test for mod to string (eg D*)
  // unmodified
  char* mod_as_text = modified_aa_to_string( mod_aa_D );
  fail_unless( strcmp(mod_as_text, "D") == 0,
               "mod as string should be 'D' but is '%s'", mod_as_text);
  free(mod_as_text);

  // modify with one
  mod_aa_D = mod_aa_D | 0x0040; // second aa mod in param list
  mod_as_text = modified_aa_to_string( mod_aa_D );
  fail_unless( strcmp(mod_as_text, "D#") == 0,
               "mod as string should be 'D#' but is '%s'", mod_as_text);
  free(mod_as_text);

  mod_aa_D = mod_aa_D | 0x0020; // first aa mod in param list
  mod_aa_D = mod_aa_D | 0x0080; // third aa mod in param list
  mod_as_text = modified_aa_to_string( mod_aa_D );
  fail_unless( strcmp(mod_as_text, "D*#@") == 0,
               "mod as string should be 'D*#@' but is '%s'", mod_as_text);
  free(mod_as_text);
}
END_TEST

START_TEST(test_mod_str_to_string){
  // init params including mods
  initialize_parameters();
  force_set_aa_mod_list(amod_list, 3);

  const char* seq = "GBBKATRM"; 
  int len = strlen(seq);
  MODIFIED_AA_T* modaa_seq = convert_to_mod_aa_seq(seq);

  char* modchar_seq = modified_aa_string_to_string(modaa_seq, len);
  fail_unless( strcmp(modchar_seq, seq) == 0,
               "Seq %s should equal %s", seq, modchar_seq);

  fail_unless( aa_mod_get_symbol(amod1) == '*', 
               "amod1 should have had symbol * but had %c",
               aa_mod_get_symbol(amod1));

  fail_unless( aa_mod_get_symbol(amod2) == '#', 
               "amod1 should have had symbol # but had %c",
               aa_mod_get_symbol(amod2));

  fail_unless( aa_mod_get_symbol(amod3) == '@', 
               "amod1 should have had symbol @ but had %c",
               aa_mod_get_symbol(amod3));

  // now modify a couple of aas
  modify_aa( &modaa_seq[1], amod2 );
  modify_aa( &modaa_seq[4], amod1 );
  modify_aa( &modaa_seq[4], amod3 );
  modify_aa( &modaa_seq[7], amod2 );
  modify_aa( &modaa_seq[7], amod3 );
  free(modchar_seq);
  modchar_seq = modified_aa_string_to_string(modaa_seq, len);

  fail_unless( strcmp(modchar_seq, "GB#BKA*@TRM#@") == 0,
               "Seq %s should equal %s", seq, modchar_seq);

  free(modchar_seq);
}
END_TEST

START_TEST(test_mod_str_to_unmod_string){
  // init params including mods
  initialize_parameters();
  force_set_aa_mod_list(amod_list, 3);

  const char* seq = "GBBKATRM"; 
  int len = strlen(seq);
  MODIFIED_AA_T* modaa_seq = convert_to_mod_aa_seq(seq);

  // change it back to a seq
  char* unmod_seq = modified_aa_to_unmodified_string(modaa_seq, len);
  fail_unless(strcmp(unmod_seq, seq) == 0,
              "The unmodified version of the mod_aa* is %s but should be %s",
              unmod_seq, seq);
  free(unmod_seq);

  // modify it and change back to a seq
  modify_aa( &modaa_seq[1], amod2 );
  modify_aa( &modaa_seq[4], amod1 );
  modify_aa( &modaa_seq[4], amod3 );
  modify_aa( &modaa_seq[7], amod2 );
  modify_aa( &modaa_seq[7], amod3 );

  unmod_seq = modified_aa_to_unmodified_string(modaa_seq, len);
  fail_unless(strcmp(unmod_seq, seq) == 0,
              "The unmodified version of the mod_aa* is %s but should be %s",
              unmod_seq, seq);
  free(unmod_seq);
  
}
END_TEST

START_TEST(test_copy_mod_seq){
  const char* seq = "GBBKATRM"; 
  int len = strlen(seq);
  MODIFIED_AA_T* mod_seq = convert_to_mod_aa_seq(seq);

  MODIFIED_AA_T* copied_seq = copy_mod_aa_seq( mod_seq, len );
  int i = 0;
  for(i = 0; i<len; i++){
    fail_unless( copied_seq[i] == mod_seq[i],
                 "Copied seq %s letter %d should be %d and is %d", 
                 seq, i, mod_seq[i], copied_seq[i]);
  }
  free(mod_seq);
  free(copied_seq);
}
END_TEST

START_TEST(test_palindrome){
  const char* seq = "GBBKATRM"; 
  int len = strlen(seq);
  MODIFIED_AA_T* mod_seq = convert_to_mod_aa_seq(seq);

  fail_unless( modified_aa_seq_is_palindrome(mod_seq, len) == FALSE,
               "The seq %s should not be a palindrome.", seq);

  // try a palindrome
  seq = "MACDDCAR";
  len = strlen(seq);
  free(mod_seq);
  mod_seq = convert_to_mod_aa_seq(seq);
  fail_unless( modified_aa_seq_is_palindrome(mod_seq, len) == TRUE,
               "The seq %s should be a palindrome.", seq);

  // try a palindrome of odd length
  seq = "MACDVDCAR";
  len = strlen(seq);
  free(mod_seq);
  mod_seq = convert_to_mod_aa_seq(seq);
  fail_unless( modified_aa_seq_is_palindrome(mod_seq, len) == TRUE,
               "The seq %s should be a palindrome.", seq);

  // modify half making it not a palindrome
  force_set_aa_mod_list(amod_list, 3);  // in setup???
  modify_aa(&mod_seq[2], amod1);
  fail_unless( modified_aa_seq_is_palindrome(mod_seq, len) == FALSE,
               "The modified seq %s should not be a palindrome.", 
               modified_aa_string_to_string(mod_seq, len));

  // modify the other half returning it to a palindrome
  modify_aa(&mod_seq[6], amod1);
  fail_unless( modified_aa_seq_is_palindrome(mod_seq, len) == TRUE,
               "The modified seq %s should be a palindrome.", 
               modified_aa_string_to_string(mod_seq, len));


  // try an almost-palindrome
  seq = "MACDVCAR";
  len = strlen(seq);
  free(mod_seq);
  mod_seq = convert_to_mod_aa_seq(seq);
  fail_unless( modified_aa_seq_is_palindrome(mod_seq, len) == FALSE,
               "The seq %s should not be a palindrome.", seq);

  // try a short seq
  seq = "MAR";
  len = strlen(seq);
  free(mod_seq);
  mod_seq = convert_to_mod_aa_seq(seq);
  fail_unless( modified_aa_seq_is_palindrome(mod_seq, len) == TRUE,
               "The seq %s should not be a palindrome.", seq);

}
END_TEST

START_TEST(test_is_modifiable){
  BOOLEAN_T* mod_us = aa_mod_get_aa_list(amod3);
  mod_us['D' - 'A'] = TRUE;
  fail_unless( is_aa_modifiable(mod_aa_D, amod3),
               "aa_D should be modifiable by amod3.");

  mod_aa_D = mod_aa_D | id_3;
  fail_unless( !is_aa_modifiable(mod_aa_D, amod3),
               "aa_D should not be modifiable AGAIN by amod3.");
}
END_TEST

START_TEST(test_modify){
  // it starts out unmodified
  fail_unless( !is_aa_modified(mod_aa_D, amod3),
               "mod_aa_D should not be modified.");
  modify_aa(&mod_aa_D, amod3);
  fail_unless( is_aa_modified(mod_aa_D, amod3),
               "mod_aa_D should be modified.");

  // and it should no longer be modifiable by amod3
  fail_unless( !is_aa_modifiable(mod_aa_D, amod3),
               "mod_aa_D should no longer be modifiable by amod3");
}
END_TEST

START_TEST(test_serialize){
  initialize_parameters();
  // create a mod with non-zero values for all fields
  aa_mod_set_mass_change(amod1, 87.9);
  aa_mod_set_max_per_peptide(amod1, 3);
  BOOLEAN_T* aa_list = aa_mod_get_aa_list(amod1);
  aa_list['S' - 'A'] = TRUE;
  aa_list['T' - 'A'] = TRUE;
  aa_list['Y' - 'A'] = TRUE;

  // write a mod to file
  FILE* file = fopen("amod.bin", "w");
  fail_unless( serialize_aa_mod(amod1, file),
               "Error in serializing amod1 to file.");
  // close the file and read it back in
  fclose(file);
  file = NULL;
  file = fopen("amod.bin", "r");
  fail_unless( parse_aa_mod(amod2, file),
               "Error in parsing amod1 from file.");
  // test that compare() calls them the same

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
  tcase_add_test(tc_core, test_is_modified);
  tcase_add_test(tc_core, test_char_to_mod);
  tcase_add_test(tc_core, test_mod_to_char);
  tcase_add_test(tc_core, test_mod_to_string);
  tcase_add_test(tc_core, test_mod_str_to_string);
  tcase_add_test(tc_core, test_mod_str_to_unmod_string);
  tcase_add_test(tc_core, test_copy_mod_seq);
  tcase_add_test(tc_core, test_palindrome);
  tcase_add_test(tc_core, test_is_modifiable);
  tcase_add_test(tc_core, test_modify);
  tcase_add_test(tc_core, test_serialize);

  tcase_add_checked_fixture(tc_core, mod_setup, mod_teardown);
  suite_add_tcase(s, tc_core);

  // Test boundry conditions
  TCase *tc_limits = tcase_create("Limits");
  tcase_add_test(tc_limits, test_too_many_mods);
  tcase_add_checked_fixture(tc_limits, mod_setup, mod_teardown);
  suite_add_tcase(s, tc_limits);

  return s;
}













