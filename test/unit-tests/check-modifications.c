#include <stdlib.h>
#include "check-modifications.h"
#include "modifications.h"
#include "parameter.h"
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

// Created three AA_MOD_T's in setup.  Test their default values via getters
START_TEST(test_create){

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

// Test the setters for AA_MOD_T
START_TEST(test_set){
  // mass
  aa_mod_set_mass_change(amod1, 45.6);
  fail_unless( aa_mod_get_mass_change(amod1) == 45.6, 
               "amod1 should have had mass change 45.6 but had %.2f",
               aa_mod_get_mass_change(amod1));
  // max
  aa_mod_set_max_per_peptide(amod1, 3);
  fail_unless( aa_mod_get_max_per_peptide(amod1) == 3, 
               "amod1 should have had max per peptide 3 but had %d",
               aa_mod_get_max_per_peptide(amod1));
  // distance
  aa_mod_set_max_distance(amod1, 1);
  fail_unless( aa_mod_get_max_distance(amod1) == 1,
               "amod1 should have had mass change 1 but had %d",
               aa_mod_get_max_distance(amod1));

  aa_mod_set_max_distance(amod1, -1);
  fail_unless( aa_mod_get_max_distance(amod1) == 40000,//MAX_PROTEIN_SEQ_LENGTH
               "amod1 should have had mass change 40000 but had %d",
               aa_mod_get_max_distance(amod1));
  // position
  aa_mod_set_position(amod1, C_TERM);
  fail_unless( aa_mod_get_position(amod1) == C_TERM, 
               "amod1 should have had position %d but had %.2f",
               C_TERM, aa_mod_get_position(amod1));
}
END_TEST

// Test the is_aa_modifieable() function
START_TEST(test_is_modifiable){
  bool* mod_us = aa_mod_get_aa_list(amod3);
  mod_us['D' - 'A'] = true;
  fail_unless( is_aa_modifiable(mod_aa_D, amod3),
               "aa_D should be modifiable by amod3.");

  mod_aa_D = mod_aa_D | id_3;
  fail_unless( !is_aa_modifiable(mod_aa_D, amod3),
               "aa_D should not be modifiable AGAIN by amod3.");
}
END_TEST

// Test the is_aa_modified() function
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

// modify and AA_T
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

// Conversion tests between MODIFIED_AA_T and char, both directions

// Convert a single char -> MODIFIED_AA_T
// Convert a char* (no mods) -> MODIFIED_AA_T*
START_TEST(test_char_to_mod){
  initialize_parameters();
  MODIFIED_AA_T mod_aa = char_aa_to_modified('C');
  fail_unless( mod_aa == 2,  
               "C should have been converted to 2, but was %d", mod_aa);

  const char* seq = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  MODIFIED_AA_T* converted = NULL;
  int len = convert_to_mod_aa_seq(seq, &converted);
  fail_unless(len == (int)strlen(seq), 
              "MOD_AA_T should be same length as seq, %d, but is %d.", 
              strlen(seq), len);

  int i = 0;
  for(i = 0; i < (int)strlen(seq); i++){ 
    fail_unless( converted[i] == i, 
                 "AA %c should have been converted to %d, instead is %d",
                 seq[i], i, converted[i]);
  }

  free(converted);
}
END_TEST

// Convert a single MODIFIED_AA_T -> char (with no mod symbols) 
START_TEST(test_mod_to_char){
  initialize_parameters();
  fail_unless( modified_aa_to_char(mod_aa_D) == 'D',
               "mod aa D should convert to D, but it is %c",
               modified_aa_to_char(mod_aa_D));

  const char* seq = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  MODIFIED_AA_T* converted = NULL;
  convert_to_mod_aa_seq(seq, &converted);

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

// Convert a single MODIFIED_AA_T to char* with symbols or masses
START_TEST(test_mod_to_string){
  // init params including mods
  initialize_parameters();
  force_set_aa_mod_list(amod_list, 3);
  // include mass shifts
  aa_mod_set_mass_change(amod1, 3.4);
  aa_mod_set_mass_change(amod2, 56.78);
  aa_mod_set_mass_change(amod3, -100);
  int precision = 2;

  // add test for mod to string (eg D*)
  // unmodified
  char* mod_as_text = modified_aa_to_string_with_symbols( mod_aa_D );
  fail_unless( strcmp(mod_as_text, "D") == 0,
               "mod as string should be 'D' but is '%s'", mod_as_text);
  free(mod_as_text);

  mod_as_text = modified_aa_to_string_with_masses(mod_aa_D, MOD_MASS_ONLY,
                                                  precision);
  fail_unless( strcmp(mod_as_text, "D") == 0,
               "mod as string should be 'D' but is '%s'", mod_as_text);
  free(mod_as_text);

  // modify with one
  mod_aa_D = mod_aa_D | 0x0040; // second aa mod in param list
  mod_as_text = modified_aa_to_string_with_symbols( mod_aa_D );
  fail_unless( strcmp(mod_as_text, "D#") == 0,
               "mod as string should be 'D#' but is '%s'", mod_as_text);
  free(mod_as_text);

  mod_as_text = modified_aa_to_string_with_masses(mod_aa_D, MOD_MASS_ONLY, 
                                                  precision);
  fail_unless( strcmp(mod_as_text, "D[56.78]") == 0,
               "mod as string should be 'D[56.78]' but is '%s'", mod_as_text);
  free(mod_as_text);

  mod_aa_D = mod_aa_D | 0x0020; // first aa mod in param list
  mod_aa_D = mod_aa_D | 0x0080; // third aa mod in param list
  mod_as_text = modified_aa_to_string_with_symbols( mod_aa_D );
  fail_unless( strcmp(mod_as_text, "D*#@") == 0,
               "mod as string should be 'D*#@' but is '%s'", mod_as_text);
  free(mod_as_text);

  mod_as_text = modified_aa_to_string_with_masses( mod_aa_D, MOD_MASS_ONLY,
                                                   precision );
  fail_unless( strcmp(mod_as_text, "D[-39.82]") == 0,
               "mod as string should be 'D[-39.82]' but is '%s'", mod_as_text);
  free(mod_as_text);

  mod_as_text = modified_aa_to_string_with_masses( mod_aa_D, 
                                                   MOD_MASSES_SEPARATE, 
                                                   precision );
  fail_unless( strcmp(mod_as_text, "D[3.40,56.78,-100.00]") == 0,
               "mod as string should be 'D[3.40,56.78,-100.00]' but is '%s'", 
               mod_as_text);
  free(mod_as_text);

  mod_as_text = modified_aa_to_string_with_masses( mod_aa_D, AA_PLUS_MOD,
                                                   precision );
  fail_unless( strcmp(mod_as_text, "D[75.27]") == 0,
               "mod as string should be 'D[75.27]' but is '%s'", 
               mod_as_text);
  free(mod_as_text);
}
END_TEST

// similar to above test but with different precision
START_TEST(test_mod_to_string_precision){
  initialize_parameters();
  force_set_aa_mod_list(amod_list, 3);
  // set mass shifts
  aa_mod_set_mass_change(amod2, 56.785);
  // modify an aa
  mod_aa_D = mod_aa_D | 0x0040; // second aa mod in param list

  char* mod_as_text = modified_aa_to_string_with_masses( mod_aa_D, MOD_MASS_ONLY, 3 );
  fail_unless( strcmp(mod_as_text, "D[56.785]") == 0,
               "precision 3 mod should be D[56.785] but is %s.", mod_as_text);
  mod_as_text = modified_aa_to_string_with_masses( mod_aa_D, MOD_MASS_ONLY, 4 );
  fail_unless( strcmp(mod_as_text, "D[56.7850]") == 0,
               "precision 3 mod should be D[56.7850] but is %s.", mod_as_text);
  mod_as_text = modified_aa_to_string_with_masses( mod_aa_D, MOD_MASS_ONLY, 2 );
  fail_unless( strcmp(mod_as_text, "D[56.78]") == 0,
               "precision 3 mod should be D[56.78] but is %s.", mod_as_text);

}
END_TEST

// Convert MODIFIED_AA_T* to char* with symbols or masses
START_TEST(test_mod_str_to_string){
  // init params including mods
  initialize_parameters();
  force_set_aa_mod_list(amod_list, 3);
  // include mass shifts
  aa_mod_set_mass_change(amod1, 3.4);
  aa_mod_set_mass_change(amod2, 56.78);
  aa_mod_set_mass_change(amod3, -100);

  const char* seq = "GGBKATRM"; 
  int len = strlen(seq);
  MODIFIED_AA_T* modaa_seq = NULL;
  convert_to_mod_aa_seq(seq, &modaa_seq);

  char* modchar_seq = modified_aa_string_to_string_with_symbols(modaa_seq, len);
  fail_unless( strcmp(modchar_seq, seq) == 0,
               "Seq %s should equal %s", seq, modchar_seq);

  free(modchar_seq);
  modchar_seq = modified_aa_string_to_string_with_masses(modaa_seq, len, 
                                                         MOD_MASS_ONLY);
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
  modchar_seq = modified_aa_string_to_string_with_symbols(modaa_seq, len);

  fail_unless( strcmp(modchar_seq, "GG#BKA*@TRM#@") == 0,
               "Seq %s should equal %s", modchar_seq, "GG#BKA*@TRM#@");

  free(modchar_seq);
  modchar_seq = modified_aa_string_to_string_with_masses(modaa_seq, len, 
                                                         MOD_MASSES_SEPARATE);
  fail_unless( strcmp(modchar_seq, 
                      "GG[56.78]BKA[3.40,-100.00]TRM[56.78,-100.00]") == 0,
               "Seq %s should equal %s", modchar_seq, 
               "GG[56.78]BKA[3.40,-100.00]TRM[56.78,-100.00]");

  free(modchar_seq);
  modchar_seq = modified_aa_string_to_string_with_masses(modaa_seq, len, 
                                                         MOD_MASS_ONLY);
  fail_unless( strcmp(modchar_seq, "GG[56.78]BKA[-96.60]TRM[-43.22]") == 0,
               "Seq %s should equal %s", modchar_seq,
               "GG[56.78]BKA[-96.60]TRM[-43.22]");

  free(modchar_seq);
  modchar_seq = modified_aa_string_to_string_with_masses(modaa_seq, len, 
                                                         AA_PLUS_MOD);
  fail_unless( strcmp(modchar_seq, "GG[113.83]BKA[-25.52]TRM[87.97]") == 0,
               "Seq %s should equal %s", modchar_seq,
               "GG[113.83]BKA[-25.52]TRM[87.97]");


  free(modchar_seq);
}
END_TEST

// Convert MODIFIED_AA_T* to a char* with no mods, even if MOD_AA_T is modified
START_TEST(test_mod_str_to_unmod_string){
  // init params including mods
  initialize_parameters();
  force_set_aa_mod_list(amod_list, 3);

  const char* seq = "GBBKATRM"; 
  int len = strlen(seq);
  MODIFIED_AA_T* modaa_seq = NULL;
  convert_to_mod_aa_seq(seq, &modaa_seq);

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

// Convert a mod symbol to an AA_MOD_T
START_TEST(test_symbol_to_aa_mod){
  initialize_parameters();
  force_set_aa_mod_list(amod_list, 3);

  const AA_MOD_T* aa = get_aa_mod_from_symbol('*');
  fail_unless( aa == amod1,
              "The symbol * is not returning the correct AA_MOD_T.  "
               "The returned AA_MOD_T has %c.", aa_mod_get_symbol(aa));
  aa = get_aa_mod_from_symbol('#');
  fail_unless( aa == amod2,
              "The symbol # is not returning the correct AA_MOD_T.  "
               "The returned AA_MOD_T has %c.", aa_mod_get_symbol(aa));
   aa = get_aa_mod_from_symbol('@');
  fail_unless( aa == amod3,
              "The symbol @ is not returning the correct AA_MOD_T.  "
               "The returned AA_MOD_T has %c.", aa_mod_get_symbol(aa));
}
END_TEST

// Convert a mod mass to an AA_MOD_T
START_TEST(test_mass_to_aa_mod){
  initialize_parameters();
  force_set_aa_mod_list(amod_list, 3);
  // include mass shifts
  aa_mod_set_mass_change(amod1, 3.4);
  aa_mod_set_mass_change(amod2, 56.78);
  aa_mod_set_mass_change(amod3, -100);

  // test exact single masses, aa_mods won't be the same, but identifiers will
  const AA_MOD_T* aa = get_aa_mod_from_mass(3.4);
  fail_unless( aa_mod_get_identifier(aa) == aa_mod_get_identifier(amod1),
               "Mass 3.4 should return id %d but returns %d.",
               aa_mod_get_identifier(amod1), aa_mod_get_identifier(aa));
  aa = get_aa_mod_from_mass(56.78);
  fail_unless( aa_mod_get_identifier(aa) == aa_mod_get_identifier(amod2),
               "Mass 56.78 should return id %d but returns %d.",
               aa_mod_get_identifier(amod2), aa_mod_get_identifier(aa));
  aa = get_aa_mod_from_mass(-100);
  fail_unless( aa_mod_get_identifier(aa) == aa_mod_get_identifier(amod3),
               "Mass -100 should return id %d but returns %d.",
               aa_mod_get_identifier(amod3), aa_mod_get_identifier(aa));

  // test single masses slightly off
  aa = get_aa_mod_from_mass(-100.00001);
  fail_unless( aa != NULL, "The mass -100.00001 did not return an AA_MOD_T.");
  fail_unless( aa_mod_get_identifier(aa) == aa_mod_get_identifier(amod3),
               "Mass -100.00001 should return id %d but returns %d.",
               aa_mod_get_identifier(amod3), aa_mod_get_identifier(aa));

  // test combinations of masses
  aa = get_aa_mod_from_mass(3.4 - 100);
  fail_unless( aa != NULL, 
               "The mass %f did not return an AA_MOD_T.", 3.4 - 100);
  fail_unless( aa_mod_get_identifier(aa) == 
               (aa_mod_get_identifier(amod1) | aa_mod_get_identifier(amod3)),
               "Combined masses %f returned id %d but should have returned %d.",
               3.4-100, aa_mod_get_identifier(aa),
               (aa_mod_get_identifier(amod1) | aa_mod_get_identifier(amod3)));
              
}
END_TEST

// Convert a char* with mod symbols -> MODIFIED_AA_T*
// Convert a char* with mod masses -> MODIFIED_AA_T*
START_TEST(test_string_to_aa_mod_str){
  // init params including mods
  initialize_parameters();
  force_set_aa_mod_list(amod_list, 3);
  // include mass shifts
  aa_mod_set_mass_change(amod1, 3.4);  // symbol *
  aa_mod_set_mass_change(amod2, 56.78);// symbol #
  aa_mod_set_mass_change(amod3, -100); // symbol @

  // create a MODIFIED_AA_T* with modifications, convert to string, convert back
  const char* seq = "GGBKATRM"; 
  int len = strlen(seq);
  MODIFIED_AA_T* modaa_seq = NULL;
  convert_to_mod_aa_seq(seq, &modaa_seq);
  // now modify a couple of aas
  modify_aa( &modaa_seq[1], amod2 );
  modify_aa( &modaa_seq[4], amod1 );
  modify_aa( &modaa_seq[4], amod3 );
  modify_aa( &modaa_seq[7], amod2 );
  modify_aa( &modaa_seq[7], amod3 );
  // turn it into a string, "GB#BKA*@TRM#@"
  seq = modified_aa_string_to_string_with_symbols(modaa_seq, len);

  // now test reverting back to MOD_AA_T*
  MODIFIED_AA_T* converted_modaa_seq = NULL;
  int convert_len = convert_to_mod_aa_seq(seq, &converted_modaa_seq);
  fail_unless( convert_len == len,
               "MOD_AA seq converted from symbols has len %d but should be %d.",
               convert_len, len);
  for(int i=0; i<len; i++){
    fail_unless(converted_modaa_seq[i] == modaa_seq[i],
                "MOD_AA[%d] converted from symbols should be %d but is %d.",
                i, modaa_seq[i], converted_modaa_seq[i]);
  }

  // repeat test from masses, GB[56.78]BKA[-96.60]TRM[-43.22]
  seq = modified_aa_string_to_string_with_masses(modaa_seq, len, MOD_MASS_ONLY);
  free(converted_modaa_seq);
  convert_len = convert_to_mod_aa_seq(seq, &converted_modaa_seq);
  fail_unless( convert_len == len,
         "MOD_AA seq converted from merged masses has len %d but should be %d.",
               convert_len, len);
  for(int i=0; i<len; i++){
    fail_unless(converted_modaa_seq[i] == modaa_seq[i],
              "MOD_AA[%d] converted from merged masses should be %d but is %d.",
                i, modaa_seq[i], converted_modaa_seq[i]);
  }

  // repeat with masses not merged, GB[56.78]BKA[3.40,-100.00]TRM[56.78,-100.00]
  seq = modified_aa_string_to_string_with_masses(modaa_seq, len, 
                                                 MOD_MASSES_SEPARATE);
  free(converted_modaa_seq);
  convert_len = convert_to_mod_aa_seq(seq, &converted_modaa_seq);
  fail_unless( convert_len == len,
    "MOD_AA seq converted from comma-separated masses has len %d but should be %d.",
               convert_len, len);
  for(int i=0; i<len; i++){
    fail_unless(converted_modaa_seq[i] == modaa_seq[i],
     "MOD_AA[%d] converted from comma-separated masses should be %d but is %d.",
                i, modaa_seq[i], converted_modaa_seq[i]);
  }

  // repeat with masses of residue plus mod, GG[113.83]BKA[-25.52]TRM[87.97]
  seq = modified_aa_string_to_string_with_masses(modaa_seq, len, 
                                                 AA_PLUS_MOD);
  free(converted_modaa_seq);
  convert_len = convert_to_mod_aa_seq(seq, &converted_modaa_seq, AA_PLUS_MOD);
  fail_unless( convert_len == len,
    "MOD_AA seq converted from comma-separated masses has len %d but should be %d.",
               convert_len, len);
  for(int i=0; i<len; i++){
    fail_unless(converted_modaa_seq[i] == modaa_seq[i],
     "MOD_AA[%d] converted from comma-separated masses should be %d but is %d.",
                i, modaa_seq[i], converted_modaa_seq[i]);
  }

}
END_TEST


// test copy function
START_TEST(test_copy_mod_seq){
  initialize_parameters();
  const char* seq = "GBBKATRM"; 
  int len = strlen(seq);
  MODIFIED_AA_T* mod_seq = NULL;
  convert_to_mod_aa_seq(seq, &mod_seq);

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

// test is_palindrome function
START_TEST(test_palindrome){
  initialize_parameters();
  const char* seq = "GBBKATRM"; 
  int len = strlen(seq);
  MODIFIED_AA_T* mod_seq = NULL;
  convert_to_mod_aa_seq(seq, &mod_seq);

  fail_unless( modified_aa_seq_is_palindrome(mod_seq, len) == false,
               "The seq %s should not be a palindrome.", seq);

  // try a palindrome
  seq = "MACDDCAR";
  len = strlen(seq);
  free(mod_seq);
  convert_to_mod_aa_seq(seq, &mod_seq);
  fail_unless( modified_aa_seq_is_palindrome(mod_seq, len) == true,
               "The seq %s should be a palindrome.", seq);

  // try a palindrome of odd length
  seq = "MACDVDCAR";
  len = strlen(seq);
  free(mod_seq);
  convert_to_mod_aa_seq(seq, &mod_seq);
  fail_unless( modified_aa_seq_is_palindrome(mod_seq, len) == true,
               "The seq %s should be a palindrome.", seq);

  // modify half making it not a palindrome
  force_set_aa_mod_list(amod_list, 3);  // in setup???
  modify_aa(&mod_seq[2], amod1);
  fail_unless( modified_aa_seq_is_palindrome(mod_seq, len) == false,
               "The modified seq %s should not be a palindrome.", 
               modified_aa_string_to_string_with_symbols(mod_seq, len));

  // modify the other half returning it to a palindrome
  modify_aa(&mod_seq[6], amod1);
  fail_unless( modified_aa_seq_is_palindrome(mod_seq, len) == true,
               "The modified seq %s should be a palindrome.", 
               modified_aa_string_to_string_with_symbols(mod_seq, len));


  // try an almost-palindrome
  seq = "MACDVCAR";
  len = strlen(seq);
  free(mod_seq);
  convert_to_mod_aa_seq(seq, &mod_seq);
  fail_unless( modified_aa_seq_is_palindrome(mod_seq, len) == false,
               "The seq %s should not be a palindrome.", seq);

  // try a short seq
  seq = "MAR";
  len = strlen(seq);
  free(mod_seq);
  convert_to_mod_aa_seq(seq, &mod_seq);
  fail_unless( modified_aa_seq_is_palindrome(mod_seq, len) == true,
               "The seq %s should not be a palindrome.", seq);

}
END_TEST

START_TEST(test_unmodify){
  const char* unmod = "MODSEQHERE";
  const char* mod1 = "MODSE*QHERE^";
  const char* mod2 = "M[8.8]ODSEQHERE";

  char* undone = unmodify_sequence(unmod);
  fail_unless( strcmp(unmod, undone) == 0,
               "The sequence '%s' should not have changed but it is '%s'",
               unmod, undone);
  free(undone);

  undone = unmodify_sequence(mod1);
  fail_unless( strcmp(unmod, undone) == 0,
               "The sequence '%s' should not have changed but it is '%s'",
               unmod, undone);
  free(undone);

  undone = unmodify_sequence(mod2);
  fail_unless( strcmp(unmod, undone) == 0,
               "The sequence '%s' should not have changed but it is '%s'",
               unmod, undone);
  free(undone);
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
  tcase_add_test(tc_core, test_mod_to_string_precision);
  tcase_add_test(tc_core, test_mod_str_to_string);
  tcase_add_test(tc_core, test_mod_str_to_unmod_string);
  tcase_add_test(tc_core, test_symbol_to_aa_mod);
  tcase_add_test(tc_core, test_mass_to_aa_mod);
  tcase_add_test(tc_core, test_string_to_aa_mod_str);
  tcase_add_test(tc_core, test_copy_mod_seq);
  tcase_add_test(tc_core, test_palindrome);
  tcase_add_test(tc_core, test_is_modifiable);
  tcase_add_test(tc_core, test_modify);
  tcase_add_test(tc_core, test_unmodify);

  tcase_add_checked_fixture(tc_core, mod_setup, mod_teardown);
  suite_add_tcase(s, tc_core);

  // Test boundry conditions
  TCase *tc_limits = tcase_create("Limits");
  tcase_add_test(tc_limits, test_too_many_mods);
  tcase_add_checked_fixture(tc_limits, mod_setup, mod_teardown);
  suite_add_tcase(s, tc_limits);

  return s;
}













