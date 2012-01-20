#include "check-mass.h"
#include "mass.h"
#include "modifications.h"
// in parameter.c
void force_set_aa_mod_list(AA_MOD_T** amod_list, int num_mods);

// declare things to set up
static AA_MOD_T *amod1, *amod2, *amod3, *amodnot;
static MODIFIED_AA_T aa1, aa2, aa3;
static AA_MOD_T* amod_list[11];
static char achar1, achar2, achar3;

int is_close_to(float a, float b){
  return (a > (b - 0.0005) && a < (b + 0.0005) );
}
void mass_setup(){
  initialize_parameters();
  //initializing the arrays happens with get_amino_acid_mass()
  amod1 = new_aa_mod(0);
  amod2 = new_aa_mod(1);
  amod3 = new_aa_mod(10);
  amodnot = new_aa_mod(5);

  amod_list[0] = amod1;
  amod_list[1] = amod2;
  amod_list[10] = amod3;
  int i=0;
  for(i=2; i<10; i++){
    amod_list[i] = amodnot;
  }
  // initialize mods in parameter.c
  force_set_aa_mod_list(amod_list, 11);

  aa_mod_set_mass_change(amod1, 500);
  aa_mod_set_mass_change(amod2, 600);
  aa_mod_set_mass_change(amod3, 700);

  achar1 = 'S';
  achar2 = 'A';
  achar3 = 'V';

  aa1 = char_aa_to_modified(achar1);
  aa2 = char_aa_to_modified(achar2);
  aa3 = char_aa_to_modified(achar3);

}

void mass_teardown(){
  free_aa_mod(amod1);
  free_aa_mod(amod2);
  free_aa_mod(amod3);
}

START_TEST(test_create){
  fail_unless( modified_aa_to_char(aa1) == achar1,
               "aa1 should be %c but is %c", achar1, modified_aa_to_char(aa1));
  fail_unless( modified_aa_to_char(aa2) == achar2,
               "aa2 should be %c but is %c", achar2, modified_aa_to_char(aa2));
  fail_unless( modified_aa_to_char(aa3) == achar3,
               "aa3 should be %c but is %c", achar3, modified_aa_to_char(aa3));

  // get masses of unmodified aas
  float avg_mass = get_mass_mod_amino_acid_average(aa1);
  fail_unless( is_close_to(avg_mass, 87.0782),
               "Mass of %c should be %.3f but is %.3f", 
               achar1, 87.0782, avg_mass);
  avg_mass = get_mass_mod_amino_acid_average(aa2);
  fail_unless( is_close_to(avg_mass, 71.0787),
               "Mass of %c should be %.3f but is %.3f", 
               achar2, 71.0787, avg_mass);
  avg_mass = get_mass_mod_amino_acid_average(aa3);
  fail_unless( is_close_to(avg_mass, 99.1326),
               "Mass of %c should be %.5f but is %.5f", 
               achar3, 99.1326, avg_mass);

  // get monoisotipic masses of unmodified aas
  float mono_mass = get_mass_mod_amino_acid_monoisotopic(aa1);
  fail_unless( is_close_to(mono_mass, 87.03203),
               "Mono mass of %c should be %.3f but is %.3f", 
               achar1, 87.03203, mono_mass);
  mono_mass = get_mass_mod_amino_acid_monoisotopic(aa2);
  fail_unless( is_close_to(mono_mass, 71.03711),
               "Mono mass of %c should be %.3f but is %.3f", 
               achar2, 71.03711, mono_mass);
  mono_mass = get_mass_mod_amino_acid_monoisotopic(aa3);
  fail_unless( is_close_to(mono_mass, 99.06841),
               "Mono mass of %c should be %.3f but is %.3f", 
               achar3, 99.06841, mono_mass);
}
END_TEST

// test that we can get the correct bitmask back from a mass shift
START_TEST(test_get_id){
  // first test single exact masses of each mod
  MODIFIED_AA_T id = get_mod_identifier(500);
  fail_unless( id == aa_mod_get_identifier(amod1),
               "Mass 500 should return identifier %d but instead returns", 
               aa_mod_get_identifier(amod1), id);
  id = get_mod_identifier(600);
  fail_unless( id == aa_mod_get_identifier(amod2),
               "Mass 600 should return identifier %d but instead returns", 
               aa_mod_get_identifier(amod2), id);
  id = get_mod_identifier(700);
  fail_unless( id == aa_mod_get_identifier(amod3),
               "Mass 700 should return identifier %d but instead returns", 
               aa_mod_get_identifier(amod3), id);

  // now test masses that are slightly different
  id = get_mod_identifier(600.000001);
  fail_unless( id == aa_mod_get_identifier(amod2),
             "Mass 600.000001 should return identifier %d but instead returns", 
               aa_mod_get_identifier(amod2), id);
  // but this is too far off
  id = get_mod_identifier(600.005);
  fail_unless( id == 0,
               "Mass 600.005 should return zero but instead returns %d", id);

  // now combine masses together
  id = get_mod_identifier(600+700);
  fail_unless( id == 
               (aa_mod_get_identifier(amod2) | aa_mod_get_identifier(amod3)),
               "Mass 600+700 should return %d but gives %d",
               (aa_mod_get_identifier(amod2) | aa_mod_get_identifier(amod3)),
               id);
}
END_TEST

START_TEST(test_avg_1mod){
  // modify each aa with each amod and test mass
  modify_aa(&aa1, amod1);
  modify_aa(&aa2, amod1);
  modify_aa(&aa3, amod1);

  float avg_mass = get_mass_mod_amino_acid_average(aa1);
  fail_unless( is_close_to(avg_mass, 87.0782 + 500),
               "Mod 1 mass of %c should be %.3f but is %.3f", 
               achar1, 87.0782 + 500, avg_mass);
  avg_mass = get_mass_mod_amino_acid_average(aa2);
  fail_unless( is_close_to(avg_mass, 71.0787 + 500),
               "Mod 1 mass of %c should be %.3f but is %.3f", 
               achar2, 71.0787 + 500, avg_mass);
  avg_mass = get_mass_mod_amino_acid_average(aa3);
  fail_unless( is_close_to(avg_mass, 99.1326 + 500),
               "Mod 1 mass of %c should be %.3f but is %.3f", 
               achar3, 99.1326 + 500, avg_mass);

  // unmodify
  aa1 = char_aa_to_modified(achar1);
  aa2 = char_aa_to_modified(achar2);
  aa3 = char_aa_to_modified(achar3);

  // mod with second aa mod
  modify_aa(&aa1, amod2);
  modify_aa(&aa2, amod2);
  modify_aa(&aa3, amod2);

  avg_mass = get_mass_mod_amino_acid_average(aa1);
  fail_unless( is_close_to(avg_mass, 87.0782 + 600),
               "Mod 2 mass of %c should be %.3f but is %.3f", 
               achar1, 87.0782 + 600, avg_mass);
  avg_mass = get_mass_mod_amino_acid_average(aa2);
  fail_unless( is_close_to(avg_mass, 71.0787 + 600),
               "Mod 2 mass of %c should be %.3f but is %.3f", 
               achar2, 71.0787 + 600, avg_mass);
  avg_mass = get_mass_mod_amino_acid_average(aa3);
  fail_unless( is_close_to(avg_mass, 99.1326 + 600),
               "Mod 2 mass of %c should be %.3f but is %.3f", 
               achar3, 99.1326 + 600, avg_mass);

  // unmodify
  aa1 = char_aa_to_modified(achar1);
  aa2 = char_aa_to_modified(achar2);
  aa3 = char_aa_to_modified(achar3);

  // mod with third aa mod
  modify_aa(&aa1, amod3);
  modify_aa(&aa2, amod3);
  modify_aa(&aa3, amod3);

  avg_mass = get_mass_mod_amino_acid_average(aa1);
  fail_unless( is_close_to(avg_mass, 87.0782 + 700),
               "Mod 3 mass of %c should be %.3f but is %.3f", 
               achar1, 87.0782 + 700, avg_mass);
  avg_mass = get_mass_mod_amino_acid_average(aa2);
  fail_unless( is_close_to(avg_mass, 71.0787 + 700),
               "Mod 3 mass of %c should be %.3f but is %.3f", 
               achar2, 71.0787 + 700, avg_mass);
  avg_mass = get_mass_mod_amino_acid_average(aa3);
  fail_unless( is_close_to(avg_mass, 99.1326 + 700),
               "Mod 3 mass of %c should be %.3f but is %.3f", 
               achar3, 99.1326 + 700, avg_mass);
}
END_TEST

START_TEST(test_avg_2mod){
  // modify each aa with each amod and test mass
  modify_aa(&aa1, amod1);
  modify_aa(&aa2, amod1);
  modify_aa(&aa3, amod1);
  // mod with second aa mod
  modify_aa(&aa1, amod2);
  modify_aa(&aa2, amod2);
  modify_aa(&aa3, amod2);


  float avg_mass = get_mass_mod_amino_acid_average(aa1);
  fail_unless( is_close_to(avg_mass, 87.0782 + 500 + 600),
               "Mod 1 mass of %c should be %.3f but is %.3f", 
               achar1, 87.0782 + 500 + 600, avg_mass);
  avg_mass = get_mass_mod_amino_acid_average(aa2);
  fail_unless( is_close_to(avg_mass, 71.0787 + 500 + 600),
               "Mod 1 mass of %c should be %.3f but is %.3f", 
               achar2, 71.0787 + 500 + 600, avg_mass);
  avg_mass = get_mass_mod_amino_acid_average(aa3);
  fail_unless( is_close_to(avg_mass, 99.1326 + 500 + 600),
               "Mod 1 mass of %c should be %.3f but is %.3f", 
               achar3, 99.1326 + 500 + 600, avg_mass);

  // unmodify
  aa1 = char_aa_to_modified(achar1);
  aa2 = char_aa_to_modified(achar2);
  aa3 = char_aa_to_modified(achar3);

  // mod with 1 and 3
  modify_aa(&aa1, amod1);
  modify_aa(&aa2, amod1);
  modify_aa(&aa3, amod1);
  modify_aa(&aa1, amod3);
  modify_aa(&aa2, amod3);
  modify_aa(&aa3, amod3);
  avg_mass = get_mass_mod_amino_acid_average(aa1);
  fail_unless( is_close_to(avg_mass, 87.0782 + 500 + 700),
               "Mod 2 mass of %c should be %.3f but is %.3f", 
               achar1, 87.0782 + 500 + 700, avg_mass);
  avg_mass = get_mass_mod_amino_acid_average(aa2);
  fail_unless( is_close_to(avg_mass, 71.0787 + 500 + 700),
               "Mod 2 mass of %c should be %.3f but is %.3f", 
               achar2, 71.0787 + 500 + 700, avg_mass);
  avg_mass = get_mass_mod_amino_acid_average(aa3);
  fail_unless( is_close_to(avg_mass, 99.1326 + 500 + 700),
               "Mod 2 mass of %c should be %.3f but is %.3f", 
               achar3, 99.1326 + 500 + 700, avg_mass);

  // unmodify
  aa1 = char_aa_to_modified(achar1);
  aa2 = char_aa_to_modified(achar2);
  aa3 = char_aa_to_modified(achar3);

  // mod with 2 and 3 
  modify_aa(&aa1, amod2);
  modify_aa(&aa2, amod2);
  modify_aa(&aa3, amod2);
  modify_aa(&aa1, amod3);
  modify_aa(&aa2, amod3);
  modify_aa(&aa3, amod3);

  avg_mass = get_mass_mod_amino_acid_average(aa1);
  fail_unless( is_close_to(avg_mass, 87.0782 + 600 + 700),
               "Mod 3 mass of %c should be %.3f but is %.3f", 
               achar1, 87.0782 + 600 + 700, avg_mass);
  avg_mass = get_mass_mod_amino_acid_average(aa2);
  fail_unless( is_close_to(avg_mass, 71.0787 + 600 + 700),
               "Mod 3 mass of %c should be %.3f but is %.3f", 
               achar2, 71.0787 + 600 + 700, avg_mass);
  avg_mass = get_mass_mod_amino_acid_average(aa3);
  fail_unless( is_close_to(avg_mass, 99.1326 + 600 + 700),
               "Mod 3 mass of %c should be %.3f but is %.3f", 
               achar3, 99.1326 + 600 + 700, avg_mass);
}
END_TEST

START_TEST(test_avg_3mod){
  // modify each aa with each amod and test mass
  modify_aa(&aa1, amod1);
  modify_aa(&aa2, amod1);
  modify_aa(&aa3, amod1);
  // mod with second aa mod
  modify_aa(&aa1, amod2);
  modify_aa(&aa2, amod2);
  modify_aa(&aa3, amod2);
  // mod with third aa mod
  modify_aa(&aa1, amod3);
  modify_aa(&aa2, amod3);
  modify_aa(&aa3, amod3);


  float avg_mass = get_mass_mod_amino_acid_average(aa1);
  fail_unless( is_close_to(avg_mass, 87.0782 + 500 + 600 + 700),
               "Mod 1 mass of %c should be %.3f but is %.3f", 
               achar1, 87.0782 + 500 + 600 + 700, avg_mass);
  avg_mass = get_mass_mod_amino_acid_average(aa2);
  fail_unless( is_close_to(avg_mass, 71.0787 + 500 + 600 + 700),
               "Mod 1 mass of %c should be %.3f but is %.3f", 
               achar2, 71.0787 + 500 + 600 + 700, avg_mass);
  avg_mass = get_mass_mod_amino_acid_average(aa3);
  fail_unless( is_close_to(avg_mass, 99.1326 + 500 + 600 + 700),
               "Mod 1 mass of %c should be %.3f but is %.3f", 
               achar3, 99.1326 + 500 + 600 + 700, avg_mass);

}
END_TEST

/* Boundry conditions test suite */
START_TEST(test_null){
}
END_TEST

Suite* mass_suite(){
  Suite* s = suite_create("Mass");
  // Test basic features
  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_create);
  tcase_add_test(tc_core, test_get_id);
  tcase_add_test(tc_core, test_avg_1mod);
  tcase_add_test(tc_core, test_avg_2mod);
  tcase_add_test(tc_core, test_avg_3mod);

  tcase_add_checked_fixture(tc_core, mass_setup, mass_teardown);
  suite_add_tcase(s, tc_core);

  // Test boundry conditions
  TCase *tc_limits = tcase_create("Limits");
  tcase_add_test(tc_limits, test_null);
  tcase_add_checked_fixture(tc_limits, mass_setup, mass_teardown);
  suite_add_tcase(s, tc_limits);

  return s;
}













