#include <math.h>
#include <stdlib.h>
#include "check-ion.h"
#include "../ion.h"
#include "../ion_series.h"
#include "../crux-utils.h"
#include "../objects.h"
// "included" from parameter.c
void force_set_aa_mod_list(AA_MOD_T** amod_list, int num_mods);

// helper functions
int is_same(float a, float b){
  float delta = .0005;
  if( b > (a + delta) ){
    return 0;
  }else if( b < (a - delta) ){
    return 0;
  }
  return 1;
}

// declare things to set up
ION_SERIES_T *is1, *is2, *is3;
char *seq, *seq2;
MODIFIED_AA_T *mod_seq, *mod_seq2;
int charge;
ION_CONSTRAINT_T *constraint, *bcnst, *ycnst;
AA_MOD_T *amod1, *amod2, *amod3;
AA_MOD_T* amod_list[3];

void ion_series_setup(){
  // init modifications
  // assigns identifiers and symbols to each aamod
  amod1 = new_aa_mod(0);
  amod2 = new_aa_mod(1);
  amod3 = new_aa_mod(2);

  aa_mod_set_mass_change(amod1, 100);
  aa_mod_set_mass_change(amod2, 30);
  aa_mod_set_mass_change(amod3, 8);

  amod_list[0] = amod1;
  amod_list[1] = amod2;
  amod_list[2] = amod3;

  initialize_parameters();
  force_set_aa_mod_list(amod_list, 3);

  // init ion series and constraints
  charge = 2;
  is1 = allocate_ion_series();

  constraint = new_ion_constraint(MONO, charge-1, BY_ION, FALSE);// no prec ion
  bcnst = new_ion_constraint(MONO, 1, B_ION, FALSE); 
  ycnst = new_ion_constraint(MONO, 1, Y_ION, FALSE); 
  seq = "ASEQ";
  seq2 = "ANOTHERSEQ";
  mod_seq = convert_to_mod_aa_seq( seq );
  mod_seq2 = convert_to_mod_aa_seq( seq2 );
  is2 = new_ion_series(seq, charge, constraint);
  is3 = new_ion_series_generic(constraint, charge);
}

void ion_series_teardown(){
  free_ion_series( is1 );
  free_ion_series( is2 );
  free_ion_series( is3 );
  free_ion_constraint( constraint );

}

START_TEST (test_new){
  fail_unless( is1 != NULL, "Failed to allocate ion series" );
  fail_unless( is2 != NULL, "Failed to create ion series" );
  char* gotten_seq = get_ion_series_peptide(is2);
  fail_unless( strcmp(gotten_seq, seq) == 0,
               "Seq of ion series is %s and should be %s", gotten_seq, seq);
  int len = get_ion_series_peptide_length(is2);
  fail_unless( len == strlen(seq), 
               "Ion series seq len is %i, should be %i", len, strlen(seq));
  int ch = get_ion_series_charge(is2);
  fail_unless( ch == charge, 
               "Ion series charge is %i, should be %i", ch, charge);
  fail_unless( get_ion_series_ion_constraint(is2) == constraint,
               "Ion series constraint is not correct");
  fail_unless( get_ion_series_num_ions(is2) == 0,
               "Ion series not yet predicted should have 0 ions.");
}
END_TEST

START_TEST (test_update){
  // update generic ion series with specific seq
  update_ion_series(is3, seq, mod_seq);

  char* gotten_seq = get_ion_series_peptide(is3);
  fail_unless( strcmp(gotten_seq, seq) == 0,
               "Seq of ion series is %s and should be %s", gotten_seq, seq);
  /*
  MODIFIED_AA_T* mod_seq = get_ion_series_mod_seq(is3);
  gotten_seq = modified_aa_string_to_string(mod_seq);
  fail_unless( strcmp(gotten_seq, seq) == 0,
  "Seq of ion series is %s and should be %s", gotten_seq, seq);
  */
  int len = get_ion_series_peptide_length(is3);
  fail_unless( len == strlen(seq), 
               "Ion series seq len is %i, should be %i", len, strlen(seq));
}
END_TEST

START_TEST (test_predict){
  // predict ions for an updated series
  update_ion_series(is3, seq, mod_seq);
  fail_unless( get_ion_series_num_ions(is2) == 0,
               "Ion series pre-prediction should have 0 ions.");

  predict_ions(is3);

  // test by getting num ions and getting ions
  int num_ions = get_ion_series_num_ions(is3);
  fail_unless( num_ions == 6,
               "Num predicted ions is %i, should be %i.", num_ions, 6);

  // look at all ions and confirm that the masses are correct
  // reuse these
  ION_T* ion = NULL;
  int i=0;

  // check b-ions
  ION_FILTERED_ITERATOR_T* iter = new_ion_filtered_iterator(is3, bcnst);
  float running_total = MASS_H_MONO;
  while( ion_filtered_iterator_has_next(iter) ){
    running_total += get_mass_amino_acid(seq[i], MONO);
    ion = ion_filtered_iterator_next(iter);
    float mz = get_ion_mass_z(ion);
    //printf("\nion mz: %f, type: %i", get_ion_mass_z(ion), get_ion_type(ion));
    fail_unless( is_same(mz, running_total),
                 "y-ion %i of seq has mz %f, should have mz %f",
                 i, mz, running_total);
    i++;
  }
  // check y-ions
  free_ion_filtered_iterator(iter);
  iter = new_ion_filtered_iterator(is3, ycnst);
  i = strlen(seq) - 1;
  running_total = MASS_H_MONO + MASS_H2O_MONO;
  while( ion_filtered_iterator_has_next(iter) ){
    running_total += get_mass_amino_acid(seq[i], MONO);
    ion = ion_filtered_iterator_next(iter);
    float mz = get_ion_mass_z(ion);
    //printf("\nion mz: %f, type: %i, calc: %f",mz,get_ion_type(ion),running_total);
    fail_unless( is_same(mz, running_total),
               "y-ion %i of seq has mz %f, should have mz %f",
                 i, mz, running_total);
    i--;
  }
  free_ion_filtered_iterator(iter);
}
END_TEST

START_TEST (test_predict2){

  // same as test_predict, but different sequence
  update_ion_series(is3, seq2, mod_seq2);
  int num_ions = get_ion_series_num_ions(is3);
  fail_unless( get_ion_series_num_ions(is2) == 0,
               "Ion series after update, pre-prediction should have 0 ions.");

  predict_ions(is3);
  num_ions = get_ion_series_num_ions(is3);
  fail_unless( num_ions == 18,
               "Num predicted ions is %i, should be %i.", num_ions, 18);

  // check b-ions
  ION_FILTERED_ITERATOR_T* iter = new_ion_filtered_iterator(is3, bcnst);
  float running_total = MASS_H_MONO;
  int i = 0;
  while( ion_filtered_iterator_has_next(iter) ){
    running_total += get_mass_amino_acid(seq2[i], MONO);
    ION_T* ion = ion_filtered_iterator_next(iter);
    float mz = get_ion_mass_z(ion);
    //printf("\nion mz: %f, type: %i", get_ion_mass_z(ion), get_ion_type(ion));
    fail_unless( is_same(mz, running_total),
                 "b-ion %i of seq 2 has mz %f, should have mz %f", 
                 i, mz, running_total);
    i++;
  }
  // check y-ions
  free_ion_filtered_iterator(iter);
  iter = new_ion_filtered_iterator(is3, ycnst);
  i = strlen(seq2) - 1;
  running_total = MASS_H_MONO + MASS_H2O_MONO;
  while( ion_filtered_iterator_has_next(iter) ){
    running_total += get_mass_amino_acid(seq2[i], MONO);
    ION_T* ion = ion_filtered_iterator_next(iter);
    float mz = get_ion_mass_z(ion);
//print("\nion mz: %f, type: %i, calc: %f",mz,get_ion_type(ion),running_total);
    fail_unless( is_same(mz, running_total),
                 "y-ion %i of seq2 has mz %f, should have mz %f",
                 i, mz, running_total);
    i--;
  }

  // modify the seq and predict again to see that masses change
}
END_TEST

START_TEST (test_predict_mod){
  // modifiy seq
  modify_aa(&mod_seq[1], amod1);
  // predict ions for an updated series
  update_ion_series(is3, seq, mod_seq);
  fail_unless( get_ion_series_num_ions(is2) == 0,
               "Ion series pre-prediction should have 0 ions.");

  predict_ions(is3);

  // test by getting num ions and getting ions
  int num_ions = get_ion_series_num_ions(is3);
  fail_unless( num_ions == 6,
               "Num predicted ions is %i, should be %i.", num_ions, 6);

  // look at all ions and confirm that the masses are correct
  // create a b- and a y-ion constraint
  ION_CONSTRAINT_T* bcnst = new_ion_constraint(MONO, 1, B_ION, FALSE); 
  ION_CONSTRAINT_T* ycnst = new_ion_constraint(MONO, 1, Y_ION, FALSE); 
  // reuse these
  ION_T* ion = NULL;
  int i=0;

  // check b-ions
  ION_FILTERED_ITERATOR_T* iter = new_ion_filtered_iterator(is3, bcnst);
  float running_total = MASS_H_MONO;
  while( ion_filtered_iterator_has_next(iter) ){
    running_total += get_mass_mod_amino_acid(mod_seq[i], MONO);
    ion = ion_filtered_iterator_next(iter);
    float mz = get_ion_mass_z(ion);
    //printf("\nion mz: %f, type: %i", get_ion_mass_z(ion), get_ion_type(ion));
    fail_unless( is_same(mz, running_total),
                 "b-ion %i of mod seq has mz %f, should have mz %f",
                 i, mz, running_total);
    i++;
  }
  // check y-ions
  free_ion_filtered_iterator(iter);
  iter = new_ion_filtered_iterator(is3, ycnst);
  i = strlen(seq) - 1;
  running_total = MASS_H_MONO + MASS_H2O_MONO;
  while( ion_filtered_iterator_has_next(iter) ){
    running_total += get_mass_mod_amino_acid(mod_seq[i], MONO);
    ion = ion_filtered_iterator_next(iter);
    float mz = get_ion_mass_z(ion);
//print("\nion mz: %f, type: %i, calc: %f",mz,get_ion_type(ion),running_total);
    fail_unless( is_same(mz, running_total),
               "y-ion %i of mod_seq has mz %f, should have mz %f",
                 i, mz, running_total);
    i--;
  }

}
END_TEST

START_TEST (test_predict_mod2){
  // modifiy seq
  modify_aa(&mod_seq[1], amod1);
  modify_aa(&mod_seq[1], amod2);
  modify_aa(&mod_seq[strlen(seq)-1], amod3);

  // predict ions for an updated series
  update_ion_series(is3, seq, mod_seq);
  fail_unless( get_ion_series_num_ions(is2) == 0,
               "Ion series pre-prediction should have 0 ions.");

  predict_ions(is3);

  // test by getting num ions and getting ions
  int num_ions = get_ion_series_num_ions(is3);
  fail_unless( num_ions == 6,
               "Num predicted ions is %i, should be %i.", num_ions, 6);

  // look at all ions and confirm that the masses are correct
  // create a b- and a y-ion constraint
  ION_CONSTRAINT_T* bcnst = new_ion_constraint(MONO, 1, B_ION, FALSE); 
  ION_CONSTRAINT_T* ycnst = new_ion_constraint(MONO, 1, Y_ION, FALSE); 
  // reuse these
  ION_T* ion = NULL;
  int i=0;

  // check b-ions
  ION_FILTERED_ITERATOR_T* iter = new_ion_filtered_iterator(is3, bcnst);
  float running_total = MASS_H_MONO;
  while( ion_filtered_iterator_has_next(iter) ){
    running_total += get_mass_mod_amino_acid(mod_seq[i], MONO);
    ion = ion_filtered_iterator_next(iter);
    float mz = get_ion_mass_z(ion);
    //printf("\nion mz: %f, type: %i", get_ion_mass_z(ion), get_ion_type(ion));
    fail_unless( is_same(mz, running_total),
                 "b-ion %i of mod seq has mz %f, should have mz %f",
                 i, mz, running_total);
    i++;
  }
  // check y-ions
  free_ion_filtered_iterator(iter);
  iter = new_ion_filtered_iterator(is3, ycnst);
  i = strlen(seq) - 1;
  running_total = MASS_H_MONO + MASS_H2O_MONO;
  while( ion_filtered_iterator_has_next(iter) ){
    running_total += get_mass_mod_amino_acid(mod_seq[i], MONO);
    ion = ion_filtered_iterator_next(iter);
    float mz = get_ion_mass_z(ion);
//print("\nion mz: %f, type: %i, calc: %f",mz,get_ion_type(ion),running_total);
    fail_unless( is_same(mz, running_total),
               "y-ion %i of mod_seq has mz %f, should have mz %f",
                 i, mz, running_total);
    i--;
  }

}
END_TEST

Suite *ion_series_suite_2(void){
  Suite *s = suite_create("Ion series");
  TCase *tc_main = tcase_create("main");
  tcase_add_test(tc_main, test_new);
  tcase_add_test(tc_main, test_update);
  tcase_add_test(tc_main, test_predict);
  tcase_add_test(tc_main, test_predict2);
  tcase_add_test(tc_main, test_predict_mod);
  tcase_add_test(tc_main, test_predict_mod2);

  tcase_add_checked_fixture(tc_main, ion_series_setup, ion_series_teardown);
  suite_add_tcase(s, tc_main);

  return s;
}

















//chris's test
START_TEST (test_create){
  
  /****************************
   * test ion_constraint
   ***************************/

  //creat new ion_constraint
  ION_CONSTRAINT_T* ion_constraint = new_ion_constraint_sequest(MONO, 3, ALL_ION, TRUE);
  ION_CONSTRAINT_T* ion_constraint2 = allocate_ion_constraint();
  ION_CONSTRAINT_T* ion_constraint3 = new_ion_constraint(MONO, 3, ALL_ION, TRUE);  
  ION_CONSTRAINT_T* ion_constraint4 = new_ion_constraint(MONO, 3, B_ION, FALSE);  

  //set ion_constraint3 modification counts
  set_ion_constraint_modification( ion_constraint3, NH3, -2);
  set_ion_constraint_modification( ion_constraint3, H2O, 1);
  set_ion_constraint_modification( ion_constraint3, ISOTOPE, 2);
  set_ion_constraint_modification( ion_constraint3, FLANK, 1);

  //check if count set correctly
  fail_unless(get_ion_constraint_modification(ion_constraint, NH3) == 1, "ion constraint not set correctly, modification NH3");
  fail_unless(get_ion_constraint_modification(ion_constraint, H2O) == 1, "ion constraint not set correctly, modification H2O");
  fail_unless(get_ion_constraint_modification(ion_constraint, ISOTOPE) == 1, "ion constraint not set correctly, modification ISOTOPE");
  fail_unless(get_ion_constraint_modification(ion_constraint, FLANK) == 1, "ion constraint not set correctly, modification FLANK");

  fail_unless(get_ion_constraint_modification(ion_constraint3, NH3) == -2, "ion constraint not set correctly, modification NH3");
  fail_unless(get_ion_constraint_modification(ion_constraint3, H2O) == 1, "ion constraint not set correctly, modification H2O");
  fail_unless(get_ion_constraint_modification(ion_constraint3, ISOTOPE) == 2, "ion constraint not set correctly, modification ISOTOPE");
  fail_unless(get_ion_constraint_modification(ion_constraint3, FLANK) == 1, "ion constraint not set correctly, modification FLANK");
  
  //try copy 
  copy_ion_constraint(ion_constraint, ion_constraint2);
  
  //check if copy is correct
  fail_unless(get_ion_constraint_modification(ion_constraint2, NH3) == 1, "ion constraint not set correctly, modification NH3");
  fail_unless(get_ion_constraint_modification(ion_constraint2, H2O) == 1, "ion constraint not set correctly, modification H2O");
  fail_unless(get_ion_constraint_modification(ion_constraint2, ISOTOPE) == 1, "ion constraint not set correctly, modification ISOTOPE");
  fail_unless(get_ion_constraint_modification(ion_constraint2, FLANK) == 1, "ion constraint not set correctly, modification FLANK");

  //set ion_constraint4 modification counts
  set_ion_constraint_modification( ion_constraint4, NH3, get_ion_constraint_modification(ion_constraint3, NH3));
  set_ion_constraint_modification( ion_constraint4, H2O, get_ion_constraint_modification(ion_constraint3, H2O));
  set_ion_constraint_modification( ion_constraint4, ISOTOPE, get_ion_constraint_modification(ion_constraint3, ISOTOPE));
  set_ion_constraint_modification( ion_constraint4, FLANK,get_ion_constraint_modification(ion_constraint3, FLANK));
  
  
  /****************************
   * test ion_series
   ***************************/
  
  //create new ion_series AKLVKNMM
  ION_SERIES_T* ion_series = new_ion_series("AKLVKNMT", 2, ion_constraint);
  ION_SERIES_T* ion_series2 = allocate_ion_series();

  //check get, set
  fail_unless(strcmp(get_ion_series_peptide(ion_series), "AKLVKNMT") == 0, "ion series not set correctly, peptide sequence"); 
  set_ion_series_peptide(ion_series, "AKLVKNMTM");
  fail_unless(strcmp(get_ion_series_peptide(ion_series), "AKLVKNMTM") == 0, "ion series not set correctly, peptide sequence"); 

  fail_unless(get_ion_series_charge(ion_series) == 2, "ion series not set correctly, peptide charge"); 
  set_ion_series_charge(ion_series, 3);
  fail_unless(get_ion_series_charge(ion_series) == 3, "ion series not set correctly, peptide charge"); 

  fail_unless(get_ion_series_ion_constraint(ion_series) == ion_constraint, "ion series not set correctly, ion constraint"); 
  set_ion_series_ion_constraint(ion_series, ion_constraint3);
  fail_unless(get_ion_series_ion_constraint(ion_series) == ion_constraint3, "ion series not set correctly, ion constraint"); 

  //now predict ions
  predict_ions(ion_series);

  //check the number of ions predicted
  fail_unless(get_ion_series_num_ions(ion_series) == 1830, "the total number of ions not predicted correctly_1");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, B_ION) == 330, "the B ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, Y_ION) == 330, "the Y ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, A_ION) == 288, "the A ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, C_ION) == 288, "the C ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, X_ION) == 288, "the X ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, Z_ION) == 288, "the Z ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, P_ION) == 18, "the P ion number of ions not predicted correctly");

  //try copy ion series
  copy_ion_series(ion_series, ion_series2);

  //check the number of ions copied
  fail_unless(get_ion_series_num_ions(ion_series) == get_ion_series_num_ions(ion_series2), "the total number of ions not predicted correctly_2");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, A_ION) == get_ion_series_num_ions_one_type(ion_series2, A_ION), "the A ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, B_ION) == get_ion_series_num_ions_one_type(ion_series2, B_ION), "the B ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, C_ION ) == get_ion_series_num_ions_one_type(ion_series2, C_ION), "the C ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, X_ION) == get_ion_series_num_ions_one_type(ion_series2, X_ION), "the X ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, Y_ION) == get_ion_series_num_ions_one_type(ion_series2, Y_ION), "the Y ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, Z_ION) == get_ion_series_num_ions_one_type(ion_series2, Z_ION), "the B ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, P_ION) == get_ion_series_num_ions_one_type(ion_series2, P_ION), "the P ion number of ions not predicted correctly");
 
  //try print ion series, should be the same
  //print_ion_series(ion_series, stdout);
  //print_ion_series(ion_series2, stdout);
  
  
  //test iterators
  ION_ITERATOR_T* ion_iterator = new_ion_iterator(ion_series);
  ION_FILTERED_ITERATOR_T* ion_filtered_iterator = new_ion_filtered_iterator(ion_series, ion_constraint4);

  int total_ions = 0;
  ION_T* temp_ion;

  //iterate over all ions
  while(ion_iterator_has_next(ion_iterator)){
    ++total_ions;
    //get next ion
    temp_ion = ion_iterator_next(ion_iterator);
    
    //check ion meets ion_constraint
    fail_unless(ion_constraint_is_satisfied(ion_constraint3, temp_ion), "ion does not meet the constraint, the ion was created under");
  }
  
  //check total ion nums
  fail_unless(total_ions == get_ion_series_num_ions(ion_series), "ion iterator failed to iterate over all ions");

  //check for filtered ion iterator
  total_ions = 0;
  while(ion_filtered_iterator_has_next(ion_filtered_iterator)){
    temp_ion = ion_filtered_iterator_next(ion_filtered_iterator);
    
    //check ion meets ion_constraint
    fail_unless(ion_constraint_is_satisfied(ion_constraint4, temp_ion), "ion does not meet the constraint, the ion was filtered under");
    ++total_ions;
  }

  //check if the filtered iterator iterated over all B_IONS
  fail_unless(total_ions == get_ion_series_num_ions_one_type(ion_series, B_ION), "ion_filtered_iterator failed to iterate over all ions that meets its constraint");
  
  //free
  free_ion_constraint(ion_constraint);
  free_ion_constraint(ion_constraint2);  
  free_ion_constraint(ion_constraint3);  
  free_ion_constraint(ion_constraint4);  
  free_ion_series(ion_series);
  free_ion_series(ion_series2);
  free_ion_iterator(ion_iterator);
  free_ion_filtered_iterator(ion_filtered_iterator);
}
END_TEST

Suite *ion_series_suite(void){
  Suite *s = suite_create("ion_series");
  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_create);
  suite_add_tcase(s, tc_core);

  return s;
}
