#include <math.h>
#include <stdlib.h>
#include "check-ion.h"
#include "Ion.h"
#include "IonSeries.h"
#include "IonConstraint.h"
#include "IonFilteredIterator.h"
#include "crux-utils.h"
#include "objects.h"
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
IonSeries *is1, *is2, *is3;
char *seq, *seq2;
MODIFIED_AA_T *mod_seq, *mod_seq2;
int charge;
IonConstraint *constraint, *bcnst, *ycnst;
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
  //  aa_mod_set_mass_change(amod3, 8);
  aa_mod_set_mass_change(amod3, 70);

  amod_list[0] = amod1;
  amod_list[1] = amod2;
  amod_list[2] = amod3;

  initialize_parameters();
  force_set_aa_mod_list(amod_list, 3);

  // init ion series and constraints
  charge = 2;
  is1 = new IonSeries();

  constraint = new IonConstraint(MONO, charge-1, BY_ION, false);// no prec ion
  bcnst = new IonConstraint(MONO, 1, B_ION, false); 
  ycnst = new IonConstraint(MONO, 1, Y_ION, false); 
  seq = my_copy_string("ASEQ");
  seq2 = my_copy_string("ANOTHERSEQ");
  convert_to_mod_aa_seq( seq, &mod_seq );
  convert_to_mod_aa_seq( seq2, &mod_seq2 );
  is2 = new IonSeries(seq, charge, constraint);
  is3 = new IonSeries(constraint, charge);
}

void ion_series_teardown(){
  delete is1;
  delete is2;
  delete is3;
  delete constraint;

}

START_TEST (test_new){
  fail_unless( is1 != NULL, "Failed to allocate ion series" );
  fail_unless( is2 != NULL, "Failed to create ion series" );
  char* gotten_seq = is2->getPeptide();
  fail_unless( strcmp(gotten_seq, seq) == 0,
               "Seq of ion series is %s and should be %s", gotten_seq, seq);
  unsigned int len = is2->getPeptideLength();
  fail_unless( len == strlen(seq), 
               "Ion series seq len is %i, should be %i", len, strlen(seq));
  int ch = is2->getCharge();
  fail_unless( ch == charge, 
               "Ion series charge is %i, should be %i", ch, charge);
  fail_unless( is2->getIonConstraint() == constraint,
               "Ion series constraint is not correct");
  fail_unless( is2->getNumIons() == 0,
               "Ion series not yet predicted should have 0 ions.");
}
END_TEST

START_TEST (test_update){
  // update generic ion series with specific seq
  is3->update(seq, mod_seq);

  char* gotten_seq = is3->getPeptide();
  fail_unless( strcmp(gotten_seq, seq) == 0,
               "Seq of ion series is %s and should be %s", gotten_seq, seq);
  /*
  MODIFIED_AA_T* mod_seq = get_ion_series_mod_seq(is3);
  gotten_seq = modified_aa_string_to_string(mod_seq);
  fail_unless( strcmp(gotten_seq, seq) == 0,
  "Seq of ion series is %s and should be %s", gotten_seq, seq);
  */
  unsigned int len = is3->getPeptideLength();
  fail_unless( len == strlen(seq), 
               "Ion series seq len is %i, should be %i", len, strlen(seq));
}
END_TEST

START_TEST (test_predict){
  // predict ions for an updated series
  is3->update(seq, mod_seq);
  fail_unless( is2->getNumIons() == 0,
               "Ion series pre-prediction should have 0 ions.");

  is3->predictIons();

  // test by getting num ions and getting ions
  int num_ions = is3->getNumIons();
  fail_unless( num_ions == 6,
               "Num predicted ions is %i, should be %i.", num_ions, 6);

  // look at all ions and confirm that the masses are correct
  // reuse these
  Ion* ion = NULL;
  int i=0;

  // check b-ions
  IonFilteredIterator* iter = new IonFilteredIterator(is3, bcnst);
  float running_total = MASS_H_MONO;
  while( iter->hasNext() ){
    running_total += get_mass_amino_acid(seq[i], MONO);
    ion = iter->next();
    float mz = ion->getMassZ();
    //printf("\nion mz: %f, type: %i", get_ion_mass_z(ion), get_ion_type(ion));
    fail_unless( is_same(mz, running_total),
                 "y-ion %i of seq has mz %f, should have mz %f",
                 i, mz, running_total);
    i++;
  }
  // check y-ions
  delete iter;
  iter = new IonFilteredIterator(is3, ycnst);
  i = strlen(seq) - 1;
  running_total = MASS_H_MONO + MASS_H2O_MONO;
  while( iter->hasNext() ){
    running_total += get_mass_amino_acid(seq[i], MONO);
    ion = iter->next();
    float mz = ion->getMassZ();
    //printf("\nion mz: %f, type: %i, calc: %f",mz,get_ion_type(ion),running_total);
    fail_unless( is_same(mz, running_total),
               "y-ion %i of seq has mz %f, should have mz %f",
                 i, mz, running_total);
    i--;
  }
  delete iter;
}
END_TEST

START_TEST (test_predict2){

  // same as test_predict, but different sequence
  is3->update(seq2, mod_seq2);
  int num_ions = is3->getNumIons();
  fail_unless( is2->getNumIons() == 0,
               "Ion series after update, pre-prediction should have 0 ions.");

  is3->predictIons();
  num_ions = is3->getNumIons();
  fail_unless( num_ions == 18,
               "Num predicted ions is %i, should be %i.", num_ions, 18);

  // check b-ions
  IonFilteredIterator* iter = new IonFilteredIterator(is3, bcnst);
  float running_total = MASS_H_MONO;
  int i = 0;
  while( iter->hasNext() ){
    running_total += get_mass_amino_acid(seq2[i], MONO);
    Ion* ion = iter->next();
    float mz = ion->getMassZ();
    //printf("\nion mz: %f, type: %i", get_ion_mass_z(ion), get_ion_type(ion));
    fail_unless( is_same(mz, running_total),
                 "b-ion %i of seq 2 has mz %f, should have mz %f", 
                 i, mz, running_total);
    i++;
  }
  // check y-ions
  delete iter;
  iter = new IonFilteredIterator(is3, ycnst);
  i = strlen(seq2) - 1;
  running_total = MASS_H_MONO + MASS_H2O_MONO;
  while( iter->hasNext() ){
    running_total += get_mass_amino_acid(seq2[i], MONO);
    Ion* ion = iter->next();
    float mz = ion->getMassZ();
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
  is3->update(seq, mod_seq);
  fail_unless( is2->getNumIons() == 0,
               "Ion series pre-prediction should have 0 ions.");

  is3->predictIons();

  // test by getting num ions and getting ions
  int num_ions = is3->getNumIons();
  fail_unless( num_ions == 6,
               "Num predicted ions is %i, should be %i.", num_ions, 6);

  // look at all ions and confirm that the masses are correct
  // reuse these
  Ion* ion = NULL;
  int i=0;

  // check b-ions
  IonFilteredIterator* iter = new IonFilteredIterator(is3, bcnst);
  float running_total = MASS_H_MONO;
  while( iter->hasNext() ){
    //running_total += get_mass_mod_amino_acid(mod_seq[i], MONO);
    running_total += get_mass_amino_acid(seq[i], MONO);
    if( i == 1 ){ running_total += aa_mod_get_mass_change(amod1); }
    ion = iter->next();
    float mz = ion->getMassZ();
    //printf("\nion mz: %f, type: %i", get_ion_mass_z(ion), get_ion_type(ion));
    fail_unless( is_same(mz, running_total),
                 "b-ion %i of mod seq has mz %f, should have mz %f",
                 i, mz, running_total);
    i++;
  }
  // check y-ions
  delete iter;
  iter = new IonFilteredIterator(is3, ycnst);
  i = strlen(seq) - 1;
  running_total = MASS_H_MONO + MASS_H2O_MONO;
  while( iter->hasNext() ){
    running_total += get_mass_mod_amino_acid(mod_seq[i], MONO);
    ion = iter->next();
    float mz = ion->getMassZ();
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
  is3->update(seq, mod_seq);
  fail_unless( is2->getNumIons() == 0,
               "Ion series pre-prediction should have 0 ions.");

  is3->predictIons();

  // test by getting num ions and getting ions
  int num_ions = is3->getNumIons();
  fail_unless( num_ions == 6,
               "Num predicted ions is %i, should be %i.", num_ions, 6);

  // look at all ions and confirm that the masses are correct
  // create a b- and a y-ion constraint
  IonConstraint* bcnst = new IonConstraint(MONO, 1, B_ION, false); 
  IonConstraint* ycnst = new IonConstraint(MONO, 1, Y_ION, false); 
  // reuse these
  Ion* ion = NULL;
  int i=0;

  // check b-ions
  IonFilteredIterator* iter = new IonFilteredIterator(is3, bcnst);
  float running_total = MASS_H_MONO;
  while( iter->hasNext() ){
    running_total += get_mass_mod_amino_acid(mod_seq[i], MONO);
    ion = iter->next();
    float mz = ion->getMassZ();
    //printf("\nion mz: %f, type: %i", get_ion_mass_z(ion), get_ion_type(ion));
    fail_unless( is_same(mz, running_total),
                 "b-ion %i of mod seq has mz %f, should have mz %f",
                 i, mz, running_total);
    i++;
  }
  // check y-ions
  delete iter;
  iter = new IonFilteredIterator(is3, ycnst);
  i = strlen(seq) - 1;
  running_total = MASS_H_MONO + MASS_H2O_MONO;
  while( iter->hasNext() ){
    running_total += get_mass_mod_amino_acid(mod_seq[i], MONO);
    ion = iter->next();
    float mz = ion->getMassZ();
//print("\nion mz: %f, type: %i, calc: %f",mz,get_ion_type(ion),running_total);
    fail_unless( is_same(mz, running_total),
               "y-ion %i of mod_seq has mz %f, should have mz %f",
                 i, mz, running_total);
    i--;
  }

}
END_TEST

START_TEST(test_masses_mod){
  // similar to the test above, but 
  // compare two ion series to eachother, one with single mod of mass m
  // and one with two mods on same aa that sum to m

  // use mods with masses that sum
  aa_mod_set_mass_change(amod1, 100);
  aa_mod_set_mass_change(amod2, 130);
  aa_mod_set_mass_change(amod3, 70);

  // modifiy seq
  modify_aa(&mod_seq[1], amod1);
  // predict ions for an updated series
  is3->update(seq, mod_seq);
  is3->predictIons();

  // get a copy of seq and modify differently
  char* seq_copy = my_copy_string(seq);
  MODIFIED_AA_T* mod_seq_copy =  NULL;
  convert_to_mod_aa_seq(seq_copy, &mod_seq_copy);
  modify_aa(&mod_seq_copy[1], amod2);
  modify_aa(&mod_seq_copy[1], amod3);

  // confirm that the modified masses of aa[1] are now the same
  float mass = get_mass_amino_acid(seq[1], MONO);
  float mod_mass = get_mass_mod_amino_acid(mod_seq[1], MONO);
  float mod_mass_copy = get_mass_mod_amino_acid(mod_seq_copy[1], MONO);
  fail_unless( is_same(mod_mass, mod_mass_copy), 
               "Modified masses of %c (%f and %f) should be %f+100=%f",
               seq[1], mod_mass, mod_mass_copy, mass, mass + 100);

  // compare is4 to is3
  IonSeries* is4 = new IonSeries(constraint, charge);
  is4->update(seq_copy, mod_seq_copy);
  is4->predictIons();

  // test that num_ions same
  int num_ions = is3->getNumIons();
  int num_ions_copy = is4->getNumIons();
  fail_unless( num_ions == num_ions_copy,
   "The number of ions in the two series should be the same but are %i and %i",
               num_ions, num_ions_copy);

  // reuse these
  Ion* ion = NULL;
  Ion* ion_copy = NULL;
  int i=0;

  IonFilteredIterator* iter = new IonFilteredIterator(is3, bcnst);
  IonFilteredIterator* iter_copy = new IonFilteredIterator(is4, bcnst);
  float predicted_mass = MASS_H_MONO;

  while( iter->hasNext() ){
    predicted_mass += get_mass_amino_acid(seq[i], MONO);
    if( i==1 ){ predicted_mass += 100; } // our modification
    
    ion = iter->next();
    ion_copy = iter_copy->next();
    // get cleavage index and ion type
    int clev = ion->getCleavageIdx();
    int clev_copy = ion->getCleavageIdx();
    fail_unless( clev == clev_copy, 
                 "Ion %i has cleavage index %i and %i.",
                 clev, clev_copy);
    fail_unless( is_same(predicted_mass, ion->getMassZ()),
                 "Predicted mass is %f and ion mass is %f.",
                 predicted_mass, ion->getMassZ());
    fail_unless( is_same(ion->getMassZ(), ion_copy->getMassZ()),
                 "Masses for ion %i, index %i, %f and %f should be %f.", i,
                 ion->getMassZ(), ion_copy->getMassZ(), predicted_mass);
    i++;
  }// next ion


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
  tcase_add_test(tc_main, test_masses_mod);

  tcase_add_checked_fixture(tc_main, ion_series_setup, ion_series_teardown);
  suite_add_tcase(s, tc_main);

  return s;
}

















//chris's test
START_TEST (test_create){
  initialize_parameters();  
  /****************************
   * test ion_constraint
   ***************************/

  //creat new ion_constraint
  IonConstraint* ion_constraint = IonConstraint::newIonConstraintSequest(MONO, 3, ALL_ION, true);
  IonConstraint* ion_constraint2 = new IonConstraint();
  IonConstraint* ion_constraint3 = new IonConstraint(MONO, 3, ALL_ION, true);  
  IonConstraint* ion_constraint4 = new IonConstraint(MONO, 3, B_ION, false);  

  //set ion_constraint3 modification counts
  ion_constraint3->setModification(NH3, -2);
  ion_constraint3->setModification(H2O, 1);
  ion_constraint3->setModification(ISOTOPE, 2);
  ion_constraint3->setModification(FLANK, 1);

  //check if count set correctly
  fail_unless(ion_constraint->getModification(NH3) == 1, "ion constraint not set correctly, modification NH3");
  fail_unless(ion_constraint->getModification(H2O) == 1, "ion constraint not set correctly, modification H2O");
  fail_unless(ion_constraint->getModification(ISOTOPE) == 1, "ion constraint not set correctly, modification ISOTOPE");
  fail_unless(ion_constraint->getModification(FLANK) == 1, "ion constraint not set correctly, modification FLANK");

  fail_unless(ion_constraint3->getModification(NH3) == -2, "ion constraint not set correctly, modification NH3");
  fail_unless(ion_constraint3->getModification(H2O) == 1, "ion constraint not set correctly, modification H2O");
  fail_unless(ion_constraint3->getModification(ISOTOPE) == 2, "ion constraint not set correctly, modification ISOTOPE");
  fail_unless(ion_constraint3->getModification(FLANK) == 1, "ion constraint not set correctly, modification FLANK");
  
  //try copy 
  IonConstraint::copy(ion_constraint, ion_constraint2);
  
  //check if copy is correct
  fail_unless(ion_constraint2->getModification(NH3) == 1, "ion constraint not set correctly, modification NH3");
  fail_unless(ion_constraint2->getModification(H2O) == 1, "ion constraint not set correctly, modification H2O");
  fail_unless(ion_constraint2->getModification(ISOTOPE) == 1, "ion constraint not set correctly, modification ISOTOPE");
  fail_unless(ion_constraint2->getModification(FLANK) == 1, "ion constraint not set correctly, modification FLANK");

  //set ion_constraint4 modification counts
  ion_constraint4->setModification(NH3, ion_constraint3->getModification(NH3));
  ion_constraint4->setModification(H2O, ion_constraint3->getModification(H2O));
  ion_constraint4->setModification(ISOTOPE, ion_constraint3->getModification(ISOTOPE));
  ion_constraint4->setModification(FLANK, ion_constraint3->getModification(FLANK));
  
  
  /****************************
   * test ion_series
   ***************************/
  
  //create new ion_series AKLVKNMM
  IonSeries* ion_series = new IonSeries("AKLVKNMT", 2, ion_constraint);
  IonSeries* ion_series2 = new IonSeries();

  //check get, set
  fail_unless(strcmp(ion_series->getPeptide(), "AKLVKNMT") == 0, "ion series not set correctly, peptide sequence"); 
  char* tmpseq = my_copy_string("AKLVKNMTM");
  ion_series->setPeptide(tmpseq);
  fail_unless(strcmp(ion_series->getPeptide(), "AKLVKNMTM") == 0, "ion series not set correctly, peptide sequence"); 
  free(tmpseq);

  fail_unless(ion_series->getCharge() == 2, "ion series not set correctly, peptide charge"); 
  ion_series->setCharge(3);
  fail_unless(ion_series->getCharge() == 3, "ion series not set correctly, peptide charge"); 

  fail_unless(ion_series->getIonConstraint() == ion_constraint, "ion series not set correctly, ion constraint"); 
  ion_series->setIonConstraint(ion_constraint3);
  fail_unless(ion_series->getIonConstraint() == ion_constraint3, "ion series not set correctly, ion constraint"); 

  //now predict ions
  ion_series->predictIons();

  //check the number of ions predicted
  fail_unless(ion_series->getNumIons() == 1578, "the total number of ions not predicted correctly_1");
  fail_unless(ion_series->getNumIonsOneType(B_ION) == 288, "the B ion number of ions not predicted correctly");
  fail_unless(ion_series->getNumIonsOneType(Y_ION) == 288, "the Y ion number of ions not predicted correctly");
  fail_unless(ion_series->getNumIonsOneType(A_ION) == 246, "the A ion number of ions not predicted correctly");
  fail_unless(ion_series->getNumIonsOneType(C_ION) == 246, "the C ion number of ions not predicted correctly");
  fail_unless(ion_series->getNumIonsOneType(X_ION) == 246, "the X ion number of ions not predicted correctly");
  fail_unless(ion_series->getNumIonsOneType(Z_ION) == 246, "the Z ion number of ions not predicted correctly");
  fail_unless(ion_series->getNumIonsOneType(P_ION) == 18, "the P ion number of ions not predicted correctly");

  //try copy ion series
  IonSeries::copy(ion_series, ion_series2);

  //check the number of ions copied
  fail_unless(ion_series->getNumIons() == ion_series2->getNumIons(), "the total number of ions not predicted correctly_2");
  fail_unless(ion_series->getNumIonsOneType(A_ION) == ion_series2->getNumIonsOneType(A_ION), "the A ion number of ions not predicted correctly");
  fail_unless(ion_series->getNumIonsOneType(B_ION) == ion_series2->getNumIonsOneType(B_ION), "the B ion number of ions not predicted correctly");
  fail_unless(ion_series->getNumIonsOneType(C_ION ) == ion_series2->getNumIonsOneType(C_ION), "the C ion number of ions not predicted correctly");
  fail_unless(ion_series->getNumIonsOneType(X_ION) == ion_series2->getNumIonsOneType(X_ION), "the X ion number of ions not predicted correctly");
  fail_unless(ion_series->getNumIonsOneType(Y_ION) == ion_series2->getNumIonsOneType(Y_ION), "the Y ion number of ions not predicted correctly");
  fail_unless(ion_series->getNumIonsOneType(Z_ION) == ion_series2->getNumIonsOneType(Z_ION), "the B ion number of ions not predicted correctly");
  fail_unless(ion_series->getNumIonsOneType(P_ION) == ion_series2->getNumIonsOneType(P_ION), "the P ion number of ions not predicted correctly");
 
  //try print ion series, should be the same
  //print_ion_series(ion_series, stdout);
  //print_ion_series(ion_series2, stdout);
  
  
  //test iterators
  IonIterator ion_iterator = ion_series->begin();
  IonFilteredIterator* ion_filtered_iterator = new IonFilteredIterator(ion_series, ion_constraint4);

  int total_ions = 0;
  Ion* temp_ion;

  //iterate over all ions
  while(ion_iterator != ion_series->end()){
    ++total_ions;
    //get next ion
    temp_ion = *ion_iterator;
    ++ion_iterator;
    
    //check ion meets ion_constraint
    fail_unless(ion_constraint3->isSatisfied(temp_ion), "ion does not meet the constraint, the ion was created under");
  }
  
  //check total ion nums
  fail_unless(total_ions == ion_series->getNumIons(), "ion iterator failed to iterate over all ions");

  //check for filtered ion iterator
  total_ions = 0;
  while(ion_filtered_iterator->hasNext()){
    temp_ion = ion_filtered_iterator->next();
    
    //check ion meets ion_constraint
    fail_unless(ion_constraint4->isSatisfied(temp_ion), "ion does not meet the constraint, the ion was filtered under");
    ++total_ions;
  }

  //check if the filtered iterator iterated over all B_IONS
  fail_unless(total_ions == ion_series->getNumIonsOneType(B_ION), "ion_filtered_iterator failed to iterate over all ions that meets its constraint");
  
  //free
  delete ion_constraint;
  delete ion_constraint2;  
  delete ion_constraint3;  
  delete ion_constraint4;  
  delete ion_series;
  delete ion_series2;
//delete ion_iterator);
  delete ion_filtered_iterator;
}
END_TEST

Suite *ion_series_suite(void){
  Suite *s = suite_create("ion_series");
  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_create);
  suite_add_tcase(s, tc_core);

  return s;
}
