#include <math.h>
#include <stdlib.h>
#include "check-ion.h"
#include "../ion.h"
#include "../ion_series.h"
#include "../crux-utils.h"
#include "../objects.h"


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
  fail_unless(get_ion_series_num_ions(ion_series) == 678, "the total number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, B_ION) == 330, "the B ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, Y_ION) == 330, "the Y ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, A_ION) == 0, "the A ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, C_ION) == 0, "the C ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, X_ION) == 0, "the X ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, Z_ION) == 0, "the Z ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, P_ION) == 18, "the P ion number of ions not predicted correctly");

  //try copy ion series
  copy_ion_series(ion_series, ion_series2);

  //check the number of ions copied
  fail_unless(get_ion_series_num_ions(ion_series) == get_ion_series_num_ions(ion_series2), "the total number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, A_ION) == get_ion_series_num_ions_one_type(ion_series2, A_ION), "the A ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, B_ION) == get_ion_series_num_ions_one_type(ion_series2, B_ION), "the B ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, C_ION ) == get_ion_series_num_ions_one_type(ion_series2, C_ION), "the C ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, X_ION) == get_ion_series_num_ions_one_type(ion_series2, X_ION), "the X ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, Y_ION) == get_ion_series_num_ions_one_type(ion_series2, Y_ION), "the Y ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, Z_ION) == get_ion_series_num_ions_one_type(ion_series2, Z_ION), "the B ion number of ions not predicted correctly");
  fail_unless(get_ion_series_num_ions_one_type(ion_series, P_ION) == get_ion_series_num_ions_one_type(ion_series2, P_ION), "the P ion number of ions not predicted correctly");
 
  //try print ion series, should be the same
  print_ion_series(ion_series, stdout);
  print_ion_series(ion_series2, stdout);
  
  
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
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  return s;
}
