#include <math.h>
#include <stdlib.h>
#include "check-ion.h"
#include "../ion.h"
#include "../crux-utils.h"
#include "../objects.h"


START_TEST (test_create){
  int modification_counts[] = {-1, 1, 0, 0};
    
  //create new ion
  ION_T* ion = new_ion(B_ION, 6, 2, "AKLVKNMM", AVERAGE);
  ION_T* ion2 = new_ion(Y_ION, 7, 3, "AKLVKNMM", MONO);
  ION_T* ion3 = new_modified_ion(B_ION, 6, 3, "AKLVKNMM", MONO, modification_counts);

  //check if ion mass
  fail_unless( compare_float(get_ion_mass_z(ion), 327.919342) == 0, "Ion mass/z not set correctly: ion");
  fail_unless( compare_float(get_ion_mass_z(ion2), 288.500122) == 0, "Ion mass/z not set correctly: ion2");
  fail_unless( compare_float(get_ion_mass_z(ion3), 219.143311) == 0, 
               "Ion3 mass/z is %f, should be %f", get_ion_mass_z(ion3), 555 );

  //try print ion
  print_ion(ion, stdout);

  //add modification
  add_modification(ion, NH3, -1, AVERAGE);
  fail_unless( compare_float(get_ion_mass_z(ion), 319.406067) == 0, "Ion modification failed: #1");
  add_modification(ion, NH3, 1, AVERAGE);
  add_modification(ion, H2O, 1, AVERAGE);
  fail_unless( compare_float(get_ion_mass_z(ion), 336.924622) == 0, "Ion modification failed: #2");
  
  //check set, get methods for ion
  set_ion_mass_z(ion,  270.87);
  fail_unless( compare_float(get_ion_mass_z(ion), 270.87) == 0, "Ion mass/z not set correctly");

  set_ion_cleavage_idx(ion, 5);
  fail_unless(get_ion_cleavage_idx(ion) == 5, "Ion cleavage idx  not set correctly");

  set_ion_charge(ion, 3);
  fail_unless(get_ion_charge(ion) == 3, "Ion charge not set correctly");

  set_ion_type(ion, A_ION);
  fail_unless(get_ion_type(ion) == A_ION, "Ion type not set correctly");

  set_ion_peptide_sequence(ion, "AKLVK");
  fail_unless(strcmp(get_ion_peptide_sequence(ion), "AKLVK") == 0, "Ion type not set correctly");

  //try print ion
  print_ion(ion, stdout);

  //free all ions
  free_ion(ion);
  free_ion(ion2);
  free_ion(ion3);
}
END_TEST

Suite *ion_suite(void){
  Suite *s = suite_create("ion");
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  return s;
}
