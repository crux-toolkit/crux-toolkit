#include "check-crux-utils.h"
#include "crux-utils.h"
#include "Peak.h"

#include <vector>

using namespace std;

// declare things to set up
int myint1, myint2, *myintptr;

// helper function for reading a file with a list of peaks
FLOAT_T read_peaks_file(const char* filename,
                        vector<Peak*>& peaks, int* num_peaks){

  FILE* file = fopen(filename, "r");
  if( file == NULL ){
    fprintf(stderr, "FATAL: couldn't read %s\n", filename);
    exit(1);
  }

  FLOAT_T mz;
  int num;
  // first line has num peaks, precursor mz
  fscanf(file, "%d %f", &num, &mz);
  *num_peaks = num;
  peaks = allocate_peak_vector(num);

  FLOAT_T peakmz, intensity;
  for(int i=0; i<num; i++){
    fscanf(file, "%f %f", &peakmz, &intensity);
    peaks[i]->setLocation(peakmz);
    peaks[i]->setIntensity(intensity);
  }

  fclose(file);
  return mz;
}

void crux_utils_setup(){
}

void crux_utils_teardown(){
}

// does has_extension() work?
START_TEST(has_extension){

  fail_unless( has_extension("", "") == true,
               "An empty string should have the empty extension.");

  fail_unless( has_extension("", "end") == false,
               "An empty string should not have any extension.");

  fail_unless( has_extension("anystring", "") == true,
               "Any string should have the empty extension.");

  fail_unless( has_extension("anystring", NULL) == true,
               "Any string should have the null extension.");

  fail_unless( has_extension(NULL, "ext") == false,
               "A NULL string should not have an extension.");

  fail_unless( has_extension("namewithend", "end") == true,
               "namewithend should end in 'end'.");

  fail_unless( has_extension("namenoend", "no") == false,
               "namenoend should not have the extension 'no'.");

}
END_TEST

// are two floating point numbers within rounding error of eachother
START_TEST(test_is_equal){
  FLOAT_T a = 4.02;
  FLOAT_T b = a;
  int prec = 2;

  // same number
  fail_unless( is_equal(a, b, prec) == true,
               "%.*f and %.*f should be equal to %d places.",
               prec+1, a, prec+1, b, prec);

  // biggest difference but still the same
  a = 4.015;
  b = 4.024;
  fail_unless( is_equal(a, b, prec) == true,
               "%.*f and %.*f should be equal to %d places.",
               prec+1, a, prec+1, b, prec);

  // not the same at higher precision
  prec = 3;
  fail_unless( is_equal(a, b, prec) == false,
               "%.*f and %.*f should NOT be equal to %d places.",
               prec+1, a, prec+1, b, prec);

  // two very close numbers not the same
  a = 100.115;
  b = 100.114;
  prec = 2;
  fail_unless( is_equal(a, b, prec) == false,
               "%.*f and %.*f should NOT be equal to %d places.",
               prec+1, a, prec+1, b, prec);

  // check integer case
  prec = 0;
  fail_unless( is_equal(a, b, prec) == true,
               "%.*f and %.*f should be equal to %d places.",
               prec+1, a, prec+1, b, prec);

  // unequal integers
  a += 0.8;
  fail_unless( is_equal(a, b, prec) == false,
               "%.*f and %.*f should NOT be equal to %d places.",
               prec+1, a, prec+1, b, prec);

}
END_TEST

// determine if a spectrum is +1 or more
START_TEST(test_choose_charge){
  // error case

  vector<Peak*> peaks;

  CHARGE_STATE_T charge = choose_charge(0, peaks);
  fail_unless( charge == INVALID_CHARGE_STATE, "Choose charge should return INVALID_CHARGE_STATE(%d) with no peaks.",
	       INVALID_CHARGE_STATE);

  // multiply charged spec
  int num_peaks;
  FLOAT_T mz = read_peaks_file("input-data/mult-charge.peaks", 
                               peaks, &num_peaks);
  charge = choose_charge(mz, peaks);

  fail_unless( charge == MULTIPLE_CHARGE_STATE, 
               "Charge of mult-charge.peaks should be  MULTIPLE_CAHRGE_STATE(%d) but is %d",
               MULTIPLE_CHARGE_STATE,charge);

  // singly charged spec
  mz = read_peaks_file("input-data/single-charge.peaks", 
                       peaks, &num_peaks);
  charge = choose_charge(mz, peaks);
  fail_unless( charge == SINGLE_CHARGE_STATE, 
               "Charge of single-charge(%d).peaks should be  but is %d",
               SINGLE_CHARGE_STATE,charge);

  // spec with no peaks above precursor
  mz = read_peaks_file("input-data/none-above-precursor.peaks", 
                       peaks, &num_peaks);
  charge = choose_charge(mz, peaks);
  fail_unless( charge == SINGLE_CHARGE_STATE, 
               "Charge of none-above-precursor.peaks should be SINGLE_CHARGE_STATE(%d) but is %d",
	       SINGLE_CHARGE_STATE,charge);

  // spec with no peaks below precursor
  mz = read_peaks_file("input-data/none-below-precursor.peaks", 
                       peaks, &num_peaks);
  charge = choose_charge(mz, peaks);
  fail_unless( charge == MULTIPLE_CHARGE_STATE, 
               "Charge of none-below-precursor.peaks should be MULTIPLE_CHARGE_STATE(%d) but is %d",
               MULTIPLE_CHARGE_STATE,charge);

}
END_TEST

Suite* crux_utils_suite(){
  Suite* s = suite_create("crux_utils");
  // Test basic features
  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_is_equal);
  tcase_add_test(tc_core, has_extension);
  tcase_add_test(tc_core, test_choose_charge);

  tcase_add_checked_fixture(tc_core, crux_utils_setup, crux_utils_teardown);
  suite_add_tcase(s, tc_core);


  return s;
}













