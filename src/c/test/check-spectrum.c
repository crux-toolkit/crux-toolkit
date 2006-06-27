#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "check-peak.h"
#include "../spectrum.h"
#include "../peak.h"

SPECTRUM_T * spectrum;
SPECTRUM_T * second_spectrum;
FILE* file;

void setup(void){
  spectrum = allocate_spectrum();
  file = fopen("test.ms2", "r" );
  second_spectrum = allocate_spectrum();
}

void teardown(void){
  free_spectrum(spectrum);
  free_spectrum(second_spectrum);
  fclose(file);
}

START_TEST (test_create){
  fail_unless(parse_spectrum(spectrum, "test.ms2"), "failed to open and create new spectrum from ms2 file");
  fail_unless( parse_spectrum_file(spectrum, file), "failed to open and create new spectrum from ms2 file");
  
  fail_unless(get_spectrum_first_scan(spectrum) == 15, "first_scan field incorrect");
  fail_unless(get_spectrum_last_scan(spectrum) == 15, "last_scan field incorrect");
  //fail_unless(get_spectrum_spectrum_type(spectrum) == ???, "spectrum_type field incorrect");
  fail_unless(get_spectrum_precursor_mz(spectrum) == 600.78 , "precursor_mz field incorrect");
  fail_unless(get_spectrum_min_peak_mz(spectrum) == 285.5, "min_peak_mz field incorrect");
  fail_unless(get_spectrum_max_peak_mz(spectrum) == 735.6, "max_peak_mz field incorrect");
  fail_unless(get_spectrum_num_peaks(spectrum) == 17, "num_peaks field incorrect");
  //fail_unless(get_spectrum_total_energy(spectrum) ==  , "total_energy field incorrect");
  printf("test file name %s\n", get_spectrum_filename(spectrum));

  copy_spectrum(spectrum, second_spectrum);
  
  fail_unless(get_spectrum_first_scan(second_spectrum) == 15, "first_scan field incorrect");
  fail_unless(get_spectrum_last_scan(second_spectrum) == 15, "last_scan field incorrect");
  //fail_unless(get_spectrum_spectrum_type(second_spectrum) == ???, "spectrum_type field incorrect");
  fail_unless(get_spectrum_precursor_mz(second_spectrum) == 600.78  , "precursor_mz field incorrect");
  fail_unless(get_spectrum_min_peak_mz(second_spectrum) == 285.5, "min_peak_mz field incorrect");
  fail_unless(get_spectrum_max_peak_mz(second_spectrum) == 735.6, "max_peak_mz field incorrect");
  fail_unless(get_spectrum_num_peaks(second_spectrum) == 17, "num_peaks field incorrect");
  //fail_unless(get_spectrum_total_energy(second_spectrum) ==  , "total_energy field incorrect");
  fail_unless(strcmp(get_spectrum_filename(second_spectrum), get_spectrum_filename(spectrum))==0 , "filename field incorrect");
  
}
END_TEST

Suite *spectrum_suite(void){
  Suite *s = suite_create("Spectrum");
  TCase *tc_core = tcase_create("Core");
  //TCase *tc_core2 = tcase_create("Core2");

  suite_add_tcase(s, tc_core);
  //suite_add_tcase(s, tc_core2);

  tcase_add_test(tc_core, test_create);
  //tcase_add_test(tc_core2, test_create2);
  tcase_add_checked_fixture(tc_core, setup, teardown);
  //tcase_add_checked_fixture(tc_core2, setup, teardown);
  return s;
}
