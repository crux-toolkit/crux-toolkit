#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "check-peak.h"
#include "../spectrum.h"
#include "../peak.h"
#include "../crux-utils.h"


SPECTRUM_T * spectrum;
SPECTRUM_T * second_spectrum;
char* file_name;
static FILE* file;


START_TEST (test_create){
  //set up
  spectrum = allocate_spectrum();
  second_spectrum = allocate_spectrum();

  file_name = my_copy_string("test.ms2"); // to avoid const warning
  fail_unless(parse_spectrum(spectrum, file_name), 
              "failed to open and create new spectrum from ms2 file");

  free_spectrum(spectrum);
  spectrum = allocate_spectrum();

  file = fopen(file_name, "r" );
  fail_unless(parse_spectrum_file(spectrum, file, file_name), 
              "failed to open and create new spectrum from ms2 file");
  
  free(file_name);

  fail_unless(get_spectrum_first_scan(spectrum) == 15, "first_scan field incorrect");
  fail_unless(get_spectrum_last_scan(spectrum) == 15, "last_scan field incorrect");
  //fail_unless(get_spectrum_id(spectrum) == 2, "id field is incorrect");
  //fail_unless(get_spectrum_spectrum_type(spectrum) == ???, "spectrum_type field incorrect");
  fail_unless(compare_float(get_spectrum_precursor_mz(spectrum), 600.78)==0 , "precursor_mz field incorrect");
  fail_unless(compare_float(get_spectrum_min_peak_mz(spectrum), 285.5)==0, "min_peak_mz field incorrect");
  fail_unless(compare_float(get_spectrum_max_peak_mz(spectrum), 735.6)==0, "max_peak_mz field incorrect");
  fail_unless(compare_float(get_spectrum_num_peaks(spectrum), 17)==0, "num_peaks field incorrect");
  fail_unless(compare_float(get_spectrum_total_energy(spectrum),64.4  ) ==0 , "total_energy field incorrect");
  file_name = get_spectrum_filename(spectrum);
  fail_unless((strcmp("test.ms2", file_name) == 0),"file name incorrect" );
  free(file_name);

  //this test copy and set_fields(copy uses set methods)
  copy_spectrum(spectrum, second_spectrum);
  
  fail_unless(get_spectrum_first_scan(second_spectrum) == 15, "first_scan field incorrect");
  fail_unless(get_spectrum_last_scan(second_spectrum) == 15, "last_scan field incorrect");
  //fail_unless(get_spectrum_spectrum_type(second_spectrum) == ???, "spectrum_type field incorrect");
  fail_unless(compare_float(get_spectrum_precursor_mz(second_spectrum), 600.78)==0 , "precursor_mz field incorrect");
  fail_unless(compare_float(get_spectrum_min_peak_mz(second_spectrum), 285.5)==0, "min_peak_mz field incorrect");
  fail_unless(compare_float(get_spectrum_max_peak_mz(second_spectrum), 735.6)==0, "max_peak_mz field incorrect");
  fail_unless(compare_float(get_spectrum_num_peaks(second_spectrum), 17)==0, "num_peaks field incorrect");
  fail_unless(compare_float(get_spectrum_total_energy(second_spectrum),64.4) ==0 , "total_energy field incorrect");
  file_name = get_spectrum_filename(second_spectrum);
  fail_unless((strcmp("test.ms2", file_name) == 0),"file name incorrect" );
  free(file_name);
  
  //free
  free_spectrum(spectrum);
  free_spectrum(second_spectrum);
  fclose(file);
}
END_TEST

Suite *spectrum_suite(void){
  Suite *s = suite_create("Spectrum");
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  return s;
}
