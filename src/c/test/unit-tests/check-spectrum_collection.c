#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "check-peak.h"
#include "spectrum.h"
#include "spectrum_collection.h"
#include "peak.h"
#include "crux-utils.h"
#include "carp.h"
#include "parameter.h"


static SPECTRUM_T* spectrum;
static SPECTRUM_T* second_spectrum;
static char* file_name;
static char* comment;

static SPECTRUM_COLLECTION_T* collection;
static SPECTRUM_COLLECTION_T* copy;
static SPECTRUM_ITERATOR_T* spectrum_iterator;
static SPECTRUM_T* tmp1;


START_TEST (test_create){
  //set up
  initialize_parameters();
  collection = new_spectrum_collection("test2.ms2");
  copy = allocate_spectrum_collection();
  spectrum = allocate_spectrum();
  second_spectrum = allocate_spectrum();
  
  fail_unless(!get_spectrum_collection_is_parsed(collection), "should not be parsed yet");

  //try parse one single spectrum
  get_spectrum_collection_spectrum(collection, 16, second_spectrum);
  fail_unless(get_spectrum_first_scan(second_spectrum) == 16 &&
              get_spectrum_last_scan(second_spectrum) == 16 /*&&
              get_spectrum_id(second_spectrum) == 2*/, "failed to get specific spectrum from .ms2");

  //parse all spectrum now...
  parse_spectrum_collection(collection);
  fail_unless(get_spectrum_collection_is_parsed(collection), "should be parsed now");
  fail_unless(get_spectrum_collection_num_spectra(collection) == 6, "failed number of specrum field");
  file_name = get_spectrum_collection_filename(collection);  
  char* filename_temp = parse_filename(file_name);
  fail_unless((strcmp("test2.ms2", filename_temp) == 0),"file name incorrect" );
  free(file_name);
  free(filename_temp);

  comment = get_spectrum_collection_comment(collection);
  fail_unless(strcmp("H       CreationDate    8/29/2004 3:06:54 AM\nH       Extractor       MakeMS2\nH       ExtractorVersion        1.0\nH       Comments        MakeMS2 written by Michael J. MacCoss, 2004\nH       ExtractorOptions        MS2/MS1\n", comment) == 0, "failed comment field");
  free(comment);

  //try copy the collection
  copy_spectrum_collection(collection, copy);
  
  fail_unless(get_spectrum_collection_is_parsed(copy), "should be parsed now");
  fail_unless(get_spectrum_collection_num_spectra(copy) == 6, "failed number of specrum field");
  file_name = get_spectrum_collection_filename(copy);
  filename_temp = parse_filename(file_name);
  fail_unless((strcmp("test2.ms2", filename_temp) == 0),"file name incorrect" );
  free(file_name);
  free(filename_temp);

  comment = get_spectrum_collection_comment(copy);
  fail_unless(strcmp("H       CreationDate    8/29/2004 3:06:54 AM\nH       Extractor       MakeMS2\nH       ExtractorVersion        1.0\nH       Comments        MakeMS2 written by Michael J. MacCoss, 2004\nH       ExtractorOptions        MS2/MS1\n", comment) == 0, "failed comment field");
  free(comment);
  
  
  //iterator
  spectrum_iterator = new_spectrum_iterator(collection);
  tmp1 = spectrum_iterator_next(spectrum_iterator);
  free_spectrum_iterator(spectrum_iterator);

  copy_spectrum(tmp1, spectrum);
  
  //remove a spectrum from collection
  remove_spectrum(collection, tmp1);

  spectrum_iterator = new_spectrum_iterator(collection);
  tmp1 = spectrum_iterator_next(spectrum_iterator);
  free_spectrum_iterator(spectrum_iterator);
  
  fail_unless((get_spectrum_collection_num_spectra(collection) == 5 &&
               get_spectrum_first_scan(tmp1) == 16), "failed to remove");
  
  
  //add the spectrum back into collection
  add_spectrum(collection, spectrum);

  spectrum_iterator = new_spectrum_iterator(collection);
  tmp1 = spectrum_iterator_next(spectrum_iterator);
  
  fail_unless((get_spectrum_collection_num_spectra(collection) == 6 &&
               get_spectrum_first_scan(spectrum) == 15), "failed to add");

  //free up
  free_spectrum(second_spectrum);
  free_spectrum_collection(copy);
  free_spectrum_collection(collection);
  free_spectrum_iterator(spectrum_iterator);

}
END_TEST

Suite *spectrum_collection_suite(void){
  Suite *s = suite_create("Spectrum_collection");
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  //tcase_add_checked_fixture(tc_core, setup, teardown);
  return s;
}
