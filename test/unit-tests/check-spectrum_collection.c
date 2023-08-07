#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "check-peak.h"
#include "Spectrum.h"
#include "spectrum_collection.h"
#include "peak.h"
#include "crux-utils.h"
#include "carp.h"
#include "parameter.h"


static Spectrum* spectrum;
static Spectrum* second_spectrum;
static char* file_name;
static char* comment;

static SPECTRUM_COLLECTION_T* collection;
static SPECTRUM_COLLECTION_T* sc_copy;
static SPECTRUM_ITERATOR_T* spectrum_iterator;
static Spectrum* tmp1;


START_TEST (test_create){
  //set up
  initialize_parameters();
  collection = new_spectrum_collection("test2.ms2");
  sc_copy = allocate_spectrum_collection();
  spectrum = NULL; //new Spectrum();
  second_spectrum = NULL;//new Spectrum();
  
  fail_unless(!get_spectrum_collection_is_parsed(collection), "should not be parsed yet");

  //try parse one single spectrum
  second_spectrum = get_spectrum_collection_spectrum(collection, 16);
  fail_unless(second_spectrum->getFirstScan() == 16 &&
              second_spectrum->getLastScan() == 16, 
              "failed to get specific spectrum from .ms2");

  //parse all spectrum now...
  parse_spectrum_collection(collection);
  fail_unless(get_spectrum_collection_is_parsed(collection), "should be parsed now");
  fail_unless(get_spectrum_collection_num_spectra(collection) == 6, "failed number of specrum field");
  file_name = get_spectrum_collection_filename(collection);  
  char* filename_temp = parse_filename(file_name);
  fail_unless((strcmp("test2.ms2", filename_temp) == 0),"file name incorrect" );
  free(file_name);
  file_name = NULL;
  free(filename_temp);
  filename_temp = NULL;

  comment = get_spectrum_collection_comment(collection);
  fail_unless(strcmp("H       CreationDate    8/29/2004 3:06:54 AM\nH       Extractor       MakeMS2\nH       ExtractorVersion        1.0\nH       Comments        MakeMS2 written by Michael J. MacCoss, 2004\nH       ExtractorOptions        MS2/MS1\n", comment) == 0, "failed comment field");
  free(comment);
  comment = NULL;

  //try copy the collection
  copy_spectrum_collection(collection, sc_copy);
  
  fail_unless(get_spectrum_collection_is_parsed(sc_copy), "should be parsed now");
  fail_unless(get_spectrum_collection_num_spectra(sc_copy) == 6, "failed number of specrum field");
  file_name = get_spectrum_collection_filename(sc_copy);
  filename_temp = parse_filename(file_name);
  fail_unless((strcmp("test2.ms2", filename_temp) == 0),"file name incorrect" );
  free(file_name);
  file_name = NULL;
  free(filename_temp);
  filename_temp = NULL;

  comment = get_spectrum_collection_comment(sc_copy);
  fail_unless(strcmp("H       CreationDate    8/29/2004 3:06:54 AM\nH       Extractor       MakeMS2\nH       ExtractorVersion        1.0\nH       Comments        MakeMS2 written by Michael J. MacCoss, 2004\nH       ExtractorOptions        MS2/MS1\n", comment) == 0, "failed comment field");
  free(comment);
  comment = NULL;
  
  //iterator
  spectrum_iterator = new_spectrum_iterator(collection);
  tmp1 = spectrum_iterator_next(spectrum_iterator);
  free_spectrum_iterator(spectrum_iterator);

  spectrum = new Spectrum(*tmp1);
  
  //remove a spectrum from collection
  remove_spectrum(collection, tmp1);

  spectrum_iterator = new_spectrum_iterator(collection);
  tmp1 = spectrum_iterator_next(spectrum_iterator);
  free_spectrum_iterator(spectrum_iterator);
  
  fail_unless((get_spectrum_collection_num_spectra(collection) == 5 &&
               tmp1->getFirstScan() == 16), "failed to remove");
  
  
  //add the spectrum back into collection
  add_spectrum(collection, spectrum);

  spectrum_iterator = new_spectrum_iterator(collection);
  tmp1 = spectrum_iterator_next(spectrum_iterator);
  
  fail_unless(get_spectrum_collection_num_spectra(collection) == 6 &&
              tmp1->getFirstScan() == 15, "failed to add");
  //free up
  delete second_spectrum;
  free_spectrum_collection(sc_copy);
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
