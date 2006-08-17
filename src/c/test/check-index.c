#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "check-peak.h"
#include "../mass.h"
#include "../objects.h"
#include "../spectrum.h"
#include "../peak.h"
#include "../peptide.h"
#include "../peptide_src.h"
#include "../protein.h"
#include "../database.h"
#include "../index.h"

INDEX_T* _index;
DATABASE_T* database;
PEPTIDE_CONSTRAINT_T* constraint;
INDEX_PEPTIDE_ITERATOR_T* iterator;
PEPTIDE_T* peptide;

void setup(void){

}

void teardown(void){

}

START_TEST (test_create){
  float mass_range = 1000.0;
  int max_size =  70;
  BOOLEAN_T ok_seq = TRUE;
  set_verbosity_level(CARP_MAX);
  
  /**************************
   *check creating index
   ***************************/

  constraint = 
    new_peptide_constraint(TRYPTIC,
                           0, 50000, 
                           0, 255, 
                           TRUE, AVERAGE);

  //create index
  _index = 
    new_index("test.fasta",
              constraint,
              mass_range,
              max_size,
              TRUE);
  

  fail_unless(compare_float(get_index_mass_range(_index),mass_range)==0, "failed to set mass_range, index");
  fail_unless(get_index_max_size(_index)== max_size, "failed to set max_size, index");
  fail_unless(get_index_is_unique(_index) == TRUE, "failed to set is_unique, index"  );

  fail_unless(create_index(_index), "failed to create a index");
  fail_unless(index_exists(_index), "index exist method failed");
  free_index(_index);

  
  /**************************
   * check read in index
   ***************************/
  constraint = 
    new_peptide_constraint(TRYPTIC,
                           0, 5000, 
                           2, 25, 
                           TRUE, AVERAGE);
  
    _index = 
    new_search_index("test.fasta",
              constraint, TRUE);

    if(_index != NULL){
      //create index peptide interator
      iterator = new_index_peptide_iterator(_index);//, ok_seq);
    
      //iterate over all peptides
      while(index_peptide_iterator_has_next(iterator)){
        peptide = index_peptide_iterator_next(iterator);
        print_peptide_in_format(peptide, ok_seq, stdout);
        free_peptide(peptide);
      }

      free_index(_index);
      free_index_peptide_iterator(iterator);
    }
}
END_TEST

Suite* index_suite(void){
  Suite *s = suite_create("index");
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  tcase_add_checked_fixture(tc_core, setup, teardown);
  return s;
}
