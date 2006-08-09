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

void setup(void){

}

void teardown(void){

}

START_TEST (test_create){
  float mass_range = 1000.0;
  int max_size =  70;
  set_verbosity_level(CARP_MAX);

  constraint = 
    new_peptide_constraint(TRYPTIC,
                           0, 5000, 
                           0, 255, 
                           TRUE, AVERAGE);
    _index = 
      new_index("test.fasta",
                constraint,
                mass_range,
                max_size);

    fail_unless(create_index(_index), "failed to create a index");
    fail_unless(index_exists(_index), "index exist method failed");
    
    free_index(_index);
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
