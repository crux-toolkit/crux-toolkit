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


START_TEST (test_create){
  float mass_range = 1000.0;
  int max_size =  70;
  //BOOLEAN_T ok_seq = TRUE;
  set_verbosity_level(CARP_MAX);
  char* name = NULL;

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
    new_index("fasta_file",
              constraint,
              mass_range,
              max_size,
              TRUE,
              FALSE);
  
  
  fail_unless(compare_float(get_index_mass_range(_index),mass_range)==0, "failed to set mass_range, index");
  fail_unless(get_index_max_size(_index)== max_size, "failed to set max_size, index");
  fail_unless(get_index_is_unique(_index) == TRUE, "failed to set is_unique, index"  );
  fail_unless(get_index_database(_index) != NULL, " failed to set database");
  fail_unless(get_index_constraint(_index) == constraint, " failed to set constraint");
  fail_unless(!get_index_on_disk(_index), "failed to set on_disk");
  fail_unless(create_index(_index), "failed to create a index");
  fail_unless(get_index_on_disk(_index), "failed to set on_disk");
  fail_unless(index_exists(_index), "index exist method failed");
  fail_unless(strcmp((name = get_index_directory(_index)), "fasta_file_crux_index") == 0, 
              "failes to set directory name");
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
    new_search_index("fasta_file",
                     constraint, TRUE);
  fail_unless(_index != NULL, " failed to re create index");

  
  
  /**** test index peptide interator ***/
  iterator = new_index_peptide_iterator(_index);
  
  //iterate over all peptides
  while(index_peptide_iterator_has_next(iterator)){
    peptide = index_peptide_iterator_next(iterator);
    //print_peptide_in_format(peptide, ok_seq, stdout);
    free_peptide(peptide);
  }
  
  free_index(_index);
  free_index_peptide_iterator(iterator);
    
  delete_dir("fasta_file_crux_index");
}
END_TEST

Suite* index_suite(void){
  Suite *s = suite_create("index");
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  //tcase_add_checked_fixture(tc_core, setup, teardown);
  return s;
}
