#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "check-peak.h"
#include "mass.h"
#include "objects.h"
#include "Spectrum.h"
#include "peak.h"
#include "peptide.h"
#include "peptide_src.h"
#include "Protein.h"
#include "database.h"
#include "index.h"

INDEX_T* _index;
DATABASE_T* database;
PEPTIDE_CONSTRAINT_T* constraint;
INDEX_PEPTIDE_ITERATOR_T* iterator;
PEPTIDE_T* peptide;


START_TEST (test_create){
  delete_dir("fasta_file_crux_index");
  float mass_range = 1000.0;
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
              ".",
              constraint,
              mass_range);
  
  
  fail_unless(compare_float(get_index_mass_range(_index),mass_range)==0, "failed to set mass_range, index");  
  fail_unless(get_index_is_unique(_index) == TRUE, "failed to set is_unique, index"  );
  fail_unless(get_index_database(_index) != NULL, " failed to set database");
  //  fail_unless(get_index_constraint(_index) == constraint, " failed to set constraint");
  fail_unless(!get_index_on_disk(_index), "failed to set on_disk");
  fail_unless(create_index(_index), "failed to create a index");
  fail_unless(get_index_on_disk(_index), "failed to set on_disk");
  fail_unless(index_exists(_index), "index exist method failed");
  fail_unless(strcmp((name = get_index_directory(_index)), "fasta_file_crux_index") == 0, 
              "failes to set directory name");
  free(name);
  free_index(_index);

  
  /**************************
   * check read in index
   ***************************/
  constraint = 
    new_peptide_constraint(TRYPTIC,
                           0, 2000, 
                           2, 25, 
                           TRUE, AVERAGE);
  
  _index = 
    //new_index_from_disk("fasta_file", TRUE);
    new_index_from_disk("fasta_file");
  /*    new_search_index("fasta_file",
        constraint, TRUE);*/
  fail_unless(_index != NULL, " failed to re create index");
  
  /**** test index peptide interator ***/
  iterator = new_index_peptide_iterator(_index);
  
  int n = 0;
  //iterate over all peptides
  while(index_peptide_iterator_has_next(iterator)){
    ++n;
    fail_unless((peptide = index_peptide_iterator_next(iterator)) != NULL, "index_peptide_iterator failed");
    print_peptide_in_format(peptide, TRUE, TRUE, stdout);
    //free_peptide_for_array(peptide);    
    free_peptide(peptide);
  }
  fail_unless(n==22, "index peptide iterator did not return expected total number of peptides");
  free_index_peptide_iterator(iterator);
  free_index(_index);
  chdir("..");
  delete_dir("fasta_file_crux_index");
  /* only if you want to check the current directory
  char* cur_dir = NULL;
  cur_dir = getcwd(cur_dir, 50);
  printf("current directory: %s\n", cur_dir);
  free(cur_dir);
  */
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
