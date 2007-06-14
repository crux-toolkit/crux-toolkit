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
#include "../carp.h"
#include "../crux-utils.h"
#include "../database.h"


PEPTIDE_T* peptide2;
PEPTIDE_T* peptide3;
PEPTIDE_T* peptide4;
PEPTIDE_T* peptide5;
PEPTIDE_T* peptide6;
PEPTIDE_T* peptide1;
PROTEIN_T* protein1;
PROTEIN_T* protein2;
PROTEIN_T* protein3;
PROTEIN_PEPTIDE_ITERATOR_T* iterator;
PEPTIDE_CONSTRAINT_T* constraint;

PEPTIDE_SRC_T* association1;
PEPTIDE_SRC_T* association2;
PEPTIDE_SRC_T* association3;

DATABASE_T* db;


START_TEST (test_create){
  char* name = NULL;
  char* name2 = NULL;
  
  //try create a new database
  db = new_database("fasta_file_binary_fasta", FALSE, TRUE);
  fail_unless(parse_database(db), "failed to parse database");
  fail_unless(strncmp((name = get_database_filename(db)), "fasta_file", 10) == 0, "database filename not set correctly");
  free(name);
  fail_unless(get_database_is_parsed(db), "database parsed field not correctly set");
  fail_unless(get_database_num_proteins(db) == 3, "database number of proteins not set correctly");
  
  //check the fasta file parsing
  FILE* file = fopen("fasta_file", "r");
  protein3 = allocate_protein();
  fail_unless(parse_protein_fasta_file(protein3, file), "failed to parse protein from fasta file");
  set_protein_database(protein3, db);
  print_protein(protein3, stdout);
  
 
  protein1 = new_protein("23 Jordan", "AADAAKAGAAKFFA", 14, "this is a test protein", 45, 3, db);
  print_protein(protein1, stdout);
 
  //try copy protein
  protein2 = allocate_protein();
  copy_protein(protein1, protein2);
  fail_unless(strcmp((name = get_protein_id(protein1)), (name2=get_protein_id(protein2))) == 0, "protein, id not correct");
  free(name); 
  free(name2);
  fail_unless(strcmp((name = get_protein_sequence(protein1)), (name2 = get_protein_sequence(protein2))) == 0, "protein, sequence not correct");
  free(name); 
  free(name2);
  fail_unless(strcmp((name = get_protein_annotation(protein1)), (name2 = get_protein_annotation(protein2))) == 0, "protein, annotation not correct");
  free(name); 
  free(name2);
  fail_unless(get_protein_length(protein1) == get_protein_length(protein2), "protein, protein length not correct");
  fail_unless(get_protein_offset(protein1) == get_protein_offset(protein2), "protein, protein offset not correct");
  fail_unless(get_protein_protein_idx(protein1) == get_protein_protein_idx(protein2), "protein, protein idx not correct");
  fail_unless(get_protein_is_light(protein1) == get_protein_is_light(protein2), "protein, protein is_light not correct");
  
  //peptide constraint
  constraint = new_peptide_constraint(TRYPTIC, 0, 1200, 1, 10, 1, AVERAGE);
  
  /** test, protein_peptide_iterator **/
  
  //create iterator
  iterator = new_protein_peptide_iterator(protein3, constraint);
  
  //iterate over all possible peptides
  while(protein_peptide_iterator_has_next(iterator)){
    peptide1 = protein_peptide_iterator_next(iterator);
    //print_peptide(peptide1, stdout);
    free_peptide(peptide1);
  }  

  //  print_peptide_in_format(peptide5, TRUE,  stdout);

  //free stuff
  free_protein(protein1);
  free_protein(protein2);
  free_protein(protein3);
  free_protein_peptide_iterator(iterator);
  free_peptide_constraint(constraint);
  free_database(db);
  fclose(file);
}
END_TEST

Suite* protein_suite(void){
  Suite *s = suite_create("protein");
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  //tcase_add_checked_fixture(tc_core, setup, teardown);
  return s;
}
