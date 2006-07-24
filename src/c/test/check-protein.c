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

PEPTIDE_T* peptide2;
PEPTIDE_T* peptide3;
PEPTIDE_T* peptide4;
PEPTIDE_T* peptide5;
PEPTIDE_T* peptide1;
PROTEIN_T* protein1;
PROTEIN_T* protein2;
PROTEIN_T* protein3;
PROTEIN_PEPTIDE_ITERATOR_T* iterator;
PEPTIDE_CONSTRAINT_T* constraint;

PEPTIDE_SRC_T* association1;
PEPTIDE_SRC_T* association2;
PEPTIDE_SRC_T* association3;

void setup(void){
  //peptide4 = allocate_peptide();
  
}

void teardown(void){
  //  free_protein(protein);
}

START_TEST (test_create){
  /*
  //check the fasta file parsing
  FILE* file = fopen("fasta_file", "r");  //"/var/noble/data/protein_database/E_coli_contam_rand.fasta", "r");
  protein3 = allocate_protein();
  if(!parse_protein_fasta_file(protein3, file)){
    fprintf(stderr, "parsing failed");
  }
  fclose(file);
  print_protein(protein3, stdout);
  
 
  protein1 = new_protein("23 Jordan", "AADAAKAGAAKFFA", 14, "this is a test protein");
  print_protein(protein1, stdout);
 
  //try copy protein
  protein2 = allocate_protein();
  copy_protein(protein1, protein2);
  print_protein(protein1, stdout);
 

  //parse protein
  */

  //peptide constraint
  constraint = new_peptide_constraint(ANY_TRYPTIC, 0, 350, 2, 3,1);
  /*
  //create iterator
  iterator = new_protein_peptide_iterator(protein3, constraint);
  
  //iterate over all possible peptides
  while(protein_peptide_iterator_has_next(iterator)){
    peptide1 = protein_peptide_iterator_next(iterator);
    print_peptide(peptide1, stdout);
    free_peptide(peptide1);
  }  

  free_protein_peptide_iterator(iterator);
  free_protein(protein1);
  free_protein(protein2);

  
  free_protein(protein3);
  */
  //create a new database
  DATABASE_T* database = new_database("small_fasta");
  
  DATABASE_PEPTIDE_ITERATOR_T* iterator =
    new_database_peptide_iterator(database, constraint);

  while(database_peptide_iterator_has_next(iterator)){
    peptide4 = database_peptide_iterator_next(iterator);
    print_peptide(peptide4, stdout);
    free_peptide(peptide4);
  }


  //free database
  free_database_peptide_iterator(iterator);
  free_peptide_constraint(constraint);
  free_database(database);

}
END_TEST

Suite* protein_suite(void){
  Suite *s = suite_create("protein");
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  tcase_add_checked_fixture(tc_core, setup, teardown);
  return s;
}
