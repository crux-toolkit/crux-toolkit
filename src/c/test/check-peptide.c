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
#include "../peptide_constraint.h"
#include "../protein.h"
#include "../database.h"

/********************************************
 * to check peptide.c & peptide_constraint.c
 *
 ********************************************/


PEPTIDE_T* peptide1;
PEPTIDE_T* peptide2;
PEPTIDE_T* peptide3;
PEPTIDE_T* peptide4;
PROTEIN_T* protein;
PROTEIN_T* protein1;
PROTEIN_T* protein2;
PROTEIN_T* protein3;
PEPTIDE_CONSTRAINT_T* constraint;
PEPTIDE_SRC_T* association1;
PEPTIDE_SRC_T* association2;
PEPTIDE_SRC_T* association3;
DATABASE_T* database;

START_TEST (test_create){
  database = new_database("test", FALSE, TRUE);
  
  //test on link_list implementaion of peptide_src
  set_peptide_src_implementation(TRUE);
  
  peptide4 = allocate_peptide();
  char* seq = NULL;
  
  protein = new_protein("test protein",

"MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDALPNISDAERIFAELLTGLAAAQPGFPLAQLKTFVDQEFAQIKHVLHGISLLGQCPDSINAALICRGEKMSIAIMAGVLEARGHNVTVIDPVEKLLAVGHYLESTVDIAESTRRIAASRIPADHMVLMAGFTAGNEKGELVVLGRNGSDYSAAVLAACLRADCCEIWTDVDGVYTCDPRQVPDARLLKSMSYQEAMELSYFGAKVLHPRTITPIAQFQIPCLIKNTGNPQAPGTLIGASRDEDELPVKGISNLNNMAMFSVSGPGMKGMVGMAARVFAAMSRARISVVLITQSSSEYSISFCVPQSDCVRAERAMQEEFYLELKEGLLEPLAVTERLAIISVVGDGMRTLRGISAKFFAALARANINIVAIAQGSSERSISVVVNNDDATTGVRVTHQMLFNTDQVIEVFVIGVGGVGGALLEQLKRQQSW"
                        , 490, "this is a my test protein", 44, 4, database); //offset and protein_idx are random
  
  protein2 = new_protein("test protein",

"MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDALPNISDAERIFAELLTGLAAAQPGFPLAQLKTFVDQEFAQIKHVLHGISLLGQCPDSINAALICRGEKMSIAIMAGVLEARGHNVTVIDPVEKLLAVGHYLESTVDIAESTRRIAASRIPADHMVLMAGFTAGNEKGELVVLGRNGSDYSAAVLAACLRADCCEIWTDVDGVYTCDPRQVPDARLLKSMSYQEAMELSYFGAKVLHPRTITPIAQFQIPCLIKNTGNPQAPGTLIGASRDEDELPVKGISNLNNMAMFSVSGPGMKGMVGMAARVFAAMSRARISVVLITQSSSEYSISFCVPQSDCVRAERAMQEEFYLELKEGLLEPLAVTERLAIISVVGDGMRTLRGISAKFFAALARANINIVAIAQGSSERSISVVVNNDDATTGVRVTHQMLFNTDQVIEVFVIGVGGVGGALLEQLKRQQSW"
                        , 490, "this is a my test protein", 44, 4, database);

  //create peptides
  peptide1 = new_peptide( 6, 684.75, protein, 239, TRYPTIC); //QVPDAR
  peptide2 = new_peptide( 6, 656.69, protein, 6, PARTIALLY_TRYPTIC);//FGGTSV
  peptide3 = new_peptide( 9, 1341.32, protein, 221, NOT_TRYPTIC); //DCCEIWTDV
  
  
  //check peptide
  fail_unless( strncmp((seq = get_peptide_sequence(peptide1)), "QVPDAR", 6) == 0, "peptide sequence no set");
  free(seq);
  //test
  printf("peptide mass %.2f\n", get_peptide_peptide_mass(peptide1));
  fail_unless( compare_float(get_peptide_peptide_mass(peptide1), 684.75) ==0, "failed mass #1");
  //debug
  printf("The peptide1 mass: %.2f\n", calc_peptide_mass(peptide1, AVERAGE));

  fail_unless(684.74 < calc_peptide_mass(peptide1, AVERAGE) &&  calc_peptide_mass(peptide1, AVERAGE) < 684.76, "failed mass #2");
  
  fail_unless( compare_float(get_peptide_neutral_mass(peptide1), 684.75) == 0, "failed mass #3");
  printf("peptide charged mass(charge2): %f\n", get_peptide_charged_mass(peptide1, 2));
  printf("peptide mz(charge2): %f\n", get_peptide_mz(peptide1, 2));


  /*************peptide constraint**************/
  constraint = new_peptide_constraint(PARTIALLY_TRYPTIC, 660, 2000, 7, 20, 1, AVERAGE);
  fail_unless(!peptide_constraint_is_satisfied(constraint, peptide1), "constraint fail1");
  fail_unless(!peptide_constraint_is_satisfied(constraint, peptide2), "constraint fail2");
  fail_unless(peptide_constraint_is_satisfied(constraint, peptide3), "constraint fail3");
  
  //get, set for peptide_constraint (TRYPTIC ,400, 3000, 4, 26)
  set_peptide_constraint_peptide_type(constraint, TRYPTIC);
  fail_unless(get_peptide_constraint_peptide_type(constraint) == TRYPTIC, "constraint peptide type not set correctly");

  set_peptide_constraint_min_mass(constraint, 400);
  fail_unless( compare_float(get_peptide_constraint_min_mass(constraint), 400) == 0 , "constraint min mass not set correctly");
  
  set_peptide_constraint_max_mass(constraint, 3000);
  fail_unless( compare_float(get_peptide_constraint_max_mass(constraint), 3000) == 0 , "constraint max mass not set correctly");
 
  set_peptide_constraint_min_length(constraint, 4);
  fail_unless( get_peptide_constraint_min_length(constraint) == 4, "constraint min length no right");

  set_peptide_constraint_max_length(constraint, 26);
  fail_unless( get_peptide_constraint_max_length(constraint) == 26, "constraint max length no right");

  set_peptide_constraint_mass_type(constraint, MONO);
  fail_unless(get_peptide_constraint_mass_type(constraint) == MONO, "constraint mass type not set correctly");
  
  //check new peptide constraint
  fail_unless(peptide_constraint_is_satisfied(constraint, peptide1), "constraint fail, peptide1");
  fail_unless(peptide_constraint_is_satisfied(constraint, peptide2), "constraint fail, peptide2");
  fail_unless(peptide_constraint_is_satisfied(constraint, peptide3), "constraint fail, peptide3");


  //check set, get for peptides again
  set_peptide_length(peptide1, 10);
  fail_unless(get_peptide_length(peptide1) == 10, "peptide length no right");
  
  set_peptide_peptide_mass(peptide1, 546.5958*2 );
  fail_unless( compare_float(get_peptide_peptide_mass(peptide1), 546.5958*2) == 0 , "peptide mass not set correctly");

  print_peptide(peptide1, stdout);
  

  //check copy peptide ADD for association
  copy_peptide(peptide1, peptide4);
  print_peptide(peptide4, stdout);


  //check peptide_src
  association1 = new_peptide_src(TRYPTIC, protein, 4);
  association2 = new_peptide_src(NOT_TRYPTIC, protein2, 8); //try to add this to end
  association3 = allocate_peptide_src();
  copy_peptide_src(association1, association3);
  
  //try adding peptide_src to a peptide
  add_peptide_peptide_src(peptide4, association1);
  add_peptide_peptide_src(peptide4, association2);
  add_peptide_peptide_src(peptide4, association3);
 
  fail_unless(get_peptide_src_peptide_type(association3) == TRYPTIC, "failed to copy | set peptide type");
  fail_unless(get_peptide_src_parent_protein(association3) == protein, "failed to copy | set parent protein");
  fail_unless(get_peptide_src_start_idx(association3) == 4, "failed to copy | set start idx");
 
  //try printing peptide in various forms..to ensure nothing blows up
  print_peptide(peptide4, stdout);
  print_peptide_in_format(peptide4, TRUE, TRUE, stdout);
  serialize_peptide(peptide4, stdout);
 
  free_database(database);
  free_peptide_constraint(constraint);
  free_peptide(peptide1);
  free_peptide(peptide2);
  free_peptide(peptide3);
  free_peptide(peptide4);
  free_protein(protein);
  free_protein(protein2);
}
END_TEST

Suite* peptide_suite(void){
  Suite *s = suite_create("peptide");
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  //tcase_add_checked_fixture(tc_core, setup, teardown);
  return s;
}
