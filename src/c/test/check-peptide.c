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

void setup(void){
  peptide4 = allocate_peptide();
  
}

void teardown(void){
  free_protein(protein);
}

START_TEST (test_create){
  protein = allocate_protein();
  
  //proteins...
  /*
  protein1 = new_protein("23 Jordan", "AADAAKAGAAKFFA", 14, "this is a test protein1");
  protein2 = new_protein("33 Magic", "AADAA", 5, "this is a test protein2");
  protein3 = new_protein("32 Shaq", "AADAAKAGAAKFFADFGTS", 19, "this is a test protein3");
  */
  protein = new_protein("yo!",

"MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDALPNISDAERI
FAELLTGLAAAQPGFPLAQLKTFVDQEFAQIKHVLHGISLLGQCPDSINAALICRGEKMSIAIMAGVLEA
RGHNVTVIDPVEKLLAVGHYLESTVDIAESTRRIAASRIPADHMVLMAGFTAGNEKGELVVLGRNGSDYS
AAVLAACLRADCCEIWTDVDGVYTCDPRQVPDARLLKSMSYQEAMELSYFGAKVLHPRTITPIAQFQIPC
LIKNTGNPQAPGTLIGASRDEDELPVKGISNLNNMAMFSVSGPGMKGMVGMAARVFAAMSRARISVVLIT
QSSSEYSISFCVPQSDCVRAERAMQEEFYLELKEGLLEPLAVTERLAIISVVGDGMRTLRGISAKFFAAL
ARANINIVAIAQGSSERSISVVVNNDDATTGVRVTHQMLFNTDQVIEVFVIGVGGVGGALLEQLKRQQSW
"
                        , 490, "this is a my test protein");
                        
                        

  //create peptides
  peptide1 = new_peptide( 6, 774.83, protein, 239, TRYPTIC); //QVPDAR
  peptide2 = new_peptide( 6, 656.69, protein, 6, PARTIALLY_TRYPTIC);//FGGTSV
  peptide3 = new_peptide( 9, 1341.32, protein, 221, NOT_TRYPTIC); //DCCEIWTDV
  
  /*
  //peptides
  peptide1 = new_peptide("CDEAK", 5, 546.5958, protein1,2, TRYPTIC);
  peptide2 = new_peptide("CDEAKCDEAKCDEAK", 5*3, 546.5958*3, protein2,2, PARTIALLY_TRYPTIC);
  peptide3 = new_peptide("CDEAKCDEAKCDEAKCDEAKCDEAKCDEAK", 5*5, 546.5958*5 , protein3,2, NOT_TRYPTIC);
  */
  /*
  print_peptide(peptide1, stdout);
  print_peptide(peptide2, stdout);
  print_peptide(peptide3, stdout);
  */
  
  //check mass
  fail_unless( compare_float(get_peptide_peptide_mass(peptide1), 774.83) ==0, "failed mass");
  fail_unless( compare_float(calc_peptide_mass(peptide1, AVERAGE), 774.83) == 0 , "failed mass");
  
  fail_unless( 546.5958 -0.1<= get_peptide_neutral_mass(peptide1) <= 546.5958 + 0.1, "failed mass");
  printf("peptide charged mass(charge2): %f\n", get_peptide_charged_mass(peptide1, 2));
  printf("peptide mz(charge2): %f\n", get_peptide_mz(peptide1, 2));


  //peptide constraint
  constraint = new_peptide_constraint(PARTIALLY_TRYPTIC, 600, 2000, 7, 20, 1);
  fail_unless(!peptide_constraint_is_satisfied(constraint, peptide1), "constraint fail1");
  fail_unless(peptide_constraint_is_satisfied(constraint, peptide2), "constraint fail1");
  fail_unless(!peptide_constraint_is_satisfied(constraint, peptide3), "constraint fail1");
  
  //get, set for constraint (TRYPTIC ,400, 3000, 4, 26)
  printf("original type: %d\n", get_peptide_constraint_peptide_type(constraint));
  set_peptide_constraint_peptide_type(constraint, TRYPTIC);
  printf("new type: %d\n", get_peptide_constraint_peptide_type(constraint));

  printf("original minmass: %f\n", get_peptide_constraint_min_mass(constraint));
  set_peptide_constraint_min_mass(constraint, 400);
  printf("new minmass: %f\n", get_peptide_constraint_min_mass(constraint));

  printf("original maxmass: %f\n", get_peptide_constraint_max_mass(constraint));
  set_peptide_constraint_max_mass(constraint, 3000);
  printf("new maxmass: %f\n", get_peptide_constraint_max_mass(constraint));

  printf("original minl: %d\n", get_peptide_constraint_min_length(constraint));
  set_peptide_constraint_min_length(constraint, 4);
  printf("new minl: %d\n", get_peptide_constraint_min_length(constraint));

  printf("original maxl: %d\n", get_peptide_constraint_max_length(constraint));
  set_peptide_constraint_max_length(constraint, 26);
  printf("new maxl: %d\n", get_peptide_constraint_max_length(constraint));

  //check new peptide constraint
  fail_unless(peptide_constraint_is_satisfied(constraint, peptide1), "constraint fail2");
  fail_unless(peptide_constraint_is_satisfied(constraint, peptide2), "constraint fail2");
  fail_unless(peptide_constraint_is_satisfied(constraint, peptide3), "constraint fail2");
  

  //check set, get
  //printf("original sequence: %s\n", get_peptide_sequence(peptide1));
  set_peptide_sequence(peptide1, "KAEDCCDEAK");
  //printf("new sequence: %s\n", get_peptide_sequence(peptide1));
  
  printf("original length: %d\n", get_peptide_length(peptide1));
  set_peptide_length(peptide1, 10);
  printf("new length: %d\n", get_peptide_length(peptide1));
  
  printf("original mass: %f\n", get_peptide_peptide_mass(peptide1));
  set_peptide_peptide_mass(peptide1,546.5958*2 );
  printf("new mass: %f\n", get_peptide_peptide_mass(peptide1));

  print_peptide(peptide1, stdout);
  

  //check copy peptide ADD for association
  copy_peptide(peptide1, peptide4);
  print_peptide(peptide4, stdout);


  //check peptide_src
  association1 = new_peptide_src(TRYPTIC, protein1, 4);
  association2 = new_peptide_src(NOT_TRYPTIC, protein2, 8); //try to add this to end
  association3 = allocate_peptide_src();
  copy_peptide_src(association1, association3);

  

  //
  add_peptide_peptide_src(peptide4, association1);
  add_peptide_peptide_src(peptide4, association2);
  add_peptide_peptide_src(peptide4, association3);
  print_peptide(peptide4, stdout);


  /*
 get_peptide_src_peptide_type
 set_peptide_src_peptide_type

 get_peptide_src_parent_protein   
 set_peptide_src_parent_protein

 get_peptide_src_start_idx
 set_peptide_src_start_idx
  */
 
  free_peptide_constraint(constraint);
  free_peptide(peptide1);
  free_peptide(peptide2);
  free_peptide(peptide3);
  free_peptide(peptide4);
  free_protein(protein1);
  free_protein(protein2);
  free_protein(protein3);
}
END_TEST

Suite* peptide_suite(void){
  Suite *s = suite_create("peptide");
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  tcase_add_checked_fixture(tc_core, setup, teardown);
  return s;
}
