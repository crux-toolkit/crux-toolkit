#include <cppunit/config/SourcePrefix.h>
#include "TestModifications.h"
#include "modifications.h"
#include "parameter.h"
// also in parameter.c
void force_set_aa_mod_list(AA_MOD_T** amod_list, int num_mods);


CPPUNIT_TEST_SUITE_REGISTRATION( TestModifications );

void TestModifications::setUp(){
  // initialize parameters
  initialize_parameters();

  // add some modifications (this normally done in parameter.c)
  // assigns identifiers and symbols to each aamod
  amod1 = new_aa_mod(0);
  amod2 = new_aa_mod(1);
  amod3 = new_aa_mod(2);

  amod_list[0] = amod1;
  amod_list[1] = amod2;
  amod_list[2] = amod3;

  force_set_aa_mod_list(amod_list, 3);

  // initialize a sequence in char* and MODIFIED_AA_T* form
  seq1 = my_copy_string("MYPEPSEQ");
  modSeq1 = new MODIFIED_AA_T[strlen(seq1)];
  for(int i=0; i < strlen(seq1); i++){
    modSeq1[i] = (MODIFIED_AA_T)seq1[i] - (MODIFIED_AA_T)'A';
  }

  // initialize a sequence in char* and MODIFIED_AA_T* form
  seq2 = my_copy_string("MYMODSEQ");
  len2 = strlen(seq2);
  modSeq2 = new MODIFIED_AA_T[len2];
  for(int i=0; i < len2; i++){
    modSeq2[i] = (MODIFIED_AA_T)seq2[i] - (MODIFIED_AA_T)'A';
  }
  // add modifications
  seq2 = my_copy_string("MYM*OD#SE*#Q");
  modSeq2[2] = modSeq2[2] | aa_mod_get_identifier(amod1); // M*
  modSeq2[4] = modSeq2[4] | aa_mod_get_identifier(amod2); // D#
  modSeq2[6] = modSeq2[6] | aa_mod_get_identifier(amod1); // E*
  modSeq2[6] = modSeq2[6] | aa_mod_get_identifier(amod2); // E#
}

void TestModifications::tearDown(){
  delete seq1;
  delete modSeq1;
}

void TestModifications::convertToModAaSeq(){
  // sequence with no modifications
  MODIFIED_AA_T* converted = NULL;
  int len = convert_to_mod_aa_seq(seq1, &converted);
  for(int i=0; i < len; i++){
    CPPUNIT_ASSERT(converted[i] == modSeq1[i]);
  }

  // sequence with modifications
  free(converted);
  len = convert_to_mod_aa_seq(seq2, &converted);
  for(int i=0; i < len; i++){
    CPPUNIT_ASSERT(converted[i] == modSeq2[i]);
  }

  // null sequence
  free(converted);
  converted = NULL;
  convert_to_mod_aa_seq(NULL, &converted);
  CPPUNIT_ASSERT(converted == NULL);

  // sequence with bad characters
  convert_to_mod_aa_seq("BADM.OD", &converted);
  CPPUNIT_ASSERT( converted == NULL );
}










