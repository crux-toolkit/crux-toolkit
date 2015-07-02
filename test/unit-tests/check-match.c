#include <stdlib.h>
#include "check-match.h"
#include "match.h"
#include "peptide.h"
// "included" from parameter.c
void force_set_aa_mod_list(AA_MOD_T** amod_list, int num_mods);

// declare things to set up
static MATCH_T *m1, *mdecoy, *mmod, *mdecoymod;
static PROTEIN_T *prot;
static PEPTIDE_T *pep, *pepmod;
static char protseq [1028] = "MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDALPNISDAERIFAELLTGLAAAQPGFPLAQLKTFVDQEFAQIKHVLHGISLLGQCPDSINAALICRGEKMSIAIMAGVLEARGHNVTVIDPVEKLLAVGHYLESTVDIAESTRRIAASRIPADHMVLMAGFTAGNEKGELVVLGRNGSDYSAAVLAACLRADCCEIWTDVDGVYTCDPRQVPDARLLKSMSYQEAMELSYFGAKVLHPRTITPIAQFQIPCLIKNTGNPQAPGTLIGASRDEDELPVKGISNLNNMAMFSVSGPGMKGMVGMAARVFAAMSRARISVVLITQSSSEYSISFCVPQSDCVRAERAMQEEFYLELKEGLLEPLAVTERLAIISVVGDGMRTLRGISAKFFAALARANINIVAIAQGSSERSISVVVNNDDATTGVRVTHQMLFNTDQVIEVFVIGVGGVGGALLEQLKRQQSW";
static AA_MOD_T *amod1, *amod2;
static AA_MOD_T* amod_list[2];
static PEPTIDE_MOD_T* pmod;

void match_setup(){
  initialize_parameters();

  // create inputs
  prot = new_protein( "Protein1", protseq, strlen(protseq), 
                      NULL, 0, 0, NULL);//description, offset, idx, dbase
  pep = new_peptide( 10, 1087.20, prot, 20);//VADILESNAR
  pepmod = new_peptide( 10, 1087.20, prot, 20);//VADILESNAR

  // create match
  m1 = new_match();
  set_match_peptide(m1, pep);
  set_match_spectrum(m1, NULL);
  set_match_charge(m1, 2);
  set_match_null_peptide(m1, FALSE);

  // create match to null (decoy) peptide
  mdecoy = new_match();
  set_match_peptide(mdecoy, copy_peptide(pep));
  set_match_spectrum(mdecoy, NULL);
  set_match_charge(mdecoy, 2);
  set_match_null_peptide(mdecoy, TRUE);

  // set up modifications
  amod1 = new_aa_mod(0);
  amod2 = new_aa_mod(1);

  //aa_mod_set_mass_change(amod1, 100);
  //BOOLEAN_T* aa_list = aa_mod_get_aa_list(amod1);
  amod_list[0] = amod1;
  amod_list[1] = amod2;
  initialize_parameters();
  force_set_aa_mod_list(amod_list, 2);

  // create a modified peptide
  MODIFIED_AA_T* mod_seq = get_peptide_modified_aa_sequence(pepmod);
  modify_aa(&mod_seq[2], amod1);
  pmod = new_peptide_mod();
  peptide_mod_add_aa_mod(pmod, 0, 1);//mod at index 0, one copy
  set_peptide_mod(pepmod, mod_seq, pmod);

  // create match to modified peptide
  mmod = new_match();
  set_match_peptide(mmod, pepmod); 
  set_match_spectrum(mmod, NULL);
  set_match_charge(mmod, 2);
  set_match_null_peptide(mmod, FALSE);

  // create match to decoy modified peptide
  mdecoymod = new_match();
  PEPTIDE_T* decoy_pep = copy_peptide(pepmod);
  transform_peptide_to_decoy(decoy_pep);
  set_match_peptide(mdecoymod, decoy_pep); 
  set_match_spectrum(mdecoymod, NULL);
  set_match_charge(mdecoymod, 2);
  set_match_null_peptide(mdecoymod, TRUE);

}

void match_teardown(){
  free_match(m1);
  free_match(mdecoy);
  free_match(mmod);
  free_match(mdecoymod);
  free_protein(prot);
  /*
    The following two lines cause problems on MacOS.  --WSN 2010-05-25
  free_peptide(pep);
  free_peptide(pepmod);
  */
}

START_TEST(test_create){
  fail_unless( m1 != NULL, "Failed to allocate a match.");

  // test getters
  char* seq = get_match_sequence(m1);
  int len = strlen(seq);
  fail_unless( strcmp(seq, "VADILESNAR") == 0,
               "Match returns %s as seq instead of %s", seq, "VADILESNAR");
  MODIFIED_AA_T* mod_seq = get_match_mod_sequence(m1);
  char* mod_str = modified_aa_string_to_string_with_symbols(mod_seq, len);
  fail_unless( strcmp(mod_str, seq) == 0,
               "MOD_AA string should be %s but is %s", seq, mod_str);
  free(mod_str);
  mod_str = modified_aa_string_to_string_with_masses(mod_seq, len, FALSE);
  fail_unless( strcmp(mod_str, seq) == 0,
               "MOD_AA string should be %s but is %s", seq, mod_str);

  // test getting seqs for null (decoy) peptide
  free(seq);
  seq = get_match_sequence(mdecoy);
  /*
  len = strlen(seq);
  fail_unless( strcmp(seq, "VASDLINEAR") == 0,
               "Match decoy seq should be %s but is %s", "VASDLINEAR", seq );
  free(mod_seq);
  mod_seq = get_match_mod_sequence(mdecoy);
  free(mod_str);
  mod_str = modified_aa_string_to_string(mod_seq, len);
  fail_unless( strcmp(mod_str, seq) != 0,
               "For peptide with seq %s, shuffled MOD_AA string should " \
               "be different but is %s", seq, mod_str);
  */
}
END_TEST

START_TEST(test_create_mod){

  // first test that the peptide returns the seqs correctly
  char* pep_seq = get_peptide_sequence(pepmod);
  int len = get_peptide_length(pepmod);
  fail_unless( strcmp(pep_seq, "VADILESNAR") == 0,
               "Modified peptide returns vanilla seq as %s instead of %s", 
               pep_seq, "VADILESNAR");

  MODIFIED_AA_T* aa_pep_seq = get_peptide_modified_aa_sequence(pepmod);
  char* mod_pep_seq = modified_aa_string_to_string_with_symbols(aa_pep_seq, len);
  fail_unless( strcmp(mod_pep_seq, "VAD*ILESNAR") == 0,
               "Modified peptide returns annotated seq as %s instead of %s",
               mod_pep_seq, "VAD*ILESNAR");

  char* mod_pep_seq2 = 
    modified_aa_string_to_string_with_masses(aa_pep_seq, len, TRUE);
  fail_unless( strcmp(mod_pep_seq2, "VAD[0.00]ILESNAR") == 0,
               "Modified peptide returns annotated seq as %s instead of %s",
               mod_pep_seq2, "VAD[0.00]ILESNAR");

  // test getters of match with modified peptides
  char* match_seq = get_match_sequence(mmod);
  len = strlen(match_seq);
  fail_unless( strcmp(match_seq, pep_seq) == 0,
               "Modified match returns %s instead of %s", match_seq, pep_seq);
  
  MODIFIED_AA_T* aa_match_seq = get_match_mod_sequence(mmod);
  char* mod_str = modified_aa_string_to_string_with_symbols(aa_match_seq, len);
  fail_unless( strcmp(mod_str, mod_pep_seq) == 0,
               "MOD_AA string from mod peptide should be %s but is %s",
               mod_pep_seq, mod_str);
  free(mod_str);
  mod_str = modified_aa_string_to_string_with_masses(aa_match_seq, len, TRUE);
  fail_unless( strcmp(mod_str, mod_pep_seq2) == 0,
               "MOD_AA string from mod peptide should be %s but is %s",
               mod_pep_seq2, mod_str);

  // repeat for decoy match
  free(match_seq);
  match_seq = get_match_sequence(mdecoymod);
  len = strlen(match_seq);
  fail_unless( strcmp(match_seq, pep_seq) != 0,
               "Sequence from decoy mod peptide is %s, same as peptide", 
               match_seq);
  free(aa_match_seq);
  /*
//This test demonstrated that modify then shuffle is broken; only shuffle then modify works
  aa_match_seq = get_match_mod_sequence(mdecoymod);
  free(mod_str);
  mod_str = modified_aa_string_to_string(aa_match_seq, len);
  fail_unless( strcmp(mod_str, mod_pep_seq) != 0,
               "MOD_AA string from mod decoy peptide is %s, same as peptide",
               mod_str);
  */

}
END_TEST

/*
START_TEST(test_create){
}
END_TEST

START_TEST(test_create){
}
END_TEST
*/

Suite* match_suite(){
  Suite* s = suite_create("Match");
  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_create);
  tcase_add_test(tc_core, test_create_mod);
  //  tcase_add_test(tc_core, ????);

  tcase_add_checked_fixture(tc_core, match_setup, match_teardown);
  suite_add_tcase(s, tc_core);

  // Test boundry conditions
  /*
  TCase *tc_limits = tcase_create("Limits");
  tcase_add_test(tc_limits, ????);
  tcase_add_checked_fixture(tc_limits, mod_setup, mod_teardown);
  suite_ad_tcase(s, tc_limits);
   */

  return s;
}

/*
#define scan_num 16
#define ms2_file "test.ms2"
#define parameter_file "test_parameter_file"

//THE parameter system will not work if set CK_FORK=no
START_TEST (test_create){
  SPECTRUM_T* spectrum = NULL;
  SPECTRUM_COLLECTION_T* collection = NULL; ///<spectrum collection
  MATCH_COLLECTION_T* match_collection = NULL;
  MATCH_ITERATOR_T* match_iterator = NULL;
  MATCH_T* match = NULL;
    
  // comment this parameter section out, when using CK_FORK=no, valgrind
  // Parameters have been already confirmed in check_scorer.
  
  //set verbbosity level
  int  verbosity = CARP_INFO;
  set_verbosity_level(verbosity);

  //parse paramter file
  //  parse_update_parameters(parameter_file);
  initialize_parameters();

  //add fasta file parameter_file fasta-file
  //add_parameter("fasta-file", "fasta_file");
  
  //parameters has been confirmed
  //parameters_confirmed();
  // ****************************************** end comment out

  //read ms2 file
  collection = new_spectrum_collection(ms2_file);
  spectrum = allocate_spectrum();

  DATABASE_T* database = new_database("fasta-file", FALSE);

  
  //search for spectrum with correct scan number
  fail_unless(get_spectrum_collection_spectrum(collection, scan_num, spectrum), "failed to find scan_num in ms3 file");
  
  //get match collection with perliminary score of SP, and main score of XCORR
  match_collection = new_match_collection_from_spectrum(spectrum, 1, 500, SP, XCORR, 0, FALSE, NULL, database);
  
  fail_unless(get_match_collection_scored_type(match_collection, SP), "failed to set match_collection scored type, SP");
  fail_unless(get_match_collection_scored_type(match_collection, XCORR), "failed to set match_collection scored type, SP");

  //LOGP_EXP_SP should not be scored yet
  fail_unless(!get_match_collection_scored_type(match_collection, LOGP_EXP_SP), "failed to set match_collection scored type, xcorr");
  fail_unless(!get_match_collection_iterator_lock(match_collection), "match_collection lock is not set correctly"); 
  
  //create match iterator
  match_iterator = new_match_iterator(match_collection, SP, TRUE);
  
  //match_collection should be locked now..
  fail_unless(get_match_collection_iterator_lock(match_collection), "match_collection lock is not set correctly"); 
  
  //iterate over all matches
  while(match_iterator_has_next(match_iterator)){
    match = match_iterator_next(match_iterator);
    print_match(match, stdout, TRUE, SP);
  }

  //free match iterator
  free_match_iterator(match_iterator);
  
  //should be unlocked
  fail_unless(!get_match_collection_iterator_lock(match_collection), "match_collection lock is not set correctly"); 

  free_match_collection(match_collection);
  free_spectrum_collection(collection);
  free_spectrum(spectrum);
  free_parameters();
}
END_TEST


Suite *match_suite(void){
  Suite *s = suite_create("match");
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  return s;
}
*/
