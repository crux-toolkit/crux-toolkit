#include "check-parameter.h"
#include "parameter.h"
#include "carp.h"
#include "modifications.h"
// also in parameter.c
void parse_custom_enzyme(const char* rule_str);

// declare things to set up
const char* ops[1];
const char* args[1];
char enzyme_rule[128];
char* fake_argv[4];
  
// set up and teardown
void param_setup(){
  initialize_parameters();
  // select command line options and arguments
  ops[0] = "parameter-file";
  select_cmd_line_options(ops, 1);
  args[0] = "protein database";
  select_cmd_line_arguments(args, 1);

  fake_argv[0] = my_copy_string("test-app");
  fake_argv[1] = my_copy_string("--parameter-file");
  fake_argv[2] = my_copy_string("params/1mod");
  fake_argv[3] = my_copy_string("prot");
}

void param_teardown(){
  // free all hashes?
  for(int i=0; i<4; i++){
    free(fake_argv[i]);
  }
}

START_TEST(test_create){
  // check getters for each type (TODO)
  fail_unless( get_boolean_parameter("overwrite") == false, 
              "Does not return the bool parameter %s correctly", "overwrite");

}
END_TEST

// TODO write setup and teardown
// break up into different tests
START_TEST(test_mod){

  // CHECK MODS SETUP
  // TODO: put this in setup
  // write param file with mods 
  // good(1mod-1max, 2mods-varmax-1c, 2n)
  // fail(bad-float, bad-letters, missing-entry, 12mods)
  FILE* paramf_1m = fopen("params/1mod", "w");
  fprintf(paramf_1m, "mod=79.9:STY:1\n");
  fprintf(paramf_1m, "mod=22.9:abc:1\n");
  fclose(paramf_1m);

  set_verbosity_level(30);
  // parse command line and param file
  fail_unless(parse_cmd_line_into_params_hash(4, fake_argv, "test-app")==true,
              "Failed to parse param file %s with mods", "1mod");

  // get standard mods
  AA_MOD_T** mod_list = NULL;
  fail_unless( get_aa_mod_list(&mod_list) == 2,
         "Got the incorrect number of mods, %d", get_aa_mod_list(&mod_list));
  // check that each field is correct
  double mass_change = aa_mod_get_mass_change(mod_list[0]);
  fail_unless( mass_change == 79.9,
               "Mod mass change should be 79.9 but is %.2f",
               mass_change);
  mass_change = aa_mod_get_mass_change(mod_list[1]);
  fail_unless( mass_change == 22.9,
               "Mod mass change should be 22.9 but is %.2f",
               mass_change);

  //clean up
  //  free(mod_list);
  free_parameters();
  initialize_parameters();

  // get c mods
  FILE* paramf_2m = fopen("params/2mod", "w");
  fprintf(paramf_2m, "overwrite=T\n");
  fprintf(paramf_2m, "mod=79.9:STY:1\n");
  fprintf(paramf_2m, "cmod=33.3:0\n");
  fclose(paramf_2m);

  free(fake_argv[2]);
  fake_argv[2] = my_copy_string("params/2mod");
  fail_unless(parse_cmd_line_into_params_hash(4, fake_argv, "test-app")==true,
              "Failed to parse param file %s with mods", "2mod");
  fail_unless( get_aa_mod_list(&mod_list) == 1,
            "Got an incorrect number of mods, %d", get_aa_mod_list(&mod_list));
  fail_unless( get_c_mod_list(&mod_list) == 1,
            "Got an incorrect number of cmods, %d", get_c_mod_list(&mod_list));
  mass_change = aa_mod_get_mass_change(mod_list[0]);
  fail_unless( mass_change == 33.3,
               "Mod mass change should be 33.3 but is %.2f",
               mass_change);
  
  //clean up
  //  free(mod_list);
  free_parameters();
  initialize_parameters();

  // get n mods
  FILE* paramf_3m = fopen("params/3mod", "w");
  fprintf(paramf_3m, "nmod=11.1:1\n");
  fprintf(paramf_3m, "nmod=22.2:0\n");
  fprintf(paramf_3m, "nmod=44.4:10\n");
  fprintf(paramf_3m, "overwrite=T\n");
  fclose(paramf_3m);

  free(fake_argv[2]);
  fake_argv[2] = my_copy_string("params/3mod");
  fail_unless(parse_cmd_line_into_params_hash(4, fake_argv, "test-app")==true,
              "Failed to parse param file %s with mods", "2mod");
  fail_unless( get_aa_mod_list(&mod_list) == 0,
            "Got an incorrect number of mods, %d", get_aa_mod_list(&mod_list));
  fail_unless( get_n_mod_list(&mod_list) == 3,
            "Got an incorrect number of nmods, %d", get_n_mod_list(&mod_list));
  mass_change = aa_mod_get_mass_change(mod_list[0]);
  fail_unless( mass_change == 11.1,
               "Mod mass change should be 11.1 but is %.2f",
               mass_change);
  



  // clean up
  free_parameters();
}
END_TEST


START_TEST(test_enzyme){

  // confirm that initialize_parameters set globals correctly
  fail_unless(pre_list_size == 0, 
              "Custom enzyme pre-cleavage site list size should be "
              "initialized to 0 but is %i", pre_list_size);
  fail_unless(post_list_size == 0, 
              "Custom enzyme post-cleavage site list size should be "
              "initialized to 0 but is %i", post_list_size);
  fail_unless(pre_for_inclusion == true,
              "Pre_for_inclusion should be initialized to true but is not.");
  fail_unless(post_for_inclusion == false,
              "Post_for_inclusion should be initialized to false but is not.");

  // try various strings and see that they work
  // can't test error cases b/c parse_ dies on error
  strcpy(enzyme_rule, "[RK]|{P}");
  parse_custom_enzyme(enzyme_rule);
  fail_unless(pre_for_inclusion == true,
              "For rule %s, pre_for_inclusion should be true but is not.",
              enzyme_rule);
  fail_unless(post_for_inclusion == false,
              "For rule %s, post_for_inclusion should be false but is not.",
              enzyme_rule);
  fail_unless(pre_list_size == 2, 
              "For rule %s pre-cleavage site list size should be "
              "2 but is %i", enzyme_rule, pre_list_size);
  fail_unless(post_list_size == 1, 
              "For rule %s, post-cleavage site list size should be "
              "1 but is %i", enzyme_rule, post_list_size);
  fail_unless(pre_cleavage_list != NULL,
              "For rule %s pre-cleavage list should not be NULL.",
              enzyme_rule);
  fail_unless(post_cleavage_list != NULL,
              "For rule %s post-cleavage list should not be NULL.",
              enzyme_rule);
  fail_unless(pre_cleavage_list[0] == 'R',
              "For rule %s, first aa in pre list should be R and is %c",
              enzyme_rule, pre_cleavage_list[0]);
  fail_unless(pre_cleavage_list[1] == 'K',
              "For rule %s, second aa in pre list should be K and is %c",
              enzyme_rule, pre_cleavage_list[0]);
  fail_unless(post_cleavage_list[0] == 'P',
              "For rule %s, first aa in post list should be P and is %c",
              enzyme_rule, pre_cleavage_list[0]);
  // clean up
  free(pre_cleavage_list);
  free(post_cleavage_list);
  pre_cleavage_list = NULL;
  post_cleavage_list = NULL;
  pre_list_size = 0;
  post_list_size = 0;


  strcpy(enzyme_rule, "[PML]|[D]");
  parse_custom_enzyme(enzyme_rule);
  fail_unless(pre_for_inclusion == true,
              "For rule %s, pre_for_inclusion should be true but is not.",
              enzyme_rule);
  fail_unless(post_for_inclusion == true,
              "For rule %s, post_for_inclusion should be true but is not.",
              enzyme_rule);
  fail_unless(pre_list_size == 3, 
              "For rule %s pre-cleavage site list size should be "
              "3 but is %i", enzyme_rule, pre_list_size);
  fail_unless(post_list_size == 1, 
              "For rule %s, post-cleavage site list size should be "
              "1 but is %i", enzyme_rule, post_list_size);
  fail_unless(pre_cleavage_list != NULL,
              "For rule %s pre-cleavage list should not be NULL.",
              enzyme_rule);
  fail_unless(post_cleavage_list != NULL,
              "For rule %s post-cleavage list should not be NULL.",
              enzyme_rule);
  fail_unless(pre_cleavage_list[0] == 'P',
              "For rule %s, first aa in pre list should be P and is %c",
              enzyme_rule, pre_cleavage_list[0]);
  fail_unless(pre_cleavage_list[1] == 'M',
              "For rule %s, second aa in pre list should be M and is %c",
              enzyme_rule, pre_cleavage_list[0]);
  fail_unless(pre_cleavage_list[2] == 'L',
              "For rule %s, third aa in pre list should be L and is %c",
              enzyme_rule, pre_cleavage_list[0]);
  fail_unless(post_cleavage_list[0] == 'D',
              "For rule %s, first aa in post list should be D and is %c",
              enzyme_rule, pre_cleavage_list[0]);
  // clean up
  free(pre_cleavage_list);
  free(post_cleavage_list);
  pre_cleavage_list = NULL;
  post_cleavage_list = NULL;
  pre_list_size = 0;
  post_list_size = 0;


  strcpy(enzyme_rule, "{ABCDEFG}|[HIGK]");
  parse_custom_enzyme(enzyme_rule);
  fail_unless(pre_for_inclusion == false,
              "For rule %s, pre_for_inclusion should be false but is not.",
              enzyme_rule);
  fail_unless(post_for_inclusion == true,
              "For rule %s, post_for_inclusion should be true but is not.",
              enzyme_rule);
  fail_unless(pre_list_size == 7, 
              "For rule %s pre-cleavage site list size should be "
              "7 but is %i", enzyme_rule, pre_list_size);
  fail_unless(post_list_size == 4, 
              "For rule %s, post-cleavage site list size should be "
              "4 but is %i", enzyme_rule, post_list_size);
  fail_unless(pre_cleavage_list != NULL,
              "For rule %s pre-cleavage list should not be NULL.",
              enzyme_rule);
  fail_unless(post_cleavage_list != NULL,
              "For rule %s post-cleavage list should not be NULL.",
              enzyme_rule);
  // clean up
  free(pre_cleavage_list);
  free(post_cleavage_list);
  pre_cleavage_list = NULL;
  post_cleavage_list = NULL;
  pre_list_size = 0;
  post_list_size = 0;

  strcpy(enzyme_rule, "[X]|[X]");
  parse_custom_enzyme(enzyme_rule);
  fail_unless(pre_for_inclusion == false,
              "For rule %s, pre_for_inclusion should be false but is not.",
              enzyme_rule);
  fail_unless(post_for_inclusion == false,
              "For rule %s, post_for_inclusion should be false but is not.",
              enzyme_rule);
  fail_unless(pre_list_size == 0, 
              "For rule %s pre-cleavage site list size should be "
              "7 but is %i", enzyme_rule, pre_list_size);
  fail_unless(post_list_size == 0, 
              "For rule %s, post-cleavage site list size should be "
              "4 but is %i", enzyme_rule, post_list_size);
  fail_unless(pre_cleavage_list == NULL,
              "For rule %s pre-cleavage list should not be NULL.",
              enzyme_rule);
  fail_unless(post_cleavage_list == NULL,
              "For rule %s post-cleavage list should not be NULL.",
              enzyme_rule);
  // clean up
  free(pre_cleavage_list);
  free(post_cleavage_list);
  pre_cleavage_list = NULL;
  post_cleavage_list = NULL;
  pre_list_size = 0;
  post_list_size = 0;

}
END_TEST

Suite* parameter_suite(){
  Suite* s = suite_create("Parameter");
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  tcase_add_test(tc_core, test_mod);
  tcase_add_test(tc_core, test_enzyme);
  tcase_add_checked_fixture(tc_core, param_setup, param_teardown);
  return s;
}













