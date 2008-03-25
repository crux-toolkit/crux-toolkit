#include "check-parameter.h"
#include "../parameter.h"
#include "../carp.h"
#include "../modifications.h"

START_TEST(test_create){
  // set up
  initialize_parameters();

  // check getters for each type (TODO)
  fail_unless( get_boolean_parameter("overwrite") == FALSE, 
              "Does not return the bool parameter %s correctly", "overwrite");

  // select command line options and arguments (try non-existant ones)
  char* ops[1] = {"parameter-file"};
  fail_unless( select_cmd_line_options(ops, 1) == TRUE,
               "Failed to select command line options");
  char* args[1] = {"protein input"};
  fail_unless( select_cmd_line_arguments(args, 1) == TRUE, 
    "Failed to select command line arguments");
  // parse command line (with and without param file)
  // look for correct option (cmd line overrides param)

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
  char* fake_argv[4] = {"test-app", "--parameter-file", "params/1mod", "prot"};
  // parse command line and param file
  fail_unless(parse_cmd_line_into_params_hash(4, fake_argv, "test-app")==TRUE,
              "Failed to parse param file %s with mods", "1mod");

  // get standard mods
  AA_MOD_T* mod_list = NULL;
  fail_unless( get_aa_mod_list(&mod_list) == 2,
               "Got the incorrect number of mods, %d", get_aa_mod_list(&mod_list));
  // check that each field is correct
  fail_unless( mod_list[0].mass_change == 79.9,
               "Mod mass change should be 79.9 but is %.2f",
               mod_list[0].mass_change);
  fail_unless( mod_list[1].mass_change == 22.9,
               "Mod mass change should be 22.9 but is %.2f",
               mod_list[1].mass_change);
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

  char* fake_argv2[4] = {"test-ap", "--parameter-file", "params/2mod", "prot"};
  fail_unless(parse_cmd_line_into_params_hash(4, fake_argv2, "test-app")==TRUE,
              "Failed to parse param file %s with mods", "2mod");
  fail_unless( get_aa_mod_list(&mod_list) == 1,
            "Got an incorrect number of mods, %d", get_aa_mod_list(&mod_list));
  fail_unless( get_c_mod_list(&mod_list) == 1,
            "Got an incorrect number of cmods, %d", get_c_mod_list(&mod_list));
  fail_unless( mod_list[0].mass_change == 33.3,
               "Mod mass change should be 33.3 but is %.2f",
               mod_list[0].mass_change);

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

  char* fake_argv3[4] = {"test-ap", "--parameter-file", "params/3mod", "prot"};
  fail_unless(parse_cmd_line_into_params_hash(4, fake_argv3, "test-app")==TRUE,
              "Failed to parse param file %s with mods", "2mod");
  fail_unless( get_aa_mod_list(&mod_list) == 0,
            "Got an incorrect number of mods, %d", get_aa_mod_list(&mod_list));
  fail_unless( get_n_mod_list(&mod_list) == 3,
            "Got an incorrect number of nmods, %d", get_n_mod_list(&mod_list));
  fail_unless( mod_list[0].mass_change == 11.1,
               "Mod mass change should be 11.1 but is %.2f",
               mod_list[0].mass_change);




  // clean up
  free_parameters();
}
END_TEST

Suite* parameter_suite(){
  Suite* s = suite_create("Parameter\n");
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  //tcase_add_check_fixture(tc_core, setup, teardown);
  return s;
}













