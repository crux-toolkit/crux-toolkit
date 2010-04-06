/*************************************************************************//**
 * \file predict_peptide_ions.cpp
 * AUTHOR: Chris Park
 * CREATE DATE: 10/05 2006
 * DESCRIPTION: Object for given a peptide and a charge state, predict
 * the ions 
 *
 * REVISION: $ $ 
 ****************************************************************************/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <unistd.h>
#include "carp.h"
#include "crux-utils.h"
#include "objects.h"
#include "parameter.h"
#include "ion.h"
#include "ion_series.h"

#define NUM_PREDICT_OPTIONS 9
#define NUM_PREDICT_ARGUMENTS 2

int main(int argc, char** argv){

  /* Define optional and required command line arguments */
  int num_options = NUM_PREDICT_OPTIONS;
  const char* option_list[NUM_PREDICT_OPTIONS] = {
    "version",
    "primary-ions",
    "precursor-ions",
    "neutral-losses",
    "isotope",
    "flanking",
    "max-ion-charge",
    "nh3",
    "h2o"
  };

  int num_arguments = NUM_PREDICT_ARGUMENTS;
  const char* argument_list[NUM_PREDICT_ARGUMENTS] = {
    "peptide sequence",
    "charge state"
  };

  /* for debugging of parameter processing */
  //set_verbosity_level( CARP_DETAILED_DEBUG );
  set_verbosity_level( CARP_ERROR );

  /* Set default values for parameters in parameter.c */
  initialize_parameters();

  /* Define optional and required command line arguments */
  select_cmd_line_options( option_list, num_options );
  select_cmd_line_arguments( argument_list, num_arguments);

  /* Parse the command line, including the optional params file */
  /* does sytnax, type, bounds checking and dies if neccessessary */
  parse_cmd_line_into_params_hash(argc, argv, "crux-predict-peptide-ions");

  /* Set verbosity */
  set_verbosity_level(get_int_parameter("verbosity"));

  /* Get Arguments */
  const char* peptide_sequence = get_string_parameter_pointer("peptide sequence");
  int charge_state = get_int_parameter("charge state");

  /* Get Options */
  ION_TYPE_T ion_type = get_ion_type_parameter("primary-ions");
  BOOLEAN_T use_precursor_ions = get_boolean_parameter("precursor-ions");
  int isotope_count = get_int_parameter("isotope");
  BOOLEAN_T is_flanking = get_boolean_parameter("flanking");
  const char* max_ion_charge = get_string_parameter_pointer("max-ion-charge");
  int nh3_count = get_int_parameter("nh3");
  int h2o_count = get_int_parameter("h2o");

  int neutral_loss_count[MAX_MODIFICATIONS];
  BOOLEAN_T is_modification = FALSE;
  int max_charge = charge_state;

  // check peptide sequence
  if(!valid_peptide_sequence(peptide_sequence)){
    carp(CARP_FATAL, "The peptide sequence '%s' is not valid", 
         peptide_sequence);
  }

   // neutral_losses
   // initialize
   int modification_idx = 0;
   for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
     neutral_loss_count[modification_idx] = 0;
   }

   is_modification = (nh3_count || 
                      h2o_count || 
                      is_flanking);

   neutral_loss_count[NH3] = nh3_count;
   neutral_loss_count[H2O] = h2o_count;
   neutral_loss_count[FLANK] = (int)is_flanking;
   neutral_loss_count[ISOTOPE] = isotope_count;

   // max_ion_charge
   if(strcmp(max_ion_charge, "1")== 0){
     max_charge = 1;
   }
   else if(strcmp(max_ion_charge, "2")== 0){
     max_charge = 2;
   }
   else if(strcmp(max_ion_charge, "3")== 0){
     max_charge = 3;
   }
   else if(strcmp(max_ion_charge, "peptide")== 0){
     max_charge = charge_state;
   }
   else{
     carp(CARP_FATAL, 
          "max-ion-charge option must be 1,2,3 or peptide. '%s' is not legal",
          max_charge);
   }
   
   // create ion_constraint
   MASS_TYPE_T frag_masses = get_mass_type_parameter("fragment-mass");
   ION_CONSTRAINT_T* ion_constraint = 
     //  new_ion_constraint(MONO, max_charge, ion_type, use_precursor_ions);
     new_ion_constraint(frag_masses, max_charge, ion_type, use_precursor_ions);

   
   // set ion_constraint3 modification counts, if modifications should occur
   if(is_modification){
     set_ion_constraint_modification( ion_constraint, NH3, 
                                      neutral_loss_count[NH3]);
     set_ion_constraint_modification( ion_constraint, H2O, 
                                      neutral_loss_count[H2O]);
     set_ion_constraint_modification( ion_constraint, ISOTOPE, 
                                      neutral_loss_count[ISOTOPE]);
     set_ion_constraint_modification( ion_constraint, FLANK, 
                                      neutral_loss_count[FLANK]);
   }

   // create ion_series
   ION_SERIES_T* ion_series = new_ion_series(peptide_sequence, 
                                             charge_state, ion_constraint);
   
   // now predict ions
   predict_ions(ion_series);
   
   // print settings
   printf("# PEPTIDE: %s\n",peptide_sequence);
   printf("# AVERAGE: %f MONO:%f",calc_sequence_mass(peptide_sequence, AVERAGE),calc_sequence_mass(peptide_sequence, MONO));
   printf("# CHARGE: %d\n", charge_state);
   printf("# MAX-ION-CHRAGE: %s\n", max_ion_charge);
   printf("# NH3 modification: %d\n", neutral_loss_count[NH3]);
   printf("# H2O modification: %d\n", neutral_loss_count[H2O] );
   printf("# ISOTOPE modification: %d\n", neutral_loss_count[ISOTOPE] );
   printf("# FLANK modification: %d\n", neutral_loss_count[FLANK]);

   // print ions
   print_ion_series(ion_series, stdout);

   // free
   free_ion_constraint(ion_constraint);
   free_ion_series(ion_series);

   carp(CARP_INFO, "crux-predict-peptide-ions finished");
 exit(0);
}
