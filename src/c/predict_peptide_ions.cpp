/*************************************************************************//**
 * \file predict_peptide_ions.cpp
 * AUTHOR: Chris Park
 * CREATE DATE: 10/05 2006
 * \brief Object for given a peptide and a charge state, predict
 * the ions 
 ****************************************************************************/
#include <algorithm>
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
#include "Ion.h"
#include "IonSeries.h"
#include "IonConstraint.h"

using namespace std;

int main(int argc, char** argv){

  /* Define optional and required command line arguments */
  const char* option_list[] = {
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
  int num_options = sizeof(option_list) / sizeof(char*);

  const char* argument_list[] = {
    "peptide sequence",
    "charge state"
  };
  int num_arguments = sizeof(argument_list) / sizeof(char*);

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
  bool use_precursor_ions = get_boolean_parameter("precursor-ions");
  int isotope_count = get_int_parameter("isotope");
  bool is_flanking = get_boolean_parameter("flanking");
  const char* max_ion_charge = get_string_parameter_pointer("max-ion-charge");
  int nh3_count = get_int_parameter("nh3");
  int h2o_count = get_int_parameter("h2o");

  int neutral_loss_count[MAX_MODIFICATIONS];
  bool is_modification = false;

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

  int max_charge = get_max_ion_charge_parameter("max-ion-charge");
  max_charge = min(max_charge, charge_state);
  // create ion_constraint
  MASS_TYPE_T frag_masses = get_mass_type_parameter("fragment-mass");
  
  IonConstraint* ion_constraint = 
  //  new_ion_constraint(MONO, max_charge, ion_type, use_precursor_ions);
    new IonConstraint(frag_masses, max_charge, ion_type, use_precursor_ions);

   
   // set ion_constraint3 modification counts, if modifications should occur
  if(is_modification){
    ion_constraint->setModification(NH3, neutral_loss_count[NH3]);
    ion_constraint->setModification(H2O, neutral_loss_count[H2O]);
    ion_constraint->setModification(ISOTOPE, neutral_loss_count[ISOTOPE]);
    ion_constraint->setModification(FLANK, neutral_loss_count[FLANK]);
  }

  // create ion_series
  IonSeries* ion_series = new IonSeries(peptide_sequence, 
                                             charge_state, ion_constraint);
   
  // now predict ions
  ion_series->predictIons();

  // print settings
  printf("# PEPTIDE: %s\n",peptide_sequence);
  printf("# AVERAGE: %f MONO:%f\n",
    Peptide::calcSequenceMass(peptide_sequence, AVERAGE),
    Peptide::calcSequenceMass(peptide_sequence, MONO));
  printf("# CHARGE: %d\n", charge_state);
  printf("# MAX-ION-CHRAGE: %s\n", max_ion_charge);
  printf("# NH3 modification: %d\n", neutral_loss_count[NH3]);
  printf("# H2O modification: %d\n", neutral_loss_count[H2O] );
  printf("# ISOTOPE modification: %d\n", neutral_loss_count[ISOTOPE] );
  printf("# FLANK modification: %d\n", neutral_loss_count[FLANK]);
  // print ions
  ion_series->print(stdout);

  // free
  IonConstraint::free(ion_constraint);

  delete ion_series;
  carp(CARP_INFO, "crux-predict-peptide-ions finished");
  exit(0);
}

