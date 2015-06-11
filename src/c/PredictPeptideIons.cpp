/**
 * \file PredictPeptideIons.cpp
 *
 * AUTHOR: Manijeh Naseri
 * CREATE DATE: January 27, 2012
 *
 * DESCRIPTION: Main method for the predict-peptide-ions.
 *              Given a peptide sequence, and a charge state, predict
 *              the fragmentation ions.
 */

#include "PredictPeptideIons.h"

#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include "carp.h"
#include "crux-utils.h"
#include "objects.h"
#include "parameter.h"
#include "Ion.h"
#include "IonSeries.h"
#include "IonConstraint.h"
#include "Peptide.h"

using namespace std;
using namespace Crux;

/**
 * \returns A blank PredictPeptideIons object.
 */
PredictPeptideIons::PredictPeptideIons() {

}

/**
 * Destructor
 */
PredictPeptideIons::~PredictPeptideIons() {
}

/**
 * Main method for PredictPeptideIons.
 */
int PredictPeptideIons::main(int argc, char** argv) {

  /* Define optional and required command line arguments */
  const char* option_list[] = {
    "primary-ions",
    "precursor-ions",
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

  initialize(argument_list,
             num_arguments,
             option_list,
             num_options,
             argc, argv);

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
                      isotope_count || 
                      is_flanking);

   neutral_loss_count[NH3] = nh3_count;
   neutral_loss_count[H2O] = h2o_count;
   neutral_loss_count[FLANK] = (int)is_flanking;
   neutral_loss_count[ISOTOPE] = isotope_count;

  // create ion_constraint
  MASS_TYPE_T frag_masses = get_mass_type_parameter("fragment-mass");
  
  IonConstraint* ion_constraint = 
  //  new_ion_constraint(MONO, max_charge, ion_type, use_precursor_ions);
    new IonConstraint(frag_masses, charge_state, ion_type, use_precursor_ions);

   
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
  printf("# MAX-ION-CHARGE: %s\n", max_ion_charge);
  printf("# NH3 modification: %d\n", neutral_loss_count[NH3]);
  printf("# H2O modification: %d\n", neutral_loss_count[H2O] );
  printf("# ISOTOPE modification: %d\n", neutral_loss_count[ISOTOPE] );
  printf("# FLANK modification: %d\n", neutral_loss_count[FLANK]);
  // print ions
  ion_series->print(stdout);

  // free
  IonConstraint::free(ion_constraint);

  delete ion_series;

  return 0;
}

/**
 * \returns The command name for PredictPeptideIons.
 */
string PredictPeptideIons::getName() {
  return "predict-peptide-ions";
}

/**
 * \returns The description for PredictPeptideIons.
 */
string PredictPeptideIons::getDescription() {
  return "Given a peptide and a charge state, predict the m/z values of the "
         "resulting fragment ions.";
}

/**
 * \returns The enum of the application, default PREDICT_PEPTIDE_IONS_COMMAND.
 */
COMMAND_T PredictPeptideIons::getCommand() {
  return PREDICT_PEPTIDE_IONS_COMMAND;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

