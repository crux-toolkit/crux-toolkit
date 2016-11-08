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
#include "io/carp.h"
#include "util/crux-utils.h"
#include "model/objects.h"
#include "parameter.h"
#include "model/Ion.h"
#include "model/IonSeries.h"
#include "model/IonConstraint.h"
#include "model/Peptide.h"
#include "util/Params.h"

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
  /* Get Arguments */
  string peptide_sequence = Params::GetString("peptide sequence");
  int charge_state = Params::GetInt("charge state");

  /* Get Options */
  ION_TYPE_T ion_type;
  string_to_ion_type(Params::GetString("primary-ions"), &ion_type);
  bool use_precursor_ions = Params::GetBool("precursor-ions");
  int isotope_count = Params::GetInt("isotope");
  bool is_flanking = Params::GetBool("flanking");
  string max_ion_charge = Params::GetString("max-ion-charge");

  int nh3_count = Params::GetInt("nh3");
  int h2o_count = Params::GetInt("h2o");

  int neutral_loss_count[MAX_MODIFICATIONS];
  bool is_modification = false;

  // check peptide sequence
  if (!valid_peptide_sequence(peptide_sequence)) {
    carp(CARP_FATAL, "The peptide sequence '%s' is not valid", peptide_sequence.c_str());
  }

  // neutral_losses
  // initialize
  int modification_idx = 0;
  for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx) {
    neutral_loss_count[modification_idx] = 0;
  }

  is_modification = (nh3_count || h2o_count || isotope_count || is_flanking);

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
  if (is_modification) {
    ion_constraint->setModification(NH3, neutral_loss_count[NH3]);
    ion_constraint->setModification(H2O, neutral_loss_count[H2O]);
    ion_constraint->setModification(ISOTOPE, neutral_loss_count[ISOTOPE]);
    ion_constraint->setModification(FLANK, neutral_loss_count[FLANK]);
  }

  // create ion_series
  IonSeries* ion_series = new IonSeries(peptide_sequence, charge_state, ion_constraint);
   
  // now predict ions
  ion_series->predictIons();

  // print settings
  printf("# PEPTIDE: %s\n", peptide_sequence.c_str());
  printf("# AVERAGE: %f MONO:%f\n",
    Peptide::calcSequenceMass(peptide_sequence, AVERAGE),
    Peptide::calcSequenceMass(peptide_sequence, MONO));
  printf("# CHARGE: %d\n", charge_state);
  printf("# MAX-ION-CHARGE: %s\n", max_ion_charge.c_str());
  printf("# NH3 modification: %d\n", neutral_loss_count[NH3]);
  printf("# H2O modification: %d\n", neutral_loss_count[H2O]);
  printf("# ISOTOPE modification: %d\n", neutral_loss_count[ISOTOPE]);
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
string PredictPeptideIons::getName() const {
  return "predict-peptide-ions";
}

/**
 * \returns The description for PredictPeptideIons.
 */
string PredictPeptideIons::getDescription() const {
  return
    "[[nohtml:Given a peptide and a charge state, predict the m/z values of "
    "the resulting fragment ions.]]"
    "[[html:<p>Given a peptide and a charge state, predict the corresponding "
    "fragment ions according to the provided options.</p>]]";
}

/**
 * \returns The command arguments
 */
vector<string> PredictPeptideIons::getArgs() const {
  string arr[] = {
    "peptide sequence",
    "charge state"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns The command options
 */
vector<string> PredictPeptideIons::getOptions() const {
  string arr[] = {
    "primary-ions",
    "precursor-ions",
    "isotope",
    "flanking",
    "max-ion-charge",
    "fragment-mass",
    "nh3",
    "h2o"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns The command outputs
 */
vector< pair<string, string> > PredictPeptideIons::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("stdout",
    "a series of lines, describing how the ions were predicted. <blockquote "
    "style=\"font-family: monospace;\"># PEPTIDE: &lt;peptide sequence&gt;<br>"
    "# CHARGE: &lt;peptide charge&gt;</blockquote>"
    "The program then prints a one-"
    "line header that labels the tab-delimited columns, followed by one line "
    "per ion:<blockquote style=\"font-family: monospace;\">m/z&nbsp;mass&nbsp;"
    "charge&nbsp;ion-series&nbsp;peptide-bond-index&nbsp;nh3&nbsp;h2o&nbsp;"
    "isotope&nbsp;flank</blockquote>"
    "<p>The columns contain the following values:"
    "<ul><li>&lt;m/z&gt; is the ion's mass-to-charge</li>"
    "<li>&lt;mass&gt; is the ion's (charged) mass</li>"
    "<li>&lt;charge&gt; is the ion's charge e.g. 1,2,3</li>"
    "<li>&lt;ion-type&gt; is a string representing the series 'a', 'b', 'c', "
    "'x', 'y', 'z', and 'p'</li>"
    "<li>&lt;peptide-bond-index&gt; is in [1...n), where n is "
    "peptide-length. Consistent with standard mass spec terminology, b-1 "
    "corresponds to a prefix ion of a single amino acid, and y-1 "
    "corresponds to a suffix ion of a single amino acid.<li>"
    "<li>&lt;nh3&gt; is the number of NH3 modifications</li>"
    "<li>&lt;h2o&gt; is the number of H2O Modifications</li>"
    "<li>&lt;isotope&gt; is the number of adjacent isotopic peaks</li>"
    "<li>&lt;flank&gt; is the number of flanking ions</li>"
    "</ul></p>"));
  return outputs;
}

/**
 * \returns The enum of the application, default PREDICT_PEPTIDE_IONS_COMMAND.
 */
COMMAND_T PredictPeptideIons::getCommand() const {
  return PREDICT_PEPTIDE_IONS_COMMAND;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

