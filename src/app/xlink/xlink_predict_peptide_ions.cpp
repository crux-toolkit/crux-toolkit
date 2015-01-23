#include "IonSeries.h"
#include "XLinkMatch.h"
#include "XLinkPeptide.h"
#include "LinearPeptide.h"
#include "SelfLoopPeptide.h"

#include <iostream>
#include <fstream>

using namespace std;

int modified_aa_string_length(
  MODIFIED_AA_T* aa_string) {

  if (aa_string == NULL) {
    return 0;
  }

  int ans = 0;

  while (aa_string[ans] != MOD_SEQ_NULL) {
    ans++;
  }

  return ans;
}


int main(int argc, char** argv) {

  char* peptideA = NULL;
  char* peptideB = NULL;
  int posA = 0;
  int posB = 0;
  FLOAT_T linker_mass = 0;
  int charge = 1; 
  bool print_spectrum = false;


  /* Verbosity level for set-up/command line reading */
  set_verbosity_level(CARP_ERROR);
  
  /* Define optional command line arguments */
  const char* option_list[] = {
    "verbosity",
    "version",
    "use-flanking-peaks",
    "print-theoretical-spectrum",
    "parameter-file"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */
  const char* argument_list[] = {"peptide A",
                                 "peptide B",
				 "pos A",
				 "pos B",
				 "charge state",
				 "link mass"};

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
  parse_cmd_line_into_params_hash(argc, argv, "xlink-predict-peptide-ions");

  /* Set verbosity */
  set_verbosity_level(get_int_parameter("verbosity"));

  /* Get Arguments */
  linker_mass = get_double_parameter("link mass");
  charge = get_int_parameter("charge state");

  peptideA = get_string_parameter("peptide A");
  peptideB = get_string_parameter("peptide B");
  
  posA = get_int_parameter("pos A");
  posB = get_int_parameter("pos B");

  print_spectrum = get_boolean_parameter("print-theoretical-spectrum");

  XLinkPeptide::setLinkerMass(linker_mass);

  XLinkMatch* linked_peptide = NULL;

  cerr <<"creating peptide"<<endl;

  cerr <<"Converting peptideA to MOD_AA"<<endl;

  MODIFIED_AA_T* mod_seqA = NULL;

  int len = strlen(peptideA);

  convert_to_mod_aa_seq(peptideA, &mod_seqA);

  cerr<<"Length of modified string:"<<modified_aa_string_length(mod_seqA)<<endl;

  //convert back.
  
  char* temp = 
    modified_aa_string_to_string_with_masses(mod_seqA, modified_aa_string_length(mod_seqA), MOD_MASS_ONLY);

  cerr <<"orig:"<<peptideA<<":"<<len<<" convert:"<<temp<<":"<<strlen(temp)<<endl;


  if (string(peptideB) == string("NULL")) {
    if (posA == -1 || posB == -1) {
      cout<<"Creating linear peptide"<<endl;
      linked_peptide = new LinearPeptide(peptideA);
    } else {
      cout<<"Creating selfloop peptide"<<endl;
      linked_peptide = new SelfLoopPeptide(peptideA, posA-1, posB-1);
    }
  } else {
    cout <<"Creating XLinkPeptide"<<endl;
    linked_peptide = new XLinkPeptide(peptideA, peptideB, posA-1, posB-1);
  }

  cerr << "Printing stuff"<<endl;
  cerr << "precursor: " << linked_peptide -> getSequenceString() << endl;
  cerr << "mass:" << linked_peptide -> getMass(MONO)<<" "<< linked_peptide -> getMass(AVERAGE) <<endl;;
  cerr << "charge:" << charge << endl;
  cerr << "link mass:"<< linker_mass << endl;
  cerr << "print_spectrum:"<< print_spectrum << endl;

  int max_charge = min(get_max_ion_charge_parameter("max-ion-charge"), charge);

  IonConstraint* ion_constraint = new IonConstraint(MONO, max_charge, BY_ION,  false);

  IonSeries* ion_series = new IonSeries(ion_constraint, charge);

  linked_peptide->predictIons(ion_series, charge);

  cout << "mz\t" 
       << "mass\t"
       << "charge\t"
       << "cleavage idx\t"
       << "sequence\t"
       << "ion type" << endl;

  for (IonIterator ion_iter = ion_series->begin();
    ion_iter != ion_series->end();
    ++ion_iter) {

    Ion* ion = *ion_iter;
    FLOAT_T mz = ion->getMassZ();
    FLOAT_T mass = ion->getMassFromMassZ();
    int charge = ion->getCharge();
    string sequence = linked_peptide->getIonSequence(ion);
    int cleavage_idx = ion->getCleavageIdx();
    ION_TYPE_T ion_type = ion->getType();

    string ion_type_string = "";
    switch (ion_type) {
      case B_ION:
        ion_type_string = "b";
        break;
      case Y_ION:
        ion_type_string = "y";
        break;
      default:
        ion_type_string = "u";
    }
  

    cout << mz << "\t"
	 << mass << "\t"
	 << charge << "\t"
         << cleavage_idx << "\t"
	 << sequence << "\t" 
         << ion_type_string << endl;
  }


  delete ion_series;
  delete linked_peptide;

  return 0;
}
#include "IonSeries.h"
#include "XLinkMatch.h"
#include "XLinkPeptide.h"
#include "LinearPeptide.h"
#include "SelfLoopPeptide.h"

#include <iostream>
#include <fstream>

using namespace std;

int modified_aa_string_length(
  MODIFIED_AA_T* aa_string) {

  if (aa_string == NULL) {
    return 0;
  }

  int ans = 0;

  while (aa_string[ans] != MOD_SEQ_NULL) {
    ans++;
  }

  return ans;
}


int main(int argc, char** argv) {

  char* peptideA = NULL;
  char* peptideB = NULL;
  int posA = 0;
  int posB = 0;
  FLOAT_T linker_mass = 0;
  int charge = 1; 
  bool print_spectrum = false;


  /* Verbosity level for set-up/command line reading */
  set_verbosity_level(CARP_ERROR);
  
  /* Define optional command line arguments */
  const char* option_list[] = {
    "verbosity",
    "version",
    "use-flanking-peaks",
    "print-theoretical-spectrum",
    "parameter-file"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */
  const char* argument_list[] = {"peptide A",
                                 "peptide B",
				 "pos A",
				 "pos B",
				 "charge state",
				 "link mass"};

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
  parse_cmd_line_into_params_hash(argc, argv, "xlink-predict-peptide-ions");

  /* Set verbosity */
  set_verbosity_level(get_int_parameter("verbosity"));

  /* Get Arguments */
  linker_mass = get_double_parameter("link mass");
  charge = get_int_parameter("charge state");

  peptideA = get_string_parameter("peptide A");
  peptideB = get_string_parameter("peptide B");
  
  posA = get_int_parameter("pos A");
  posB = get_int_parameter("pos B");

  print_spectrum = get_boolean_parameter("print-theoretical-spectrum");

  XLinkPeptide::setLinkerMass(linker_mass);

  XLinkMatch* linked_peptide = NULL;

  cerr <<"creating peptide"<<endl;

  cerr <<"Converting peptideA to MOD_AA"<<endl;

  MODIFIED_AA_T* mod_seqA = NULL;

  int len = strlen(peptideA);

  convert_to_mod_aa_seq(peptideA, &mod_seqA);

  cerr<<"Length of modified string:"<<modified_aa_string_length(mod_seqA)<<endl;

  //convert back.
  
  char* temp = 
    modified_aa_string_to_string_with_masses(mod_seqA, modified_aa_string_length(mod_seqA), MOD_MASS_ONLY);

  cerr <<"orig:"<<peptideA<<":"<<len<<" convert:"<<temp<<":"<<strlen(temp)<<endl;


  if (string(peptideB) == string("NULL")) {
    if (posA == -1 || posB == -1) {
      cout<<"Creating linear peptide"<<endl;
      linked_peptide = new LinearPeptide(peptideA);
    } else {
      cout<<"Creating selfloop peptide"<<endl;
      linked_peptide = new SelfLoopPeptide(peptideA, posA-1, posB-1);
    }
  } else {
    cout <<"Creating XLinkPeptide"<<endl;
    linked_peptide = new XLinkPeptide(peptideA, peptideB, posA-1, posB-1);
  }

  cerr << "Printing stuff"<<endl;
  cerr << "precursor: " << linked_peptide -> getSequenceString() << endl;
  cerr << "mass:" << linked_peptide -> getMass(MONO)<<" "<< linked_peptide -> getMass(AVERAGE) <<endl;;
  cerr << "charge:" << charge << endl;
  cerr << "link mass:"<< linker_mass << endl;
  cerr << "print_spectrum:"<< print_spectrum << endl;

  int max_charge = min(get_max_ion_charge_parameter("max-ion-charge"), charge);

  IonConstraint* ion_constraint = new IonConstraint(MONO, max_charge, BY_ION,  false);

  IonSeries* ion_series = new IonSeries(ion_constraint, charge);

  linked_peptide->predictIons(ion_series, charge);

  cout << "mz\t" 
       << "mass\t"
       << "charge\t"
       << "cleavage idx\t"
       << "sequence\t"
       << "ion type" << endl;

  for (IonIterator ion_iter = ion_series->begin();
    ion_iter != ion_series->end();
    ++ion_iter) {

    Ion* ion = *ion_iter;
    FLOAT_T mz = ion->getMassZ();
    FLOAT_T mass = ion->getMassFromMassZ();
    int charge = ion->getCharge();
    string sequence = linked_peptide->getIonSequence(ion);
    int cleavage_idx = ion->getCleavageIdx();
    ION_TYPE_T ion_type = ion->getType();

    string ion_type_string = "";
    switch (ion_type) {
      case B_ION:
        ion_type_string = "b";
        break;
      case Y_ION:
        ion_type_string = "y";
        break;
      default:
        ion_type_string = "u";
    }
  

    cout << mz << "\t"
	 << mass << "\t"
	 << charge << "\t"
         << cleavage_idx << "\t"
	 << sequence << "\t" 
         << ion_type_string << endl;
  }


  delete ion_series;
  delete linked_peptide;

  return 0;
}
#include "IonSeries.h"
#include "XLinkMatch.h"
#include "XLinkPeptide.h"
#include "LinearPeptide.h"
#include "SelfLoopPeptide.h"

#include <iostream>
#include <fstream>

using namespace std;

int modified_aa_string_length(
  MODIFIED_AA_T* aa_string) {

  if (aa_string == NULL) {
    return 0;
  }

  int ans = 0;

  while (aa_string[ans] != MOD_SEQ_NULL) {
    ans++;
  }

  return ans;
}


int main(int argc, char** argv) {

  char* peptideA = NULL;
  char* peptideB = NULL;
  int posA = 0;
  int posB = 0;
  FLOAT_T linker_mass = 0;
  int charge = 1; 
  bool print_spectrum = false;


  /* Verbosity level for set-up/command line reading */
  set_verbosity_level(CARP_ERROR);
  
  /* Define optional command line arguments */
  const char* option_list[] = {
    "verbosity",
    "version",
    "use-flanking-peaks",
    "print-theoretical-spectrum",
    "parameter-file"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */
  const char* argument_list[] = {"peptide A",
                                 "peptide B",
				 "pos A",
				 "pos B",
				 "charge state",
				 "link mass"};

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
  parse_cmd_line_into_params_hash(argc, argv, "xlink-predict-peptide-ions");

  /* Set verbosity */
  set_verbosity_level(get_int_parameter("verbosity"));

  /* Get Arguments */
  linker_mass = get_double_parameter("link mass");
  charge = get_int_parameter("charge state");

  peptideA = get_string_parameter("peptide A");
  peptideB = get_string_parameter("peptide B");
  
  posA = get_int_parameter("pos A");
  posB = get_int_parameter("pos B");

  print_spectrum = get_boolean_parameter("print-theoretical-spectrum");

  XLinkPeptide::setLinkerMass(linker_mass);

  XLinkMatch* linked_peptide = NULL;

  cerr <<"creating peptide"<<endl;

  cerr <<"Converting peptideA to MOD_AA"<<endl;

  MODIFIED_AA_T* mod_seqA = NULL;

  int len = strlen(peptideA);

  convert_to_mod_aa_seq(peptideA, &mod_seqA);

  cerr<<"Length of modified string:"<<modified_aa_string_length(mod_seqA)<<endl;

  //convert back.
  
  char* temp = 
    modified_aa_string_to_string_with_masses(mod_seqA, modified_aa_string_length(mod_seqA), MOD_MASS_ONLY);

  cerr <<"orig:"<<peptideA<<":"<<len<<" convert:"<<temp<<":"<<strlen(temp)<<endl;


  if (string(peptideB) == string("NULL")) {
    if (posA == -1 || posB == -1) {
      cout<<"Creating linear peptide"<<endl;
      linked_peptide = new LinearPeptide(peptideA);
    } else {
      cout<<"Creating selfloop peptide"<<endl;
      linked_peptide = new SelfLoopPeptide(peptideA, posA-1, posB-1);
    }
  } else {
    cout <<"Creating XLinkPeptide"<<endl;
    linked_peptide = new XLinkPeptide(peptideA, peptideB, posA-1, posB-1);
  }

  cerr << "Printing stuff"<<endl;
  cerr << "precursor: " << linked_peptide -> getSequenceString() << endl;
  cerr << "mass:" << linked_peptide -> getMass(MONO)<<" "<< linked_peptide -> getMass(AVERAGE) <<endl;;
  cerr << "charge:" << charge << endl;
  cerr << "link mass:"<< linker_mass << endl;
  cerr << "print_spectrum:"<< print_spectrum << endl;

  int max_charge = min(get_max_ion_charge_parameter("max-ion-charge"), charge);

  IonConstraint* ion_constraint = new IonConstraint(MONO, max_charge, BY_ION,  false);

  IonSeries* ion_series = new IonSeries(ion_constraint, charge);

  linked_peptide->predictIons(ion_series, charge);

  cout << "mz\t" 
       << "mass\t"
       << "charge\t"
       << "cleavage idx\t"
       << "sequence\t"
       << "ion type" << endl;

  for (IonIterator ion_iter = ion_series->begin();
    ion_iter != ion_series->end();
    ++ion_iter) {

    Ion* ion = *ion_iter;
    FLOAT_T mz = ion->getMassZ();
    FLOAT_T mass = ion->getMassFromMassZ();
    int charge = ion->getCharge();
    string sequence = linked_peptide->getIonSequence(ion);
    int cleavage_idx = ion->getCleavageIdx();
    ION_TYPE_T ion_type = ion->getType();

    string ion_type_string = "";
    switch (ion_type) {
      case B_ION:
        ion_type_string = "b";
        break;
      case Y_ION:
        ion_type_string = "y";
        break;
      default:
        ion_type_string = "u";
    }
  

    cout << mz << "\t"
	 << mass << "\t"
	 << charge << "\t"
         << cleavage_idx << "\t"
	 << sequence << "\t" 
         << ion_type_string << endl;
  }


  delete ion_series;
  delete linked_peptide;

  return 0;
}
