#include "LinkedIonSeries.h"
#include "xhhc_scorer.h"
#include "LinkedPeptide.h"
#include "XHHC_Peptide.h"

#include <iostream>
#include <fstream>

#define NUM_ARGUMENTS 6
#define NUM_OPTIONS 4
using namespace std;

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
  int num_options = NUM_OPTIONS;
  const char* option_list[NUM_OPTIONS] = {
    "verbosity",
    "version",
    "use-flanking-peaks",
    "print-theoretical-spectrum"
  };

  /* Define required command line arguments */
  int num_arguments = NUM_ARGUMENTS;
  const char* argument_list[NUM_ARGUMENTS] = {"peptide A",
                                              "peptide B",
					      "pos A",
					      "pos B",
					      "charge state",
					      "link mass"};



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

  LinkedPeptide::setLinkerMass(linker_mass); 
  LinkedPeptide linked_peptide;

  if (string(peptideB) == "NULL") {
    linked_peptide = LinkedPeptide(peptideA, NULL, posA-1, posB-1, charge);
  }
  else {
    linked_peptide = LinkedPeptide( peptideA, peptideB, posA-1, posB-1, charge);  
  }
 
  ostringstream oss;
  oss << linked_peptide;
  string precursor_sequence = oss.str();
  carp(CARP_INFO, "XLink Product:%s", precursor_sequence.c_str());
  carp(CARP_INFO, "Charge:%i", charge);
  carp(CARP_INFO, "Linker Mass:%g", linker_mass);
  carp(CARP_INFO, "Product Mass: %g (MONO) %g (AVERAGE)", linked_peptide.getMass(MONO), linked_peptide.getMass(AVERAGE));
  carp(CARP_INFO, "====================");

  LinkedIonSeries ion_series(charge);
  ion_series.addLinkedIons(linked_peptide);
  ion_series.print(); 

  if (print_spectrum) {
    carp(CARP_INFO, "Writing XCORR theoretical spectrum to theoretical.out");
    map<int, FLOAT_T> theoretical;
    XHHC_Scorer::xlinkCreateMapTheoretical(ion_series,
					 theoretical);

    ofstream fout("theoretical.out");
    fout <<"> "<< linked_peptide<<"\t"<<linked_peptide.getMZ(MONO)<<endl;
    
    map<int, FLOAT_T>::iterator iter;

    for (iter = theoretical.begin();
	 iter != theoretical.end();
	 ++iter) {
      fout << iter -> first << "\t";
      fout << iter -> second << "\t" << endl;
    }

    fout.close();
  }
  return 0;
}
