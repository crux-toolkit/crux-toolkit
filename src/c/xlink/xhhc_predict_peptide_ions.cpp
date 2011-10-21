#include "xhhc_ion_series.h"
#include "xhhc_scorer.h"

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

  LinkedPeptide::linker_mass = linker_mass; 
  LinkedPeptide linked_peptide;

  if (string(peptideB) == "NULL")
    linked_peptide = LinkedPeptide(peptideA, NULL, posA, posB, charge);
  else  
    linked_peptide = LinkedPeptide( peptideA, peptideB, posA, posB, charge);  

  cout << "precursor: " << linked_peptide << endl;
  cout << "charge:" << charge << endl;
  cout << "link mass:"<< linker_mass << endl;
  cout << "print_spectrum:"<< print_spectrum << endl;


  LinkedIonSeries ion_series;
  cout <<"Adding linked ions"<<endl;
  ion_series.add_linked_ions(linked_peptide);
  cout <<"Printing ions"<<endl;
  ion_series.print(); 

  if (print_spectrum) {
    carp(CARP_INFO,"writing theoretical spectrum");
    map<int, FLOAT_T> theoretical;

    XHHC_Scorer::xlinkCreateMapTheoretical(ion_series,
					 theoretical);


    ofstream fout("theoretical.out");
    fout <<"> "<< linked_peptide<<"\t"<<linked_peptide.get_mz(MONO)<<endl;
    
    map<int, FLOAT_T>::iterator iter;

    for (iter = theoretical.begin();
	 iter != theoretical.end();
	 ++iter) {
      fout << iter -> first << "\t";
      fout << iter -> second << "\t";
      fout << "nolabel" <<endl;// << "\t";
      //fout << "red" << endl;
    }

    fout.close();
  }



  return 0;
}
