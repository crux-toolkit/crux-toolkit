/**
 * \file GeneratePeptides.cpp
 *
 * AUTHOR: Barbara Frewen
 * CREATE DATE: January 20, 2012
 *
 * DESCRIPTION: Main method for the generate-peptides command.
 *              Output all peptide sequences in the given fasta file
 *              that fall within all peptide constraints.
 */

#include "GeneratePeptides.h"
#include "ModifiedPeptidesIterator.h"

using namespace std;
using namespace Crux; 

/**
 * \returns A blank GeneratePeptides object.
 */
GeneratePeptides::GeneratePeptides() {

}

/**
 * Destructor
 */
GeneratePeptides::~GeneratePeptides() {
}

/**
 * Main method for GeneratePeptides.
 */
int GeneratePeptides::main(int argc, char** argv) {

  // Define optional command line arguments
  const char* option_list[] = { 
    "min-length",
    "max-length",
    "min-mass",
    "max-mass",
    "isotopic-mass",
    "enzyme",
    "custom-enzyme",
    "digestion",
    "missed-cleavages",
    "unique-peptides",
    "output-sequence",
    "verbosity",
    "parameter-file" 
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  // Define required command line arguments
  const char* argument_list[] = { "protein-database" }; 
  int num_arguments = sizeof(argument_list) / sizeof(char*);
  
  initialize(argument_list, 
             num_arguments,
             option_list, 
             num_options,
             argc, argv);

  // Get arguments and options
  bool output_sequence = get_boolean_parameter("output-sequence");
  string filename = get_string_parameter("protein-database");
  bool use_index = is_directory(filename.c_str());
  Index* index = NULL;
  Database* database = NULL;

  if( use_index == true ){
    index = new Index(filename.c_str()); 
  }else{
    database = new Database(filename.c_str(), false); // not memmapped
  }

  // get list of mods
  PEPTIDE_MOD_T** peptide_mods = NULL;
  int num_peptide_mods = generate_peptide_mod_list( &peptide_mods );
  carp(CARP_DEBUG, "Got %d peptide mods", num_peptide_mods);

  printHeader();

  // one iterator for each peptide mod
  long total_peptides = 0;
  int mod_idx = 0;
  for(mod_idx = 0; mod_idx < num_peptide_mods; mod_idx++){
    carp(CARP_DETAILED_DEBUG, "Using peptide mod %d with %d aa mods", 
         mod_idx, peptide_mod_get_num_aa_mods(peptide_mods[mod_idx]));
    // create peptide iterator
    ModifiedPeptidesIterator* peptide_iterator = 
      new ModifiedPeptidesIterator(peptide_mods[mod_idx],
                                   index,
                                   database);

    // iterate over all peptides
    int orders_of_magnitude = 1000; // for counting when to print
    while(peptide_iterator->hasNext()){
      ++total_peptides;
      Peptide* peptide = peptide_iterator->next();
      peptide->printInFormat(output_sequence, 
                             stdout);
    
      // free peptide
      delete peptide;
    
      if(total_peptides % orders_of_magnitude == 0){
        if( (total_peptides)/10 == orders_of_magnitude){
          orders_of_magnitude *= 10;
        }
        carp(CARP_INFO, "Reached peptide %d", total_peptides);
      }
    }// last peptide
    delete peptide_iterator;

  }// last peptide modification

  // debug purpose
  carp(CARP_INFO, "total peptides: %d", total_peptides);

  delete index;
  delete database;

  return(0);

}

void GeneratePeptides::printHeader(){
  bool bool_val;
  int missed_cleavages;

  char* database_name = get_string_parameter("protein-database");
  printf("# PROTEIN DATABASE: %s\n", database_name);

  printf("# OPTIONS:\n");
  printf("#\tmin-mass: %.2f\n", get_double_parameter("min-mass"));
  printf("#\tmax-mass: %.2f\n", get_double_parameter("max-mass"));
  printf("#\tmin-length: %d\n", get_int_parameter("min-length"));
  printf("#\tmax-length: %d\n", get_int_parameter("max-length"));
  printf("#\tenzyme: %s\n", get_string_parameter_pointer("enzyme"));
  printf("#\tdigestion: %s\n", get_string_parameter_pointer("digestion"));
  missed_cleavages = get_int_parameter("missed-cleavages");
  printf("#\tnumber of allowed missed-cleavages: %d\n", missed_cleavages);
  printf("#\tisotopic mass type: %s\n", 
         get_string_parameter_pointer("isotopic-mass"));
  printf("#\tverbosity: %d\n", get_verbosity_level());

  bool_val = is_directory(database_name);
  printf("#\tuse index: %s\n", boolean_to_string(bool_val));
  
  AA_MOD_T** aa_mod_list = NULL;
  int num_aa_mods = get_all_aa_mod_list(&aa_mod_list);
  int mod_idx = 0;
  for(mod_idx=0; mod_idx < num_aa_mods; mod_idx++){
    printf("#\tmodification: ");
    char* mod_str = aa_mod_to_string(aa_mod_list[mod_idx]);
    printf("%s\n", mod_str);
    free(mod_str);
  }
  free(database_name);
}


/**
 * \returns The command name for GeneratePeptides.
 */
string GeneratePeptides::getName() {
  return "generate-peptides";
}

/**
 * \returns The description for GeneratePeptides.
 */
string GeneratePeptides::getDescription() {
  return "Extract from a given set of protein sequences a list of target and "
         "decoy peptides fitting the specified criteria.";
}

/**
 * \returns The file stem of the application.
 */
string GeneratePeptides::getFileStem() {
  return "peptides";
}

COMMAND_T GeneratePeptides::getCommand() {
  return GENERATE_PEPTIDES_COMMAND;
}

bool GeneratePeptides::needsOutputDirectory() {
  return false;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

