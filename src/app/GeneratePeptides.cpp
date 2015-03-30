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
#include "model/ModifiedPeptidesIterator.h"

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
  initialize(argc, argv);

  // Get arguments and options
  bool output_sequence = get_boolean_parameter("output-sequence");
  string filename = get_string_parameter("protein fasta file");
  Database* database = new Database(filename.c_str(), false); // not memmapped

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
      new ModifiedPeptidesIterator(peptide_mods[mod_idx], database);

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

  delete database;

  return(0);

}

void GeneratePeptides::printHeader(){
  bool bool_val;
  int missed_cleavages;

  string database_name = get_string_parameter("protein fasta file");
  printf("# PROTEIN DATABASE: %s\n", database_name.c_str());

  printf("# OPTIONS:\n");
  printf("#\tmin-mass: %.2f\n", get_double_parameter("min-mass"));
  printf("#\tmax-mass: %.2f\n", get_double_parameter("max-mass"));
  printf("#\tmin-length: %d\n", get_int_parameter("min-length"));
  printf("#\tmax-length: %d\n", get_int_parameter("max-length"));
  printf("#\tenzyme: %s\n", get_string_parameter("enzyme").c_str());
  printf("#\tdigestion: %s\n", get_string_parameter("digestion").c_str());
  missed_cleavages = get_int_parameter("missed-cleavages");
  printf("#\tnumber of allowed missed-cleavages: %d\n", missed_cleavages);
  printf("#\tisotopic mass type: %s\n", 
         get_string_parameter("isotopic-mass").c_str());
  printf("#\tverbosity: %d\n", get_verbosity_level());

  AA_MOD_T** aa_mod_list = NULL;
  int num_aa_mods = get_all_aa_mod_list(&aa_mod_list);
  int mod_idx = 0;
  for(mod_idx=0; mod_idx < num_aa_mods; mod_idx++){
    printf("#\tmodification: ");
    char* mod_str = aa_mod_to_string(aa_mod_list[mod_idx]);
    printf("%s\n", mod_str);
    free(mod_str);
  }
}


/**
 * \returns The command name for GeneratePeptides.
 */
string GeneratePeptides::getName() const {
  return "generate-peptides";
}

/**
 * \returns The description for GeneratePeptides.
 */
string GeneratePeptides::getDescription() const {
  return
    "[[nohtml:Extract from a given set of protein sequences a list of target "
    "and decoy peptides fitting the specified criteria.]]"
    "[[html:<p>This command takes as input a protein FASTA file and outputs "
    "the corresponding list of peptides, as well as a matched list of decoy "
    "peptides and decoy proteins. Decoys are generated either by reversing or "
    "shuffling the non-terminal amino acids of each peptide. The program will "
    "shuffle each peptide multiple times to attempt to ensure that there is no "
    "overlap between the target and decoy peptides. For homopolymers, this is "
    "not possible. In this case, the occurrence of these target/decoy overlaps "
    "is recorded in the log file.</p><p>The program considers only the "
    "standard set of 20 amino acids. Peptides containing non-amino acid "
    "alphanumeric characters (BJOUXZ) are skipped. Non-alphanumeric characters "
    "are ignored completely.</p>]]";
}

/**
 * \returns The command arguments
 */
vector<string> GeneratePeptides::getArgs() const {
  string arr[] = {
    "protein fasta file"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns The command options
 */
vector<string> GeneratePeptides::getOptions() const {
  string arr[] = {
    "min-length",
    "max-length",
    "min-mass",
    "max-mass",
    "isotopic-mass",
    "decoys",
    "enzyme",
    "custom-enzyme",
    "digestion",
    "missed-cleavages",
    "unique-peptides",
    "output-sequence",
    "verbosity",
    "parameter-file" 
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns The command outputs
 */
map<string, string> GeneratePeptides::getOutputs() const {
  map<string, string> outputs;
  outputs["peptides.target.txt"] =
    "a text file containing the target peptides, one per line.";
  outputs["peptides.decoy.txt"] =
    "a text file containing the decoy peptides, one per line. There is a "
    "one-to-one correspondence between targets and decoys.";
  outputs["proteins.decoy.txt"] =
    "a FASTA format file containing decoy proteins, in which all of the "
    "peptides have been replaced with their shuffled or reversed counterparts. "
    "Note that this file will only be created if the enzyme specificity is "
    "\"full-digest\" and no missed cleavages are allowed.";
  outputs["generate-peptides.params.txt"] =
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs.";
  outputs["generate-peptides.log.txt"] =
    "a log file containing a copy of all messages that were printed to the "
    "screen during execution.";
  return outputs;
}

/**
 * \returns The file stem of the application.
 */
string GeneratePeptides::getFileStem() const {
  return "peptides";
}

COMMAND_T GeneratePeptides::getCommand() const {
  return GENERATE_PEPTIDES_COMMAND;
}

bool GeneratePeptides::needsOutputDirectory() const {
  return false;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

