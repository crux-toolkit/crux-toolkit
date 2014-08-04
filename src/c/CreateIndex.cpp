/**
 * \file CreateIndex.cpp 
 * ORIGINAL AUTHOR: Chris Park
 * CRUX APPLICATION CONVERSION: Sean McIlwain
 * Missed-cleavage Conversion: Kha Nguyen
 * \brief Given a protein fasta sequence database as input, generate
 * crux_index files that contain list of peptides in the database that
 * meet certain criteria (e.g. mass, length, trypticity) as output.
 ****************************************************************************/

#include "CreateIndex.h"

#include "WinCrux.h"

using namespace std;

/**
 * \returns a blank CreateIndex object
 */
CreateIndex::CreateIndex() {

}

/**
 * Destructor
 */
CreateIndex::~CreateIndex() {
}

/**
 * main method for CreateIndex
 */
int CreateIndex::main(int argc, char** argv) {

  /* Declarations */
  int min_length;
  int max_length;
  double min_mass;
  double max_mass;
  MASS_TYPE_T mass_type;
  ENZYME_T enzyme;
  DIGEST_T digest;
  int missed_cleavages; 

  double mass_range;
  PeptideConstraint* constraint;
  char* in_file = NULL;
  Index* crux_index;
  char* binary_fasta_file = NULL;

  /* Define optional command line arguments */ 
  const char* option_list[] = { 
    "verbosity",
    "parameter-file", 
    "overwrite",
    "min-length", 
    "max-length", 
    "min-mass", 
    "max-mass", 
    "isotopic-mass",
    "enzyme", 
    "custom-enzyme", 
    "digestion", 
    "missed-cleavages",
    "peptide-list",
    "decoys"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */ 
  const char* argument_list[] = { "protein fasta file", 
                                  "index name"}; 
  int num_arguments = sizeof(argument_list) / sizeof(char*);


  /* For debugging of parameter processing */
  set_verbosity_level(CARP_ERROR);  
  carp(CARP_DETAILED_DEBUG, "Starting create_index");

  /* initialize the application */
  initialize(argument_list, num_arguments,
    option_list, num_options, argc, argv);

  /* Get parameter values */
  min_mass = get_double_parameter("min-mass");
  max_mass = get_double_parameter("max-mass");
  mass_range = (max_mass - min_mass)/MAX_INDEX_FILES;

  min_length = get_int_parameter("min-length");
  max_length = get_int_parameter("max-length");

  missed_cleavages = get_int_parameter("missed-cleavages");
  enzyme = get_enzyme_type_parameter("enzyme");
  digest = get_digest_type_parameter("digestion");
  mass_type = get_mass_type_parameter("isotopic-mass");
  DECOY_TYPE_T decoys = get_decoy_type_parameter("decoys");

  /* create peptide constraint */
  constraint = new PeptideConstraint(enzyme, digest, min_mass, max_mass, 
                                      min_length, max_length, 
                                      missed_cleavages, mass_type);
  
  /* check if input file exist */
  in_file = get_string_parameter("protein fasta file");
  if(access(in_file, F_OK)){
    carp(CARP_FATAL, "The file \"%s\" does not exist " 
         "(or is not readable or is empty).", in_file);
  }
  carp(CARP_INFO,"Creating index from fasta file '%s'", in_file);
  
  /* check if output name already exists
     fail if --overwrite is false */
  char* out_dir = get_string_parameter("index name");
  carp(CARP_DEBUG, "New index name is '%s'", out_dir);
  bool overwrite = get_boolean_parameter("overwrite");
  if( (!overwrite) && (chdir(out_dir) == 0)){
      carp(CARP_FATAL, "Index '%s' already exists. Use " \
           "--overwrite T to replace.", out_dir);
  }

  /* create new index object */
  crux_index = new Index(in_file,
                         out_dir,
                         constraint,
                         mass_range,
                         decoys
                         );
  
  /* create crux_index files */
  if(!crux_index->create(get_boolean_parameter("peptide-list"))){
    carp(CARP_FATAL, "Failed to create index");
  }
  
  /* free index(frees constraint together) */
  Index::free(crux_index);
     
  free(binary_fasta_file);
  free(out_dir);
  free(in_file);
  free_parameters();

  return 0;

}

/**
 * \returns the command name for CreateIndex
 */
string CreateIndex::getName() {
  return "create-index";
}

/**
 * \returns the description for CreateIndex
 */
string CreateIndex::getDescription() {
  return "Create an index for all peptides in a fasta file.";

}

/**
 * \returns the file stem of the application, default getName.
 */
string CreateIndex::getFileStem() {
  return "index";
}

/**
 * \returns the enum of the application, default MISC_COMMAND
 */
COMMAND_T CreateIndex::getCommand() {
  return INDEX_COMMAND;
}

bool CreateIndex::hidden() {
  return true;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
