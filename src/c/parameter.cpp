/***********************************************************************//**
 * \file parameter.cpp
 * FILE: parameter.cpp
 * AUTHOR: written by Tobias Mann, CRUXified by Chris Park
 * Missed-cleavage conversion: Kha Nguyen
 * \brief General parameter handling utilities. MUST declare ALL
 * optional command parameters here inside initalialize_parameters.
 ****************************************************************************/

#include "crux-utils.h"
#include "LineFileReader.h"
#include "parameter.h"
#include "WinCrux.h"
#include "Peptide.h"
#include <iostream>

using namespace std;

//TODO:  in all set, change result=add_... to result= result && add_...

/**
 * Starting location for zeroth m/z bin.
 */
static const FLOAT_T SMART_MZ_OFFSET = 0.40;

/*
 * Global variables
 */

static const char* parameter_type_strings[NUMBER_PARAMETER_TYPES] = { 
  "INT_ARG", "DOUBLE_ARG", "STRING_ARG", "MASS_TYPE_T", "DIGEST_T", 
  "ENZYME_T", 
  "bool", "SCORER_TYPE_T", "ION_TYPE_T",
  "ALGORITHM_T", "HARDKLOR_ALGORITHM_TYPE_T", "SPECTRUM_PARSER_T" ,
  "WINDOW_TYPE_T", "MEASURE_TYPE_T", "THRESHOLD_T", 
  "PARSIMONY_TYPE_T", "QUANT_LEVEL_TYPE_T", "DECOY_TYPE_T", "MASS_FORMAT_T"};

//one hash for parameter values, one for usage statements, one for types
// all hashes keyed on parameter/option name
HASH_T* parameters; // values of parameters
HASH_T* usages;     // usage statments
HASH_T* types;      // PARAMETER_TYPE_T
HASH_T* file_notes; // additional notes for param file
HASH_T* for_users;  // false to hide from param file (args and ops for
                    // dev/research only)
HASH_T* min_values; // for numeric parameters
HASH_T* max_values; // for numeric parameters

AA_MOD_T* list_of_mods[MAX_AA_MODS]; // list containing all aa mods
                                    // in param file, c-,n-term mods at end
AA_MOD_T** list_of_variable_mods = NULL; // pointer to first non-fixed mod
AA_MOD_T** list_of_c_mods = NULL; // pointer to first c_term mod in list
AA_MOD_T** list_of_n_mods = NULL; //pointer to first n_term mod in list
int fixed_c_mod = -1; // position in list_of_mods or -1 if not present
int fixed_n_mod = -1; // position in list_of_mods or -1 if not present
int num_fixed_mods = 0;
int num_mods = 0;     // ANY_POSITION mods
int num_c_mods = 0; // variable c-term mods
int num_n_mods = 0; // variable n-term mods
//require num_mods + num_c_mods + num_n_mods + 
//(fixed_c_mod > -1) + (fixed_n_mod > -1) <= MAX_AA_MODS

bool parameter_initialized = false; //have param values been initialized
bool usage_initialized = false; // have the usages been initialized?
bool type_initialized = false; // have the types been initialized?

bool parameter_plasticity = true; // can the parameters be changed?

char* pre_cleavage_list;
char* post_cleavage_list;
int pre_list_size;
int post_list_size;
bool pre_for_inclusion;
bool post_for_inclusion;

vector<string> comet_enzyme_info_lines_;


/**
 * \returns the comet enzyme info lines parsed from the file
 * or generated defaults
 */
const std::vector<std::string>& get_comet_enzyme_info_lines() {
  return comet_enzyme_info_lines_;
}

/************************************
 * Private function declarations
 ************************************ 
 */

// parse the parameter file given the filename
// called by parse command line
void parse_parameter_file(
  char* parameter_filename ///< the parameter file to be parsed -in
  );

/**
 * Examine each option in list to determine if the values
 * are within the proper range and of the correct type
 * Requires (or at least only makes sense after) 
 * parse_cmd_line_into_params_hash() has been run.
 */
bool check_option_type_and_bounds(const char* name);

void check_parameter_consistency();
void parse_custom_enzyme(const char* rule_str);

/**
 *
 */
bool string_to_param_type(const char*, PARAMETER_TYPE_T* );

bool set_boolean_parameter(
 const char* name,       ///< the name of the parameter looking for -in
 bool   set_value,  ///< the value to be set -in
 const char* usage,      ///< message for the usage statement
 const char* filenotes,  ///< additional information for the params file
 const char* foruser     ///< "true" if should be revealed to user
 );

bool set_int_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 int set_value,  ///< the value to be set -in
 int min_value,  ///< the value to be set -in
 int max_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,  ///< additional info for param file
 const char* foruser     ///< "true" if should be revealed to user
 );

bool set_double_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 double set_value,  ///< the value to be set -in
 double min_value,  ///< the value to be set -in
 double max_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,  ///< additional info for param file
 const char* foruser     ///< "true" if should be revealed to user
  );

bool set_string_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 const char* set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser
  );

bool set_mass_type_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 MASS_TYPE_T set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser
  );

bool set_digest_type_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 DIGEST_T set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser
  );

bool set_enzyme_type_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 ENZYME_T set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser
  );

bool set_digest_type_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 DIGEST_T set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser
  );

bool set_enzyme_type_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 ENZYME_T set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser
  );

bool set_window_type_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 WINDOW_TYPE_T set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser
  );

bool set_threshold_type_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 THRESHOLD_T set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser
  );

bool set_algorithm_type_parameter(
 const char* name,
 ALGORITHM_TYPE_T set_value,
 const char* usage,
 const char* filenotes,
 const char* foruser);

bool set_hardklor_algorithm_type_parameter(
  const char* name,
  HARDKLOR_ALGORITHM_T set_value,
  const char* usage,
  const char* filenotes,
  const char* foruser);

bool set_spectrum_parser_parameter(
  const char* name, ///< the name of the parameter looking for -in
  SPECTRUM_PARSER_T set_value, ///< the value to be set -in
  const char* usage, ///< string to print in usage statement
  const char* filenotes, ///< additional info for param file
  const char* foruser);

bool set_scorer_type_parameter(
 const char* name,
 SCORER_TYPE_T set_value,
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser);

bool set_ion_type_parameter(
 const char* name,
 ION_TYPE_T set_value,
 const char* usage,
 const char* filenotes,
 const char* foruser);

bool set_parsimony_type_parameter(
  const char* name, ///< the name of the parameter looking for -in
  PARSIMONY_TYPE_T set_value, ///< the value to be set -in
  const char* usage, ///< string to print in usage statement
  const char* filenotes, ///<additional infor for param file
  const char* foruser); 

bool set_quant_level_parameter(
  const char* name, ///< the name of the parameter looking for -in
  QUANT_LEVEL_TYPE_T set_value, ///< the value to be set -in
  const char* usage, ///< string to print in usage statement
  const char* filenotes, ///<additional infor for param file
  const char* foruser); 

bool set_measure_type_parameter(
  const char* name, ///< the name of the parameter looking for -in
  MEASURE_TYPE_T set_value, ///< the value to be set -in
  const char* usage, ///< string to print in usage statement
  const char* filenotes, ///<additional infor for param file
  const char* foruser);
 
bool set_decoy_type_parameter(
  const char* name, ///< the name of the parameter looking for -in
  DECOY_TYPE_T set_value, ///< the value to be set -in
  const char* usage, ///< string to print in usage statement
  const char* filenotes, ///<additional infor for param file
  const char* foruser);

bool set_mass_format_type_parameter(
  const char* name, ///< the name of the parameter looking for -in
  MASS_FORMAT_T set_value, ///< the value to be set -in
  const char* usage, ///< string to print in usage statement
  const char* filenotes, ///<additional infor for param file
  const char* foruser);
 
bool select_cmd_line(  
  const char** option_names, ///< list of options to be allowed for main -in
  int    num_options,  ///< number of optons in that list -in
  int (*parse_argument_set)(const char*, const char*, void*, enum argument_type, bool print) ///< function point to choose arguments or options 
  );

bool update_aa_masses();
void read_mods_from_file(char* param_file);

/************************************
 * Function definitions
 ************************************
 */

/**
 * initialize parameters
 * ONLY add optional parameters here!!!
 * MUST declare ALL optional parameters in array to be used
 * Every option and its default value for every executable 
 * must be declared here
 */
void initialize_parameters(void){
  carp(CARP_DETAILED_DEBUG, "Initializing parameters in parameter.c");

  // check if parameters been initialized
  if(parameter_initialized){
    carp(CARP_ERROR, "parameters have already been initialized");
    return;
  }
  
  /* allocate the hash tables */
  parameters = new_hash(NUM_PARAMS);
  usages = new_hash(NUM_PARAMS);
  file_notes = new_hash(NUM_PARAMS);
  for_users = new_hash(NUM_PARAMS);
  types = new_hash(NUM_PARAMS);
  min_values = new_hash(NUM_PARAMS);
  max_values = new_hash(NUM_PARAMS);

  /* set number of parameters to zero */

  /* initialize the list of mods */                           
  int mod_idx = 0;                                            
  for(mod_idx = 0; mod_idx < MAX_AA_MODS; mod_idx++){         
    //initialize_aa_mod(&list_of_mods[mod_idx], mod_idx);     
    list_of_mods[mod_idx] = new_aa_mod(mod_idx);              
  }                                                           

  /* initialize custom enzyme variables */
  pre_list_size = 0;
  post_list_size = 0;
  pre_cleavage_list = NULL;
  post_cleavage_list = NULL;
  pre_for_inclusion = true;
  post_for_inclusion = false;

  /* *** Initialize Arguments *** */

  // set with name, default value, [max, min], usage, notes, for param file
  // all arguments are left out of param file

  /* generate_peptide arguments */
  set_string_parameter("protein-database", NULL, 
      "Fasta file of proteins or directory containing an index.",
      "Argument for generate, index, search. Optional for analyze and spectral-counts.", 
      "false");

  set_string_parameter("search results directory", NULL, 
      "Directory containing the results of one search.",
      "Argument for q-ranker, percolator, compute-q-values.", "false");

  /* create_index arguments */
  set_string_parameter("protein fasta file", NULL,
                       "File containing protein sequences in fasta format.",
                       "Argument for crux-create-index, tide-index, and "
                       "generate-decoys.", "false");
  set_string_parameter("index name", NULL,
    "Name to give the new directory containing index files.",
    "Argument for create index.", "false");

  /* search-for-matches arguments */
  set_string_parameter("ms2 file", NULL,
                       "File containing spectra to be searched.",
    "Argument, not option, for create-psm-files, get-ms2-spec, and search",
    "false");

  /* get-ms2-spectrum */
  set_int_parameter("scan number", 0, 1, BILLION, 
                    "Scan number identifying the spectrum.",
                    "Argument for get-ms2-spectrum", "false");
  set_string_parameter("output file", NULL, 
                       "File where spectrum will be written.",
                       "Argument for get-ms2-spectrum.", "false");

  /* predict-peptide-ions */
  set_string_parameter("peptide sequence", NULL, 
      "The sequence of the peptide.",
      "Argument for predict-peptide-ions.", "false");
  set_int_parameter("charge state", 0, 0, 10, 
      "The charge state of the peptide.",
      "Argument for predict-peptide-ions", "false");

  /* hardklor arguments */
  set_string_parameter("spectra", NULL,
                       "The name of a file from which to parse "
                       "high-resolution spectra. The file may be "
                       "in MS1 (.ms1), binary MS1 (.bms1), compressed MS1 (.cms1), or "
                       "mzXML (.mzXML) format.",
                       "Argument, not option, for hardklor",
                       "false");

  /*Percolator arguments*/
  set_string_parameter(
    "pin", NULL,
    "PIN files are tab-delimited files for PIN format. "
    "Also, this argument can be \"-\" which indicates the pin file will come from standard input. "
    "Alternately, a SQT, PepXML, or tab-delimited file may be given (a corresponding decoy"
    "file must also exist in the same directory), in which case a pin file will be "
    "generated in the output directory prior to execution.",
    "Argument, not option for percolator",
    "false"
  );
  /*make-pin arguments*/
  set_string_parameter(
    "target input", NULL,
    "search results file in sqt, tab-delimited or pep.xml format.  "
    "Also, this argument can be - which indicates the result file will come from standard input",
    "Argument for make-pin and calibrate-scores",
    "false"
  );
  set_string_parameter(
    "decoy input", NULL,
    "make-pin can convert any file format in sqt, tab-delimited and pep.xml file "
    "to pin file ",
    "Argument, not option for make-pin",
    "false"
  );
  set_string_parameter(
    "output-file", NULL,
    "Path where pin file will be written",
    "It is optional for make-pin",
    "false"
  );
  set_boolean_parameter(
    "filestem-prefixes", false,
    "Prefix PSM IDs with filestems instead of target or decoy and file index.",
    "Available for make-pin",
    "false");
  /* *** Initialize Options (command line and param file) *** */

  /* options for all executables */
  set_boolean_parameter("version", false, "Print version number and quit.",
      "Available for all crux programs.  On command line use '--version T'.",
      "true");
  set_int_parameter("verbosity", CARP_INFO, CARP_FATAL, CARP_MAX,
      "Set level of output to stderr (0-100).  Default=30.",
      "Available for all crux programs.  Each level prints the following "
      "messages, including all those at lower verbosity levels: 0-fatal "
      "errors, 10-non-fatal errors, 20-warnings, 30-information on the "
      "progress of execution, 40-more progress information, 50-debug info, "
      "60-detailed debug info.", "true");
  set_string_parameter("parameter-file", NULL, 
      "Set additional options with values in the given file.",
      "Available for all crux programs. Any options specified on the "
      "command line will override values in the parameter file.", "true");
  set_boolean_parameter("overwrite", false, 
      "Replace existing files (T) or exit if attempting to "
      "overwrite (F). Default=F.",
      "Available for all crux programs.  Applies to parameter file "
      "as well as index, search, and analysis output files.", "true");
    
  /* generate_peptide, create_index parameters  */
  set_int_parameter("min-length", 6, 1, MAX_PEPTIDE_LENGTH,
      "The minimum length of peptides to consider. Default=6.",
      "Used from the command line or parameter file by "
      "crux-create-index, crux-generate-peptides, crux tide-index, and crux "
      "generate-decoys. Parameter file only for crux-search-for-matches.",
      "true");
  set_int_parameter("max-length", 50, 1, MAX_PEPTIDE_LENGTH,
      "The maximum length of peptides to consider. Default=50.",
      "Available from command line or parameter file for crux-create-index, "
      "crux-generate-peptides, crux tide-index, and crux generate-decoys. "
      "Parameter file only for crux-search-for-matches.", "true");
  set_double_parameter("min-mass", 200, 0, BILLION,
      "The minimum mass of peptides to consider. Default=200.",
      "Available from command line or parameter file for crux-create-index, "
      "crux-generate-peptides, crux tide-index, and crux generate-decoys. "
      "Parameter file only for crux-search-for-matches.", "true");
  set_double_parameter("max-mass", 7200, 1, BILLION, 
      "The maximum mass of peptides to consider. Default=7200.",
      "Available from command line or parameter file for crux-create-index, "
      "crux-generate-peptides, crux tide-index, and crux generate-decoys. "
      "Parameter file only for crux-search-for-matches.", "true");
  set_mass_type_parameter("isotopic-mass", AVERAGE, 
      "Which isotopes to use in calcuating peptide mass. "
      "<string>=average|mono. Default=average.", 
      "Used from command line or parameter file by crux-create-index, "
      "crux-generate-peptides, and crux generate-decoys.  Parameter file only "
      "for crux-search-for-matches.", "true");
  set_int_parameter("min-peaks", 20, 0, BILLION,
      "The minimum number of peaks a spectrum must have for it to be searched."
      " Default=20.", 
      "Available from command line for tide-search or parameter file for "
      "search-for-matches and sequest-search.", "true");
  set_digest_type_parameter("digestion", FULL_DIGEST,
      "Degree of digestion used to generate peptides. "
      "<string>=full-digest|partial-digest. Either both ends or one end "
      "of a peptide must conform to enzyme specificity rules. "
      "Default=full-digest.",
      "Used in conjunction with enzyme option when enzyme is not set to "
      "to 'no-enzyme'.  Available from command line or parameter file for "
      "crux-generate-peptides, crux create-index, crux tide-index, and crux "
      "generate-decoys. Available from parameter file for crux "
      "search-for-matches. Digestion rules are as "
      "follows: enzyme name [cuts after one of these residues][but not before "
      "one of these residues].  trypsin [RK][P], elastase [ALIV][P], "
      "chymotrypsin [FWYL][P].",
      "true");
  set_enzyme_type_parameter("enzyme", TRYPSIN,
      "Enzyme to use for in silico digestion of proteins. "
      "<string>=no-enzyme|trypsin|trypsin/p|chymotrypsin| " 
      "elastase|clostripain|cyanogen-bromide|iodosobenzoate| " 
      "proline-endopeptidase|staph-protease|asp-n|lys-c "
      "lys-n|arg-c|glu-c|pepsin-a| "
      "|elastase-trypsin-chymotrypsin|custom-enzyme. "
      "Default=trypsin.", 
      "Used in conjunction with the options digestion and missed-cleavages. "
      "Use 'no-enzyme' for non-specific digestion.  Available "
      "from command line or parameter file for crux-generate-peptides, "
      "crux create-index, crux tide-index, and crux generate-decoys.  "
      "Available from parameter file for crux search-for-matches. "
      "Digestion rules: enzyme name [cuts after one of these residues]|{but "
      "not before one of these residues}. trypsin [RK]|{P}, trypsin/p [RK]|[], "
      "elastase [ALIV]|{P}, chymotrypsin [FWYL]|{P}, clostripain [R]|[], "
      "cyanogen-bromide [M]|[], iodosobenzoate [W]|[], proline-endopeptidase "
      "[P]|[], staph-protease [E]|[], elastase-trypsin-chymotrypsin "
      "[ALIVKRWFY]|{P},asp-n []|[D] (cuts before D), lys-c [K]|{P}, lys-n "
      "[]|[K] (cuts before K), arg-c [R]|{P}, glu-c [DE]|{P}, pepsin-a "
      "[FL]|{P}.", "true");

  set_window_type_parameter("precursor-window-type", WINDOW_MASS,
      "Window type to use for selecting candidate "
      "peptides.  <string>=mass|mz|ppm. Default=mass.",
      "Available for search-for-matches, search-for-xlinks, "
      "and tide-search.",
      "true");

  set_spectrum_parser_parameter("spectrum-parser", PROTEOWIZARD_SPECTRUM_PARSER,
    "Parser to use for reading in spectra "
    "<string>=pwiz|mstoolkit. Default=pwiz.",
    "Available for search-for-matches, search-for-xlinks.",
    "true");

  set_string_parameter("custom-enzyme", NULL, 
      "Specify rules for in silico digestion of proteins. "
      "See HTML documentation for syntax. Default is trypsin.",
      "Overrides the enzyme option.  Two lists of residues are given enclosed "
      "in square brackets or curly braces and separated by a |. The first list "
      "contains residues required/prohibited before the cleavage site and the "
      "second list is residues after the cleavage site.  If the residues are "
      "required for digestion, they are in square brackets, '[' and ']'.  "
      "If the residues prevent digestion, then they are enclosed in curly "
      "braces, '{' and '}'.  Use X to indicate all residues.  For example, "
      "trypsin cuts after R or K but not before P which is represented as "
      "[RK]|{P}.  AspN cuts after any residue but only before D which is "
      "represented as [X]|[D].",
                       "true");
  
  set_int_parameter("missed-cleavages",
                    0, 0, 500,
      "Include peptides with up to n missed cleavage sites. Default=0.",
      "Available from command line or parameter file for crux-create-index, "
      "crux-generate-peptides, and crux generate-decoys.  Parameter file only "
      "for crux-search-for-matches.  When used with enzyme=<trypsin|elastase|"
      "chymotrpysin> includes peptides containing one or more potential "
      "cleavage sites.", "true");   
  set_string_parameter("keep-terminal-aminos", "NC", 
      "When creating decoy peptides using decoy-format=shuffle or decoy-format="
      "peptide-reverse, this option specifies whether the N-terminal and "
      "C-terminal amino acids are kept in place or allowed to be shuffled or "
      "reversed. Default = NC.",
      "Available for tide-index and generate-decoys.", "true");

  set_boolean_parameter("unique-peptides", true,
      "Generate peptides only once, even if they appear in more "
      "than one protein (T,F).  Default=T.",
      "Available from command line or parameter file for "
      "crux-genereate-peptides. Returns one line per peptide "
      "when true or one line per peptide per protein occurence when false.  ",
      "true");
  set_boolean_parameter("peptide-list", false,
                        "Create an ASCII version of the peptide list.  "
                        "Default=F.",
                        "Creates an ASCII file in the output directory "
                        "containing one peptide per line.",
                        "true");
  
  /* more generate_peptide parameters */
  set_boolean_parameter("output-sequence", false, 
      "Print peptide sequence (T,F). Default=F.",
      "Available only for crux-generate-peptides.", "true");

  /* search-for-matches command line options */
  set_boolean_parameter("sqt-output", false,
      "Output SQT in the output directory.  Default=F",
      "Available for search-for-matches.", "true");
  set_boolean_parameter("mzid-output", false,
      "Output MZID in the output directory.  Default=F",
      "Available for search-for-matches, percolator.", "true");
  set_boolean_parameter("pin-output", false,
      "Output PIN XML in the output directory.  Default=F",
      "Available for search-for-matches.", "true");
  set_boolean_parameter("pout-output", false,
      "Output POUT XML in the output directory.  Default=F",
      "Available for percolator.", "true");
  set_boolean_parameter("pepxml-output", false,
      "Output pepXML in the output directory.  Default=F",
      "Available for search-for-matches, q-ranker, barista, percolator.",
      "true");
  set_boolean_parameter("txt-output", true,
      "Output tab-delimited text in the output directory.  Default=T",
      "Available for search-for-matches, percolator, q-ranker, barista.",
      "true");
  set_scorer_type_parameter("prelim-score-type", SP, 
      "Initial scoring (sp, xcorr). Default=sp,", 
      "Available for crux-search-for-matches.  The score applied to all "
      "possible psms for a given spectrum.  Typically used to filter out "
      "the most plausible for further scoring. See max-rank-preliminary and "
      "score-type.", "false");
  set_scorer_type_parameter("score-type", XCORR, 
      "The primary scoring method to use (xcorr, sp, xcorr-pvalue, sp-pvalue)."
      " Default=xcorr.", 
      "Only available for crux-search-for-matches.  Primary scoring is "
      "typically done on a subset (see max-rank-preliminary) of all "
      "possible psms for each spectrum. Default is the SEQUEST-style xcorr."
      " Crux also offers a p-value calculation for each psm based on xcorr "
      "or sp (xcorr-pvalue, sp-pvalue).", "false"); 
  set_boolean_parameter("compute-sp", false,
      "Compute the Sp score for all candidate peptides.  Default=F",
      "Available for search-for-matches and tide-search.  Sp scoring is always "
      "done for sequest-search.", "true");
  set_boolean_parameter("compute-p-values", false, 
      "Compute p-values for the main score type. Default=F.",
      "Currently only implemented for XCORR.", "true");
  set_string_parameter("scan-number", NULL,
      "Search only select spectra specified as a single "
      "scan number or as a range as in x-y.  Default=search all.",
      "The search range x-y is inclusive of x and y.", "true");
  /* N.B. Use NaN to indicate that no user preference was specified.
   * In this case, the default value depends on the mass type.
   * S.M. Also prevent a width of 0.                                */
  set_double_parameter("mz-bin-width", 1.0005079, 1e-4, BILLION,
      "Specify the width of the bins used to "
      "discretize the m/z axis.  Also used as tolerance for assigning "
      "ions.  Default=1.0005079 for monoisotopic mass "
      "or 1.0011413 for average mass.",
      "Available for crux-search-for-matches, tide-search, and xlink-assign-ions.", "true");
  set_double_parameter("mz-bin-offset", SMART_MZ_OFFSET, 0.0, 1.0,
      "Specify the location of the left edge of the "
      "first bin used to discretize the m/z axis. Default=0.40",
      "Available for crux-search-for-matches and tide-search.", "true");
  // initialize as "unset", then set as bool after cmdline parsed
  set_boolean_parameter("use-flanking-peaks", false,
      "Include peaks +/- 1da around b/y ions in theoretical spectrum.  "
      "default=F.",
      "Available for the tide-search and search-for-xlinks commands.",
      "true");
  set_double_parameter("spectrum-min-mz", 0.0, 0, BILLION, 
      "The lowest spectrum m/z to search. Default=0.0.",
      "Available for crux-search-for-matches and tide-search.", "true");
  set_double_parameter("spectrum-max-mz", BILLION, 1, BILLION, 
      "The highest spectrum m/z to search. Default=no maximum.",
      "Available for crux-search-for-matches and tide-search.", "true");
  set_string_parameter("spectrum-charge", "all", 
      "Spectrum charge states to search. "
      "<string>=1|2|3|all. Default=all.",
      "Used by tide-search and crux-search-for-matches to limit the charge "
      "states considered in the search.  With 'all' every spectrum will be "
      "searched and spectra with multiple charge states will be searched "
      "once at each charge state.  With 1, 2 ,or 3 only spectra with that "
      "that charge will be searched.", "true");
  set_string_parameter("fileroot", NULL, 
      "Prefix added to output file names. Default=none. ",
      "Used by crux search-for-matches, crux sequest-search, crux percolator "
      "crux compute-q-values, and crux q-ranker.", "true");
  set_string_parameter("output-dir", "crux-output", 
      "Folder to which results will be written. "
      "Default='crux-output'.",
      "Used by crux create-index, crux search-for-matches, "
      "crux compute-q-values, and crux percolator.", "true");
  set_string_parameter("search-decoy-pvalue-file", "search.decoy.p.txt", 
      "Output filename for complete list of decoy p-values.  Default="
      "'search.decoy.p.txt'",
      "Only available for crux search-for-matches. The location of this "
      "file is controlled by --output-dir.", "true");

  // user options regarding decoys
  set_string_parameter("decoys", "peptide-shuffle",
      "Include a decoy version of every peptide by shuffling or reversing the "
      "target sequence.  <string>=none|reverse|protein-shuffle|peptide-shuffle."
      " Use 'none' for no decoys.  Default=peptide-shuffle.",
      "For create-index, store the decoys in the index.  For search, either "
      "use decoys in the index or generate them from the fasta file.", "true");
  set_int_parameter("num-decoys-per-target", 1, 0, 10,
      "Number of decoy peptides to search for every target peptide searched."
      "Only valid for fasta searches when --decoys is not none. Default=1.",
      "Use --decoy-location to control where they are returned (which "
      "file(s)) and --decoys to control how targets are randomized.  Available "
      "for search-for-matches when searching a fasta file. ",
      "true");
  set_string_parameter("decoy-location", "separate-decoy-files",
      "Specify location of decoy search results. "
      "<string>=target-file|one-decoy-file|separate-decoy-files. "
      "Default=separate-decoy-files.",
      "Applies when decoys is not none.  Use 'target-file' to mix "
      "target and decoy search results in one file. 'one-decoy-file' will "
      "return target results in one file and all decoys in another. "
      "'separate-decoy-files' will create as many decoy files as "
      "num-decoys-per-target.",
      "true");

  // coder options regarding decoys
  set_int_parameter("num-decoy-files", 2, 0, 10,
                    "Replaces number-decoy-set.  Determined by decoy-location"
                    " and num-decoys-per-target",
                    "", "false");
  set_boolean_parameter("tdc", false,
      "Target-decoy competition. puts decoy psms in target file. ",
      "Now hidden from the user", "false");
  set_boolean_parameter("decoy-p-values", false,
                        "Store all decoy p-values in a file",
                        "", "false");
  set_int_parameter("max-rank-preliminary", 500, 0, BILLION, 
      "Number of psms per spectrum to score with xcorr after preliminary "
      "scoring with Sp. "
      "Set to 0 to score all psms with xcorr. Default=500.",
      "Used by crux-search-for-matches.  For positive values, the Sp "
      "(preliminary) score acts as a filter; only high scoring psms go "
      "on to be scored with xcorr.  This saves some time.  If set to 0, "
      "all psms are scored with both scores. ", "true");
  set_int_parameter("top-match", 5, 1, BILLION, 
      "The number of PSMs per spectrum writen to the output " 
      " file(s).  Default=5.",
      "Available for tide-search, crux percolator, and from parameter file for "
      "crux-search-for-matches.",
      "true");
  set_int_parameter("psms-per-spectrum-reported", 0, 0, BILLION,
                   "place holder", "this may be replaced by top-match","false");
  set_string_parameter("seed", "1",
      "When given a unsigned integer value seeds the random number generator with that value. "
      "When given the string \"time\" seeds the random number generator with the system time. "
      "Default = 1.",
      "Available for all Crux commands.",
      "true");
  set_double_parameter("precursor-window", 3.0, 0, 100, 
      "Search peptides within +/- 'precursor-window' "
      "of the spectrum mass.  Definition of precursor window depends "
      "upon precursor-window-type. Default=3.0.",
      "Available for tide-search and from the parameter file only for "
      "crux-search-for-matches, crux-create-index, crux-generate-peptides.",
      "true");
  set_mass_type_parameter("fragment-mass", MONO, 
      "Which isotopes to use in calculating fragment ion mass. "
      "<string>=average|mono. Default=mono.", 
      "Parameter file only.  "
      "Used by crux-search-for-matches and crux-predict-peptide-ions.",
                          "true");
  set_string_parameter("mod", "NO MODS",
      "Specify a variable modification to apply to peptides.  " 
      "<mass change>:<aa list>:<max per peptide>:<prevents cleavage>:<prevents cross-link>."
      "  Sub-parameters prevents cleavage and prevents cross-link are optional (T|F)."
      "Default=no mods.",
      "Available from parameter file for crux-generate-peptides and "
      "crux-search-for-matches and the "
      "the same must be used for crux compute-q-value.", "true");
  set_string_parameter("cmod", "NO MODS",
      "Specify a variable modification to apply to C-terminus of peptides. " 
      "<mass change>:<max distance from protein c-term (-1 for no max)>. " 
      "Default=no mods.",       
      "Available from parameter file for crux-generate-peptides and "
      "crux-search-for-matches and the "
      "the same must be used for crux compute-q-value.", "true");
  set_string_parameter("nmod", "NO MODS",
      "Specify a variable modification to apply to N-terminus of peptides.  " 
      "<mass change>:<max distance from protein n-term (-1 for no max)>",
      "Available from parameter file for crux-generate-peptides and "
      "crux-search-for-matches and the "
      "the same must be used for crux compute-q-value.", "true");
  set_string_parameter("cmod-fixed", "NO MODS",
      "Specify a fixed modification to apply to the C-terminus of peptides.",
      "Available from parameter file for crux sequest-search and "
      "search-for-matches.", "true");
  set_string_parameter("nmod-fixed", "NO MODS",
      "Specify a fixed modification to apply to the N-terminus of peptides.",
      "Available from parameter file for crux sequest-search and "
      "search-for-matches.", "true");
  set_int_parameter("min-mods", 0, 0, MAX_PEPTIDE_LENGTH,
      "The minimum number of modifications that can be applied to a single " 
      "peptide.  Default=0.",
      "Available for tide-index.", "true");
  set_int_parameter("max-mods", MAX_PEPTIDE_LENGTH, 0, MAX_PEPTIDE_LENGTH,
      "The maximum number of modifications that can be applied to a single " 
      "peptide.  Default=no limit.",
      "Available for tide-index.", "true");
  set_int_parameter("max-aas-modified", MAX_PEPTIDE_LENGTH, 0,
      MAX_PEPTIDE_LENGTH,
      "The maximum number of modified amino acids that can appear in one "
      "peptide.  Each aa can be modified multiple times.  Default=no limit.",
      "Available from parameter file for search-for-matches.", "true");
  set_mass_format_type_parameter("mod-mass-format", MOD_MASS_ONLY,
      "Print the masses of modified sequences in one of three ways 'mod-only', "
      "'total' (residue mass plus modification), or 'separate' (for multiple "
      "mods to one residue): Default 'mod-only'.",
      "Available in the parameter file for search-for-matches, sequest-search "
      "and generate-peptides.",
      "true");
  set_int_parameter("mod-precision", MOD_MASS_PRECISION, 0, 20,//arbitrary
      "Set the precision for modifications as written to .txt files.",
      "Also changes mods written to parameter file. Set internally based on "
      "the max mod precision in the param file.",
      "false");
  set_int_parameter("precision", 8, 1, 100, //max is arbitrary
      "Set the precision for scores written to sqt and text files. "
      "Default=8.",
      "Available from parameter file for crux search-for-matches, percolator, "
      "and compute-q-values.", "true");
  set_int_parameter("mass-precision", 4, 1, 100, // max is arbitrary
      "Set the precision for masses and m/z written to sqt and .txt files.  "
      "Default=4",
      "Available from parameter file for all commands.", "true");
  set_int_parameter("print-search-progress", 1000, 0, BILLION,
      "Show search progress by printing every n spectra searched.  Default="
      "1000.", "Set to 0 to show no search progress.  Available for crux "
      "search-for-matches and tide-search from parameter file.",
      "true");

  // Sp scoring params
  set_double_parameter("beta", 0.075, 0, 1, "Not for general users.",
      "Only used to set scorer->sp_beta which is used to score sp.", 
      "false"); 
  set_double_parameter("max-mz", 4000, 0, BILLION, 
      "Used in scoring sp.",
      "Hide from users", "false");
  set_double_parameter("fraction-top-scores-to-fit", 0.55, 0, 1, 
      "The fraction of psms per spectrum to use for estimating the "
      "score distribution for calculating p-values. "
      "Not compatible with 'number-top-scores-to-fig'. Default=0.55.",
      "For developers/research only.", "false");

  /* analyze-matches options */
  set_algorithm_type_parameter("algorithm", PERCOLATOR_ALGORITHM, 
      "The analysis algorithm to use (percolator, curve-fit, none)."
      " Default=percolator.",
      "Available only for crux-analyze-matches.  Using 'percolator' will "
      "assign a q-value to the top-ranking psm for each spectrum based on "
      "the decoy searches.  Using 'curve-fit' will assign a q-value to same "
      "using the p-values calculated with score-type=<xcorr-pvalue|"
      "sq-pvalue>.  Incorrect combinations of score-type and algorithm cause"
      " undefined behavior. Using 'none' will turn the binary .csm files "
      "into text.", "false");

  // **** percolator options. ****
  set_boolean_parameter("feature-file", false,
     "Optional file into which psm features are printed. Default=F.",
     "Available for percolator and q-ranker.  File will be named "
     "<fileroot>.percolator.features.txt or <fileroot>.qranker.features.txt.",
     "true");
  
  set_boolean_parameter(
    "protein",
    false,
    "output protein level probability. Default=F",
    "Available for crux percolator",
    "true"
  );
 
  set_boolean_parameter(
    "decoy-xml-output",
    false,
    "Include decoys (PSMs, peptides, and/or proteins) in the "
    "xml-output. Only available if -X is used. Default=F",
    "Available for crux percolator",
    "true"
  );
  set_string_parameter(
    "decoy-prefix",
    "decoy_",
    "Option for single SQT file mode defining the name pattern "
    "used for decoy database. Default=decoy_.",
    "Available for percolator",
    "true"
  );
 
  set_double_parameter(
    "c-pos",
    0.01,-BILLION,BILLION,
    "Penalty for mistakes made on positive examples. Set by "
    "cross validation if not specified. Default=cross-validate ",
    "Available for crux percolator",
    "true"
  );
  set_double_parameter(
    "c-neg",
    0.0,0.0,0.90,
    "Penalty for mistake made on negative examples. Set by cross "
    "validation if not specified or --c-pos not specified.",
    "Available for crux percolator",
    "true"
  );
 
  set_double_parameter(
    "test-fdr",
    0.01,0.0,1.0,
    "False discovery rate threshold for evaluating best cross validation result "
    "and the reported end result. Default is 0.01.",
    "Availble for crux percolator.",
    "true"
  );
 
  set_int_parameter(
    "maxiter",
    10,0,100000000,
    "Maximum number of iterations for training (default 10).",
    "Available for crux percolator",
    "false"
  );
  set_double_parameter(
    "train-ratio",
    0.6,0.0,1.0,
    "Fraction of the negative data set to be used as train set when only providing"
    " one negative set, remaining examples will be used as test set.Default 0.6",
    "Available for crux percolator.",
    "true"
  );
  set_boolean_parameter(
    "output-weights",
    false,
    "Output final weights to percolator.target.weights.txt.Default=T.",
    "Available for crux percolator",
    "true"
  );
  set_string_parameter(
    "input-weights",
    NULL,
    " Read initial weights from the given file (one per line). Default do not read " 
    "initial weights. Default=F",
    "Available for crux percolator ",
    "true"
  );
  set_string_parameter(
    "default-direction",
    NULL,
    "The most informative feature given as the feature name.  The name can be "
    "preceded by a hyphen (e.g., \"-XCorr\") to indicate that a lower value is better."
    " By default, Percolator will select the feature that produces"
    " the largest set of target PSMs at a specified FDR threshold"
    " (cf. --train-fdr).",
    "Available for crux percolator",
    "true"
  );
  set_boolean_parameter(
    "unitnorm",
    false,
    "Use unit normalization [0-1] instead of standard deviation normalization",
    "Available for crux percolator.",
    "true"
  );
 
  set_double_parameter(
    "train-fdr",
    0.01,0,BILLION,
    "False discovery rate threshold to define positive examples in training. "
    "Set by cross validation if 0. Default is 0.01.",
    "Available for crux percolator",
    "true"
  );

  set_double_parameter(
    "alpha",
    0.0,0.0,1.0,
    "Probability with which a present protein emits an associated peptide (--protein T "
    "must be set). Set by grid search if not specified.",
    "Available for crux percolator if --protein T is set.",
    "true"
  );
  set_double_parameter(
    "beta",
    0.0,0.0,10.0,
    "Probability of the creation of a peptide from noise (--protein T "
    "must be set). Set by grid search if not specified.",
    "Available for crux percolator if --protein T is set.",
    "true"
  );
 
  set_double_parameter(
    "gamma",
    0.0,0.0,10.0,
    "Prior probability of that a protein is present in the sample (--protein T "
    "must be set). Set by grid search if not specified.",
    "Available for crux percolator if --protein T is set.",
    "true"
  );
  set_boolean_parameter(
    "test-each-iteration",
    false,
    "Measure performance on test set each iteration",
    "Available for crux percolator.",
    "true"
  );
 
  set_boolean_parameter(
    "static-override",
    false,
    "Override error check and do not fall back on default score vector in case of suspect score vector.",
    "Available for crux percolator.",
    "true"
  );
 
  set_boolean_parameter(
    "klammer",
    false,
    "Using retention time features calculated as in Klammer et al.",
    "Available for crux percolator",
    "true"
  );


  set_int_parameter(
    "doc",
    -1,-1,15,
    "Include description of correct features.",
    "Avilable for crux percolator",
    "true"
  );

  set_boolean_parameter(
    "only-psms",
    false,
    "Do not remove redundant peptides, keep all PSMs and exclude peptide level probability.",
    "Available for crux percolator",
    "true"
  );
 
  set_boolean_parameter(
    "allow-protein-group",
    false,
    "Treat ties as if it were one protein ",
    "Available for crux percolator.",
    "true"
  );

  set_boolean_parameter(
    "protein-level-pi0",
    false,
    "Use pi_0 value when calculating empirical q-values (--protein T must be set).",
    "Available for crux percolator if --protein T is set.",
    "true"
  );

  set_boolean_parameter(
    "group-proteins",
    false,
    "Proteins with same probabilities will be grouped (--protein T must be set).",
    "Available for crux percolator if --protein T is set.",
    "true"
  );

  set_boolean_parameter(
    "empirical-protein-q",
    false,
    "Output empirical q-values from target-decoy analysis (--protein T must be set).",
    "Available for crux percolator if --protein T is set.",
    "true"
  );
 
  set_boolean_parameter(
    "no-prune-proteins",
    false,
    "Peptides with low score will not be pruned before calculating protein probabilities "
    "(--protein T must be set).",
    "Available for crux percolator if --protein T is set.",
    "true"
  );

  set_int_parameter(
    "deepness",
    0,0,2,
    "Setting deepness 0, 1, or 2 from low depth to high depth (less computational time) "
    "of the grid search for estimation Alpha,Beta and Gamma parameters for fido "
    "(--protein T must be set). Default value is 0.",
    "Available for crux percolator if --protein T is set.",
    "true"
  );

  set_boolean_parameter(
    "original-output", false,
    "Output the standalone Percolator tab-delimited output.",
    "Available for crux percolator.",
    "false"
  );

  // **** Tide arguments ****
  set_string_parameter("spectrum records file", NULL,
    "A spectrum records file generated by a previous run of crux tide-search "
    "using the store-spectra parameter.",
    "Available for read-spectrumrecords",
    "true"
  );

  set_string_parameter("tide spectra file", NULL,
    "The name of the file from which to parse the fragmentation spectra, in any "
    "of the file formats supported by ProteoWizard. Alternatively, the argument "
    "may be a binary spectrum file produced by a previous run of crux "
    "tide-search using the store-spectra parameter.",
    "Available for tide-search",
    "true"
  );

  set_string_parameter("tide database index", NULL,
    "A directory containing a database index created by a previous run of crux "
    "tide-index.",
    "Available for tide-search",
    "true"
  );

  // **** Tide options ****
  set_string_parameter("decoy-format", "shuffle",
    "Include a decoy version of every peptide by shuffling or reversing the "
    "target sequence or protein. <string>=none|shuffle|peptide-reverse|"
    "protein-reverse. Default=shuffle.",
    "Available for tide-index",
    "true"
  );
  set_boolean_parameter("monoisotopic-precursor", true,
    "Use monoisotopic precursor masses rather than average mass for precursor. "
    "Default=T.",
    "Available for tide-index",
    "true"
  );
  set_string_parameter("mods-spec", "C+57.02146",
    "Expression for static and variable mass modifications to include. "
    "Specify a comma-separated list of modification sequences of the form: "
    "C+57.02146,2M+15.9949,1STY+79.966331,...",
    "Available for tide-index",
    "true"
  );
  set_string_parameter("cterm-peptide-mods-spec", "",
    "Specifies C-terminal static and variable mass modifications on peptides. "
    "Specify a comma-separated list of C-terminal modification sequences of the form: "
    "X+21.9819",
    "Available for tide-index",
    "true"
  );
  set_string_parameter("nterm-peptide-mods-spec", "",
    "Specifies N-terminal static and variable mass modifications on peptides. "
    "Specify a comma-separated list of N-terminal modification sequences of the form: "
    "1E-18.0106,C-17.0265",
    "Available for tide-index",
    "true"
  );
  set_string_parameter("cterm-protein-mods-spec", "",
    "Specifies C-terminal static and variable mass modifications on proteins. "
    "Specify a comma-separated list of C-terminal protein modification sequences of the form: "
    ",...",
    "Available for tide-index",
    "true"
  );
  set_string_parameter("nterm-protein-mods-spec", "",
    "Specifies N-terminal static and variable mass modifications on proteins. "
    "Specify a comma-separated list of N-terminal protein modification sequences of the form: "
    ",...",
    "Available for tide-index",
    "true"
  );
  set_string_parameter("store-spectra", "",
    "Specify the name of the file where the binarized fragmentation spectra "
    "will be stored.",
    "Available for tide-search",
    "true"
  );
  set_boolean_parameter("exact-p-value", false,
    "Uses exact P-value calculation for peptide-spectrum-matching. "
    "Default=F.",
    "Available for tide-search",
    "true"
  );
  set_boolean_parameter("concat", false,
    "Output target and decoy PSMs into a single file.",
    "Available for tide-search",
    "true"
  );
  // Same as remove_precursor_peak and remove_precursor tolerance in Comet
  set_boolean_parameter("remove-precursor-peak", false,
    "Remove peaks around the precursor m/z.",
    "Available for tide-search.",
    "true"
  );
  set_double_parameter("remove-precursor-tolerance", 1.5, 0, BILLION,
    "+- m/z tolerance for precursor peak removal. Default = 1.5.",
    "Available for tide-search.",
    "true"
  );
  set_boolean_parameter("clip-nterm-methionine", false,
    "This parameter controls whether Tide will automatically "
    "remove the N-terminal methionine from a sequence entry.",
    "Available for tide-index.",
    "true"
  );
  set_boolean_parameter("use-neutral-loss-peaks", false,
    "Controls whether neutral loss ions are considered in the search. "
    "Two types of neutral losses are included and are applied only to "
    "singly charged b- and y-ions: loss of ammonia (NH3, 17.0086343 Da) "
    "and H2O (18.0091422). Each neutral loss peak has intensity 1/5 of "
    "the primary peak",
    "Available for tide-search.",
    "true"
  );
  /*
   * Comet parameters
   */
  set_string_parameter("input spectra", NULL,
     "Tha name of file (in MS2 format) from which to parse the spectra.",
     "Available for comet.",
     "false");

  set_string_parameter("database_name", NULL,
                      "A full or relative path to the sequence database, "
                      "in FASTA format, to search. Example databases include "
                      "RefSeq or UniProt.  Database can contain amino acid "
                      "sequences or nucleic acid sequences. If sequences are "
                      "amino acid sequences, set the parameter \"nucleotide_reading_frame = 0\". "
                      "If the sequences are nucleic acid sequences, you must instruct Comet to "
                      "translate these to amino acid sequences. Do this by setting "
                      "nucleotide_reading_frame\" to a value between 1 and 9. ",
                      "Comet only", "true");

  set_int_parameter("decoy_search", 0, 0, 2,
                    "0=no (default), 1=concatenated search, 2=separate search",
                    "option for Comet only", "true");

  set_int_parameter("num_threads",0,0,32, 
    "0=poll CPU to set num threads; else specify num threads directly (max 32)",
    "option for Comet only",
    "true"
  );

  set_string_parameter("output_suffix","",
          "specifies the suffix string that is appended to the base output name "
          "for the pep.xml, pin.xml, txt and sqt output files."
          "Default = \"\"",
          "Available for comet.","true");

  set_double_parameter("peptide_mass_tolerance", 3.0, 0, BILLION,
                       "Controls the mass tolerance value.  The mass tolerance "
                       "is set at +/- the specified number i.e. an entered value "
                       "of \"1.0\" applies a -1.0 to +1.0 tolerance. "
                       "The units of the mass tolerance is controlled by the parameter "
                       "\"peptide_mass_units\". ", 
                       "option for Comet only","true");

  set_int_parameter("peptide_mass_units", 0,0,2,
                    "0=amu, 1=mmu, 2=ppm",
                    "option for Comet only", "true");

  set_int_parameter("mass_type_parent", 1,0,1,
                    "0=average masses, 1=monoisotopic masses","option for Comet only", "true");
  set_int_parameter("mass_type_fragment", 1,0,1,
                    "0=average masses, 1=monoisotopic masses","option for Comet only", "true");
  
  set_int_parameter("precursor_tolerance_type", 0, 0, 1,
                    "0=MH+ (default), 1=precursor m/z","option for Comet only", "true");

  set_int_parameter("isotope_error",0,0,2, 
    "0=off, 1=on -1/0/1/2/3 (standard C13 error), 2= -8/-4/0/4/8 (for +4/+8 labeling)",
    "option for Comet only",
    "true"
  );

  set_int_parameter("search_enzyme_number", 1, 0, BILLION,
                    "choose from list at end of this params file",
                    "option for Comet only",
                    "true");

  set_int_parameter("num_enzyme_termini", 2,1,2,
                    "valid values are 1 (semi-digested), "
                    "2 (fully digested, default), 8 N-term, 9 C-term",
                    "option for Comet only",
                    "true");

  set_int_parameter("allowed_missed_cleavage", 2, 0, 5,
                    "maximum value is 5; for enzyme search",
                    "option for Comet only",
                    "true");

  set_double_parameter("fragment_bin_tol", 1.000507, 0, BILLION,
                       "binning to use on fragment ions",
                       "option for Comet only",
                       "true");

  set_double_parameter("fragment_bin_offset", SMART_MZ_OFFSET, 0, 1.0,
                       "offset position to start the binning (0.0 to 1.0)",
                       "option for Comet only",
                       "true");

  set_int_parameter("theoretical_fragment_ions", 1, 0, 1,
                    "0=default peak shape, 1=M peak only",
                    "option for Comet only",
                    "true");
  
  set_int_parameter("use_A_ions",
                    0, 0, 1, 
                    "Controls whether or not A-ions are considered in the "
                    "search (0 - no, 1 - yes)",
                    "option for Comet only",
                    "true");

  set_int_parameter("use_B_ions",
                    1, 0, 1, 
                    "Controls whether or not B-ions are considered in the "
                    "search (0 - no, 1 - yes)",
                    "option for Comet only",
                    "true");

  set_int_parameter("use_C_ions",
                    0, 0, 1, 
                    "Controls whether or not C-ions are considered in the "
                    "search (0 - no, 1 - yes)",
                    "option for Comet only",
                    "true");

  set_int_parameter("use_X_ions",
                    0, 0, 1, 
                    "Controls whether or not X-ions are considered in the "
                    "search (0 - no, 1 - yes)",
                    "option for Comet only",
                    "true");

  set_int_parameter("use_Y_ions",
                    1, 0, 1,
                    "Controls whether or not Y-ions are considered in the "
                    "search (0 - no, 1 - yes)",
                    "option for Comet only",
                    "true");

  set_int_parameter("use_Z_ions",
                    0, 0, 1, 
                    "Controls whether or not Z-ions are considered in the "
                    "search (0 - no, 1 - yes)",
                    "option for Comet only",
                    "true");

  set_int_parameter("use_NL_ions",
                    1, 0, 1,
                    "0=no, 1= yes to consider NH3/H2O neutral loss peak",
                    "option for Comet only",
                    "true");

  set_int_parameter("use_sparse_matrix",
                    0, 0, 1,
                    "Controls whether or not internal sparse matrix data "
                    "representation is used.",
                    "option for Comet only",
                    "true");

  set_int_parameter("output_sqtfile",
                    0, 0, 1,
                    "0=no, 1=yes  write sqt file",
                    "option for Comet only",
                    "true");

  set_int_parameter("output_pepxmlfile",
        1, 0, 1,
        "0=no, 1=yes  write pep.xml file",
        "option for Comet only",
        "true");

  set_int_parameter("output_txtfile",
        1, 0, 1,
        "0=no, 1=yes  write tab-delimited text file",
        "option for Comet only (default 1)",
        "true");
                    
  set_int_parameter("output_outfiles",
        0, 0, 1,
        "0=no, 1=yes  write .out files",
        "option for Comet only",
        "true");

  set_int_parameter("print_expect_score",
        1, 0, 1,
        "0=no, 1=yes to replace Sp with expect in out & sqt",
        "option for Comet.",
        "true"
  );

  set_int_parameter("num_output_lines",
        5, 1, BILLION,
        "num peptide results to show",
        "option for Comet.",
        "true"
  );

  set_int_parameter("show_fragment_ions",
        0, 0, 1,
        "0=no, 1=yes for out files only",
        "option for Comet.",
        "true"
  );

  set_int_parameter("sample_enzyme_number",
    1,0,10, 
    "Sample enzyme which is possibly different than the one applied to the search."
    "Used to calculate NTT & NMC in pepXML output (default=1 for trypsin).",
    "option for Comet. ",
    "true"
  );

  set_string_parameter("scan_range", "0 0",
           "start and scan scan range to search; 0 as 1st entry "
           "ignores parameter",
           "option for Comet",
           "true"
           );
  

  set_string_parameter("precursor_charge", "0 0",
           "precursor charge range to analyze; does not override "
           "mzXML charge; 0 as 1st entry ignores parameter",
           "option for Comet.",
           "true"
           );
  
  set_int_parameter("ms_level",
    2,2,3, 
    "MS level to analyze, valid are levels 2 (default) or 3",
    "option for Comet. ",
    "true"
  );

  set_string_parameter("activation_method",
    "ALL" ,
    "<string>= ALL|CID|ECD|ETD|PQD|HCD|IRMPD. Default=All",
    "option for Comet. ",
    "true"
  );

  set_string_parameter("digest_mass_range", "600.0 5000.0",
           "MH+ peptide mass range to analyze",
           "option for Comet.",
           "true"
           );
  set_int_parameter("num_results", 50,0,BILLION,
        "number of search hits to store internally",
        "option for Comet.",
        "true");

  set_int_parameter("skip_researching", 1, 0, 1,
        "for '.out' file output only, 0=search everything again "
        "(default), 1=don't search if .out exists",
        "option for Comet",
        "true");

  set_int_parameter("max_fragment_charge", 3, 1, 5,
        "set maximum fragment charge state to analyze (allowed max 5)",
        "option for Comet",
        "true");
  
  set_int_parameter("max_precursor_charge", 6, 1, 9,
        "set maximum precursor charge state to analyze (allowed max 9)",
        "option for Comet",
        "true");
  
  set_int_parameter("nucleotide_reading_frame", 0, 0, 9,
        "0=proteinDB, 1-6, 7=forward three, 8=reverse three, 9=all six",
        "option for Comet",
        "true");

  set_int_parameter("clip_nterm_methionine", 0, 0, 1,
        "0=leave sequences as-is; 1=also consider sequence w/o N-term methionine",
        "option for Comet",
        "true");

  set_int_parameter("spectrum_batch_size", 0, 0, BILLION,
        "max. # of spectra to search at a time; 0 to search the "
        "entire scan range in one loop",
        "option for Comet",
        "true");

  set_int_parameter("minimum_peaks", 10, 1, BILLION,
        "minimum num. of peaks in spectrum to search (default 10)",
        "option for Comet",
        "true");

  set_double_parameter("minimum_intensity", 0, 0, BILLION,
    "minimum intensity value to read in",
    "option for comet. ",
    "true"
  );

  set_int_parameter("remove_precursor_peak", 0, 0, 2, 
    "0=no, 1=yes, 2=all charge reduced precursor peaks (for ETD)",
    "option for Comet. ",
    "true"
  );

  set_double_parameter("remove_precursor_tolerance", 1.5, -BILLION, BILLION, 
    "+- Da tolerance for precursor removal",
    "option for Comet. ",
    "true"
  );

  set_string_parameter("clear_mz_range", "0.0 0.0",
           "for iTRAQ/TMT type data; will clear out all peaks in the specified m/z range",
           "option for Comet",
           "true"
           );

  set_string_parameter("variable_mod1", "__NULL_STR",
           "Up to 6 variable modifications are supported\n"
           "format:  <mass> <residues> <0=variable/1=binary> <max mods per a peptide>\n"
           "    e.g. 79.966331 STY 0 3",
           "option for Comet",
           "true"
           );
  
  set_string_parameter("variable_mod2", "__NULL_STR",
           "Up to 6 variable modifications are supported\n"
           "format:  <mass> <residues> <0=variable/1=binary> <max mods per a peptide>\n"
           "    e.g. 79.966331 STY 0 3",
           "option for Comet",
           "true"
           );
  
  set_string_parameter("variable_mod3",  "__NULL_STR",
           "Up to 6 variable modifications are supported\n"
           "format:  <mass> <residues> <0=variable/1=binary> <max mods per a peptide>\n"
           "    e.g. 79.966331 STY 0 3",
           "option for Comet",
           "true"
           );
  
  set_string_parameter("variable_mod4",  "__NULL_STR",
           "Up to 6 variable modifications are supported\n"
           "format:  <mass> <residues> <0=variable/1=binary> <max mods per a peptide>\n"
           "    e.g. 79.966331 STY 0 3",
           "option for Comet",
           "true"
           );
  
  set_string_parameter("variable_mod5",  "__NULL_STR",
           "Up to 6 variable modifications are supported\n"
           "format:  <mass> <residues> <0=variable/1=binary> <max mods per a peptide>\n"
           "    e.g. 79.966331 STY 0 3",
           "option for Comet",
           "true"
           );
  
  set_string_parameter("variable_mod6",  "__NULL_STR",
           "Up to 6 variable modifications are supported\n"
           "format:  <mass> <residues> <0=variable/1=binary> <max mods per a peptide>\n"
           "    e.g. 79.966331 STY 0 3",
           "option for Comet",
           "true"
           );
  
  set_int_parameter("max_variable_mods_in_peptide", 5, 0, BILLION,
        "Specifies the total/maximum number of residues that can "
        "be modified in a peptide",
        "option for Comet",
        "true"
        );

  set_double_parameter("variable_C_terminus", 0, 0, BILLION,
           "Specifiy a variable modification to peptide's c-terminus"
           "Works in conjunction with variable_c_terminus_distance",
           "option for Comet",
           "true"
           );

  set_double_parameter("variable_N_terminus", 0, 0, BILLION,
           "Specifiy a variable modification to peptide's c-terminus"
           "Works in conjunction with variable_c_terminus_distance",
           "option for Comet",
           "true");

  set_int_parameter("variable_C_terminus_distance", -1, -1, BILLION,
        "-1=all peptides, 0=protein terminus, 1-N = maximum offset from C-terminus",
        "option for Comet",
        "true");

  set_int_parameter("variable_N_terminus_distance", -1, -1, BILLION,
        "-1=all peptides, 0=protein terminus, 1-N = maximum offset from N-terminus",
        "option for Comet",
        "true");

  set_double_parameter("add_Cterm_peptide", 0, 0, BILLION,
           "Specifiy a static modification to the c-terminus of all peptides",
           "option for Comet",
           "true");

  set_double_parameter("add_Nterm_peptide", 0, 0, BILLION,
           "Specify a static modification to the n-terminus of all peptides",
           "option for Comet",
           "true");
  
  set_double_parameter("add_Cterm_protein", 0, 0, BILLION,
           "Specify a static modification to the c-terminal peptide of each protein",
           "option for Comet",
           "true");

  set_double_parameter("add_Nterm_protein", 0, 0, BILLION,
           "Specify a static modification to the n-terminal peptide of each protein",
           "option for Comet",
           "true");

  set_double_parameter("add_G_glycine", 0, 0, BILLION,
           "added to G - avg.  57.0513, mono.  57.02146",
           "option for Comet",
           "true");

  set_double_parameter("add_A_alanine", 0, 0, BILLION,
           "added to A - avg.  71.0779, mono.  71.03711",
           "option for Comet",
           "true");

  set_double_parameter("add_S_serine", 0, 0, BILLION,
           "added to S - avg.  87.0773, mono.  87.03203",
           "option for Comet",
           "true");

  set_double_parameter("add_P_proline", 0, 0, BILLION,
           "added to P - avg.  97.1152, mono.  97.05276",
           "option for Comet",
           "true");
  
  set_double_parameter("add_V_valine", 0, 0, BILLION,
           "added to V - avg.  99.1311, mono.  99.06841",
           "option for Comet",
           "true");

  set_double_parameter("add_T_threonine", 0, 0, BILLION,
           "added to T - avg. 101.1038, mono. 101.04768",
           "option for Comet",
           "true");

  set_double_parameter("add_C_cysteine", 57.021464, 0, BILLION,
           "added to C - avg. 103.1429, mono. 103.00918",
           "option for Comet",
           "true");

  set_double_parameter("add_L_leucine", 0, 0, BILLION,
           "added to L - avg. 113.1576, mono. 113.08406",
           "option for Comet",
           "true");

  set_double_parameter("add_I_isoleucine", 0, 0, BILLION,
           "added to I - avg. 113.1576, mono. 113.08406",
           "option for Comet",
           "true");

  set_double_parameter("add_N_asparagine", 0, 0, BILLION,
           "added to N - avg. 114.1026, mono. 114.04293",
           "option for Comet",
           "true");

  set_double_parameter("add_D_aspartic_acid", 0, 0, BILLION,
           "added to D - avg. 115.0874, mono. 115.02694",
           "option for Comet",
           "true");

  set_double_parameter("add_Q_glutamine", 0, 0, BILLION,
           "added to Q - avg. 128.1292, mono. 128.05858",
           "option for Comet",
           "true");

  set_double_parameter("add_K_lysine", 0, 0, BILLION,
           "added to K - avg. 128.1723, mono. 128.09496",
           "option for Comet",
           "true");

  set_double_parameter("add_E_glutamic_acid", 0, 0, BILLION,
           "added to E - avg. 129.1140, mono. 129.04259",
           "option for Comet",
           "true");

  set_double_parameter("add_M_methionine", 0, 0, BILLION,
           "added to M - avg. 131.1961, mono. 131.04048",
           "option for Comet",
           "true");

  set_double_parameter("add_O_ornithine", 0, 0, BILLION,
           "added to O - avg. 132.1610, mono  132.08988",
           "option for Comet",
           "true");

  set_double_parameter("add_H_histidine", 0, 0, BILLION,
           "added to H - avg. 137.1393, mono. 137.05891",
           "option for Comet",
           "true");

  set_double_parameter("add_F_phenylalanine", 0, 0, BILLION,
           "added to F - avg. 147.1739, mono. 147.06841",
           "option for Comet",
           "true");

  set_double_parameter("add_R_arginine", 0, 0, BILLION,
           "added to R - avg. 156.1857, mono. 156.10111",
           "option for Comet",
           "true");

  set_double_parameter("add_Y_tyrosine", 0, 0, BILLION,
           "added to Y - avg. 163.0633, mono. 163.06333",
           "option for Comet",
           "true");

  set_double_parameter("add_W_tryptophan", 0, 0, BILLION,
           "added to W - avg. 186.0793, mono. 186.07931",
           "option for Comet",
           "true");

  set_double_parameter("add_B_user_amino_acid", 0, 0, BILLION,
           "added to B - avg.   0.0000, mono.   0.00000",
           "option for Comet",
           "true");

  set_double_parameter("add_J_user_amino_acid", 0, 0, BILLION,
           "added to J - avg.   0.0000, mono.   0.00000",
           "option for Comet",
           "true");

  set_double_parameter("add_U_user_amino_acid", 0, 0, BILLION,
           "added to U - avg.   0.0000, mono.   0.00000",
           "option for Comet",
           "true");

  set_double_parameter("add_X_user_amino_acid", 0, 0, BILLION,
           "added to X - avg.   0.0000, mono.   0.00000",
           "option for Comet",
           "true");

  set_double_parameter("add_Z_user_amino_acid", 0, 0, BILLION,
           "added to Z - avg.   0.0000, mono.   0.00000",
           "option for Comet",
           "true");

  // **** q-ranker-barista arguments ****
  set_string_parameter("database", NULL,
     "The program requires the FASTA format protein database files against "
     "which the search was performed. The protein database input may be a "
     "concatenated database or separate target and decoy databases; the "
     "latter is supported with the --separate-searches option, described "
     "below. In either case, Barista distinguishes between target and decoy "
     "proteins based on the presence of a decoy prefix on the sequence "
     "identifiers (see the --decoy-prefix option, below). The database can "
     "be provided in three different ways: (1) as a a single FASTA file "
     "with suffix \".fa\", \".fsa\" or \".fasta\", (2) as a text file "
     "containing a list of FASTA files, one per line, or (3) as a directory "
     "containing multiple FASTA files (identified via the filename suffixes "
     "\".fa\", \".fsa\" or \".fasta\").", 
     "argument for barista", "false"); 

  set_string_parameter("spectra", NULL,
     "The fragmentation spectra must be provided in MS2 format. Like the "
     "database, the spectra can be specified in three different ways: (1) "
     "as a single file with suffix \".ms2\", (2) as a text file containing a "
     "list of MS2 files or (3) as a directory in which all the MS2 files can "
     "be found.", 
     "argument for q-ranker and barista", "false");
  
  set_string_parameter("search results", NULL,
     "Q-ranker recognizes search results in tab-delimited format. Like the spectra, the "
     "search results can be provided as a single file, a list of files or a "
     "directory of files. Note, however, that the input mode for spectra and "
     "for search results must be the same; i.e., if you provide a list of "
     "files for the spectra, then you must also provide a list of files "
     "containing your search results. When the MS2 files and tab-delimited text files are "
     "provided via a file listing, Q-ranker assumes that the order of the MS2 "
     "files matches the order of the tab-delimited files. Alternatively, when the MS2 "
     "files and tab-delimited files are provided via directories, Q-ranker will search "
     "for pairs of files with the same root name but different extensions "
     "(\".ms2\" and \".txt\").", 
     "argument for q-ranker and barista", "false");
  

  // **** q-ranker options. ****
  set_boolean_parameter("skip-cleanup", false, 
     "Q-ranker analysis begins with a pre-processsing step that creates a "
     "set of lookup tables which are then used during training. Normally, "
     "these lookup tables are deleted at the end of the Q-ranker analysis, "
     "but setting this option to T prevents the deletion of these tables. "
     "Subsequently, the Q-ranker analysis can be repeated more efficiently "
     "by specifying the --re-run option. Default = F.", 
     "Available for q-ranker and barista.", "true");

  set_boolean_parameter("use-spec-features", true, 
     "Q-ranker uses enriched feature set derived from the spectra in ms2 "
     "files. It can be forced to use minimal feature set by setting the "
     "--use-spec-features option to F. Default T.", 
     "Available for q-ranker and barista.", "true");

  set_string_parameter("decoy_prefix", "decoy_",
     "Specifies the prefix of the protein names that indicates a decoy. "
     "Default = decoy_.",
     " Available for q-ranker and barista.", "true");

  set_string_parameter("re-run", "__NULL_STR",
     "Re-run a previous Q-ranker analysis using a previously computed set of"
     " lookup tables.",
     " Available for q-ranker and barista.", "true");

  set_string_parameter("separate-searches", "__NULL_STR",
     "If the target and decoy searches were run separately, rather than" 
     " using a concatenated database, then Q-ranker will assume that the"
     " database search results provided as a required argument are from the"
     " target database search. This option then allows the user to specify"
     " the location of the decoy search results. Like the required arguments,"
     " these search results can be provided as a single file, a list of files"
     " or a directory. However, the choice (file, list or directory) must be"
     " consistent for the MS2 files and the target and decoy tab-delimited files. Also,"
     " if the MS2 and tab-delimited files are provided in directories, then Q-ranker"
     " will use the MS2 filename (foo.ms2) to identify corresponding target"
     " and decoy tab-delimited files with names like foo*.target.txt and"
     " foo*.decoy.txt. This naming convention allows the target and decoy txt"
     " files to reside in the same directory.",
     " Available for q-ranker and barista.", "true");
 //**** Barista and QRanker options. ******
 set_boolean_parameter("list-of-files",false, 
    "Search result can be as a file or a list of files. This option"
    " allows users to specify the search results are provided as a list of files by " 
    "setting the --list-of-files option to T."
    " Default= false.", 
    "Available for barista.","true");

  set_string_parameter("optimization", "protein",
     "Specifies whether to do optimization at the protein, peptide or psm level. "
     "Default = protein.",
     "Available for barista.", "true");


  /* analyze-matches parameter options */
  set_double_parameter("pi-zero", 1.0, 0, 1, 
      "The estimated percent of target scores that are drawn from the "
      "null distribution.",
      "Used by compute-q-values, percolator and q-ranker", "true");
  set_string_parameter("percolator-intraset-features", "F",
      "Set a feature for percolator that in later versions is not an option.",
      "Shouldn't be variable; hide from user.", "false");

  // **** predict-peptide-ions options. ****
  set_ion_type_parameter("primary-ions", BY_ION,
      "The ion series to predict (b,y,by,bya). Default='by' (both b and y ions).",
      "Only available for crux-predict-peptide-ions.  Set automatically to "
      "'by' for searching.", "true");
  set_boolean_parameter("precursor-ions", false,
      "Predict the precursor ions, and all associated ions "
      "(neutral-losses, multiple charge states) consistent with the "
      "other specified options. (T,F) Default=F.",
      "Only available for crux-predict-peptide-ions.", "true");
  set_int_parameter("isotope", 0, 0, 2,
      "Predict the given number of isotope peaks (0|1|2). Default=0.",
      "Only available for crux-predict-peptide-ion.  Automatically set to "
      "0 for Sp scoring and 1 for xcorr scoring.", "true");
  set_boolean_parameter("flanking", false, 
      "Predict flanking peaks for b and y ions (T,F). Default=F.",
      "Only available for crux-predict-peptide-ion.", "true");
  set_string_parameter("max-ion-charge", "peptide",
      "Predict theoretical ions up to max charge state (1,2,...,6) or up to the charge state "
      "of the peptide (peptide).  If the max-ion-charge is greater than the "
      "charge state of the peptide, then the max is the peptide charge. "
      "Default='peptide'.",
      "Available for predict-peptide-ions and search-for-xlinks. "
      "Set to 'peptide' for search.",
      "true");
  set_int_parameter("nh3",0, -100, BILLION, 
      "Predict peaks with the given maximum number of nh3 neutral loss "
      "modifications. Default=0.",
      "Only available for crux-predict-peptide-ions.", "true");
  set_int_parameter("h2o",0, -100, BILLION,
      "Predict peaks with the given maximum number of h2o neutral loss "
      "modifications. Default=0.",
      "Only available for crux-predict-peptide-ions.", "true");

  // ***** spectral-counts aguments *****
  set_string_parameter("input PSMs", NULL,
       "Name of file in text format which holds match results.",
       "For quantify to retrieve scores for protein and peptides.",
       "false");
  // also uses "protein-database"

  // ***** spectral-counts options *****
   set_string_parameter("input-ms2", NULL,
       "MS2 file corresponding to the psm file. Required for SIN.",
       "Available for spectral-counts with measure=SIN.",
       "true");

  set_threshold_type_parameter("threshold-type", THRESHOLD_QVALUE,
    "What type of threshold to use when parsing matches "
    "none|qvalue|custom. Default=qvalue.",
    "used for crux spectral-counts",
    "true");

  set_double_parameter("threshold", 0.01, -BILLION, BILLION, 
       "The threshold to use for filtering matches. " 
       "Default=0.01.",
       "Available for spectral-counts.  All PSMs with higher (or lower) than "
       "this will be ignored.",
       "true");

  set_boolean_parameter("custom-threshold-min", true,
    "Direction of threshold for matches.  If true, then all matches "
     "whose value is <= threshold will be accepted.  If false, "
     "then all matches >= threshold will be accepted.  Used with "
     "custom-threshold parameter (default True).",
     "Available for spectral-counts.",
     "true");

  set_string_parameter("custom-threshold-name", "__NULL_STR",
     "Use a name of custom threshold rather than (default NULL)",
     "Available for spectral-counts.", "true");


  set_measure_type_parameter("measure", MEASURE_NSAF,
       "Type of analysis to make on the match results: "
       "(RAW|NSAF|dNSAF|SIN|EMPAI). "
       "Default=NSAF. ", 
       "Available for spectral-counts.  RAW is raw counts, "
       "NSAF is Normalized Spectral Abundance Factor, "
       "dNSAF is Distributed Spectral Abundance Factor, "
       "SIN is Spectral Index Normalized and EMPAI is "
       "Exponentially Modified Protein Abundance Index",
       "true");
  set_boolean_parameter("unique-mapping", false,
       "Ignore peptides with multiple mappings to proteins (T,F). Default=F.",
       "Available for spectral-counts.",
       "true");
  set_quant_level_parameter("quant-level", PROTEIN_QUANT_LEVEL,
       "Quantification at protein or peptide level (PROTEIN,PEPTIDE). "
       "Default=PROTEIN.",
       "Available for spectral-counts and either NSAF and SIN.",
       "true"); 
  set_parsimony_type_parameter("parsimony", PARSIMONY_NONE,
       "Perform parsimony analysis on the proteins and report "
       "a parsimony rank column in output file. "
       "Default=none. Can be <string>=none|simple|greedy",
       "Available for spectral-counts.",
       "true");

  set_boolean_parameter("mzid-use-pass-threshold", false,
    "Use mzid's passThreshold attribute to filter matches. Default false.",
    "Used when parsing mzIdentML files.",
    "true");


  // ***** static mods *****
  set_double_parameter("A", 0.0, -100, BILLION, 
      "Change the mass of all amino acids 'A' by the given amount.",
      "For parameter file only.  Default=no mass change.", "true");
  set_double_parameter("B", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'B' by the given amount.",
      "For parameter file only.  Default=no mass change.", "false");
  set_double_parameter("C", 57.0214637206, -100, BILLION,
      "Change the mass of all amino acids 'C' by the given amount.",
      "For parameter file only.  Default=+57.0214637206.", "true");
  set_double_parameter("D", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'D' by the given amount.",
      "For parameter file only.  Default=no mass change.", "true");
  set_double_parameter("E", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'E' by the given amount.",
      "For parameter file only.  Default=no mass change.", "true");
  set_double_parameter("F", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'F' by the given amount.",
      "For parameter file only.  Default=no mass change.", "true");
  set_double_parameter("G", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'G' by the given amount.",
      "For parameter file only.  Default=no mass change.", "true");
  set_double_parameter("H", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'H' by the given amount.",
      "For parameter file only.  Default=no mass change.", "true");
  set_double_parameter("I", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'I' by the given amount.",
      "For parameter file only.  Default=no mass change.", "true");
  set_double_parameter("J", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'J' by the given amount.",
      "For parameter file only.  Default=no mass change.", "false");
  set_double_parameter("K", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'K' by the given amount.",
      "For parameter file only.  Default=no mass change.", "true");
  set_double_parameter("L", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'L' by the given amount.",
      "For parameter file only.  Default=no mass change.", "true");
  set_double_parameter("M", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'M' by the given amount.",
      "For parameter file only.  Default=no mass change.", "true");
  set_double_parameter("N", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'N' by the given amount.",
      "For parameter file only.  Default=no mass change.", "true");
  set_double_parameter("O", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'O' by the given amount.",
      "For parameter file only.  Default=no mass change.", "false");
  set_double_parameter("P", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'P' by the given amount.",
      "For parameter file only.  Default=no mass change.", "true");
  set_double_parameter("Q", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'Q' by the given amount.",
      "For parameter file only.  Default=no mass change.", "true");
  set_double_parameter("R", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'R' by the given amount.",
      "For parameter file only.  Default=no mass change.", "true");
  set_double_parameter("S", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'S' by the given amount.",
      "For parameter file only.  Default=no mass change.", "true");
  set_double_parameter("T", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'T' by the given amount.",
      "For parameter file only.  Default=no mass change.", "true");
  set_double_parameter("U", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'U' by the given amount.",
      "For parameter file only.  Default=no mass change.", "false");
  set_double_parameter("V", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'V' by the given amount.",
      "For parameter file only.  Default=no mass change.", "true");
  set_double_parameter("W", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'W' by the given amount.",
      "For parameter file only.  Default=no mass change.", "true");
  set_double_parameter("X", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'X' by the given amount.",
      "For parameter file only.  Default=no mass change.", "false");
  set_double_parameter("Y", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'Y' by the given amount.",
      "For parameter file only.  Default=no mass change.", "true");
  set_double_parameter("Z", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'Z' by the given amount.",
      "For parameter file only.  Default=no mass change.", "false");

  /* get-ms2-spectrum options */
  set_boolean_parameter("stats", false, 
      "Print to stdout additional information about the spectrum.",
      "Avaliable only for crux-get-ms2-spectrum.  Does not affect contents "
      "of the output file.", "true");

  // **** xlink-predict-peptide-ions options ****
  set_string_parameter("peptide A", NULL, 
      "The sequence of peptide A.",
      "Argument for xlink-predict-peptide-ions.", "false");

  set_string_parameter("peptide B", NULL, 
      "The sequence of peptide B.",
      "Argument for xlink-predict-peptide-ions.", "false");
  
  set_int_parameter("pos A", 0 , 0, BILLION, 
      "Position of xlink on peptide A",
      "Available for xlink-predict-peptide-ions.", "false");

  set_int_parameter("pos B", 0 , 0, BILLION, 
      "Position of xlink on peptide B",
      "Available for xlink-predict-peptide-ions.", "false");

  set_boolean_parameter("print-theoretical-spectrum", false,
      "Print the theoretical spectrum",
      "Available for xlink-predict-peptide-ions (Default=F).",
      "true");

  set_boolean_parameter("use-old-xlink", true /* Turn to false later */,
      "Use old xlink searching algorithm",
      "Available for search-for-xlinks program (Default=T).",
      "false");

  // **** xlink-score-spectrum options ****
  set_string_parameter("xlink-score-method", "composite", 
      "Score method for xlink {composite, modification, concatenated}. Default=composite.",
      "Argument for xlink-score-spectrum.", "false");

  // **** search-xlink options ****
  set_string_parameter("isotope-windows", "0",
    "List of integers of isotopic masses to search",
    "Used for crux search-for-xlinks", "true");

  set_boolean_parameter("xlink-print-db", false,
    "Print the database in tab delimited format to xlink_peptides.txt",
    "Used for testing the candidate generatation (Default=F).",
    "false");

  set_boolean_parameter("xlink-include-linears", true, 
      "Include linear peptides in the "
      "database.  Default=T.",
      "Available for crux search-for-xlinks program (Default=T).",
      "true");
  set_boolean_parameter("xlink-include-deadends", true, 
      "Include dead-end peptides in the "
      "database.  Default=T.",
      "Available for crux search-for-xlinks program.",
      "true");
  set_boolean_parameter("xlink-include-selfloops", true, 
      "Include self-loop peptides in the "
      "database.  Default=T.",
      "Available for crux search-for-xlinks program.",
      "true");

  
  set_string_parameter("xlink-prevents-cleavage", "K",
      "List of amino acids that xlinker can prevent cleavage",
      "Available for search-for-xlinks program (Default=K).",
      "false" /*TODO - turn this to true after new
                                      xlink code is implemented */
                        );

  set_double_parameter("precursor-window-weibull", 20.0, 0, 1e6, 
      "Search decoy-peptides within +/- "
      " 'mass-window-decoy' of the spectrum mass.  Default=20.0.",
      "Available for crux search-for-xlinks. ",
      "true");

  set_window_type_parameter("precursor-window-type-weibull", WINDOW_MASS,
      "Window type to use for selecting "
      "decoy peptides from precursor mz. <string>=mass|mz|ppm. "
      "Default=mass.",
      "Available for crux search-for-xlinks",
      "true");

  set_string_parameter("link sites", NULL, 
      "Comma delimited pair of amino acid link sites, e.g., A:K,A:D.",
      "Argument for crux search-for-xlinks.", "false");

  set_double_parameter("link mass", 0.0, -100, BILLION,
      "Mass modification of a cross link between two amino acids.",
      "Argument for crux search-for-xlinks.","false");

  set_int_parameter("min-weibull-points", 4000, 1, BILLION, 
      "Minimum number of points for estimating the "
      "Weibull parameters.  Default=4000.",
      "Available for crux search-for-xlinks", "true");

  set_int_parameter("max-xlink-mods", 0, 0, BILLION,
    "Maximum number of modifications allowed on a crosslinked peptide "
    " Default=0.",
    "Available for crux search-for-xlinks", "true");

  /* hardklor parameters */
  set_hardklor_algorithm_type_parameter(
    "hardklor-algorithm", FAST_FEWEST_PEPTIDES_HK_ALGORITHM, 
    "Choose the algorithm for analyzing combinations of "
    "multiple peptide or protein isotope distributions. "
    "(basic | fewest-peptides | fast-fewest-peptides | "
    "fewest-peptides-choice | fast-fewest-peptides-choice) "
    "Default=fast-fewest-peptides.", "Available for crux hardklor", "true");

  set_string_parameter("cdm", "Q",
    "Choose the charge state determination method. (B|F|P|Q|S). "
    "Default=Q.",
    "Available for crux hardklor", "true");

  set_int_parameter("min-charge", 1, 1, BILLION,
    "Set the minimum charge state to look for when analyzing a spectrum. "
    "Default=1.",
    "Available for crux hardklor", "true");

  set_int_parameter("max-charge", 5, 1, BILLION,
    "Set the maximum charge state to look for when analyzing a spectrum. "
    "Default=5.",
    "Available for crux hardklor", "true");

  set_double_parameter("corr", 0.85, 0,1.0, 
    "Set the correlation threshold [0,1.0] to accept a predicted "
    "isotope distribution.  Default=0.85",
    "Available for crux hardklor", "true");

  set_int_parameter("depth", 3, 1, BILLION,
    "Set the depth of combinatorial analysis. Default 3.",
    "Available for crux hardklor", "true");

  set_boolean_parameter("distribution-area", false,
    "Reports peptide intensities as the distribution area. Default false.",
    "Available for crux hardklor",
    "true");

  set_string_parameter("averagine-mod", "__NULL_STR",
    "Include alternative averagine models in the analysis that  "
    "incorporate additional atoms or isotopic enrichments.",
    "Available for crux hardklor",
    "true");

  set_string_parameter("mzxml-filter", "none",
    "Set a filter for mzXML files. Default=none",
    "Available for crux hardklor",
    "true");

  set_boolean_parameter("no-base", false,
    "Specify \"no base\" averagine. Only modified averagine models "
    "will be used in the analysis. Default = F ",
    "Available for crux hardklor",
    "true");

  set_int_parameter("max-p", 10, 1, BILLION,
    "Set the maximum number of peptides or proteins that are "
    "estimated from the peaks found in a spectrum segment. The "
    "default value is 10.",
    "Available for crux hardklor", "true");  

  set_double_parameter("resolution", 100000, 1, BILLION,
    "Set the resolution of the observed spectra at m/z 400. "
    "Used in conjunction with --instrument The default is 100000.",
    "Available for crux hardklor",
    "true");

  set_boolean_parameter("centroided", false,
    "Are spectra centroided?  Default false.",
    "Available for crux hardklor",
    "true");

  set_string_parameter("instrument", "fticr",
    "Type of instrument (fticr|orbi|tof|qit) on which the data was "
    "collected. Used in conjuction with --resolution. The default is fticr.",
    "Available for crux hardklor",
    "true");

  //scan-number already defined.

  set_int_parameter("sensitivity", 2, 0, 3,
    "Set the sensitivity level. There are four levels, 0 (low), 1 (moderate), "
    "2 (high), and 3 (max). The default value is 2.",
    "Available for crux hardklor", "true");

  set_double_parameter("signal-to-noise", 1.0, 0.0, BILLION,
    "Set the signal-to-noise threshold. Any integer or decimal "
    "value greater than or equal to 0.0 is valid. The default value is 1.0.",
    "Available for crux hardklor", "true");

  set_double_parameter("sn-window", 250.0, 0.0, BILLION,
    "Set the signal-to-noise window length (in m/z). Because noise may "
    "be non-uniform across a spectra, this value adjusts the segment size "
    "considered when calculating a signal-over-noise ratio. The default "
    "value is 250.0.",
    "Available for crux hardklor", "true");

  set_boolean_parameter("static-sn", true, 
    "If true, Hardklor will calculate the local noise levels across the "
    "spectrum using --sn-window, then select a floor of this set of noise "
    "levels to apply to the whole spectrum.",
    "Available for crux hardklor", "true");

  set_string_parameter("mz-window", "__NULL_STR",
    "Restrict analysis to only a small window in each segment ( (min-max) in m/z). "
    "The user must specify the starting and ending m/z values between which "
    "the analysis will be performed. By default the whole spectrum is analyzed.",
    "Available for crux hardklor", "true");

  set_double_parameter("max-width", 4.0, 0.0, BILLION,
    "Set the maximum width of any set of peaks in a spectrum when computing the "
    "results (in m/z). Thus, if the value was 5.0, then sets of peaks greater "
    "than 5 m/z are divided into smaller sets prior to analysis. The default value is 4.0.",
    "Available for crux hardklor", "true");

  set_string_parameter("hardklor-options", "__NULL_STR", 
    "Directly set hardklor options",
    "Available for crux hardklor", "false");

  /* bullseye parameters */
  set_string_parameter("MS1 spectra", NULL, 
    "The name of a file from which to parse high-resolution spectra of intact peptides.\n"
    " The file may be in MS1 (.ms1), binary MS1 (.bms1), compressed MS1 (.cms1), or "
    "mzXML (.mzXML) format. ",
    "Argument for crux bullseye.", "false");

  set_string_parameter("MS2 spectra", NULL, 
    "The name of a file from which to parse peptide fragmentation spectra.\n The file may "
    "be in MS2 (.ms2), binary MS2 (.bms2), compressed MS2 (.cms2) or mzXML (.mzXML) format. ",
    "Argument for crux bullseye.", "false");

  set_string_parameter("hardklor-file", "__NULL_STR",
    "Input hardklor file into bullseye",
    "Hidden option for crux bullseye.", "false");

  set_double_parameter("max-persist", 2.0, 0, BILLION,
    "Ignore peptides that persist for this length. The unit of time is whatever unit is "
    "used in your data file (usually minutes). These peptides are considered contaminants." 
    " Default = 2.0.",
    "Available for crux bullseye", "true");

  set_boolean_parameter("exact-match", false, 
    "Require an exact match to the precursor ion. Rather than use wide precursor boundaries, "
    "this flag forces Bullseye to match precursors to the base isotope peak identified in "
    "Hardklor. The tolerance is set with the --persist-tolerance flag. Default = F.",
    "Available for crux bullseye", "true");

  set_int_parameter("gap-tolerance", 1, 0, BILLION,
    "Gap size tolerance when checking for peptides across consecutive MS1 scans. Used in "
    "conjunction with --scan-tolerance. Default = 1.",
    "Available for crux bullseye", "true");

  set_double_parameter("bullseye-min-mass", 600, 0, BILLION,
      "The minimum mass of peptides to consider. Default=600.",
      "Available from command line or parameter file for crux bullseye",
      "true");

  set_double_parameter("bullseye-max-mass", 8000, 1, BILLION, 
      "The maximum mass of peptides to consider. Default=8000.",
      "Available from command line or parameter file for crux bullseye",
      "true");
  
  set_double_parameter("exact-tolerance", 10.0, 0, BILLION,
    "Set the tolerance (+/-ppm) for exact match searches. Default = 10.0.",
    "Available for crux bullseye", "true");

  set_double_parameter("persist-tolerance", 10.0, 0, BILLION,
    "Set the tolerance (+/-ppm) for finding persistent peptides. Default = 10.0.",
    "Available for crux bullseye", "true");

  set_int_parameter("scan-tolerance", 3, 0, BILLION,
    "Number of consecutive MS1 scans over which a peptide must be observed to "
    "be considered real. Gaps in persistence are allowed when setting --gap-tolerance. "
    "Default = 3.",
    "Available for crux bullseye", "true");

  set_double_parameter("retention-tolerance", 0.5, 0, BILLION,
    "Set the tolerance (+/-units) around the retention time over which a peptide "
    "can be matched to the MS/MS spectrum. The unit of time is whatever unit is "
    "used in your data file (usually minutes). Default = 0.5.",
    "Available for crux bullseye", "true");

  set_string_parameter("spectrum-format", "__NULL_STR",
    "The format to write the output spectra to. By default, the spectra will be "
    "output in the same format as the MS/MS input.",
    "Available for crux bullseye", "true");

  /* crux-util parameters */

  set_boolean_parameter("ascending", true,
    "Sort in ascending order.  Otherwise, descending. "
    "Default: True.",
    "Available for sort-by-column", "true");

  set_string_parameter("tsv file", NULL,
    "Path to a delimited file (-) for standard input",
    "Available for the delimited utility programs", "false");

  set_string_parameter("delimiter", "tab",
    "Character delimiter to use when parsing a delimited file (Default tab).",
    "Available for the delimited utility programs.", "false");

  set_string_parameter("column names", NULL,
    "List of column names separated by a comma",
    "Available for the delimited utility", "false");

  set_string_parameter("column name", NULL,
    "Name of the column to do the operation on",
    "Available for the delimited utility programs", "false");

  set_string_parameter("column value", NULL,
    "value of the column",
    "Available for the delimited utility programs", "false");

  set_boolean_parameter("header", true,
    "Print the header line of the tsv file. Default=T.",
    "Available for crux extract-columns and extract-rows",
    "true");

  set_string_parameter("column-type","string",
    "Specifies the data type the column contains "
    "(int|real|string) Default: string",
    "Available for crux extract-rows",
    "true");

  set_string_parameter("comparison", "eq",
    "Specifies the operator that is used to compare an "
    "entry in the specified column to the value given "
    "on the command line.  (eq|gt|gte|lt|lte|neq). "
    "Default: eq.",
    "Available for crux extract-rows",
    "true");

  // Default comet enzyme lines
  comet_enzyme_info_lines_.push_back("0.  No_enzyme\t\t\t\t0       -           -");
  comet_enzyme_info_lines_.push_back("1.  Trypsin\t\t\t\t1      KR           P");
  comet_enzyme_info_lines_.push_back("2.  Trypsin/P\t\t\t\t1      KR           -");
  comet_enzyme_info_lines_.push_back("3.  Lys_C\t\t\t\t1      K            P");
  comet_enzyme_info_lines_.push_back("4.  Lys_N\t\t\t\t0      K            -");
  comet_enzyme_info_lines_.push_back("5.  Arg_C\t\t\t\t1      R            P");
  comet_enzyme_info_lines_.push_back("6.  Asp_N\t\t\t\t0      D            -");
  comet_enzyme_info_lines_.push_back("7.  CNBr\t\t\t\t1      M            -");
  comet_enzyme_info_lines_.push_back("8.  Glu_C\t\t\t\t1      DE           P");
  comet_enzyme_info_lines_.push_back("9.  PepsinA\t\t\t\t1      FL           P");
  comet_enzyme_info_lines_.push_back("10. Chymotrypsin\t\t\t1      FWYL         P");
 

 
  // now we have initialized the parameters
  parameter_initialized = true;
  usage_initialized = true;
  type_initialized = true;

}


/*
 * Main calls this to determine which required arguments
 * must be used.  Must be called after initialize_parameters
 * This is the interface for parse_argument_set_req()
 */
bool select_cmd_line_arguments(  //remove options from name
  const char** option_names,
  int    num_options 
  ){
  select_cmd_line( option_names, num_options, 
                   &parse_arguments_set_req);
  return true;
}

/*
 * Main calls this to determine which options
 * can be used.  Must be called after initialize_parameters
 * This is the interface for parse_argument_set_opt()
 */
bool select_cmd_line_options(  //remove options from name
  const char** option_names,
  int    num_options 
  ){

  vector<bool> option_found(num_options, false);
  bool success = false;  

  HASH_ITERATOR_T* hash_iter = new_hash_iterator(parameters);
  while(hash_iterator_has_next(hash_iter)) {
    char* option_name = hash_iterator_next(hash_iter);
    if (strstr(option_name, " ") == NULL) {
      void* value_ptr = get_hash_value(parameters, option_name);
      void* usage_ptr = get_hash_value(usages, option_name);
      void* type_ptr =  get_hash_value(types, option_name);

      /* check that the option is in the params hash */
      if( value_ptr == NULL || usage_ptr == NULL || type_ptr == NULL ){
        carp(CARP_FATAL, 
           "Cannot select parameter '%s'. Value, usage or type not found. "
           "Found value: %s, usage: %s, type: %s", 
           option_name,
           value_ptr,
           usage_ptr,
           type_ptr);
      
      }

      if( //strcmp(type_ptr, "PEPTIDE_TYPE_T") == 0 ||
          strcmp((char*)type_ptr, "MASS_TYPE_T") == 0 ||
          strcmp((char*)type_ptr, "bool") == 0 ||
          strcmp((char*)type_ptr, "SCORER_TYPE_T") == 0 ){
        type_ptr = (void*)"STRING_ARG";
      }
      carp(CARP_DETAILED_DEBUG, 
           "Found value: %s, usage: %s, type(to be passed to parse_args): %s", 
           (char*)value_ptr, (char*)usage_ptr, (char*)type_ptr);
      /* is it in the list of options to print?
       */
      bool print = false;
      for (size_t idx=0;idx < num_options;idx++) {
        if (strcmp(option_name, option_names[idx]) == 0) {
          print = true;
          option_found[idx] = true;
          break;
        }
      }
    
      /* add the option via parse_arguments.c. pointer decides opt or req */
      success &= parse_arguments_set_opt(option_name,
                                      (const char*)usage_ptr,
                                      value_ptr, 
                                      string_to_argument_type((char*)type_ptr),
                                      print); 
    }
  }
  
  if (success) {
    /* Check to see if all options to print are in the entire list of options*/
    for (size_t idx = 0;idx < num_options;idx++) {
      if (!option_found[idx]) {
        void* value_ptr = get_hash_value(parameters, option_names[idx]);
        void* usage_ptr = get_hash_value(usages, option_names[idx]);
        void* type_ptr =  get_hash_value(types, option_names[idx]);
        carp(CARP_ERROR, 
           "Cannot select parameter '%s'. Value, usage or type not found. "
           "Found value: %s, usage: %s, type: %s", 
           option_names[idx],
           value_ptr,
           usage_ptr,
           type_ptr);
      }
      success = false;
    }
  }
  
  carp(CARP_DETAILED_DEBUG, "Did setting the arguments work? %i", success);
  return success;
}
/*
 * Private function for doing the work of select_cmd_line_options
 * and select_cmd_line_arguments which is all the same except for
 * the last function call which is now set with a function pointer
 * 
 */
bool select_cmd_line(  //remove options from name
  const char** option_names,
  int    num_options, 
  int (*parse_arguments_set_ptr)(const char*, const char*, void*, enum argument_type, bool print) 
  ){

  carp(CARP_DETAILED_DEBUG, "Selecting options");
  bool success = true;

  if( (num_options < 1) || (option_names == NULL) ){
    success = false; //?
    return success;
  }

  /* for each option name in list */
  int i;
  for( i=0; i< num_options; i++){
    //carp(CARP_INFO, "%i Option is: %s", i, option_names[i]);
    /* get value, usage, types */
    void* value_ptr = get_hash_value(parameters, option_names[i]);
    void* usage_ptr = get_hash_value(usages, option_names[i]);
    void* type_ptr =  get_hash_value(types, option_names[i]);

    /* check that the option is in the params hash */
    if( value_ptr == NULL || usage_ptr == NULL || type_ptr == NULL ){
      carp(CARP_FATAL, 
           "Cannot select parameter '%s'. Value, usage or type not found. "
           "Found value: %s, usage: %s, type: %s", 
           option_names[i],
           value_ptr,
           usage_ptr,
           type_ptr);
      
    }

    if( //strcmp(type_ptr, "PEPTIDE_TYPE_T") == 0 ||
        strcmp((char*)type_ptr, "MASS_TYPE_T") == 0 ||
        strcmp((char*)type_ptr, "bool") == 0 ||
        strcmp((char*)type_ptr, "SCORER_TYPE_T") == 0 ){
      type_ptr = (void*)"STRING_ARG";
    }
    carp(CARP_DETAILED_DEBUG, 
         "Found value: %s, usage: %s, type(to be passed to parse_args): %s", 
         (char*)value_ptr, (char*)usage_ptr, (char*)type_ptr);
    
    /* add the option via parse_arguments.c. pointer decides opt or req */
    success = parse_arguments_set_ptr(option_names[i],
                                      (const char*)usage_ptr,
                                      value_ptr, 
                                      string_to_argument_type((char*)type_ptr), true); 
  }

  carp(CARP_DETAILED_DEBUG, "Did setting the arguments work? %i", success);
  return success;
}

/**
 * helper used below.  look for param file name, die if error
 * return null if not found
 */
bool find_param_filename(int argc, 
                              char** argv, 
                              char* filename_buffer, 
                              int buffer_size){
  bool success = true;
  int i;
  int param_file_index = -1;
  for( i=0; i< argc; i++){
    if( strcmp(argv[i], "--parameter-file") == 0){
      param_file_index = i+1;
      break;
    }
  }
  // check for error
  if( param_file_index >= argc ){
    carp(CARP_FATAL, "Option '--parameter-file' requires argument");
  }
  //return the filename
  else if( param_file_index > 0 ){  
    char* param_filename = argv[param_file_index];
    carp(CARP_DETAILED_DEBUG, "Parameter file name is %s", param_filename);

    if( strlen(param_filename) < (unsigned)buffer_size ){
      strcpy(filename_buffer, param_filename);
      success = true;
    }
    else{
      carp(CARP_FATAL, "Parameter filename is too long");
    }
  }
  else{ //parameter_file_index < 0, i.e. no paramter file option
    success = false;
  }

  return success;
}

/**
 * \brief Take the user options related to decoys to set the number of
 * decoy files and to adjust the number of top Sp-ranked PSMs to
 * score with xcorr.  Number of decoy files is 0 if decoy location is
 * "target", 1 if location is "one-decoy-file", or
 * num-decoys-per-target.  Max rank prelimiary is adjusted if location
 * is target and num-decoys-per-target > 1.
 */
void translate_decoy_options(){

  // get user values
  int num_decoy_per_target = get_int_parameter("num-decoys-per-target");
  string decoy_type_str = get_string_parameter_pointer("decoys");
  DECOY_TYPE_T decoy_type = string_to_decoy_type(decoy_type_str.c_str());

  char* location = get_string_parameter("decoy-location");
  int max_rank_preliminary = get_int_parameter("max-rank-preliminary");

  // store new values here
  bool tdc = false;  // target-decoy competitition
  int new_num_decoy_files = -1;
  int new_max_rank_preliminary = max_rank_preliminary; 

  // user may not have set num-decoys-per-target and default is 1
  if( decoy_type == NO_DECOYS ){
    num_decoy_per_target = 0;
  } 

  // user may not have set target-location if no decoys requested
  if( num_decoy_per_target == 0 ){
    free(location);
    location = my_copy_string("separate-decoy-files");
    new_num_decoy_files = 0;
  }

  // set new values
  if( strcmp(location, "target-file") == 0 ){
    tdc = true;
    new_num_decoy_files = 0;

    if( max_rank_preliminary > 0 ){  // scale to num decoys
      new_max_rank_preliminary = max_rank_preliminary * 
                                (1 + num_decoy_per_target);
    }
  }else if( strcmp(location, "one-decoy-file") == 0 ){
    tdc = false;
    new_num_decoy_files = 1;
  }else if( strcmp(location, "separate-decoy-files") == 0 ){
    tdc = false;
    new_num_decoy_files = num_decoy_per_target;
  }else{
    carp(CARP_FATAL, "Unrecoginzed decoy location '%s'."
         "Must be target-file, one-decoy-file, or separate-decoy-files",
         location);
  }

  free(location);

  // now update all values
  char buffer[PARAMETER_LENGTH];
  sprintf(buffer, "%i", new_num_decoy_files);
  update_hash_value(parameters, "num-decoy-files", buffer);

  sprintf(buffer, "%i", new_max_rank_preliminary);
  update_hash_value(parameters, "max-rank-preliminary", buffer);

  if( tdc == true ){
    update_hash_value(parameters, "tdc", (void*)"true");
  }else{
    update_hash_value(parameters, "tdc", (void*)"false");
  }
}

/**
 * Set the m/z bin width.  If the user does not request a specific
 * width, then use pre-defined values depending on the fragment mass
 * type.
 */
static void set_mz_bin_width()
{
  double new_value = get_double_parameter("mz-bin-width");

#ifdef _MSC_VER
  // Peculiarities of Windows floating point handling 
  // results in us getting 0.0 here rather than Nan
  // FIXME: is there a more portable way of checking
  // that a floating point value has not been set?
  if (new_value == 0.0) {
#else
  if (isnan(new_value)) {
#endif
    // If no width specified, choose based on mass type.
    if (get_mass_type_parameter("fragment-mass") == MONO) {
      new_value = BIN_WIDTH_MONO;
    } else {
      new_value = BIN_WIDTH_AVERAGE;
    }

    // Update the parameter hash.
    char buffer[PARAMETER_LENGTH];
    snprintf(buffer, PARAMETER_LENGTH, "%f", new_value);
    add_or_update_hash(parameters, "mz-bin-width", buffer);
  }
}

/**
 * Take the command line string from main, find the parameter file 
 * option (if present), parse it's values into the hash, and parse
 * the command line options and arguments into the hash.
 * Main then retrieves the values through get_<type>_parameter.
 * \returns true is command line is successfully parsed.
 */
bool parse_cmd_line_into_params_hash(int argc, 
                                     char** argv, 
                                     const char* exe_name){
  carp(CARP_DETAILED_DEBUG, "Parameter.c is parsing the command line");
  assert(parameter_initialized && usage_initialized && type_initialized);
  bool success = true;
  int i;
  /* first look for parameter-file option and parse values in file before
     command line values.  Checks types and bounds, exiting if invalid */

  char param_filename[SMALL_BUFFER];
  if(find_param_filename(argc, argv, param_filename, SMALL_BUFFER)){
    parse_parameter_file(param_filename);  
    read_mods_from_file(param_filename);
  }
  else{ 
    carp(CARP_INFO, 
      "No parameter file specified.  Using defaults and command line values");
  }

  /* now parse the command line using parse_arguments.c, 
     check options for legal values, and put values in hash 
     overwriting file parameters */ 

  success = parse_arguments_into_hash(argc, argv, parameters, 0); 

  if( success ){
    // check each option value
    for(i=1; i<argc; i++){
      char* word = argv[i];
      if( word[0] == '-' && word[1] == '-' ){   //if word starts with --
        word = word + 2;      //ignore the --
        check_option_type_and_bounds(word);
      }//else skip this word
    }
  }
  else{  // parse_arguments encountered an error
    char* error_message = NULL;
    char* usage = parse_arguments_get_usage(exe_name);
    int error_code = parse_arguments_get_error(&error_message);
    carp(
      CARP_FATAL, 
      "Error in command line. Error # %d\n%s\n%s", 
      error_code,
      error_message,
      usage
    );
  }
  
  // do global checks on parameters
  check_parameter_consistency();

  // do global parameter-specific tasks: set static aa modifications
  //   set verbosity, make adjustments for tdc, print param file (if
  //   requested), set custom enzyme rules, set max psms per spec to file

  update_aa_masses();

  translate_decoy_options();
  

  // for compute-q-values, set algorithm to q-value (perc by default)
  if( strcmp(argv[0], "compute-q-values") == 0 ){
    char value_str[SMALL_BUFFER];
    algorithm_type_to_string(QVALUE_ALGORITHM, value_str);
    update_hash_value(parameters, "algorithm", value_str);
  }
  // for q-ranker, set algorithm to qranker (perc by default)
  else if( strcmp(argv[0], "q-ranker") == 0 ){
    char value_str[SMALL_BUFFER];
    algorithm_type_to_string(QRANKER_ALGORITHM, value_str);
    update_hash_value(parameters, "algorithm", value_str);
  }

  // if custom-enzyme used, set values
  char* enzyme_rule_str = get_string_parameter("custom-enzyme");
  if( enzyme_rule_str != NULL ){
    parse_custom_enzyme(enzyme_rule_str);
    free(enzyme_rule_str);

    char* value_str = enzyme_type_to_string(CUSTOM_ENZYME);
    update_hash_value(parameters, "enzyme", value_str);
    free(value_str);
  }

  // user may set --enzyme OR --digestion, update the other
  // if no-enzyme is used, set digestion as non-specific
  if( get_enzyme_type_parameter("enzyme") == NO_ENZYME ){
    char* value_str = digest_type_to_string(NON_SPECIFIC_DIGEST);
    update_hash_value(parameters, "digestion", value_str);
    free (value_str);

  }else if(get_digest_type_parameter("digestion") == NON_SPECIFIC_DIGEST ){
    char* value_str = enzyme_type_to_string(NO_ENZYME);
    update_hash_value(parameters, "enzyme", value_str);
    free (value_str);
  }

  set_verbosity_level(get_int_parameter("verbosity"));

  // what is the max psms per spec to keep for printing
  int max = get_int_parameter("top-match");
  char value_str[SMALL_BUFFER];
  sprintf(value_str, "%i", max);
  update_hash_value(parameters, "psms-per-spectrum-reported", value_str);

  parameter_plasticity = false;

  // Set m/z bin width based on mass type.
  set_mz_bin_width();

  return success;
}

/**
 * Read the value given for custom-enzyme and enter values into global
 * params.  Correct syntax is [A-Z]|[A-Z] or {A-Z}|{A-Z}.  An X
 * indicates that any residue is legal. Sets pre/post_list size and
 * allocates memory for pre/post_cleavage_list.  Sets
 * pre/post_for_inclusion as true if [] encloses list or false if {}
 * encloses list.  For special case of [X], set p_cleavage_list as
 * empty and inclusion as false.
 */
// NOTE (BF mar-11-09): for testing would be nice if this returned
// error code instead of dying
void parse_custom_enzyme(const char* rule_str){

  bool success = true;
  int len = strlen(rule_str);
  int idx = 0;
  int pipe_idx = 0;

  // 1. find the |
  for(idx = 0; idx < len; idx++){
    if( rule_str[idx] == '|' ){
      pipe_idx = idx;
      break;
    }
  }
  // check that there isn't a second
  for(idx = idx+1; idx < len; idx++){
    if( rule_str[idx] == '|' ){
      success = false;      
      break;
    }
  }


  // 2. set beginning and end of strings relative to pipe, start, end
  //    0 1    p-1 p p+1 p+2     len-1 len
  //    [ X    ]   | [   X       ]     '0'
  int pre_first_idx = 1;
  int pre_end_idx = pipe_idx - 1;
  int post_first_idx = pipe_idx + 2;
  int post_end_idx = len -1;

  // 3. check that braces match and set inclusion
  // pre-list
  if(pipe_idx < 1){
    success = false;
  }else if(rule_str[pre_first_idx-1] == '[' && 
           rule_str[pre_end_idx] == ']'){
    pre_for_inclusion = true;
  }else if(rule_str[pre_first_idx-1] == '{' && 
           rule_str[pre_end_idx] == '}'){
    pre_for_inclusion = false;
  }else{
    success = false;
  }

  // post list
  if(pipe_idx + 2 >= len ){
    success = false;
  }else if(rule_str[post_first_idx-1] == '[' && 
           rule_str[post_end_idx] == ']'){
    post_for_inclusion = true;
  }else if(rule_str[post_first_idx-1] == '{' && 
           rule_str[post_end_idx] == '}'){
    post_for_inclusion = false;
  }else{
    success = false;
  }

  // check that braces aren't empty 
  if(pre_first_idx >= pre_end_idx || post_first_idx >= post_end_idx ){
    success = false;
  }

  if( success == false ){
    carp(CARP_FATAL, "Custom enzyme syntax '%s' is incorrect.  "
         "Must be of the form [AZ]|[AZ] or with [] replaced by {}. "
         "AZ is a list of residues (letters A-Z) required [] or prohibited {}. "
         "Use [X] to indicate any reside is legal.",
         rule_str);
  }

  // 4. allocate lists and fill
  pre_list_size = pre_end_idx - pre_first_idx;
  pre_cleavage_list = (char*)mycalloc(pre_list_size, sizeof(char));
  for(idx = 0; idx < pre_list_size; idx++){
    pre_cleavage_list[idx] = rule_str[pre_first_idx+idx];
  }

  post_list_size = post_end_idx - post_first_idx;
  post_cleavage_list = (char*)mycalloc(post_list_size, sizeof(char));
  for(idx = 0; idx < post_list_size; idx++){
    post_cleavage_list[idx] = rule_str[post_first_idx+idx];
  }


  // 5. check special case of [X]
  if(strncmp( rule_str, "[X]", pre_list_size+2) == 0){
    free(pre_cleavage_list);
    pre_cleavage_list = NULL;
    pre_list_size = 0;
    pre_for_inclusion = false;
  }

  if(strncmp( rule_str+post_first_idx-1, "[X]", post_list_size+2) == 0){
    free(post_cleavage_list);
    post_cleavage_list = NULL;
    post_list_size = 0;
    post_for_inclusion = false;
  }

}

void check_parameter_consistency(){

  /* Min length/mass is less than max */
  int min_length = get_int_parameter("min-length");
  int max_length = get_int_parameter("max-length");

  if( min_length > max_length){
    carp(CARP_FATAL, "Parameter inconsistency.  Minimum peptide length (%i)"
         " must be less than max (%i).", min_length, max_length);
  }

  double min_mass = get_double_parameter("min-mass");
  double max_mass = get_double_parameter("max-mass");

  if( min_mass > max_mass){
    carp(CARP_FATAL, "Parameter inconsistency.  Minimum peptide mass (%.2f)"
         " must be less than max (%.2f).", min_mass, max_mass);
  }

  double min_spec_mz = get_double_parameter("spectrum-min-mz");
  double max_spec_mz = get_double_parameter("spectrum-max-mz");

  if( min_spec_mz > max_spec_mz){
    carp(CARP_FATAL, "Parameter inconsistency. Minimum spectrum mz (%.2f)"
         " must be less than max (%.2f).", min_spec_mz, max_spec_mz);
  }

  /* If no-enzyme, set digestion to non-specific and missed to true */
  if( get_enzyme_type_parameter("enzyme") == NO_ENZYME ){
    char* val_str = digest_type_to_string(NON_SPECIFIC_DIGEST);
    update_hash_value(parameters, "digestion", val_str);
    free(val_str);
    update_hash_value(parameters, "missed-cleavages", (void*)"500");
  }

  // spectral counts SIN requires an MS2 file
  if( get_measure_type_parameter("measure") == MEASURE_SIN ){
    if( strcmp(get_string_parameter_pointer("input-ms2"), "__NULL_STR") == 0 ){
      carp(CARP_FATAL, "The SIN computation for spectral-counts requires "
           "that the --input-ms2 option specify a file.");
    }
  }

  // decoys must be one of "none", "protein-shuffle", "peptide-shuffle",
  // "reverse"
  const char* decoys = get_string_parameter_pointer("decoys");
  DECOY_TYPE_T decoy_type = string_to_decoy_type(decoys);
  if( decoy_type == INVALID_DECOY_TYPE ){
    carp(CARP_FATAL, "The 'decoys' option must be 'none', 'reverse', "
         "'protein-shuffle', or 'peptide-shuffle'.  '%s' is invalid.", decoys);
  }
}


/*
 * Checks the current value of the named option
 *   as stored in the parameters hash and checks 
 *   that it is a legal value (within min/max for
 *   numeric, correct word for specialized type 
 */
bool check_option_type_and_bounds(const char* name){

  bool success = true;
  char die_str[SMALL_BUFFER];
  char* type_str = (char*)get_hash_value(types, name);
  char* value_str = (char*)get_hash_value(parameters, name);
  char* min_str = (char*)get_hash_value(min_values, name);
  char* max_str = (char*)get_hash_value(max_values, name);

  MASS_TYPE_T mass_type;
  SCORER_TYPE_T scorer_type;
  ALGORITHM_TYPE_T algorithm_type;
  ION_TYPE_T ion_type;

  PARAMETER_TYPE_T param_type;
  string_to_param_type( type_str, &param_type ); 

  carp(CARP_DETAILED_DEBUG, 
       "Checking option '%s' of type '%s' for type and bounds", 
       name, type_str);

  switch( param_type ){
  case INT_P:
  case DOUBLE_P:
    if( atof(value_str) < atof(min_str) || 
        atof(value_str) > atof(max_str) ){
      success = false;
      sprintf(die_str, 
              "The option '%s' must be between %s and %s.  %s is out of bounds",
              name, min_str, max_str, value_str);
    }
    break;
  case STRING_P:
    carp(CARP_DETAILED_DEBUG, "found string opt with value %s ", value_str);
    //check list of legal values?
    break;
  case MASS_TYPE_P:
    carp(CARP_DETAILED_DEBUG, "found mass_type opt with value %s ", 
         value_str);
    if( ! string_to_mass_type( value_str, &mass_type )){
      success = false;
      sprintf(die_str, "Illegal mass-type.  Must be 'mono' or 'average'");
    }
    break;
  case DIGEST_TYPE_P:
      carp(CARP_DETAILED_DEBUG, "found digest_type param, value '%s' ", 
           value_str);
    if( string_to_digest_type(value_str) == INVALID_DIGEST){
      success = false;
      sprintf(die_str, "Illegal digest value. "
              "Must be full-digest or partial-digest.");
    }
    break;
  case ENZYME_TYPE_P:
      carp(CARP_DETAILED_DEBUG, "found enzyme_type param, value '%s' ", 
           value_str);
    if( string_to_enzyme_type(value_str) == INVALID_ENZYME){
      success = false;
      sprintf(die_str, "Illegal enzyme. Must be no-enzyme, trypsin,trypsin/p, chymotrypsin, " 
      "elastase,clostripain, cyanogen-bromide, iodosobenzoate, " 
      "proline-endopeptidase, staph-protease, aspn, lys-c, "
      "lys-n, arg-c , glue-c, pepsin-a, modified-chymotrypsin, elastase-trypsin-chymotrypsin, "
      "custom-enzyme trypsin, chymotrypsin, elastase, or no-enzyme.");
    }
    break;
  case BOOLEAN_P:
    carp(CARP_DETAILED_DEBUG, "found boolean_type param, value '%s'", 
         value_str);
    if( (value_str[0] != 'T') && (value_str[0] != 'F') && 
        (value_str[0] != 't') && (value_str[0] != 'f') ){
      success =  false;
      sprintf(die_str, 
              "Illegal boolean value '%s' for option '%s'.  Must be T/t or F/f",
              value_str, name);
    }
    break;
  case SCORER_TYPE_P:
    carp(CARP_DETAILED_DEBUG, "found scorer_type param, value '%s'",
         value_str);
    //check for legal type
    if(! string_to_scorer_type( value_str, &scorer_type)){
      success = false;
      sprintf(die_str, "Illegal score value '%s' for option '%s'.  "
      "Must be sp, xcorr or xcorr-pvalue.", value_str, name);
    }else if((scorer_type != SP ) &&   //check for one of the accepted types
             (scorer_type != XCORR ) &&
             (scorer_type != LOGP_BONF_WEIBULL_XCORR )){
      success = false;
      sprintf(die_str, "Illegal score value '%s' for option '%s'.  "
      "Must be sp, xcorr or xcorr-pvalue.", value_str, name);
    }
    break;
  case ALGORITHM_TYPE_P:
    carp(CARP_DETAILED_DEBUG, "found algorithm_type param, value '%s'",
         value_str);
    if(! string_to_algorithm_type( value_str, &algorithm_type)){
      success = false;
      sprintf(die_str, "Illegal score value '%s' for option '%s'.  "
              "Must be percolator, curve-fit, or none.", value_str, name);
    }
    break;
  case HARDKLOR_ALGORITHM_TYPE_P:
    if (string_to_hardklor_algorithm_type(value_str) == INVALID_HK_ALGORITHM) {
      success = false;
      sprintf(die_str, "Illegal value '%s' for option '%s'.   "
                       "Must be basic, fewest-peptides, fast-fewest-peptides, "
                       "fewest-peptides-choice, or fast-fewest-peptides-choice",
                        value_str, name);
    }
    break;
  case SPECTRUM_PARSER_P:
    if (string_to_spectrum_parser_type(value_str) == INVALID_SPECTRUM_PARSER) {
      success = false;
      sprintf(die_str, "Illegal value '%s' for option '%s'.   "
                       "Must be pwiz or mstoolkit",
                       value_str, name);
    }
    break;
  case ION_TYPE_P:
    carp(CARP_DETAILED_DEBUG, "found ion_type param, value '%s'",
         value_str);
    if( !string_to_ion_type(value_str, &ion_type)){
      success = false;
      sprintf(die_str, "Illegal ion type '%s' for option '%s'.  "
              "Must be b,y,by,bya.", value_str, name);
    }
    break;
  case WINDOW_TYPE_P:
    carp(CARP_DETAILED_DEBUG, "found window type param, value '%s'",
         value_str);
    if(string_to_window_type(value_str) == WINDOW_INVALID) {
      success = false;
      sprintf(die_str, "Illegal window type '%s' for option '%s'.  "
              "Must be (mass, mz, ppm)", value_str, name);
    }
    break;
  case MEASURE_TYPE_P:
    carp(CARP_DETAILED_DEBUG, "found measure type param, value '%s'",
         value_str);
    if (string_to_measure_type(value_str) == MEASURE_INVALID){
      success = false;
      sprintf(die_str, "Illegal measure type '%s' for option '%s'. "
              "Must be (NSAF, SIN, EMPAI)", value_str, name);
    }
    break;
  case THRESHOLD_P:
    carp(CARP_DETAILED_DEBUG, "found threshold type param, value '%s'",
       value_str);
    if (string_to_threshold_type(value_str) == THRESHOLD_INVALID) {
      success = false;
      sprintf(die_str, "Illegal threshold type '%s' for option '%s'."
                       "Must be (none, qvalue, custom).", value_str, name);
    }
    break;
  case PARSIMONY_TYPE_P:
    carp(CARP_DETAILED_DEBUG, "found parsimony type param, value '%s'",
     value_str);
    if (string_to_parsimony_type(value_str) == PARSIMONY_INVALID){
      success = false;
      sprintf(die_str, "Illegal parsimony type '%s' for option '%s'. "
              "Must be (none, simple, greedy)", value_str, name);
    }
    break;
  case QUANT_LEVEL_TYPE_P:
    carp(CARP_DETAILED_DEBUG, "found quant level type param, value");// '%s'",
    //value_str);
    if (string_to_quant_level_type(value_str) == QUANT_LEVEL_INVALID){
      success = false;
      sprintf(die_str, "Illegal quantitation level type '%s' for option '%s'. "
              "Must be (peptide, protein)", value_str, name);
    }
    break;
  case DECOY_TYPE_P:
    carp(CARP_DETAILED_DEBUG, "found decoy type param, value");
    if (string_to_decoy_type(value_str) == INVALID_DECOY_TYPE){
      success = false;
      sprintf(die_str, "Illegal decoy type '%s' for option '%s'. "
              "Must be one of none, reverse, protein-shuffle, peptide-shuffle",
              value_str, name);
    }
    break;
  case MASS_FORMAT_P:
    carp(CARP_DEBUG, "found mass format param, value %s", value_str);
    if (string_to_mass_format(value_str) == INVALID_MASS_FORMAT){
      success = false;
      sprintf(die_str, "Illegal mass format '%s' for option '%s'. "
              "Must be one of mod-only, total, separate",
              value_str, name);
    }
    break;

  case NUMBER_PARAMETER_TYPES:
    carp(CARP_FATAL, "Your param type '%s' wasn't found (code %i)", 
        type_str, (int)param_type);
  } // end switch

  if( ! success ){
    carp(CARP_FATAL, die_str);
  }
  return success;
}

/**
 * Maximum size of the description of a parameter.
 */
static const int PARAMETER_BUFFER = 10000;

/**
 *  Print modifications to the given file in the format used in the
 *  parameter file.  Expected names are "mod", "cmod", "nmod".  The
 *  function to get the list of modifications should correspond to the
 *  mod name.  (i.e. for "mod", use get_aa_mod_list(), for "cmod" use
 *  get_c_mod_list(), for "nmod" use get_n_mod_list())
 */
void print_mods_parameter_file(FILE* param_file, 
                               const char* name,
                               int (*mod_getter)(AA_MOD_T***)){
  // get mod description
  char comments[PARAMETER_BUFFER] = "";
  strcat_formatted(comments, "# ", (char*)get_hash_value(usages, name));
  strcat_formatted(comments, "# ", (char*)get_hash_value(file_notes, name));
  int precision = get_int_parameter("mod-precision");

  // get list of mods to print
  AA_MOD_T** mod_list = NULL;
  int total_mods = (*mod_getter)(&mod_list);
  for( int mod_idx = 0 ; mod_idx < total_mods; mod_idx++){
    float mass = aa_mod_get_mass_change(mod_list[mod_idx]);

    // standard mods have the format mass:aa_list:max
    if( strcmp(name, "mod") == 0 ){
      int max = aa_mod_get_max_per_peptide(mod_list[mod_idx]);
      bool* aas_modified = aa_mod_get_aa_list(mod_list[mod_idx]);
      char aa_str[PARAMETER_BUFFER] = "";
      char* aa_str_ptr = aa_str;
      for(int aa_idx = 0; aa_idx < AA_LIST_LENGTH; aa_idx++){
        if( aas_modified[aa_idx] == true ){
          sprintf(aa_str_ptr, "%c", (aa_idx + 'A'));
          aa_str_ptr++;
        }
      }
      fprintf(param_file, "%s%s=%.*f:%s:%i\n\n", comments, name, precision, 
              mass, aa_str, max); 
    } else { // nmod, cmod have the format mass:end_distance
      int distance = aa_mod_get_max_distance(mod_list[mod_idx]);
      fprintf(param_file, "%s%s=%.*f:%i\n\n", comments, name, precision, 
              mass, distance);
    }
  }

  // if there were no mods, print placeholder
  if( total_mods == 0 ){
    fprintf(param_file, "%s%s=NO MODS\n\n", comments, name);
  }
}


/**
 * \brief prints all parameters except mods into the output stream
 * in xml format. 
 *
 * Each parameter has a self closing tag and has attributes name for 
 * parameter name and value for parameter value
 */
void print_parameters_xml(FILE* output){
  if (output == NULL){
    return;
  }
  carp(CARP_DEBUG, "Printing parameters to xml output");
  HASH_ITERATOR_T* iterator = new_hash_iterator(parameters);
  while (hash_iterator_has_next(iterator)){
    char * key = hash_iterator_next(iterator);
    char * show_users = (char *)get_hash_value(for_users, key);
    if ( strcmp(show_users, "true") == 0 ){
      if (strcmp(key, "mod") == 0 || strcmp(key, "cmod") == 0
          || strcmp(key, "nmod") == 0){
        continue;
      }
      fprintf(output, "<parameter name=\"%s\" value=\"%s\"/>\n",
              key, (char*) get_hash_value(parameters, key));
    }
  }
  free_hash_iterator(iterator);
}


/**
 * \brief Creates a file containing all parameters and their current
 * values in the parameter file format. Created in the output directory
 * named by the parameter "output-dir".
 */
void print_parameter_file(char** filename){

  carp(CARP_DEBUG, "Printing parameter file");
  prefix_fileroot_to_name(filename);
  char* output_dir = get_string_parameter("output-dir");
  bool overwrite = get_boolean_parameter("overwrite");
  FILE* param_file = create_file_in_path(*filename, 
                                         output_dir, 
                                         overwrite);

  // Add header to file for comet parsing
  fprintf(param_file, "# comet_version 2014.01 rev. 0"
          "\n# Comet MS/MS search engine parameters file."
          "\n# Everything following the \'#\' symbol is treated as a comment.\n");

  // iterate over all parameters and print to file
  HASH_ITERATOR_T* iterator = new_hash_iterator(parameters);
  while(hash_iterator_has_next(iterator)){
    string key = hash_iterator_next(iterator);
    string show_users = (char*)get_hash_value(for_users, key.c_str());
    if( show_users == "true") {
      // print mods separately at the end
      if( key == "mod"  || key == "cmod" || key == "nmod"){ 
        continue;
      }
      // print comet parameters after these.
      if (key.find("_") == string::npos) {
        char buffer[PARAMETER_BUFFER] = "";
        strcat_formatted(buffer, "# ", (char*)get_hash_value(usages, key.c_str()));
        strcat_formatted(buffer, "# ", (char*)get_hash_value(file_notes, key.c_str()));
        fprintf(param_file, "%s%s=%s\n\n", buffer, key.c_str(), 
          (char*)get_hash_value(parameters, key.c_str()));
      }
    }
  }

  // now print all mods
  print_mods_parameter_file(param_file, "mod", get_aa_mod_list);
  print_mods_parameter_file(param_file, "nmod", get_n_mod_list);
  print_mods_parameter_file(param_file, "cmod", get_c_mod_list);
  
  free_hash_iterator(iterator);

  //now print out Comet parameters
  fprintf(param_file, "#################\n");
  fprintf(param_file, "#Comet Parameters\n");
  fprintf(param_file, "#################\n");

  iterator = new_hash_iterator(parameters);
  while(hash_iterator_has_next(iterator)){
    string key = hash_iterator_next(iterator);
    string show_users = (char*)get_hash_value(for_users, key.c_str());
    if( show_users == "true") {
      // print mods separately at the end                                                                                                                                                                                                                                       
      if( key == "mod"  || key == "cmod" || key == "nmod"){
        continue;
      }
      // print comet parameters after these.                                                                                                                                                                                                                                    
      if (key.find("_") != string::npos) {
        char buffer[PARAMETER_BUFFER] = "";
        strcat_formatted(buffer, "# ", (char*)get_hash_value(usages, key.c_str()));
        strcat_formatted(buffer, "# ", (char*)get_hash_value(file_notes, key.c_str()));
        fprintf(param_file, "%s%s=%s\n\n", buffer, key.c_str(),
          (char*)get_hash_value(parameters, key.c_str()));
      } 
    }
  }
  free_hash_iterator(iterator);


  // now print out Comet enzyme information
      
  fprintf(param_file, "#\n");
  fprintf(param_file, "# COMET_ENZYME_INFO _must_ be at the end of this parameters file\n");
  fprintf(param_file, "#\n");
  fprintf(param_file, "[COMET_ENZYME_INFO]\n");
  
  if (comet_enzyme_info_lines_.size() == 0) {
    fprintf(param_file, "0.  No_enzyme\t\t\t\t");
    fprintf(param_file, "0");
    fprintf(param_file, "       -           -\n");
  
    fprintf(param_file, "1.  Trypsin\t\t\t\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      KR           P\n");
  
    fprintf(param_file, "2.  Trypsin/P\t\t\t\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      KR           -\n");
  
    fprintf(param_file, "3.  Lys_C\t\t\t\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      K            P\n");
  
    fprintf(param_file, "4.  Lys_N\t\t\t\t");
    fprintf(param_file, "0");
    fprintf(param_file, "      K            -\n");
    
    fprintf(param_file, "5.  Arg_C\t\t\t\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      R            P\n");
    
    fprintf(param_file, "6.  Asp_N\t\t\t\t");
    fprintf(param_file, "0");
    fprintf(param_file, "      D            -\n");
    
    fprintf(param_file, "7.  CNBr\t\t\t\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      M            -\n");
    
    fprintf(param_file, "8.  Glu_C\t\t\t\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      DE           P\n");
    
    fprintf(param_file, "9.  PepsinA\t\t\t\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      FL           P\n");
    
    fprintf(param_file, "10. Chymotrypsin\t\t\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      FWYL         P\n");
  
  /*
    TODO: Put these back in after we figure out what to do
    with enzyme info.
    fprintf(param_file, "11. Elastase \t\t\t\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      ALIV         P\n");
    
    fprintf(param_file, "12. Clostripai\t\t\t\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      R            -\n");
    
    fprintf(param_file, "13. Iodosobenzoate\t\t\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      W            -\n");
    
    fprintf(param_file, "14. Proline_Endopeptidase\t\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      P            -\n");
    
    fprintf(param_file, "15. Staph_Protease\t\t\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      E            -\n");
    
    fprintf(param_file, "16. Modified_Chymotrypsin\t\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      FWYL         P\n");
    
    
    fprintf(param_file, "17. Elastase_Trypisn_Chymotrypsin\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      ALIVKRWFY    P\n");
    
  */

  } else {
    for (int idx=0;idx < comet_enzyme_info_lines_.size();idx++) {
      fprintf(param_file, "%s\n", comet_enzyme_info_lines_[idx].c_str());
    }
  }

  fclose(param_file);
  free(output_dir);
}

/**
 * free heap allocated parameters
 */
void free_parameters(void){
  if(parameter_initialized){
    free_hash(parameters);
    free_hash(usages);
    free_hash(file_notes);
    free_hash(for_users);
    free_hash(types);
    free_hash(min_values);
    free_hash(max_values);

    for (int mod_idx=0; mod_idx < MAX_AA_MODS; mod_idx++) {
      free_aa_mod(list_of_mods[mod_idx]);
      list_of_mods[mod_idx] = NULL;
    }

  }
  // this is mostly so I can test repeatedly
  num_mods = 0;
  num_c_mods = 0;
  num_n_mods = 0;
  list_of_c_mods = NULL;
  list_of_n_mods = NULL;
  parameter_initialized = false;
  parameter_plasticity = true;
}

/**
 *
 * parse the parameter file given the filename
 */
void parse_parameter_file(
  char* parameter_filename ///< the parameter file to be parsed -in
  )
{
  FILE *file;
  char *line;
  int idx;

  carp(CARP_DETAILED_DEBUG, "Parsing parameter file '%s'",parameter_filename);

  /* check if parameters can be changed */
  if(!parameter_plasticity){
    carp(CARP_FATAL, "Can't change parameters once they are confirmed");
  }

  /* check if parameter file exists, if not die */
  if(access(parameter_filename, F_OK)){
    carp(CARP_FATAL, "Could not open parameter file '%s'", parameter_filename);
  }

  LineFileReader line_reader(parameter_filename);
  
  bool found_comet = false;

  while(line_reader.hasNext()) {
    string line = line_reader.next();
    string option_name = "";
    string option_value = "";
    bool found_equal = false;
    if (found_comet) {
      comet_enzyme_info_lines_.push_back(line);
    } else {
      if (line.find("[COMET_ENZYME_INFO]") != string::npos) {
        comet_enzyme_info_lines_.clear();
        found_comet = true;
      } else {
        for (string::iterator iter = line.begin(); iter != line.end(); ++iter) {

          //ignore everything after commet
          if (*iter == '#') {
            break;
          } else if (!found_equal) {
            if (*iter != '=') {
              option_name += *iter;
            } else {
              //okay we found an equal sign, rest of text should be the parameter value.
              found_equal = true;
            }
          } else {
            option_value += *iter;
          }
        }

        if (found_equal) {
          //trim the spaces on both sides of the name and parameter and update the hash.
          option_name = trim(option_name);
          option_value = trim(option_value);
          if(! update_hash_value(parameters, option_name.c_str(), option_value.c_str()) ){
            carp(CARP_WARNING, "Unexpected parameter file option '%s'", option_name.c_str());
          } else {
            check_option_type_and_bounds(option_name.c_str());
          }
        } else {
          carp(CARP_DEBUG, "equal not found:", line.c_str());
        }
      }
    }
  }
  if (comet_enzyme_info_lines_.size() == 0) {
    carp(CARP_WARNING, "putting in default comet enzyme lines");
    comet_enzyme_info_lines_.push_back("0.  No_enzyme\t\t\t\t0       -           -");
    comet_enzyme_info_lines_.push_back("1.  Trypsin\t\t\t\t1      KR           P");
    comet_enzyme_info_lines_.push_back("2.  Trypsin/P\t\t\t\t1      KR           -");
    comet_enzyme_info_lines_.push_back("3.  Lys_C\t\t\t\t1      K            P");
    comet_enzyme_info_lines_.push_back("4.  Lys_N\t\t\t\t0      K            -");
    comet_enzyme_info_lines_.push_back("5.  Arg_C\t\t\t\t1      R            P");
    comet_enzyme_info_lines_.push_back("6.  Asp_N\t\t\t\t0      D            -");
    comet_enzyme_info_lines_.push_back("7.  CNBr\t\t\t\t1      M            -");
    comet_enzyme_info_lines_.push_back("8.  Glu_C\t\t\t\t1      DE           P");
    comet_enzyme_info_lines_.push_back("9.  PepsinA\t\t\t\t1      FL           P");
    comet_enzyme_info_lines_.push_back("10. Chymotrypsin\t\t\t1      FWYL         P");
  }
  /*
    TODO: Put these back in after we figure out what to do
    with enzyme info.
    fprintf(param_file, "11. Elastase \t\t\t\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      ALIV         P\n");
    
    fprintf(param_file, "12. Clostripai\t\t\t\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      R            -\n");
    
    fprintf(param_file, "13. Iodosobenzoate\t\t\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      W            -\n");
    
    fprintf(param_file, "14. Proline_Endopeptidase\t\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      P            -\n");
    
    fprintf(param_file, "15. Staph_Protease\t\t\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      E            -\n");
    
    fprintf(param_file, "16. Modified_Chymotrypsin\t\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      FWYL         P\n");
    
    
    fprintf(param_file, "17. Elastase_Trypisn_Chymotrypsin\t");
    fprintf(param_file, "1");
    fprintf(param_file, "      ALIVKRWFY    P\n");
    
  */

}


/**************************************************
 *   GETTERS (public)
 **************************************************
 */

/**
 * Each of the following functions searches through the hash table of
 * parameters, looking for one whose name matches the string.  The
 * function returns the corresponding value.
 * \returns true if paramater value is true, else false
 */ 
bool get_boolean_parameter(
 const char*     name  ///< the name of the parameter looking for -in
 )
{
  static char buffer[PARAMETER_LENGTH];
  
  char* value = (char*)get_hash_value(parameters, name);
 
  // can't find parameter
  if(value == NULL){
    carp(CARP_FATAL, "Parameter name '%s' doesn't exist", name);
  }
  
  //check type
  char* type_str = (char*)get_hash_value(types, name);
  PARAMETER_TYPE_T type;
  bool found = string_to_param_type(type_str, &type);
 
  if(found == false || type != BOOLEAN_P){
    carp(CARP_ERROR, "Request for boolean parameter '%s' which is of type %s",
         name, type_str);
  }

 // make sure that there is enough storage allocated in the string
  if((int)strlen(value) 
     > PARAMETER_LENGTH) {
    carp(CARP_FATAL, "parameter %s with value %s was too long to copy to string ",
        name,
        value);
  }
  strncpy(buffer,
          value,
          PARAMETER_LENGTH);

  
  if ((strncmp(buffer,"TRUE", 4) == 0) || 
      (strncmp(buffer, "true", 4) == 0) || 
      (strncmp(buffer, "T", 1) == 0)){
    return(true);
  } 
  else if ((strncmp(buffer,"FALSE", 5) == 0) || 
           (strncmp(buffer, "false", 5) == 0) || 
           (strncmp(buffer, "F", 1) == 0)){

    return(false);
  } 
  else {
    carp(CARP_FATAL, "Invalid Boolean parameter %s. ", buffer);
  }
  
  carp(CARP_FATAL, "parameter name: %s, doesn't exist", name);
  return false; // Return value to avoid compiler warning
}

/**
 * Searches through the list of parameters, looking for one whose
 * name matches the string.  This function returns the parameter value if the
 * parameter is in the parameter hash table.  This
 * function exits if there is a conversion error. 
 *\returns the int value of the parameter
 */
int get_int_parameter(
  const char* name  ///< the name of the parameter looking for -in
  )
{
  //char *endptr;
  //long int value;
  int value;

  char* int_value = (char*)get_hash_value(parameters, name);

  //  carp(CARP_DETAILED_DEBUG, "int value string is %s", int_value);

  // can't find parameter
  if(int_value == NULL){
    carp(CARP_FATAL, "parameter name: %s, doesn't exist", name);
  }
  //check type
  char* type_str = (char*)get_hash_value(types, name);
  PARAMETER_TYPE_T type;
  bool found = string_to_param_type(type_str, &type);

  if(found==false || type != INT_P){
    carp(CARP_ERROR, "Request for int parameter '%s' which is of type %s",
         name, type_str);
  }

  /* there is a parameter with the right name.  Now 
     try to convert it to a base 10 integer*/
  value = atoi(int_value);
  return value;
}

vector<int> get_int_vector_parameter(
  const char* name
  ) {

  vector<int> ans;
  
  vector<string> sans = get_string_vector_parameter(name);

  for (unsigned int idx=0;idx<sans.size();idx++) {

    int ival;
    from_string(ival,sans[idx]);
    ans.push_back(ival);
  }


  return ans;
}

/**
 * Searches through the list of parameters, looking for one whose
 * name matches the string.  This function returns the parameter value if the
 * parameter is in the parameter hash table.  This
 * function exits if there is a conversion error. 
 *\returns the double value of the parameter
 */
double get_double_parameter(
  const char* name   ///< the name of the parameter looking for -in
  )
{
  char *endptr;
  double value;
  
  // check if parameter file has been parsed
  if(!parameter_initialized){
    carp(CARP_FATAL, "parameters have not been set yet");
  }

  char* double_value = (char*)get_hash_value(parameters, name);
 
  // can't find parameter
  if(double_value == NULL){
    carp(CARP_FATAL, "parameter name '%s', doesn't exit", name);
  }
 
  //check type
  char* type_str = (char*)get_hash_value(types, name);
  PARAMETER_TYPE_T type;
  bool found = string_to_param_type(type_str, &type);

  if(found==false || type != DOUBLE_P){
    carp(CARP_ERROR, "Request for double parameter '%s' which is of type %s",
         name, type_str);
  }
  
  /* there is a parameter with the right name.  Now 
     try to convert it to a double*/
  value = strtod(double_value, NULL);
 
  return(value);
  
  carp(CARP_FATAL, "parameter name: %s, doesn't exist", name);
}

vector<double> get_double_vector_parameter(
  const char* name
  ) {

  vector<double> ans;
  
  vector<string> sans = get_string_vector_parameter(name);

  for (unsigned int idx=0;idx<sans.size();idx++) {

    double dval;
    from_string(dval,sans[idx]);
    ans.push_back(dval);
  }


  return ans;
}

/**
 * \brief Get the value of a parameter whose type is char*
 *
 * Searches through the list of parameters, looking for one whose
 * parameter_name matches the given name string. 
 * The return value is allocated here and must be freed by the caller.
 * If the value is not found, abort.
 * \returns The hash-alllocated string value of the given parameter name
 */
char* get_string_parameter(
  const char* name  ///< the name of the parameter looking for -in
  )
{
  
  char* string_value = (char*)get_hash_value(parameters, name);
  
  // can't find parameter
  if(string_value == NULL){
    carp(CARP_FATAL, "Parameter name: %s, doesn't exist", name);
  }

  //change "__NULL_STR" to NULL
  if( (strcmp(string_value, "__NULL_STR"))==0){
    string_value = NULL;
  }
  //check type
  char* type_str = (char*)get_hash_value(types, name);
  PARAMETER_TYPE_T type;
  string_to_param_type(type_str, &type);

  return my_copy_string(string_value);
}

vector<string> get_string_vector_parameter(
  const char* name
  ) {

  vector<string> ans;

  string value(get_string_parameter_pointer(name));

  tokenize(value, ans, ',');

  return ans;
}

// TODO (BF 04-Feb-08): Should we delete this since it allows caller
//      to change the value of a parameter?
/**
 * \brief Get the value of a parameter whose type is char*
 * 
 * Searches through the list of parameters, looking for one whose
 * parameter_name matches the given name string. 
 * The return value is a pointer to the original string
 * Thus, user should not free, good for printing
 * \returns the string value to which matches the parameter name, else aborts
 */
const char* get_string_parameter_pointer(
  const char* name  ///< the name of the parameter looking for -in
  )
{
  
  char* string_value = (char*)get_hash_value(parameters, name);

  // can't find parameter
  if(string_value == NULL){
    carp(CARP_FATAL, "parameter name: %s, doesn't exist", name);
  }
  //check type
  char* type_str = (char*)get_hash_value(types, name);
  PARAMETER_TYPE_T type;
  string_to_param_type(type_str, &type);

  return string_value;

}

DIGEST_T get_digest_type_parameter( const char* name ){

  char* param = (char*)get_hash_value(parameters, name);

  DIGEST_T digest_type = string_to_digest_type(param);
  if( digest_type == INVALID_DIGEST ){
    carp(CARP_FATAL, "Digest_type parameter %s has the value %s " 
         "which is not of the correct type.", name, param);
  }
  return digest_type;
}

ENZYME_T get_enzyme_type_parameter( const char* name ){

  char* param = (char*)get_hash_value(parameters, name);

  ENZYME_T enzyme_type = string_to_enzyme_type(param);
  if( enzyme_type == INVALID_ENZYME ){
    carp(CARP_FATAL, "Enzyme_type parameter %s has the value %s " 
         "which is not of the correct type.", name, param);
  }
  return enzyme_type;
}

HARDKLOR_ALGORITHM_T get_hardklor_algorithm( const char* name ){

  char* param = (char*)get_hash_value(parameters, name);
  
  HARDKLOR_ALGORITHM_T hk_algorithm = 
    string_to_hardklor_algorithm_type(param);

  if ( hk_algorithm == INVALID_HK_ALGORITHM) {
    carp(CARP_FATAL, "Hardklor-algorithm parameter %s has "
      "the value %s which is not of the correct type.",name,param);
  }
  return hk_algorithm;
}

SPECTRUM_PARSER_T get_spectrum_parser_parameter( const char* name ) {

  char* param = (char*)get_hash_value(parameters, name);
  SPECTRUM_PARSER_T spectrum_parser = 
    string_to_spectrum_parser_type(param);
  
  if ( spectrum_parser == INVALID_SPECTRUM_PARSER) {
    carp(CARP_FATAL, "spectrum parser parameter %s has "
      "the value of %s which is not of the correct type.", name, param);
  }

  return spectrum_parser;


}



MASS_TYPE_T get_mass_type_parameter(
   const char* name
   ){
  char* param_value_str = (char*)get_hash_value(parameters, name);
  MASS_TYPE_T param_value;
  bool success = string_to_mass_type(param_value_str, &param_value);

  if( ! success ){
    carp(CARP_FATAL, 
         "Mass_type parameter %s has the value %s which is not of "
          "the correct type", name, param_value_str);
  }
  return param_value;
}

WINDOW_TYPE_T get_window_type_parameter(
  const char* name
  ){
  char* param_value_str = (char*)get_hash_value(parameters, name);
  WINDOW_TYPE_T param_value =  
    string_to_window_type(param_value_str);

  return param_value;
}

THRESHOLD_T get_threshold_type_parameter(
  const char* name
  ){
  char* param_value_str = (char*)get_hash_value(parameters, name);
  THRESHOLD_T param_value =  
    string_to_threshold_type(param_value_str);

  return param_value;
}


PARSIMONY_TYPE_T get_parsimony_type_parameter(
  const char* name
  ){
  char * param_value_str = (char*)get_hash_value(parameters, name);
  PARSIMONY_TYPE_T param_value =
    string_to_parsimony_type(param_value_str);
  
  return param_value;
}

QUANT_LEVEL_TYPE_T get_quant_level_type_parameter(
  const char* name
  ){
  char * param_value_str = (char*)get_hash_value(parameters, name);
  QUANT_LEVEL_TYPE_T param_value =
    string_to_quant_level_type(param_value_str);
  
  return param_value;
}


MEASURE_TYPE_T get_measure_type_parameter(
  const char* name
  ){
  char * param_value_str = (char*)get_hash_value(parameters, name);
  MEASURE_TYPE_T param_value =
    string_to_measure_type(param_value_str);
  
  return param_value;
}

DECOY_TYPE_T get_decoy_type_parameter(
  const char* name
  ){
  char * param_value_str = (char*)get_hash_value(parameters, name);
  DECOY_TYPE_T param_value =
    string_to_decoy_type(param_value_str);
  
  return param_value;
}

DECOY_TYPE_T get_tide_decoy_type_parameter(
  const char* name
) {
  char* param_value_str = (char*)get_hash_value(parameters, name);
  return string_to_tide_decoy_type(param_value_str);
}

MASS_FORMAT_T get_mass_format_type_parameter(
  const char* name
  ){
  char * param_value_str = (char*)get_hash_value(parameters, name);
  MASS_FORMAT_T param_value =
    string_to_mass_format(param_value_str);
  
  return param_value;
}

int get_max_ion_charge_parameter(
  const char* name
  ){
  char* param_value_str = (char*)get_hash_value(parameters, name);
  if (strcmp(param_value_str,"peptide") == 0) {
    return BILLION; //using this with min function on peptide charge.
  } else {
    int ans = atoi(param_value_str);
    if (ans <= 0 || ans > 6) {
      carp(CARP_FATAL,
        "Max_ion_charge parameter %s has the value %s which is not a "
        "legal value", name, param_value_str);
    }
    return ans;
  }
}

char get_delimiter_parameter(
  const char* name
  ) {
  
  char* param_value_str = (char*)get_hash_value(parameters, name);
  if (strcmp(param_value_str,"tab") == 0) {
    return '\t'; //using this with min function on peptide charge.
  } else {
    if (strlen(param_value_str) != 1) {
      carp(CARP_FATAL,
        "delimiter parameter %s with value %s is not a single character or "
        "'tab'", name, param_value_str);
    }
    return param_value_str[0];
  }
}

ALGORITHM_TYPE_T get_algorithm_type_parameter(const char* name){
  char* param_value_str = (char*)get_hash_value(parameters, name);
  ALGORITHM_TYPE_T param_value;
  bool success = string_to_algorithm_type(param_value_str, &param_value);

  if(!success){
    carp(CARP_FATAL, "Algorithm_type parameter %s has the value %s "
         "which is not of the correct type.", name, param_value_str);
  }
  return param_value;
}


SCORER_TYPE_T get_scorer_type_parameter(const char* name){
  char* param_value_str = (char*)get_hash_value(parameters, name);
  SCORER_TYPE_T param_value;
  bool success = string_to_scorer_type(param_value_str, &param_value);

  if(!success){
    carp(CARP_FATAL, "Scorer_type parameter %s has the value %s " 
         "which is not of the correct type.", name, param_value_str);
  }
  return param_value;
}

ION_TYPE_T get_ion_type_parameter(const char* name){
  char* param_value_str = (char*)get_hash_value(parameters, name);
  ION_TYPE_T param_value;
  bool success = string_to_ion_type(param_value_str, &param_value);

  if(!success){
    carp(CARP_FATAL, 
   "Ion_type parameter %s has the value %s which is not of the correct type.",
         name, param_value_str);
  }
  return param_value;
}

COLTYPE_T get_column_type_parameter(const char* name){
  char* param_value_str = (char*)get_hash_value(parameters, name);
  COLTYPE_T param_value = string_to_column_type(param_value_str);

  return param_value;
}

COMPARISON_T get_comparison_parameter(const char* name) {
  char* param_value_str = (char*)get_hash_value(parameters, name);
  COMPARISON_T param_value = string_to_comparison(param_value_str);

  return param_value;
}


/**************************************************
 *   SETTERS (private)
 **************************************************
 */

bool reset_parameter(const char* name, const char* value){
  return add_or_update_hash(parameters, name, value);
}

bool set_boolean_parameter(
 const char* name,       ///< the name of the parameter looking for -in
 bool   set_value,  ///< the value to be set -in
 const char* usage,      ///< message for the usage statement
 const char* filenotes,  ///< additional information for the params file
 const char* foruser     ///< "true" if should be revealed to user
)
{
  bool result;
    
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return false;
  }

  const char* bool_str;
  if(set_value){
    bool_str = "true";
  }
  else{
    bool_str = "false";
  }
  result = add_or_update_hash(parameters, name, bool_str);
  result = add_or_update_hash(usages, name, usage);
  result = add_or_update_hash(file_notes, name, filenotes);
  result = add_or_update_hash(for_users, name, foruser);
  result = add_or_update_hash(types, name, (void*)"bool");
  return result;
}

bool set_int_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 int set_value,  ///< the value to be set -in
 int min_value,  ///< the minimum accepted value -in
 int max_value,  ///< the maximum accepted value -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,  ///< additional info for param file
 const char* foruser     ///< true if should be revealed to user
  )
{
  bool result;
  char buffer[PARAMETER_LENGTH];
  
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return false;
  }
  
  //stringify default, min, and max values and set
  snprintf(buffer, PARAMETER_LENGTH, "%i", set_value);
  result = add_or_update_hash(parameters, name, buffer);

  snprintf(buffer, PARAMETER_LENGTH, "%i", min_value);
  result = add_or_update_hash(min_values, name, buffer);

  snprintf(buffer, PARAMETER_LENGTH, "%i", max_value);
  result = add_or_update_hash(max_values, name, buffer);

  result = add_or_update_hash(usages, name, usage);
  result = add_or_update_hash(file_notes, name, filenotes);
  result = add_or_update_hash(for_users, name, foruser);
  result = add_or_update_hash(types, name, (void*)"INT_ARG");
  return result;
}

bool set_double_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 double set_value,  ///< the value to be set -in
 double min_value,  ///< the value to be set -in
 double max_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,  ///< additional info for param file
 const char* foruser     ///< "true" if should be revealed to user
  )
{
  bool result;
  char buffer[PARAMETER_LENGTH];
  
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return false;
  }
  
  // convert to string
  snprintf(buffer, PARAMETER_LENGTH, "%f", set_value);
  result = add_or_update_hash(parameters, name, buffer);    

  snprintf(buffer, PARAMETER_LENGTH, "%f", min_value);
  result = add_or_update_hash(min_values, name, buffer);    

  snprintf(buffer, PARAMETER_LENGTH, "%f", max_value);
  result = add_or_update_hash(max_values, name, buffer);    

  result = add_or_update_hash(usages, name, usage);    
  result = add_or_update_hash(file_notes, name, filenotes);
  result = add_or_update_hash(for_users, name, foruser);
  result = add_or_update_hash(types, name, (void*)"DOUBLE_ARG");    
  return result;
}

/**
 * temporary replacement for function, return name once all exe's are fixed
 * \returns true if paramater value is set, else false
 */ 
bool set_string_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 const char* set_value,  ///< the value to be set -in
 const char* usage,
 const char* filenotes,  ///< additional comments for parameter file
 const char* foruser
  )
{
  bool result;

  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return false;
  }
  
  if( set_value == NULL ){
    set_value = "__NULL_STR";
  }

  result = add_or_update_hash(parameters, name, set_value);
  result = add_or_update_hash(usages, name, usage);
  result = add_or_update_hash(file_notes, name, filenotes);
  result = add_or_update_hash(for_users, name, foruser);
  result = add_or_update_hash(types, name, (void*)"STRING_ARG");

  return result;
}

bool set_mass_type_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 MASS_TYPE_T set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,  ///< additional info for param file
 const char* foruser
 )
{
  bool result;
  char value_str[265] ;

  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return false;
  }
  
  /* stringify the value */
  mass_type_to_string(set_value, value_str);
  
  result = add_or_update_hash(parameters, name, value_str);
  result = add_or_update_hash(usages, name, usage);
  result = add_or_update_hash(file_notes, name, filenotes);
  result = add_or_update_hash(for_users, name, foruser);
  result = add_or_update_hash(types, name, (void*)"MASS_TYPE_T");
  return result;

}

bool set_digest_type_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 DIGEST_T set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser
  )
{
  bool result = true;
  
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return false;
  }
  
  /* stringify the value */
  char* value_str = digest_type_to_string(set_value);
  carp(CARP_DETAILED_DEBUG, "Setting digest param '%s' to value '%s'.", name, value_str);

  result = add_or_update_hash(parameters, name, value_str);
  result = add_or_update_hash(usages, name, usage);
  result = add_or_update_hash(file_notes, name, filenotes);
  result = add_or_update_hash(for_users, name, foruser);
  result = add_or_update_hash(types, name, (void*)"DIGEST_T");
  free(value_str);
  return result;

}

bool set_enzyme_type_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 ENZYME_T set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser
  )
{
  bool result = true;
  
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return false;
  }
  
  /* stringify the value */
  char* value_str = enzyme_type_to_string(set_value);

  result = add_or_update_hash(parameters, name, value_str);
  result = add_or_update_hash(usages, name, usage);
  result = add_or_update_hash(file_notes, name, filenotes);
  result = add_or_update_hash(for_users, name, foruser);
  result = add_or_update_hash(types, name, (void*)"ENZYME_T");
  free(value_str);
  return result;

}

bool set_window_type_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 WINDOW_TYPE_T set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser
  ) {
  bool result = true;
  
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return false;
  }
  
  /* stringify the value */
  char* value_str = window_type_to_string(set_value);

  result = add_or_update_hash(parameters, name, value_str);
  result = add_or_update_hash(usages, name, usage);
  result = add_or_update_hash(file_notes, name, filenotes);
  result = add_or_update_hash(for_users, name, foruser);
  result = add_or_update_hash(types, name, "WINDOW_TYPE_T");
  free(value_str);
  return result;

}

bool set_threshold_type_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 THRESHOLD_T set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser
  ) {
  bool result = true;
  
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return false;
  }
  
  /* stringify the value */
  char* value_str = threshold_type_to_string(set_value);

  result = add_or_update_hash(parameters, name, value_str);
  result = add_or_update_hash(usages, name, usage);
  result = add_or_update_hash(file_notes, name, filenotes);
  result = add_or_update_hash(for_users, name, foruser);
  result = add_or_update_hash(types, name, "THRESHOLD_T");
  free(value_str);
  return result;

}


bool set_measure_type_parameter(
  const char* name, ///< the name of the parameter looking for -in
  MEASURE_TYPE_T set_value, ///< the value to be set -in
  const char* usage, ///< string to print in usage statement
  const char* filenotes, ///<additional infor for param file
  const char* foruser 
  ){
  bool result = true;

  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return false;
  }
  
  /* stringify the value */
  char* value_str = measure_type_to_string(set_value);
  
  result = add_or_update_hash(parameters, name, value_str);
  result = add_or_update_hash(usages, name, usage);
  result = add_or_update_hash(file_notes, name, filenotes);
  result = add_or_update_hash(for_users, name, foruser);
  result = add_or_update_hash(types, name, "MEASURE_TYPE_T");
  free(value_str);
  return result;
  
}

bool set_decoy_type_parameter(
  const char* name, ///< the name of the parameter looking for -in
  DECOY_TYPE_T set_value, ///< the value to be set -in
  const char* usage, ///< string to print in usage statement
  const char* filenotes, ///<additional infor for param file
  const char* foruser 
  ){
  bool result = true;

  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return false;
  }
  
  /* stringify the value */
  char* value_str = decoy_type_to_string(set_value);
  
  result = add_or_update_hash(parameters, name, value_str);
  result = add_or_update_hash(usages, name, usage);
  result = add_or_update_hash(file_notes, name, filenotes);
  result = add_or_update_hash(for_users, name, foruser);
  result = add_or_update_hash(types, name, "DECOY_TYPE_T");
  free(value_str);
  return result;
  
}

bool set_mass_format_type_parameter(
  const char* name, ///< the name of the parameter looking for -in
  MASS_FORMAT_T set_value, ///< the value to be set -in
  const char* usage, ///< string to print in usage statement
  const char* filenotes, ///<additional infor for param file
  const char* foruser 
  ){
  bool result = true;

  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return false;
  }
  
  /* stringify the value */
  char* value_str = mass_format_type_to_string(set_value);
  
  result = add_or_update_hash(parameters, name, value_str);
  result = add_or_update_hash(usages, name, usage);
  result = add_or_update_hash(file_notes, name, filenotes);
  result = add_or_update_hash(for_users, name, foruser);
  result = add_or_update_hash(types, name, "MASS_FORMAT_T");
  free(value_str);
  return result;
  
}

bool set_parsimony_type_parameter(
  const char* name, ///< the name of the parameter looking for -in
  PARSIMONY_TYPE_T set_value, ///< the value to be set -in
  const char* usage, ///< string to print in usage statement
  const char* filenotes, ///<additional infor for param file
  const char* foruser 
  ){
  bool result = true;

  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return false;
  }
  
  /* stringify the value */
  char* value_str = parsimony_type_to_string(set_value);
  
  result = add_or_update_hash(parameters, name, value_str);
  result = add_or_update_hash(usages, name, usage);
  result = add_or_update_hash(file_notes, name, filenotes);
  result = add_or_update_hash(for_users, name, foruser);
  result = add_or_update_hash(types, name, "PARSIMONY_TYPE_T");
  free(value_str);
  return result;
  
}

bool set_quant_level_parameter(
  const char* name, ///< the name of the parameter looking for -in
  QUANT_LEVEL_TYPE_T set_value, ///< the value to be set -in
  const char* usage, ///< string to print in usage statement
  const char* filenotes, ///<additional infor for param file
  const char* foruser 
  ){
  bool result = true;

  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return false;
  }
  
  /* stringify the value */
  char* value_str = quant_level_type_to_string(set_value);
  
  result = add_or_update_hash(parameters, name, value_str);
  result = add_or_update_hash(usages, name, usage);
  result = add_or_update_hash(file_notes, name, filenotes);
  result = add_or_update_hash(for_users, name, foruser);
  result = add_or_update_hash(types, name, "QUANT_LEVEL_TYPE_T");
  free(value_str);
  return result;
  
}

bool set_algorithm_type_parameter(
 const char* name,
 ALGORITHM_TYPE_T set_value,
 const char* usage,
 const char* filenotes,
 const char* foruser)
{
  bool result = true;
  char value_str[SMALL_BUFFER];
  
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return false;
  }
  /* stringify value */
  algorithm_type_to_string(set_value, value_str);
  carp(CARP_DETAILED_DEBUG, "setting algorithm type to %s", value_str);  

  result = add_or_update_hash(parameters, name, value_str);
  result = add_or_update_hash(usages, name, usage);
  result = add_or_update_hash(file_notes, name, filenotes);
  result = add_or_update_hash(for_users, name, foruser);
  result = add_or_update_hash(types, name, (void*)"ALGORITHM_TYPE_T");
  return result;
}

bool set_hardklor_algorithm_type_parameter(
  const char* name,
  HARDKLOR_ALGORITHM_T set_value,
  const char* usage,
  const char* filenotes,
  const char* foruser) {

  bool result = true;
  
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return false;
  }
  /* stringify value */
  char* value_str = hardklor_algorithm_type_to_string(set_value);
  carp(CARP_DETAILED_DEBUG, "setting algorithm type to %s", value_str);  

  result = add_or_update_hash(parameters, name, value_str);
  result = add_or_update_hash(usages, name, usage);
  result = add_or_update_hash(file_notes, name, filenotes);
  result = add_or_update_hash(for_users, name, foruser);
  result = add_or_update_hash(types, name, (void*)"HARDKLOR_ALGORITHM_TYPE_T");
  free(value_str);
  return result;
  
}

bool set_spectrum_parser_parameter(
  const char* name,
  SPECTRUM_PARSER_T set_value,
  const char* usage,
  const char* filenotes,
  const char* foruser) {

  bool result = true;
  char value_str[SMALL_BUFFER];

  // check if parameters can be changed
  if (!parameter_plasticity) {
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return false;
  }
  /* stringify value */
  strcpy(value_str, spectrum_parser_type_to_string(set_value));
  carp(CARP_DETAILED_DEBUG, "setting spectrum_parser type to %s", value_str);

  result = add_or_update_hash(parameters, name, value_str);
  result &= add_or_update_hash(usages, name, usage);
  result &= add_or_update_hash(file_notes, name, filenotes);
  result &= add_or_update_hash(for_users, name, foruser);
  result &= add_or_update_hash(types, name, (void*)"SPECTRUM_PARSER_T");

  return result;

}


bool set_scorer_type_parameter(
 const char* name,
 SCORER_TYPE_T set_value,
 const char* usage, 
 const char* filenotes,
 const char* foruser)
{
  bool result = true;
  char value_str[SMALL_BUFFER];
  
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return false;
  }
  /* stringify value */
  strcpy(value_str, scorer_type_to_string(set_value));
  carp(CARP_DETAILED_DEBUG, "setting score type to %s", value_str);  

  result = add_or_update_hash(parameters, name, value_str);
  result = add_or_update_hash(usages, name, usage);
  result = add_or_update_hash(file_notes, name, filenotes);
  result = add_or_update_hash(for_users, name, foruser);
  result = add_or_update_hash(types, name, (void*)"SCORER_TYPE_T");


  return result;
}

bool set_ion_type_parameter(
 const char* name,
 ION_TYPE_T set_value,
 const char* usage,
 const char* filenotes,
 const char* foruser)
{
  bool result = true;
  char value_str[SMALL_BUFFER];

  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return false;
  }
  
  /* stringify value */
  ion_type_to_string(set_value, value_str);

  result = add_or_update_hash(parameters, name, value_str);
  result = add_or_update_hash(usages, name, usage);
  result = add_or_update_hash(file_notes, name, filenotes);
  result = add_or_update_hash(for_users, name, foruser);
  result = add_or_update_hash(types, name, (void*)"ION_TYPE_T");
  return result;
}
/**
 * Routines that return crux enumerated types. 
 */

bool string_to_param_type(const char* name, PARAMETER_TYPE_T* result ){
  bool success = true;
  if( name == NULL ){
    return false;
  }

  int param_type = convert_enum_type_str(
                   name, -10, parameter_type_strings, NUMBER_PARAMETER_TYPES);
  (*result) = (PARAMETER_TYPE_T)param_type;

  if( param_type == -10 ){
    success = false;
  }
  return success;
}

/*
 * Applies any static mods to the aa masses
 */
bool update_aa_masses(){
  bool success = true;
  int aa;
  char aa_str[2];
  aa_str[1] = '\0';
  int alphabet_size = (int)'A' + ((int)'Z'-(int)'A'); //get_alpha_size(ALL_SIZE);
  carp(CARP_DETAILED_DEBUG, "updating masses, last is %d", alphabet_size);

  for(aa=(int)'A'; aa< alphabet_size -1; aa++){
    aa_str[0] = (char)aa;
    double delta_mass = get_double_parameter(aa_str);
    carp(CARP_DETAILED_DEBUG, "aa: %i, aa_str: %s, mass: %f", aa, aa_str, delta_mass);
    increase_amino_acid_mass((char)aa, delta_mass);
  }
  
  return success;
}

/**
 * \brief Get the pointer to the list of AA_MODs requested by the
 * user.  Does NOT include the c- and n-term mods
 * \returns The number of items pointed to by mods
 */
int get_aa_mod_list
  (AA_MOD_T*** mods) ///< the address of an array of pointers
{
  *mods = list_of_variable_mods;
  return num_mods;

}

/**
 * \brief Get the pointer to the list of AA_MODs for the peptide
 * c-terminus.  Return 0 and set mods==NULL if there are no c-term
 * mods.
 *
 * \returns The number of items pointed to by mods
 */
int get_c_mod_list
  (AA_MOD_T*** mods) ///< the address of an array of pointers
{
  *mods = list_of_c_mods;
  return num_c_mods;
}
/**
 * \brief Get the pointer to the list of AA_MODs for the peptide
 * n-terminus.  Return 0 and set mods==NULL if there are no n-term
 * mods.
 *
 * \returns The number of items pointed to by mods
 */
int get_n_mod_list
  (AA_MOD_T*** mods) ///< the address of an array of pointers
{
  *mods = list_of_n_mods;
  return num_n_mods;
}

/**
 * \brief Get the pointer to the list of AA_MODs requested by the
 * user.  Includes aa_mods, c- and n-term mods.
 * \returns The number of items pointed to by mods
 */
int get_all_aa_mod_list
  (AA_MOD_T*** mods) ///< the address of an array of pointers
{
  *mods = list_of_mods;
  return num_mods + num_c_mods + num_n_mods + num_fixed_mods;
}


 
/**
 * \returns The index of the C_TERM or N_TERM fixed modification in
 * the global list of modifications.
 */
int get_fixed_mod_index(MOD_POSITION_T position){
  int index = -1;
  switch(position){
  case N_TERM:
    index = fixed_n_mod;
    break;
  case C_TERM:
    index = fixed_c_mod;
    break;
  case ANY_POSITION:
    carp(CARP_ERROR, "Getting non-terminal fixed mods not implemented.");
    break;
  }
  return index;
}

/**
 * \returns the number of fixed terminal modifications: 0, 1, or 2.
 */
int get_num_fixed_mods(){
  return num_fixed_mods;
}

/* Helper functions for read_mods_from_file */
/* TODO: these reads could be made more general as in
   read_double(&double, char*, char)   */

/**
 * \brief Set the mass_change field in an AA_MOD based on a line from
 * a parameter file.
 *
 * Assumes that the line points to a FLOAT_T followed by separator.
 * Converts the number and sets the appropriate field in the mod.
 * Dies with error if line does not point to a number.  Returns a
 * pointer to the character after the first instance of separator or
 * to the end of the line if separator does not appear.
 *
 * \returns A pointer to the next token in the line.
 */
char* read_mass_change(AA_MOD_T* mod, char* line, char separator,
                       int& max_precision){
  //carp(CARP_DEBUG, "token points to %s", line);

  aa_mod_set_mass_change(mod, atof(line));
  if( aa_mod_get_mass_change(mod) == 0){
    carp(CARP_FATAL, "The mass change is not valid for mod %s", line);
  }
  char* next = line;
  char* decimal = NULL;
  // read to end of line, colon, or whitespace
  while(*next != '\0' && *next != separator && *next != ' ' && *next != '\t'){
    if(*next == '.'){
      decimal = next;
    }
    next++;
  }
  // distance from decimal to separator is number of digits
  if( decimal != NULL ){
    int distance = next - decimal -1;
    if( distance > max_precision ){
      max_precision = distance;
    }
  }

  next++;  // point past the separator

  return next;
}
/**
 * \brief Read the line from the parameter file and set the bool vaues
 * in the mod's aa_list appropriately.
 *
 * Assumes that line points to a list of letters, a-Z, followed by
 * separator.  Fails if characters before separator are not a-Z.
 * Returns a pointer to the character after separator, or to the end
 * of the line if separator is not found.
 *
 * \returns A pointer to the next token in the line.
 */
char* set_aa_list(AA_MOD_T* mod, char* line, char separator){
  carp(CARP_DETAILED_DEBUG, "token points to %s", line);

  bool* aa_list = aa_mod_get_aa_list(mod);
  while( *line != '\0' && *line != ':'){
    char aa = toupper( *line );
    carp(CARP_DETAILED_DEBUG, "aa is %c", aa);

    if( aa < 'A' || aa > 'Z' ){
      carp(CARP_FATAL, "The letter '%c' in the aa list is invalid.", aa);
    }
    carp(CARP_DETAILED_DEBUG, "aa index is %d", aa - 'A');
    aa_list[aa - 'A'] = true;
    //mod->aa_list[aa - 'A'] = true;
    carp(CARP_DETAILED_DEBUG, "Set %c to true index %d", aa, (int)(aa-'A'));
    line++;
  }

  if( *line == separator ){
    line++;
  }
  return line;
}

/**
 * \brief Set the max_per_peptide field in the mod with the value
 * pointed to by line.
 *
 * Fails if line does not point to a valid integer.
 * \returns void
 */
char* read_max_per_peptide(AA_MOD_T* mod, char* line, char separator){
  //carp(CARP_DETAILED_DEBUG, "token points to %s", line);
  if( *line == '\0' ){
    carp(CARP_FATAL, "Missing maximum mods per peptide for mod %s", line);
  }

  aa_mod_set_max_per_peptide(mod, atoi(line));
  if( aa_mod_get_max_per_peptide(mod) == 0 ){
    carp(CARP_FATAL, "Maximum mods per peptide is invalid for mod %s", line);
  }

  char* next = line;
  while(*next != '\0' && *next != separator){
    next++;
  }
  if (*next != '\0') {
    next++;  // point past the separator
  }
  return next;
}

char* read_prevents_cleavage(AA_MOD_T* mod, char* line, char separator) {
  char* next = line;
  switch(*line) {
    case ':':
      next++;
    case '\0':
      carp(CARP_DEBUG, "No prevents_cleavage property found %s",line);
      break;
    case 'T':
    case 't':
      carp(CARP_DEBUG, "prevents_cleavage set to true %s",line);
      aa_mod_set_prevents_cleavage(mod, true);
      break;
    case 'F':
    case 'f':
      carp(CARP_DEBUG, "prevents_cleavage set to false %s", line);
      aa_mod_set_prevents_cleavage(mod, false);
      break;
  }

  while (*next != '\0' && *next != separator) {
    next++;
  }
  if (*next != '\0') {
    next++; //points past the separator
  }

  return next;
}

char* read_prevents_xlink(AA_MOD_T* mod, char* line, char separator) {
  char* next = line;
  switch(*line) {
    case ':':
      next++;
    case '\0':
      carp(CARP_DEBUG, "No prevents_xlink property found %s",line);
      break;
    case 'T':
    case 't':
      carp(CARP_DEBUG, "prevents_xlink set to true %s",line);
      aa_mod_set_prevents_xlink(mod, true);
      break;
    case 'F':
    case 'f':
      carp(CARP_DEBUG, "prevents_xlink set to false %s",line);
      aa_mod_set_prevents_xlink(mod, false);
      break;
  }

  while (*next != '\0' && *next != separator) {
    next++;
  }
  if (*next != '\0') {
    next++; //points past the separator
  }

  return next;
}



/**
 * \brief Set the max_distance field in the mod with the value
 * pointed to by line.
 *
 * Fails if line does not point to a valid integer.
 * \returns void
 */
void read_max_distance(AA_MOD_T* mod, char* line){
  int max_distance = atoi(line);

  // make sure it's 0 because the file says 0; if not, default to max
  if( max_distance == 0 && *line != '0' ){
    max_distance = MAX_PROTEIN_SEQ_LENGTH;
  }
  aa_mod_set_max_distance(mod, max_distance);

  return;
}

/**
 * \brief This is the private function, detailed reader used by
 * read_mods_from_file().
 *
 * Reads through whole file looking for lines that start with
 * line_tag.  Parses information into each of those lines putting it
 * into the list_of_aa_mods, beginning with index cur_index and not
 * exceeding MAX_AA_MODS.  Distinguishes between regular
 * aa_mods and c-term or n-term mods for which struct fields it
 * fills.
 * \returns Returns the index of the next mod in the list.
 */
int read_mods(
 FILE* param_file, ///< file from which to read mod info
 int cur_index,    ///< index of next mod to be entered
 const char* line_tag,///< text at beginning of mod line (eg mod=)
 MOD_POSITION_T position,///< type of mod (any, c-, n-term)
 int& max_precision///< most digits in mass change
){

  carp(CARP_DEBUG, "Reading mods for %d position", (int)position);
  char* line = (char*)mycalloc(MAX_LINE_LENGTH, sizeof(char));


  // read the whole file looking for mods
  while(fgets(line, MAX_LINE_LENGTH, param_file)==line){
    // read line until one starts with tag (mod=, cmod=, nmod=)
    if( 0 != strncmp(line, line_tag, strlen(line_tag)) ){
      continue;
    }

    // check bounds on index
    if( cur_index == MAX_AA_MODS ){
      carp(CARP_FATAL, "Too many modifications in parameter file, " \
           "%d maximum", MAX_AA_MODS);
    }
    AA_MOD_T* cur_mod = list_of_mods[cur_index];

    // prepare for reading line
    carp(CARP_DEBUG, "mod line: %s", line);
    char* token = line + strlen(line_tag);

    // check for default value "NO MODS" written to default.parameter
    if( strncmp(token, "NO MODS", strlen("NO MODS")) == 0 ){
      return cur_index;
    }

    // get the FLOAT_T and check for ok-ness
    token = read_mass_change(cur_mod, token, ':', max_precision);

    // fill in values for standard mods
    if( position == ANY_POSITION ){
      // read the aa list and set the values in mod
      token = set_aa_list(cur_mod, token, ':');

      // get max per peptide
      token = read_max_per_peptide(cur_mod, token, ':');

      // read whether this modification prevents cleavage (OPTIONAL)
      token = read_prevents_cleavage(cur_mod, token, ':');

      // read whether this modification prevents xlink (OPTIONAL)
      token = read_prevents_xlink(cur_mod, token, ':');


    }// fill in values for c- or n-mod
    else{
      // get the max distance
      read_max_distance(cur_mod, token);

      // set all bools to true
      int i = 0;
      bool* aa_list = aa_mod_get_aa_list(cur_mod);
      for(i=0; i<AA_LIST_LENGTH; i++){
        aa_list[i] = true;
      }
      // set type to c-/n-term and max to 1
      aa_mod_set_position(cur_mod, position);
      aa_mod_set_max_per_peptide(cur_mod, 1);
    }

    //  increment counter and get next mod
    cur_index++;

  }// repeat until end of file

  free(line);
  return cur_index;
}

/**
 * \brief Read the paramter file and populate the static parameter
 * list of AA_MODS, inlcuding the list of position mods.
 *
 * Also updates the array of amino_masses.  Dies with an error if the
 * number of mods in the parameter file is greater than MAX_AA_MODS.
 * \returns void
 */
void read_mods_from_file(char* param_filename){
  carp(CARP_DEBUG, "Reading mods from parameter file '%s'", param_filename);

  // open file
  FILE* param_file = fopen(param_filename, "rb");
  if( param_file == NULL ){
    carp(CARP_FATAL, "Could not open parameter file '%s'", param_filename);
  }

  // get first mod
  int total_num_mods = 0;
  int max_precision = 0;  // gets updated with max seen in param file

  // start with fixed terminal mods
  total_num_mods = read_mods(param_file, 0, "nmod-fixed=", 
                             N_TERM, max_precision);
  // keep track of where we stored the fixed mod
  if( total_num_mods == 1 ){
    fixed_n_mod = 0;  // first in list
  } else if (total_num_mods > 1){
    carp(CARP_FATAL, 
         "Cannot specify more than one fixed n-terminal modification.");
  }
  rewind( param_file );

  total_num_mods = read_mods(param_file, total_num_mods, "cmod-fixed=", 
                             C_TERM, max_precision);
  // keep track of where we stored the fixed mod
  if( total_num_mods == 1 && fixed_n_mod == -1 ){
    fixed_c_mod = 0;  // first in list
  } else if( total_num_mods == 2 ){
    fixed_c_mod = 1;  // second in list
  } else if (total_num_mods > 2){
    carp(CARP_FATAL, 
         "Cannot specify more than one fixed n-terminal modification.");
  }
  rewind( param_file );
  num_fixed_mods = total_num_mods;

  // now get the variable mods
  list_of_variable_mods = &list_of_mods[total_num_mods];
  total_num_mods = read_mods(param_file, total_num_mods,
                             "mod=", ANY_POSITION, max_precision);
  num_mods = total_num_mods - num_fixed_mods;  // set global var

  // Read the file again to get the cmods
  rewind( param_file );

  // set cmod pointer to next in array
  list_of_c_mods = &list_of_mods[total_num_mods];

  total_num_mods = read_mods(param_file, total_num_mods, "cmod=", 
                             C_TERM, max_precision);
  num_c_mods = total_num_mods - num_mods - num_fixed_mods;

  // if no cmods present, don't point to the list of mods
  if( num_c_mods == 0){
    list_of_c_mods = NULL;
  }

  // Read the file again to get the nmods
  rewind( param_file );

  // set nmod pointer to next in array
  list_of_n_mods = &list_of_mods[total_num_mods];

  total_num_mods = read_mods(param_file, total_num_mods, "nmod=", 
                             N_TERM, max_precision);
  num_n_mods = total_num_mods - num_mods - num_c_mods - num_fixed_mods;

  // if no nmods present, don't point to the list of mods
  if( num_n_mods == 0){
    list_of_n_mods = NULL;
  }

  // set the mod-precision option
  char val_str[16];
  sprintf(val_str, "%i", max_precision);
  update_hash_value(parameters, "mod-precision", val_str);


  // close file
  fclose(param_file);
  carp(CARP_DEBUG, "Finished reading mods file");
}

void incrementNumMods() {
    num_mods++;
}

// Secret functions used in testing

void force_set_aa_mod_list(AA_MOD_T** amod_list, int new_num_mods){
  int i=0;
  for(i=0; i<new_num_mods; i++){
    list_of_mods[i] = amod_list[i];
  }
  num_mods = new_num_mods;
}
void force_set_c_mod_list(AA_MOD_T** cmod_list, int new_num_mods){
  list_of_c_mods = cmod_list;
  num_c_mods = new_num_mods;
}

void force_set_n_mod_list(AA_MOD_T** nmod_list, int new_num_mods){
  list_of_n_mods = nmod_list;
  num_n_mods = new_num_mods;
}



/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

