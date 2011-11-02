/***********************************************************************//**
 * \file parameter.cpp
 * FILE: parameter.cpp
 * AUTHOR: written by Tobias Mann, CRUXified by Chris Park
 * Missed-cleavage conversion: Kha Nguyen
 * \brief General parameter handling utilities. MUST declare ALL
 * optional command parameters here inside initalialize_parameters.
 ****************************************************************************/

#include "crux-utils.h"
#include "parameter.h"


//TODO:  in all set, change result=add_... to result= result && add_...

/**
 * Starting location for zeroth m/z bin.
 */
static const FLOAT_T SMART_MZ_OFFSET = 0.68;

/*
 * Global variables
 */

static const char* parameter_type_strings[NUMBER_PARAMETER_TYPES] = { 
  "INT_ARG", "DOUBLE_ARG", "STRING_ARG", "MASS_TYPE_T", "DIGEST_T", 
  "ENZYME_T", 
  "bool", "SCORER_TYPE_T", "ION_TYPE_T",
  "ALGORITHM_TYPE_T", "WINDOW_TYPE_T", "MEASURE_TYPE_T", 
  "PARSIMONY_TYPE_T", "QUANT_LEVEL_TYPE_T", "DECOY_TYPE_T"};

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

bool set_algorithm_type_parameter(
 const char* name,
 ALGORITHM_TYPE_T set_value,
 const char* usage,
 const char* filenotes,
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
 
bool select_cmd_line(  
  const char** option_names, ///< list of options to be allowed for main -in
  int    num_options,  ///< number of optons in that list -in
  int (*parse_argument_set)(const char*, const char*, void*, enum argument_type) ///< function point to choose arguments or options 
  );

bool update_aa_masses();
void read_mods_from_file(char* param_file);

/************************************
 * Function definitions
 ************************************
 */


/**
 * The size of the bins for discretizing the m/z axis of the
 * observed spectrum.  For use with monoisotopic mass.
 */
static const FLOAT_T BIN_WIDTH_MONO = 1.0005079;

/**
 * The size of the bins for discretizing the m/z axis of the
 * observed spectrum.  For use with average mass.
 */
static const FLOAT_T  BIN_WIDTH_AVERAGE = 1.0011413;

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
  set_string_parameter("protein database", NULL, 
      "Fasta file of proteins or directory containing an index.",
      "Argument for generate, index, search, analyze.", "false");

  set_string_parameter("search results directory", NULL, 
      "Directory containing the results of one search.",
      "Argument for q-ranker, percolator, compute-q-values.", "false");

  /* create_index arguments */
  set_string_parameter("protein fasta file", NULL,
                       "File containing protein sequences in fasta format.",
                       "Argument for crux-create-index.", "false");
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
      "crux-create-index and crux-generate-peptides.  Parameter file "
      "only for crux-search-for-matches.", "true");
  set_int_parameter("max-length", 50, 1, MAX_PEPTIDE_LENGTH,
      "The maximum length of peptides to consider. Default=50.",
      "Available from command line or parameter file for crux-create-index "
      " and crux-generate-peptides. Parameter file only for crux-search-"
      "for-matches.", "true");
  set_double_parameter("min-mass", 200, 0, BILLION,
      "The minimum mass of peptides to consider. Default=200.",
      "Available from command line or parameter file for crux-create-index "
      "and crux-generate-peptides. Parameter file only for crux-search-"
      "for-matches.", "true");
  set_double_parameter("max-mass", 7200, 1, BILLION, 
      "The maximum mass of peptides to consider. Default=7200.",
      "Available from command line or parameter file for crux-create-index "
      "and crux-generate-peptides. Parameter file only for crux-search-"
      "for-matches.", "true");
  set_mass_type_parameter("isotopic-mass", AVERAGE, 
      "Which isotopes to use in calcuating peptide mass. "
      "<string>=average|mono. Default=average.", 
      "Used from command line or parameter file by crux-create-index and "
      "crux-generate-peptides.  Parameter file only for "
      "crux-search-for-matches.", "true");
  set_int_parameter("min-peaks", 20, 0, BILLION,
      "The minimum number of peaks a spectrum must have for it to be searched."
      " Default=20.", 
      "Parameter file only for search-for-matches and sequest-search.", "true");
  set_digest_type_parameter("digestion", FULL_DIGEST,
      "Degree of digestion used to generate peptides. "
      "<string>=full-digest|partial-digest. Either both ends or one end "
      "of a peptide must conform to enzyme specificity rules. "
      "Default=full-digest.",
      "Used in conjunction with enzyme option when enzyme is not set to "
      "to 'no-enzyme'.  Available from command line or parameter file for "
      "crux-generate-peptides and crux create-index.  Available from parameter"
      " file for crux search-for-matches.  Digestion rules are as "
      "follows: enzyme name [cuts after one of these residues][but not before "
      "one of these residues].  trypsin [RK][P], elastase [ALIV][P], "
      "chymotrypsin [FWY][P].",
      "true");
  set_enzyme_type_parameter("enzyme", TRYPSIN,
      "Enzyme to use for in silico digestion of proteins. "
      "<string>=trypsin|chymotrypsin|elastase|clostripain| "
      "cyanogen-bromide|iodosobenzoate|proline-endopeptidase| "
      "staph-protease|aspn|modified-chymotrypsin|no-enzyme. "
      "Default=trypsin.", 
      "Used in conjunction with the options digestion and missed-cleavages. "
      "Use 'no-enzyme' for non-specific digestion.  Available "
      "from command line or parameter file for crux-generate-peptides "
      "and crux create-index.  Available from parameter file for crux "
      "search-for-matches.   Digestion rules: enzyme name [cuts after one "
      "of these residues]|{but not before one of these residues}.  "
      "trypsin [RK]|{P}, elastase [ALIV]|{P}, chymotrypsin [FWY]|{P}, "
      "clostripain [R]|[], cyanogen-bromide [M]|[], "
      "iodosobenzoate [W]|[], proline-endopeptidase [P]|[], staph-protease "
      "[E]|[], modified-chymotrypsin [FWYL]|{P}, elastase-trypsin-chymotrypsin "
      "[ALIVKRWFY]|{P},aspn []|[D] (cuts before D).",
      "true");

  set_window_type_parameter("precursor-window-type", WINDOW_MASS,
      "Window type to use for selecting candidate "
      "peptides.  <string>=mass|mz|ppm. Default=mass.",
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
		    "Available from command line or parameter file for crux-create-index "
		    "and crux-generate-peptides.  Parameter file only for crux-search-"
		    "for-matches.  When used with enzyme=<trypsin|elastase|chymotrpysin> "
		    " includes peptides containing one or more potential cleavage sites.",
		    "true");	    

  set_boolean_parameter("unique-peptides", true,
      "Generate peptides only once, even if they appear in more "
      "than one protein (T,F).  Default=F.",
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
      "Available for search-for-matches.  Sp scoring is always done for "
      "sequest-search.", "true");
  set_boolean_parameter("compute-p-values", false, 
      "Compute p-values for the main score type. Default=F.",
      "Currently only implemented for XCORR.", "true");
  set_boolean_parameter("use-mstoolkit", false,
      "Use MSToolkit to parse spectra. Default=F.",
      "Available for crux-search-for-matches", "true");
  set_string_parameter("scan-number", NULL,
      "Search only select spectra specified as a single "
      "scan number or as a range as in x-y.  Default=search all.",
      "The search range x-y is inclusive of x and y.", "true");
  /* N.B. Use NaN to indicate that no user preference was specified.
   * In this case, the default value depends on the mass type. */
  set_double_parameter("mz-bin-width", NaN(), 0.0, BILLION,
      "Specify the width of the bins used to "
      "discretize the m/z axis.  Also used as tolerance for assigning "
      "ions.  Default=1.0005079 for monoisotopic mass "
      "or 1.0011413 for average mass.",
      "Available for crux-search-for-matches and xlink-assign-ions.", "true");
  set_double_parameter("mz-bin-offset", SMART_MZ_OFFSET, -1.0, 1.0,
      "Specify the location of the left edge of the "
      "first bin used to discretize the m/z axis. Default=0.68",
      "Available for crux-search-for-matches.", "true");
  // initialize as "unset", then set as bool after cmdline parsed
  set_string_parameter("use-flanking-peaks", "unset",
      "Include peaks +/- 1da around b/y ions in theoretical spectrum.  "
      "sequest-search and search-for-xlinks default=T. search-for-matches "
      "default=F.",
      "Available in the paramter file for all search commands.",
      "true");
  set_double_parameter("spectrum-min-mass", 0.0, 0, BILLION, 
      "Minimum mass of spectra to be searched. Default=0.",
      "Available for crux-search-for-matches.", "true");
  set_double_parameter("spectrum-max-mass", BILLION, 1, BILLION, 
      "Maximum mass of spectra to search. Default no maximum.",
      "Available for crux-search-for-matches.", "true");
  set_string_parameter("spectrum-charge", "all", 
      "Spectrum charge states to search. "
      "<string>=1|2|3|all. Default=all.",
      "Used by crux-search-for-matches to limit the charge states "
      "considered in the search.  With 'all' every spectrum will be "
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
      " Use 'none' for no decoys.  Default=protein-shuffle.",
      "For create-index, store the decoys in the index.  For search, either "
      "use decoys in the index or generate them from the fasta file.", "true");
  set_int_parameter("num-decoys-per-target", 1, 0, 10,
      "Number of decoy peptides to search for every target peptide searched."
      "Only valid for fasta searches when --decoys is not none. Default=0.",
      "Use --decoy-location to control where they are returned (which "
      "file(s)) and --decoys to control how targets are randomized.  Available "
      "for search-for-matches and sequest-search when searching a fasta file. ",
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
      "Available from parameter file for crux-search-for-matches.",
      "true");
  set_int_parameter("psms-per-spectrum-reported", 0, 0, BILLION,
                   "place holder", "this may be replaced by top-match","false");
  set_string_parameter("seed", "time", "HIDE ME FROM USER",
      "Given a real-number value, will always produce the same decoy seqs",
      "false");
  set_double_parameter("precursor-window", 3.0, 0, 100, 
      "Search peptides within +/- 'precursor-window' "
      "of the spectrum mass.  Definition of precursor window depends "
      "upon precursor-window-type. Default=3.0.",
      "Available from the parameter file only for crux-search-for-matches, "
      "crux-create-index, and crux-generate-peptides.",
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

  set_int_parameter("max-mods", MAX_PEPTIDE_LENGTH, 0, MAX_PEPTIDE_LENGTH,
      "The maximum number of modifications that can be applied to a single " 
      "peptide.  Default=no limit.",
      "Available from parameter file for crux-search-for-matches.", "true");
  set_int_parameter("max-aas-modified", MAX_PEPTIDE_LENGTH, 0,
      MAX_PEPTIDE_LENGTH,
      "The maximum number of modified amino acids that can appear in one "
      "peptide.  Each aa can be modified multiple times.  Default=no limit.",
      "Available from parameter file for search-for-matches.", "true");
  set_boolean_parameter("display-summed-mod-masses", true,
      "When a residue has multiple modifications, print the sum of those "
      "modifications rather than listing each in a comma-separated list.  "
      "Default=T.",
      "Available in the parameter file for any command that prints peptides "
      "sequences.  Example: true is SE[12.40]Q and false is SE[10.00,2.40]Q",
      "true" );
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
  set_int_parameter("print-search-progress", 10, 0, BILLION,
      "Show search progress by printing every n spectra searched.  Default="
      "10.", "Set to 0 to show no search progress.  Available for crux "
      "search-for-matches from parameter file.",
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

  // **** q-ranker options. ****
  set_boolean_parameter("no-xval", false, 
      "Turn off cross-validation to select hyperparameters.",
      "Available for q-ranker.", "true");

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
      "The ion series to predict (b,y,by). Default='by' (both b and y ions).",
      "Only available for crux-predict-peptide-ions.  Set automatically to "
      "'by' for searching.", "true");
  set_boolean_parameter("precursor-ions", false,
      "Predict the precursor ions, and all associated ions "
      "(neutral-losses, multiple charge states) consistent with the "
      "other specified options. (T,F) Default=F.",
      "Only available for crux-predict-peptide-ions.", "true");
  set_string_parameter("neutral-losses", "all", 
      "Predict neutral loss ions (none, h20, nh3, all). Default='all'.",
      "Only available for crux-predict-peptide-ions. Set to 'all' for "
      "sp and xcorr scoring.", "true");
  set_int_parameter("isotope", 0, 0, 2,
      "Predict the given number of isotope peaks (0|1|2). Default=0.",
      "Only available for crux-predict-peptide-ion.  Automatically set to "
      "0 for Sp scoring and 1 for xcorr scoring.", "true");
  set_boolean_parameter("flanking", false, 
      "Predict flanking peaks for b and y ions (T,F). Default=F.",
      "Only available for crux-predict-peptide-ion.", "true");
  set_string_parameter("max-ion-charge", "peptide",
      "Predict ions up to max charge state (1,2,...,6) or up to the charge state "
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
  // also uses "protein database"

  // ***** spectral-counts options *****
   set_string_parameter("input-ms2", NULL,
       "MS2 file corresponding to the psm file. Required for SIN.",
       "Available for spectral-counts with measure=SIN.",
       "true");
  set_double_parameter("threshold", 0.01, 0, 1, 
       "The p-value or q-value threshold. Default=0.01.",
       "Available for spectral-counts.  All PSMs with q-value higher than "
       "this will be ignored.",
       "true");
  set_measure_type_parameter("measure", MEASURE_NSAF,
       "Type of analysis to make on the match results: (NSAF|SIN|EMPAI). "
       "Default=NSAF. ", 
       "Available for spectral-counts.  NSAF is Normalized Spectral "
       "Abundance Factor, SIN is Spectral Index Normalized and EMPAI is "
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

  set_boolean_parameter("use-mgf", false,
      "Use MGF file format for parsing files",
      "Available for search-for-xlinks program (Default=F).",
      "true");

  set_boolean_parameter("use-old-xlink", true /* Turn to false later */,
      "Use old xlink searching algorihtm",
      "Available for search-for-xlinks program (Default=F).",
      "false");

  // **** xlink-score-spectrum options ****
  set_string_parameter("xlink-score-method", "composite", 
      "Score method for xlink {composite, modification, concatenated}. Default=composite.",
      "Argument for xlink-score-spectrum.", "false");

  // **** search-xlink options ****

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

  set_double_parameter("precursor-window-decoy", 20.0, 0, 1e6, 
      "Search decoy-peptides within +/- "
      " 'mass-window-decoy' of the spectrum mass.  Default=20.0.",
      "Available for crux search-for-xlinks. ",
      "true");

  set_window_type_parameter("precursor-window-type-decoy", WINDOW_MASS,
      "Window type to use for selecting "
      "decoy peptides from precursor mz. <string>=mass|mz|ppm. "
      "Default=mass.",
      "Available for crux search-for-matches",
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
  select_cmd_line( option_names, num_options, 
                   &parse_arguments_set_opt);
  return true;
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
  int (*parse_arguments_set_ptr)(const char*, const char*, void*, enum argument_type) 
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
    carp(CARP_DETAILED_DEBUG, "Option is: %s", option_names[i]);

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
                                      string_to_argument_type((char*)type_ptr)); 
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

  if (isnan(new_value)) {

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
 * Set the use-flanking-peaks parameter.  If the value was set by
 * user, set to that value.  Otherwise, set according to what command
 * is being run.  Defaults are F for search-for-matches and T for
 * others (sequest-search, search-for-xlinks).
 */
void set_flanking_peaks(const char* exe_name){

  const char* value = get_string_parameter_pointer("use-flanking-peaks");
  // if it is the default value, it was not set by the user
  if( strcmp(value, "unset") == 0 ){
    if( strcmp(exe_name, "search-for-matches") == 0 ){
      value = "false";
    } else {
      value = "true";
    }
  } else { // use the value set by the user
    if( value[0] == 'T' || value[0] == 't' ){
      value = "true";
    } else if( value[0] == 'F' || value[0] == 'f' ){
      value = "false";
    } //else don't change, let check find error
  }
  // set the new value and change type to bool
  update_hash_value(parameters, "use-flanking-peaks", value);
  update_hash_value(types, "use-flanking-peaks", (void*)"bool");
  check_option_type_and_bounds("use-flanking-peaks");
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
    free(usage);
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

  set_flanking_peaks(exe_name);


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

  double min_spec_mass = get_double_parameter("spectrum-min-mass");
  double max_spec_mass = get_double_parameter("spectrum-max-mass");

  if( min_spec_mass > max_spec_mass){
    carp(CARP_FATAL, "Parameter inconsistency. Minimum spectrum mass (%.2f)"
         " must be less than max (%.2f).", min_spec_mass, max_spec_mass);
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
      sprintf(die_str, "Illegal enzyme. Must be trypsin, chymotrypsin, "
              ", elastase, or no-enzyme.");
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
  case ION_TYPE_P:
    carp(CARP_DETAILED_DEBUG, "found ion_type param, value '%s'",
         value_str);
    if( !string_to_ion_type(value_str, &ion_type)){
      success = false;
      sprintf(die_str, "Illegal ion type '%s' for option '%s'.  "
              "Must be b,y,by.", value_str, name);
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
    //	 value_str);
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

  // TODO (BF Nov-12-08): could add header to file

  // iterate over all parameters and print to file
  HASH_ITERATOR_T* iterator = new_hash_iterator(parameters);
  while(hash_iterator_has_next(iterator)){
    char* key = hash_iterator_next(iterator);
    char* show_users = (char*)get_hash_value(for_users, key);
    if( strcmp(show_users, "true") == 0 ){
      // print mods separately at the end
      if( strcmp(key, "mod") == 0 || strcmp(key, "cmod") == 0
          || strcmp(key, "nmod") == 0 ){ 
        continue;
      }
      char buffer[PARAMETER_BUFFER] = "";
      strcat_formatted(buffer, "# ", (char*)get_hash_value(usages, key));
      strcat_formatted(buffer, "# ", (char*)get_hash_value(file_notes, key));
      fprintf(param_file, "%s%s=%s\n\n", buffer, key, 
              (char*)get_hash_value(parameters, key));
    }
  }

  // now print all mods
  print_mods_parameter_file(param_file, "mod", get_aa_mod_list);
  print_mods_parameter_file(param_file, "nmod", get_n_mod_list);
  print_mods_parameter_file(param_file, "cmod", get_c_mod_list);
  
  free_hash_iterator(iterator);
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
    carp(CARP_FATAL, "Could not open parameter file.");
  }

  file = fopen(parameter_filename, "r");
  if(file == NULL){
    carp(CARP_FATAL, "Couldn't open parameter file '%s'", parameter_filename);
  }

  line = (char*)mycalloc(MAX_LINE_LENGTH, sizeof(char));

  while(fgets(line, MAX_LINE_LENGTH, file)==line){

    idx = 0;
    
    /* Change the newline to a '\0' ignoring trailing whitespace */
    for(idx = MAX_LINE_LENGTH - 1; idx >= 0; idx--){
      if(line[idx] == '\n' || line[idx] == '\r' || 
         line[idx] == '\f' || line[idx] == ' ' || line[idx] == '\t')
        line[idx] = '\0';
    }
    /* empty lines and those beginning with '#' are ignored */
    if(line[0] != '#' && line[0] != '\0'){

      /* find the '=' in the line.  Exit with error if the line 
         has no equals sign. */
      idx = 0;
      while(idx < (int)strlen(line) && line[idx] != '='){
        idx++;
      }
      if(idx == 0 || idx >= (int)(strlen(line)-1)){
        carp(CARP_FATAL, "Lines in a parameter file must have the form: "
             "\n\tname=value\n "
             "In file %s, the line '%s' does not have this format",
             parameter_filename, line);
      }

      line[idx] = '\0';
      char* option_name = line;
      char* option_value = &(line[idx+1]);
      carp(CARP_DETAILED_DEBUG, "Found option '%s' and value '%s'", 
           option_name, option_value);

      if(! update_hash_value(parameters, option_name, option_value) ){
        carp(CARP_FATAL, "Unexpected parameter file option '%s'", option_name);
      }

      check_option_type_and_bounds(option_name);

    }
  }

  fclose(file);
  myfree(line);

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

  
  if ((strcmp(buffer,"TRUE") == 0) || 
      (strcmp(buffer, "true") == 0) || 
      (strcmp(buffer, "T") == 0)){
    return(true);
  } 
  else if ((strcmp(buffer,"FALSE") == 0) || 
           (strcmp(buffer, "false") == 0) || 
            (strcmp(buffer, "F") == 0)){

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
  value = strtod(double_value, &endptr);
 
  return(value);
  
  carp(CARP_FATAL, "parameter name: %s, doesn't exist", name);
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
  scorer_type_to_string(set_value, value_str);
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
  FILE* param_file = fopen(param_filename, "r");
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

