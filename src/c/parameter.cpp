/***********************************************************************//**
 * \file parameter.cpp
 * FILE: parameter.cpp
 * AUTHOR: written by Tobias Mann, CRUXified by Chris Park
 * CREATE DATE: 2006 Oct 09
 * DESCRIPTION: General parameter handling utilities. MUST declare ALL
 * optional command parameters here inside initalialize_parameters.
 ****************************************************************************/

#include "crux-utils.h"
#include "parameter.h"

//TODO:  in all set, change result=add_... to result= result && add_...

/*
 * Global variables
 */

static const char* parameter_type_strings[NUMBER_PARAMETER_TYPES] = { 
  "INT_ARG", "DOUBLE_ARG", "STRING_ARG", "MASS_TYPE_T", "DIGEST_T", 
  "ENZYME_T", //"PEPTIDE_TYPE_T", 
  "BOOLEAN_T", "SORT_TYPE_T", "SCORER_TYPE_T", "ION_TYPE_T",
  "ALGORITHM_TYPE_T"};

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
AA_MOD_T** list_of_c_mods = NULL; // pointer to first c_term mod in list
AA_MOD_T** list_of_n_mods = NULL; //pointer to first n_term mod in list
int num_mods = 0;
int num_c_mods = 0;
int num_n_mods = 0;//require num_mods + num_c_mods + num_n_mods <= MAX_AA_MODS

BOOLEAN_T parameter_initialized = FALSE; //have param values been initialized
BOOLEAN_T usage_initialized = FALSE; // have the usages been initialized?
BOOLEAN_T type_initialized = FALSE; // have the types been initialized?

BOOLEAN_T parameter_plasticity = TRUE; // can the parameters be changed?

char* pre_cleavage_list;
char* post_cleavage_list;
int pre_list_size;
int post_list_size;
BOOLEAN_T pre_for_inclusion;
BOOLEAN_T post_for_inclusion;

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
BOOLEAN_T check_option_type_and_bounds(const char* name);

void check_parameter_consistency();
void parse_custom_enzyme(char* rule_str);

/**
 *
 */
BOOLEAN_T string_to_param_type(const char*, PARAMETER_TYPE_T* );

BOOLEAN_T set_boolean_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 BOOLEAN_T set_value,  ///< the value to be set -in
 const char* usage,          ///< message for the usage statement
 const char* filenotes,      ///< additional information for the params file
 const char* foruser         ///< "true" if should be revealed to user
 );

BOOLEAN_T set_int_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 int set_value,  ///< the value to be set -in
 int min_value,  ///< the value to be set -in
 int max_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,  ///< additional info for param file
 const char* foruser     ///< "true" if should be revealed to user
 );

BOOLEAN_T set_double_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 double set_value,  ///< the value to be set -in
 double min_value,  ///< the value to be set -in
 double max_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,  ///< additional info for param file
 const char* foruser     ///< "true" if should be revealed to user
  );

BOOLEAN_T set_string_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 const char* set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser
  );

BOOLEAN_T set_mass_type_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 MASS_TYPE_T set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser
  );

BOOLEAN_T set_digest_type_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 DIGEST_T set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser
  );

BOOLEAN_T set_enzyme_type_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 ENZYME_T set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser
  );

BOOLEAN_T set_digest_type_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 DIGEST_T set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser
  );

BOOLEAN_T set_enzyme_type_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 ENZYME_T set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser
  );

BOOLEAN_T set_sort_type_parameter(
 const char* name,
 SORT_TYPE_T set_value,
 const char* usage,
 const char* filenotes,
 const char* foruser);

BOOLEAN_T set_algorithm_type_parameter(
 const char* name,
 ALGORITHM_TYPE_T set_value,
 const char* usage,
 const char* filenotes,
 const char* foruser);

BOOLEAN_T set_scorer_type_parameter(
 const char* name,
 SCORER_TYPE_T set_value,
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser);

BOOLEAN_T set_ion_type_parameter(
 const char* name,
 ION_TYPE_T set_value,
 const char* usage,
 const char* filenotes,
 const char* foruser);

BOOLEAN_T select_cmd_line(  
  const char** option_names, ///< list of options to be allowed for main -in
  int    num_options,  ///< number of optons in that list -in
  int (*parse_argument_set)(const char*, const char*, void*, enum argument_type) ///< function point to choose arguments or options 
  );

BOOLEAN_T update_aa_masses();
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
  pre_for_inclusion = TRUE;
  post_for_inclusion = FALSE;

  /* *** Initialize Arguments *** */

  // set with name, default value, [max, min], usage, notes, for param file
  // all arguments are left out of param file

  /* generate_peptide arguments */
  set_string_parameter("protein input", NULL, 
      "Fasta file of protein sequences or directory containing an index.",
      "Argument for generate, index, search, analyze.", "false");

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
  //and uses 'protein input'

  /* analyze-matches arguments */
  set_string_parameter("psm-folder", NULL, 
     "Directory containing the binary psm files created by "
     "crux-search-for-matches.",
     "Argument for analyze-matches.", "false");
  //and uses protein input

  /* get-ms2-spectrum */
  set_int_parameter("scan number", 0, 1, BILLION, 
                    "Scan number identifying the spectrum.",
                    "Argument for get-ms2-spectrum", "false");
  //uses ms2 file
  set_string_parameter("output file", NULL, 
                       "File where spectrum will be written.",
                       "Argument for get-ms2-spectrum.", "false");

  /* predict-peptide-ions */
  set_string_parameter("peptide sequence", NULL, 
      "The sequence of the peptide.",
      "Argument for predict-peptide-ions.", "false");
  set_int_parameter("charge state", 0, 0, 3, 
      "The charge state of the peptide.",
      "Argument for predict-peptide-ions", "false");

  /* create-psm-files */
  set_string_parameter("peptide-file-name", NULL,
      "A file containing peptides for which to create ion files",
      "Only for create-psm-files, which is not being distributed", "false"); 
  set_string_parameter("output-dir", NULL,
      "A directory in which to place the ion files",
      "Argument for create-psm-files", "false");
  set_string_parameter("model-type", NULL,
      "The kind of model (paired or single)",
      "Argument for create-psm-files", "false");

  /* *** Initialize Options (command line and param file) *** */

  /* options for all executables */
  set_boolean_parameter("version", FALSE, "Print version number and quit.",
      "Available for all crux programs.  On command line use '--version T'.",
      "true");
  set_int_parameter("verbosity", CARP_INFO, CARP_FATAL, CARP_MAX,
      "Set level of output to stderr (0-100).  Default 30.",
      "Available for all crux programs.  Each level prints the following "
      "messages, including all those at lower verbosity levels: 0-fatal "
      "errors, 10-non-fatal errors, 20-warnings, 30-information on the "
      "progress of execution, 40-more progress information, 50-debug info, "
      "60-detailed debug info.", "true");
  set_string_parameter("parameter-file", NULL, 
      "Set additional options with values in the given file. Default "
      "to use only command line options and default values.",
      "Available for all crux programs. Any options specified on the "
      "command line will override values in the parameter file.", "true");
  set_boolean_parameter("overwrite", FALSE, 
      "Replace existing files (T) or exit if attempting to overwrite "
      "(F). Default F.",
      "Available for all crux programs.  Applies to parameter file "
      "as well as index, search, and analysis output files.", "true");
    
  /* create-psm-files */
  set_int_parameter("starting-sentence-idx", 0, 0, BILLION, 
      "Starting sentence idx",
      "Only for create-psm-file, not distributed", "false"); 
  set_int_parameter("charge", 2, 1, 4,
      "Charge for peptide for which to predict ions.",
      "for create-psm-files and score-peptide-spectrum (neither in distro",
      "false"); 

  /* generate_peptide, create_index parameters  */
  set_int_parameter("min-length", 6, 1, MAX_PEPTIDE_LENGTH,
      "The minimum length of peptides to consider. Default 6.",
      "Used from the command line or parameter file by "
      "crux-create-index and crux-generate-peptides.  Parameter file "
      "only for crux-search-for-matches.", "true");
  set_int_parameter("max-length", 50, 1, MAX_PEPTIDE_LENGTH,
      "The maximum length of peptides to consider. Default 50.",
      "Available from command line or parameter file for crux-create-index "
      " and crux-generate-peptides. Parameter file only for crux-search-"
      "for-matches.", "true");
  set_double_parameter("min-mass", 200, 0, BILLION,
      "The minimum mass of peptides to consider. Default 200.",
      "Available from command line or parameter file for crux-create-index "
      "and crux-generate-peptides. Parameter file only for crux-search-"
      "for-matches.", "true");
  set_double_parameter("max-mass", 7200, 1, BILLION, 
      "The maximum mass of peptides to consider. Default 7200.",
      "Available from command line or parameter file for crux-create-index "
      "and crux-generate-peptides. Parameter file only for crux-search-"
      "for-matches.", "true");
  set_mass_type_parameter("isotopic-mass", AVERAGE, 
      "Which isotopes to use in calcuating peptide mass (average, mono). "
      "Default average.", 
      "Used from command line or parameter file by crux-create-index and "
      "crux-generate-peptides.  Parameter file only for "
      "crux-search-for-matches.", "true");
  set_digest_type_parameter("digestion", FULL_DIGEST,
      "Degree of digestion used to generate peptides (full-digest, "
      "partial-digest). Either both ends or one end of a peptide "
      "must conform to enzyme specificity rules. Default full-digest.",
      "Used in conjunction with enzyme option when enzyme is not set to "
      "to 'no-enzyme'.  Available from command line or parameter file for "
      "crux-generate-peptides and crux create-index.  Available from parameter"
      " file for crux search-for-matches.  Digestion rules are as "
      "follows: enzyme name [cuts after one of these residues][but not before "
      "one of these residues].  trypsin [RK][P], elastase [ALIV][P], "
      "chymotrypsin [FWY][P].",
      "true");
  set_enzyme_type_parameter("enzyme", TRYPSIN,
      "Enzyme to use for in silico digestion of protein sequences "
      "(trypsin, chymotrypsin, elastase, clostripain, cyanogen-bromide, "
      "iodosobenzoate, proline-endopeptidase, staph-protease, aspn, "
      "modified-chymotrypsin, no-enzyme).  Default trypsin.", 
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
  set_string_parameter("custom-enzyme", NULL, 
      "Specify rules for in silico digestion of protein sequences. See html "
      "docs for syntax. Default to use pre-defined enzyme trypsin.",
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
  set_boolean_parameter("missed-cleavages", FALSE, 
      "Include peptides with missed cleavage sites (T,F). Default FALSE.",
      "Available from command line or parameter file for crux-create-index "
      "and crux-generate-peptides.  Parameter file only for crux-search-"
      "for-matches.  When used with enzyme=<trypsin|elastase|chymotrpysin> "
      " includes peptides containing one or more potential cleavage sites.",
      "true");
  //  set_boolean_parameter("unique-peptides", FALSE,
  set_boolean_parameter("unique-peptides", TRUE,
      "Generate peptides only once, even if they appear in more "
      "than one protein (T,F).  Default FALSE.",
      "Available from command line or parameter file for "
      "crux-genereate-peptides. Returns one line per peptide "
      "when true or one line per peptide per protein occurence when false.  ",
      "true");
  
  /* more generate_peptide parameters */
  set_boolean_parameter("output-sequence", FALSE, 
      "Print peptide sequence (T,F). Default FALSE.",
      "Available only for crux-generate-peptides.", "true");
  set_sort_type_parameter("sort", SORT_NONE, 
      "Sort peptides according to which property "
      "(mass, length, lexical, none).  Default none.",
      "Only available for crux-generate-peptides.", "true");

  /* search-for-matches command line options */
  set_scorer_type_parameter("prelim-score-type", SP, 
      "Initial scoring (sp, xcorr). Default sp,", 
      "Available for crux-search-for-matches.  The score applied to all "
      "possible psms for a given spectrum.  Typically used to filter out "
      "the most plausible for further scoring. See max-rank-preliminary and "
      "score-type.", "false");
  set_scorer_type_parameter("score-type", XCORR, 
      "The primary scoring method to use (xcorr, sp, xcorr-pvalue, sp-pvalue)."
      " Default xcorr.", 
      "Only available for crux-search-for-matches.  Primary scoring is "
      "typically done on a subset (see max-rank-preliminary) of all "
      "possible psms for each spectrum. Default is the SEQUEST-style xcorr."
      " Crux also offers a p-value calculation for each psm based on xcorr "
      "or sp (xcorr-pvalue, sp-pvalue).", "false"); 
  set_boolean_parameter("compute-p-values", FALSE, 
      "Compute p-values for the main score type. Default F.",
      "Currently only implemented for XCORR.", "true");
  set_boolean_parameter("use-mz-window", FALSE,
      "Use mass-to-charge rather than mass for finding the window of peptides. Default F.",
      "Available for crux-search-for-matches", "true");

  set_boolean_parameter("use-mstoolkit", FALSE,
      "Use MSToolkit to parse spectra. Default F.",
      "Available for crux-search-for-matches", "false");
  set_string_parameter("scan-number", NULL,
      "Search only select spectra specified as a single scan "
      "number or as a range as in x-y.  Default search all.",
      "The search range x-y is inclusive of x and y.", "true");

  set_double_parameter("spectrum-min-mass", 0.0, 0, BILLION, 
      "Minimum mass of spectra to be searched.  Default 0.",
      "Available for crux-search-for-matches.", "true");
  set_double_parameter("spectrum-max-mass", BILLION, 1, BILLION, 
      "Maximum mass of spectra to search.  Default no maximum.",
      "Available for crux-search-for-matches.", "true");
  set_string_parameter("spectrum-charge", "all", 
      "Spectrum charge states to search (1,2,3,all). Default all.",
      "Used by crux-search-for-matches to limit the charge states "
      "considered in the search.  With 'all' every spectrum will be "
      "searched and spectra with multiple charge states will be searched "
      "once at each charge state.  With 1, 2 ,or 3 only spectra with that "
      "that charge will be searched.", "true");
  set_string_parameter("fileroot", NULL, 
      "Prefix added to output file names. Default None. ",
      "Used by crux create-index, crux search-for-matches, "
      "crux compute-q-values, and crux percolator.", "true");
  set_string_parameter("output-dir", "crux-output", 
      "Folder to which results will be written. Default 'crux-output'. ",
      "Used by crux create-index, crux search-for-matches, "
      "crux compute-q-values, and crux percolator.", "true");
  set_string_parameter("search-decoy-pvalue-file", "search.decoy.p.txt", 
      "Output filename for complete list of decoy p-values.  Default "
      "'search.decoy.p.txt'",
      "Only available for crux search-for-matches. The location of this "
      "file is controlled by --output-dir.", "true");
  // user options regarding decoys
  set_int_parameter("num-decoys-per-target", 2, 0, 10,
      "Number of decoy peptides to search for every target peptide searched."
      "  Default 2.",
      "Use --decoy-location to control where they are returned (which "
      "file(s)).  Available only from the command line for search-for-matches.",
      "true");
  set_string_parameter("decoy-location", "separate-decoy-files",
      "Where the decoy search results are returned.  (target-file|"
      "one-decoy-file|separate-decoy-files'. Default separate-decoy-files.",
      "Applies when num-decoys-per-target > 0.  Use 'target-file' to mix "
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
  set_boolean_parameter("tdc", FALSE,
      "Target-decoy competition. puts decoy psms in target file. ",
      "Now hidden from the user", "false");
  set_boolean_parameter("decoy-p-values", FALSE,
			"Store all decoy p-values in a file",
			"", "false");

  set_boolean_parameter("reverse-sequence", FALSE,
      "Generate decoys by reversing the peptide string rather than shuffling."
      " Default F.",
      "The first and last residues of the sequence are not changed.  If the "
      "reversed decoy sequence is the same as the target then the decoy is "
      "generated by shuffling.",
      "true");
  set_int_parameter("max-rank-preliminary", 500, 0, BILLION, 
      "Number of psms per spectrum to score with xcorr after preliminary "
      "scoring with Sp. "
      "Set to 0 to score all psms with xcorr. Default 500.",
      "Used by crux-search-for-matches.  For positive values, the Sp "
      "(preliminary) score acts as a filter; only high scoring psms go "
      "on to be scored with xcorr.  This saves some time.  If set to 0, "
      "all psms are scored with both scores. ", "true");
  set_int_parameter("top-match", 5, 1, BILLION, 
      "The number of psms per spectrum writen to the output file(s)." 
      "Default 5.",
      "Available from parameter file for crux-search-for-matches.",
      "true");
  set_int_parameter("psms-per-spectrum-reported", 0, 0, BILLION,
                   "place holder", "this may be replaced by top-match","false");
  set_double_parameter("mass-offset", 0.0, 0, 0, "obsolete",
      "Was used in search.", "false");
  set_string_parameter("seed", "time", "HIDE ME FROM USER",
      "Given a real-number value, will always produce the same decoy seqs",
      "false");
  set_double_parameter("mass-window", 3.0, 0, 100, 
      "Search peptides within +/- 'mass-window' of the "
      "spectrum mass.  Default 3.0.",
      "Available from the parameter file only for crux-search-for-matches, "
      "crux-create-index, and crux-generate-peptides.",
      "true");
  set_mass_type_parameter("fragment-mass", MONO, 
      "Which isotopes to use in calcuating fragment ion mass "
      "(average, mono). Default mono.", 
      "Parameter file only.  "
      "Used by crux-search-for-matches and crux-predict-peptide-ions.",
                          "true");
  set_double_parameter("ion-tolerance", 0.5, 0, BILLION,
      "Tolerance used for matching observed peaks to predicted "
      "fragment ions.  Default 0.5.",
      "Available from parameter-file for crux-search-for-matches.", "true");
  set_string_parameter("mod", "NO MODS",
      "Specify a variable modification to apply to peptides.  " 
      "<mass change>:<aa list>:<max per peptide>. Default no mods.",
      "Available from parameter file for crux-generate-peptides and "
      "crux-search-for-matches and the "
      "the same must be used for crux compute-q-value.", "true");
  set_string_parameter("cmod", "NO MODS",
      "Specify a variable modification to apply to C-terminus of peptides. " 
      "<mass change>:<max distance from protein c-term (-1 for no max)>. " 
      "Default no mods.",       
      "Available from parameter file for crux-generate-peptides and "
      "crux-search-for-matches and the "
      "the same must be used for crux compute-q-value.", "true");
  set_string_parameter("nmod", "NO MODS",
      "Specify a variable modification to apply to N-terminus of peptides.  " 
      "<mass change>:<max distance from protein n-term (-1 for no max)>",
      "Available from parameter file for crux-generate-peptides and "
      "crux-search-for-matches and the "
      "the same must be used for crux compute-q-value.", "true");
  set_int_parameter("max-mods", MAX_PEPTIDE_LENGTH, 0, MAX_PEPTIDE_LENGTH,
      "The maximum number of modifications that can be applied to a single " 
      "peptide.  Default no limit.",
      "Available from parameter file for crux-search-for-matches.", "true");
  set_int_parameter("max-aas-modified", MAX_PEPTIDE_LENGTH, 0,
      MAX_PEPTIDE_LENGTH,
      "The maximum number of modified amino acids that can appear in one "
      "peptide.  Each aa can be modified multiple times.  Default no limit.",
      "Available from parameter file for search-for-matches.", "true");
  set_int_parameter("precision", 8, 1, 100, //max is arbitrary
      "Set the precision for masses and scores written to sqt and text files. "
      "Default 8.",
      "Available from parameter file for crux search-for-matches, percolator, "
      "and compute-q-values.", "true");
  set_int_parameter("print-search-progress", 10, 0, BILLION,
      "Show search progress by printing every n spectra searched.  Default "
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
  set_int_parameter("top-fit-sp", 1000, 1, BILLION, 
       "Hide from user",
       "used in estimating parameters for exp sp", "false");
  // how many peptides to sample for EVD parameter estimation
  set_int_parameter("sample-count", 500, 0, BILLION, "NOT FOR USER",
       "Number of psms to use for weibul estimation.", "false");

  //in estimate_weibull_parameters
  // no longer used
  set_int_parameter("number-top-scores-to-fit", -1, -10, BILLION, 
      "Not for general users", 
      "The number of psms per spectrum to use for estimating the "
      "score distribution for calculating p-values. 0 to use all. "
      "Not compatible with 'fraction-top-scores-to-fit'. Default 0 (all).",
      "false");
  set_double_parameter("fraction-top-scores-to-fit", 0.55, 0, 1, 
      "The fraction of psms per spectrum to use for estimating the "
      "score distribution for calculating p-values. "
      "Not compatible with 'number-top-scores-to-fig'. Default 0.55.",
      "For developers/research only.", "false");

  /* analyze-matches options */
  set_algorithm_type_parameter("algorithm", PERCOLATOR_ALGORITHM, 
      "The analysis algorithm to use (percolator, curve-fit, none)."
      " Default percolator.",
      "Available only for crux-analyze-matches.  Using 'percolator' will "
      "assign a q-value to the top-ranking psm for each spectrum based on "
      "the decoy searches.  Using 'curve-fit' will assign a q-value to same "
      "using the p-values calculated with score-type=<xcorr-pvalue|"
      "sq-pvalue>.  Incorrect combinations of score-type and algorithm cause"
      " undefined behavior. Using 'none' will turn the binary .csm files "
      "into text.", "false");
  set_string_parameter("feature-file", NULL,//"match_analysis.features"
     "Optional file into which psm features are printed.",
     "Available only for crux-analyze-matches.  File will contain features "
     "used by percolator.", "true");

  /* analyze-matches parameter options */
  set_double_parameter("pi0", 0.9, 0, 1, "Hide from user",
      "Used in curve fitting for assigning q-values from p-values and "
      "used by percolator", "false");
  set_string_parameter("percolator-intraset-features", "F",
      "Set a feature for percolator that in later versions is not an option.",
      "Shouldn't be variable; hide from user.", "false");

  //parameter to parse the tab-delimited result files rather than
  //the binary csm files when analyzing matches.
  set_boolean_parameter("parse-tab-files", FALSE,
    "Parse the tab delimited .txt result files rather "
    "than the .csm files. Default F.",
    "Available for compute-q-values, percolator, and q-ranker.",
    "false");


  /* predict-peptide-ions */
  set_ion_type_parameter("primary-ions", BY_ION,
      "The ion series to predict (b,y,by). Default 'by' (both b and y ions).",
      "Only available for crux-predict-peptide-ions.  Set automatically to "
                         "'by' for searching.", "true");
  set_boolean_parameter("precursor-ions", FALSE,
      "Predict the precursor ions, and all associated ions "
      "(neutral-losses, multiple charge states) consistent with the "
      "other specified options. (T,F) Default F.",
      "Only available for crux-predict-peptide-ions.", "true");
  set_string_parameter("neutral-losses", "all", 
      "Predict neutral loss ions (none, h20, nh3, all). Default 'all'.",
      "Only available for crux-predict-peptide-ions. Set to 'all' for "
      "sp and xcorr scoring.", "true");
  set_int_parameter("isotope", 0, 0, 2,
      "Predict the given number of isotope peaks (0|1|2). Default 0.",
      "Only available for crux-predict-peptide-ion.  Automatically set to "
      "0 for Sp scoring and 1 for xcorr scoring.", "true");
  set_boolean_parameter("flanking", FALSE, 
      "Predict flanking peaks for b and y ions (T,F). Default F.",
      "Only available for crux-predict-peptide-ion.", "true");
  set_string_parameter("max-ion-charge", "peptide",
      "Predict ions up to this charge state (1,2,3) or to the charge state "
      "of the peptide (peptide).  Default 'peptide'.",
      "Available only for predict-peptide-ions.  Set to 'peptide' for search.",
      "true");
  set_int_parameter("nh3",0, -100, BILLION, 
      "Predict peaks with the given maximum number of nh3 neutral loss "
      "modifications. Default 0.",
      "Only available for crux-predict-peptide-ions.", "true");
  set_int_parameter("h2o",0, -100, BILLION,
      "Predict peaks with the given maximum number of h2o neutral loss "
      "modifications. Default 0.",
      "Only available for crux-predict-peptide-ions.", "true");

  /* static mods */
  set_double_parameter("A", 0.0, -100, BILLION, 
      "Change the mass of all amino acids 'A' by the given amount.",
      "For parameter file only.  Default no mass change.", "true");
  set_double_parameter("B", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'B' by the given amount.",
      "For parameter file only.  Default no mass change.", "false");
  set_double_parameter("C", 57.0, -100, BILLION,
      "Change the mass of all amino acids 'C' by the given amount.",
      "For parameter file only.  Default +57.0.", "true");
  set_double_parameter("D", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'D' by the given amount.",
      "For parameter file only.  Default no mass change.", "true");
  set_double_parameter("E", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'E' by the given amount.",
      "For parameter file only.  Default no mass change.", "true");
  set_double_parameter("F", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'F' by the given amount.",
      "For parameter file only.  Default no mass change.", "true");
  set_double_parameter("G", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'G' by the given amount.",
      "For parameter file only.  Default no mass change.", "true");
  set_double_parameter("H", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'H' by the given amount.",
      "For parameter file only.  Default no mass change.", "true");
  set_double_parameter("I", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'I' by the given amount.",
      "For parameter file only.  Default no mass change.", "true");
  set_double_parameter("J", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'J' by the given amount.",
      "For parameter file only.  Default no mass change.", "false");
  set_double_parameter("K", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'K' by the given amount.",
      "For parameter file only.  Default no mass change.", "true");
  set_double_parameter("L", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'L' by the given amount.",
      "For parameter file only.  Default no mass change.", "true");
  set_double_parameter("M", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'M' by the given amount.",
      "For parameter file only.  Default no mass change.", "true");
  set_double_parameter("N", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'N' by the given amount.",
      "For parameter file only.  Default no mass change.", "true");
  set_double_parameter("O", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'O' by the given amount.",
      "For parameter file only.  Default no mass change.", "false");
  set_double_parameter("P", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'P' by the given amount.",
      "For parameter file only.  Default no mass change.", "true");
  set_double_parameter("Q", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'Q' by the given amount.",
      "For parameter file only.  Default no mass change.", "true");
  set_double_parameter("R", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'R' by the given amount.",
      "For parameter file only.  Default no mass change.", "true");
  set_double_parameter("S", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'S' by the given amount.",
      "For parameter file only.  Default no mass change.", "true");
  set_double_parameter("T", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'T' by the given amount.",
      "For parameter file only.  Default no mass change.", "true");
  set_double_parameter("U", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'U' by the given amount.",
      "For parameter file only.  Default no mass change.", "false");
  set_double_parameter("V", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'V' by the given amount.",
      "For parameter file only.  Default no mass change.", "true");
  set_double_parameter("W", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'W' by the given amount.",
      "For parameter file only.  Default no mass change.", "true");
  set_double_parameter("X", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'X' by the given amount.",
      "For parameter file only.  Default no mass change.", "false");
  set_double_parameter("Y", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'Y' by the given amount.",
      "For parameter file only.  Default no mass change.", "true");
  set_double_parameter("Z", 0.0, -100, BILLION,
      "Change the mass of all amino acids 'Z' by the given amount.",
      "For parameter file only.  Default no mass change.", "false");

  /* get-ms2-spectrum options */
  set_boolean_parameter("stats", FALSE, 
      "Print to stdout additional information about the spectrum.",
      "Avaliable only for crux-get-ms2-spectrum.  Does not affect contents "
      "of the output file.", "true");

  // now we have initialized the parameters
  parameter_initialized = TRUE;
  usage_initialized = TRUE;
  type_initialized = TRUE;

}


/*
 * Main calls this to determine which required arguments
 * must be used.  Must be called after initialize_parameters
 * This is the interface for parse_argument_set_req()
 */
BOOLEAN_T select_cmd_line_arguments(  //remove options from name
  const char** option_names,
  int    num_options 
  ){
  select_cmd_line( option_names, num_options, 
                   &parse_arguments_set_req);
  return TRUE;
}

/*
 * Main calls this to determine which options
 * can be used.  Must be called after initialize_parameters
 * This is the interface for parse_argument_set_opt()
 */
BOOLEAN_T select_cmd_line_options(  //remove options from name
  const char** option_names,
  int    num_options 
  ){
  select_cmd_line( option_names, num_options, 
                   &parse_arguments_set_opt);
  return TRUE;
}
/*
 * Private function for doing the work of select_cmd_line_options
 * and select_cmd_line_arguments which is all the same except for
 * the last function call which is now set with a function pointer
 * 
 */
BOOLEAN_T select_cmd_line(  //remove options from name
  const char** option_names,
  int    num_options, 
  int (*parse_arguments_set_ptr)(const char*, const char*, void*, enum argument_type) 
  ){

  carp(CARP_DETAILED_DEBUG, "Selecting options");
  BOOLEAN_T success = TRUE;

  if( (num_options < 1) || (option_names == NULL) ){
    success = FALSE; //?
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
           "Cannot select parameter '%s'. Value, usage or type not found.\n"
           "Found value: %s, usage: %s, type: %s", 
           option_names[i],
           value_ptr,
           usage_ptr,
           type_ptr);
      
    }

    if( //strcmp(type_ptr, "PEPTIDE_TYPE_T") == 0 ||
        strcmp((char*)type_ptr, "MASS_TYPE_T") == 0 ||
        strcmp((char*)type_ptr, "BOOLEAN_T") == 0 ||
        strcmp((char*)type_ptr, "SORT_TYPE_T") == 0 ||
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
BOOLEAN_T find_param_filename(int argc, 
                              char** argv, 
                              char* filename_buffer, 
                              int buffer_size){
  BOOLEAN_T success = TRUE;
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
      success = TRUE;
    }
    else{
      carp(CARP_FATAL, "Parameter filename is too long");
    }
  }
  else{ //parameter_file_index < 0, i.e. no paramter file option
    success = FALSE;
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
  char* location = get_string_parameter("decoy-location");
  int max_rank_preliminary = get_int_parameter("max-rank-preliminary");

  // store new values here
  BOOLEAN_T tdc = FALSE;  // target-decoy competitition
  int new_num_decoy_files = -1;
  int new_max_rank_preliminary = max_rank_preliminary; 

  // user may not have set target-location if no decoys requested
  if( num_decoy_per_target == 0 ){
    free(location);
    location = my_copy_string("separate-decoy-files");
  }

  // set new values
  if( strcmp(location, "target-file") == 0 ){
    tdc = TRUE;
    new_num_decoy_files = 0;

    if( max_rank_preliminary > 0 ){  // scale to num decoys
      new_max_rank_preliminary = max_rank_preliminary * 
                                (1 + num_decoy_per_target);
    }
  }else if( strcmp(location, "one-decoy-file") == 0 ){
    tdc = FALSE;
    new_num_decoy_files = 1;
  }else if( strcmp(location, "separate-decoy-files") == 0 ){
    tdc = FALSE;
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

  if( tdc == TRUE ){
    update_hash_value(parameters, "tdc", (void*)"TRUE");
  }else{
    update_hash_value(parameters, "tdc", (void*)"FALSE");
  }
}

/**
 * Take the command line string from main, find the parameter file 
 * option (if present), parse it's values into the hash, and parse
 * the command line options and arguments into the hash.
 * Main then retrieves the values through get_<type>_parameter.
 * \returns TRUE is command line is successfully parsed.
 */
BOOLEAN_T parse_cmd_line_into_params_hash(int argc, 
                                          char** argv, 
                                          const char* exe_name){
  carp(CARP_DETAILED_DEBUG, "Parameter.c is parsing the command line");
  assert(parameter_initialized && usage_initialized && type_initialized);
  BOOLEAN_T success = TRUE;
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

  // For version option, print version and quit
  if( get_boolean_parameter("version") ){
    printf("Crux version %s\n", VERSION);
    exit(0);
  }

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

  parameter_plasticity = FALSE;

  return success;
}

/**
 * Read the value given for custom-enzyme and enter values into global
 * params.   Correct syntax is <brace>A-Z<brace>|<brace>A-Z<brace> 
 * where <brace> can be [] or {}.  An X indicates that any residue is
 * legal. Sets pre/post_list size and allocates memory for
 * pre/post_cleavage_list.  Sets pre/post_for_inclusion as true if []
 * encloses list or false if {} encloses list. 
 * For special case of [X], set p_cleavage_list as empty and inclusion
 * as false.
 */
// NOTE (BF mar-11-09): for testing would be nice if this returned
// error code instead of dying
void parse_custom_enzyme(char* rule_str){

  BOOLEAN_T success = TRUE;
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
      success = FALSE;      
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
    success = FALSE;
  }else if(rule_str[pre_first_idx-1] == '[' && 
           rule_str[pre_end_idx] == ']'){
    pre_for_inclusion = TRUE;
  }else if(rule_str[pre_first_idx-1] == '{' && 
           rule_str[pre_end_idx] == '}'){
    pre_for_inclusion = FALSE;
  }else{
    success = FALSE;
  }

  // post list
  if(pipe_idx + 2 >= len ){
    success = FALSE;
  }else if(rule_str[post_first_idx-1] == '[' && 
           rule_str[post_end_idx] == ']'){
    post_for_inclusion = TRUE;
  }else if(rule_str[post_first_idx-1] == '{' && 
           rule_str[post_end_idx] == '}'){
    post_for_inclusion = FALSE;
  }else{
    success = FALSE;
  }

  // check that braces aren't empty 
  if(pre_first_idx >= pre_end_idx || post_first_idx >= post_end_idx ){
    success = FALSE;
  }

  if( success == FALSE ){
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
    pre_for_inclusion = FALSE;
  }

  if(strncmp( rule_str+post_first_idx-1, "[X]", post_list_size+2) == 0){
    free(post_cleavage_list);
    post_cleavage_list = NULL;
    post_list_size = 0;
    post_for_inclusion = FALSE;
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
    update_hash_value(parameters, "missed-cleavages", (void*)"TRUE");
  }
}


/*
 * Checks the current value of the named option
 *   as stored in the parameters hash and checks 
 *   that it is a legal value (within min/max for
 *   numeric, correct word for specialized type 
 */
BOOLEAN_T check_option_type_and_bounds(const char* name){

  BOOLEAN_T success = TRUE;
  char die_str[SMALL_BUFFER];
  char* type_str = (char*)get_hash_value(types, name);
  char* value_str = (char*)get_hash_value(parameters, name);
  char* min_str = (char*)get_hash_value(min_values, name);
  char* max_str = (char*)get_hash_value(max_values, name);

  MASS_TYPE_T mass_type;
  //PEPTIDE_TYPE_T pep_type;
  SORT_TYPE_T sort_type;
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
      success = FALSE;
      sprintf(die_str, 
              "The option '%s' must be between %s and %s.  %s is out of bounds",
              name, min_str, max_str, value_str);
    }
    break;
  case STRING_P:
    carp(CARP_DETAILED_DEBUG, "found string opt with value %s\n", value_str);
    //check list of legal values?
    break;
  case MASS_TYPE_P:
    carp(CARP_DETAILED_DEBUG, "found mass_type opt with value %s\n", 
         value_str);
    if( ! string_to_mass_type( value_str, &mass_type )){
      success = FALSE;
      sprintf(die_str, "Illegal mass-type.  Must be 'mono' or 'average'");
    }
    break;
  case DIGEST_TYPE_P:
      carp(CARP_DETAILED_DEBUG, "found digest_type param, value '%s'\n", 
           value_str);
    if( string_to_digest_type(value_str) == INVALID_DIGEST){
      success = FALSE;
      sprintf(die_str, "Illegal digest value. "
              "Must be full-digest or partial-digest.");
    }
    break;
  case ENZYME_TYPE_P:
      carp(CARP_DETAILED_DEBUG, "found enzyme_type param, value '%s'\n", 
           value_str);
    if( string_to_enzyme_type(value_str) == INVALID_ENZYME){
      success = FALSE;
      sprintf(die_str, "Illegal enzyme. Must be trypsin, chymotrypsin, "
              ", elastase, or no-enzyme.");
    }
    break;
  case BOOLEAN_P:
    carp(CARP_DETAILED_DEBUG, "found boolean_type param, value '%s'", 
         value_str);
    if( value_str[0] != 'T' && value_str[0] != 'F'){
      success =  FALSE;
      sprintf(die_str, 
              "Illegal boolean value '%s' for option '%s'.  Must be T or F",
              value_str, name);
    }
    break;
  case SORT_TYPE_P:
    carp(CARP_DETAILED_DEBUG, "found sort_type param, value '%s'",
         value_str);
    if( ! string_to_sort_type( value_str, &sort_type)){
      success = FALSE;
      sprintf(die_str, "Illegal sort value '%s' for option '%s'. "
              "Must be mass, length, lexical, or none.", 
              value_str, name);
    }
    break;
  case SCORER_TYPE_P:
    carp(CARP_DETAILED_DEBUG, "found scorer_type param, value '%s'",
         value_str);
    //check for legal type
    if(! string_to_scorer_type( value_str, &scorer_type)){
      success = FALSE;
      sprintf(die_str, "Illegal score value '%s' for option '%s'.  "
      "Must be sp, xcorr, sp-pvalue, or xcorr-pvalue.", value_str, name);
    }else if((scorer_type != SP ) &&   //check for one of the accepted types
             (scorer_type != XCORR ) &&
             (scorer_type != DOTP ) &&
             (scorer_type != LOGP_BONF_WEIBULL_SP ) &&
             (scorer_type != LOGP_BONF_WEIBULL_XCORR )){
      success = FALSE;
      sprintf(die_str, "Illegal score value '%s' for option '%s'.  "
      "Must be sp, xcorr, sp-pvalue, or xcorr-pvalue.", value_str, name);
    }
    break;
  case ALGORITHM_TYPE_P:
    carp(CARP_DETAILED_DEBUG, "found algorithm_type param, value '%s'",
         value_str);
    if(! string_to_algorithm_type( value_str, &algorithm_type)){
      success = FALSE;
      sprintf(die_str, "Illegal score value '%s' for option '%s'.  "
              "Must be percolator, curve-fit, or none.", value_str, name);
    }
    break;
  case ION_TYPE_P:
    carp(CARP_DETAILED_DEBUG, "found ion_type param, value '%s'",
         value_str);
    if( !string_to_ion_type(value_str, &ion_type)){
      success = FALSE;
      sprintf(die_str, "Illegal ion type '%s' for option '%s'.  "
              "Must be b,y,by.", value_str, name);
    }
    break;
  default:
    carp(CARP_FATAL, "Your param type '%s' wasn't found (code %i)", 
        type_str, (int)param_type);
  }

  if( ! success ){
    carp(CARP_FATAL, die_str);
  }
  return success;
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
  BOOLEAN_T overwrite = get_boolean_parameter("overwrite");
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
      fprintf(param_file, "# %s\n# %s\n%s=%s\n\n",
              (char*)get_hash_value(usages, key),
              (char*)get_hash_value(file_notes, key),
              key,
              (char*)get_hash_value(parameters, key));
    }
  }

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
  }
  // this is mostly so I can test repeatedly
  num_mods = 0;
  num_c_mods = 0;
  num_n_mods = 0;
  list_of_c_mods = NULL;
  list_of_n_mods = NULL;
  parameter_initialized = FALSE;
  parameter_plasticity = TRUE;
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
    /* why does this segfault?  only with break, not without
    if(line[0] == '#' || line[0] == '\0'){
      printf("comment or blank line");
      break;
    }
    */
    /* empty lines and those beginning with '#' are ignored */
    if(line[0] != '#' && line[0] != '\0'){

      /* find the '=' in the line.  Exit with error if the line 
         has no equals sign. */
      idx = 0;
      while(idx < (int)strlen(line) && line[idx] != '='){
        idx++;
      }
      if(idx == 0 || idx >= (int)(strlen(line)-1)){
        carp(CARP_FATAL, "Lines in a parameter file must have the form:\n"
             "\n\tname=value\n\n"
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
 * \returns TRUE if paramater value is TRUE, else FALSE
 */ 
BOOLEAN_T get_boolean_parameter(
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
  BOOLEAN_T found = string_to_param_type(type_str, &type);
 
  if(found == FALSE || type != BOOLEAN_P){
    carp(CARP_ERROR, "Request for boolean parameter '%s' which is of type %s",
         name, type_str);
  }

 // make sure that there is enough storage allocated in the string
  if((int)strlen(value) 
     > PARAMETER_LENGTH) {
    carp(CARP_FATAL, "parameter %s with value %s was too long to copy to string\n",
        name,
        value);
  }
  strncpy(buffer,
          value,
          PARAMETER_LENGTH);

  if ((strcmp(buffer, "TRUE") == 0) || (strcmp(buffer, "T") == 0)){
    return(TRUE);
  } 
  else if ((strcmp(buffer, "FALSE") == 0) || (strcmp(buffer, "F") == 0)){
    return(FALSE);
  } 
  else {
    carp(CARP_FATAL, "Invalid Boolean parameter %s.\n", buffer);
  }
  
  carp(CARP_FATAL, "parameter name: %s, doesn't exist", name);
  return FALSE; // Return value to avoid compiler warning
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
  BOOLEAN_T found = string_to_param_type(type_str, &type);

  if(found==FALSE || type != INT_P){
    carp(CARP_ERROR, "Request for int parameter '%s' which is of type %s",
         name, type_str);
  }

  /* there is a parameter with the right name.  Now 
     try to convert it to a base 10 integer*/
  value = atoi(int_value);
  /*  value = strtol(int_value, &endptr, 10);
  if ((value == LONG_MIN) || 
      (value == LONG_MAX) || 
      (endptr == int_value)) {
    carp(CARP_FATAL, "Conversion error when trying to convert parameter %s with value %s to an int\n",
        name, 
        int_value);
        exit(1);
  } 
  return((int)value);
  */
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
  BOOLEAN_T found = string_to_param_type(type_str, &type);

  if(found==FALSE || type != DOUBLE_P){
    carp(CARP_ERROR, "Request for double parameter '%s' which is of type %s",
         name, type_str);
  }
  
  /* there is a parameter with the right name.  Now 
     try to convert it to a double*/
  value = strtod(double_value, &endptr);
  /*if((value == HUGE_VALF) ||  // AAK removed //BF: why?
    (value == -HUGE_VALF) || 
    (endptr == double_value)) {
    capr(CARP_FATAL, "Conversion error when trying to convert parameter %s with value %s to an double\n",
    name,
    double_value);
    exit(1); */
  // } else {  
  return(value);
  // }
  
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

  /*  Let any type be retrieved as string
  if(found==FALSE || type != STRING_P){
    carp(CARP_ERROR, "Request for string parameter '%s' which is of type %s",
         name, type_str);
  }
  */

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

  /*if(found==FALSE || type != STRING_P){
    carp(CARP_ERROR, "Request for string parameter '%s' which is of type %s",
         name, type_str);
         }*/

  return string_value;

}

DIGEST_T get_digest_type_parameter( const char* name ){

  char* param = (char*)get_hash_value(parameters, name);

  DIGEST_T digest_type = string_to_digest_type(param);
  if( digest_type == INVALID_DIGEST ){
    carp(CARP_FATAL, "Digest_type parameter %s has the value %s " 
         "which is not of the correct type\n", name, param);
  }
  return digest_type;
}

ENZYME_T get_enzyme_type_parameter( const char* name ){

  char* param = (char*)get_hash_value(parameters, name);

  ENZYME_T enzyme_type = string_to_enzyme_type(param);
  if( enzyme_type == INVALID_ENZYME ){
    carp(CARP_FATAL, "Enzyme_type parameter %s has the value %s " 
         "which is not of the correct type\n", name, param);
  }
  return enzyme_type;
}

MASS_TYPE_T get_mass_type_parameter(
   const char* name
   ){
  char* param_value_str = (char*)get_hash_value(parameters, name);
  MASS_TYPE_T param_value;
  BOOLEAN_T success = string_to_mass_type(param_value_str, &param_value);

  if( ! success ){
    carp(CARP_FATAL, 
         "Mass_type parameter %s has the value %s which is not of "
          "the correct type", name, param_value_str);
  }
  return param_value;
}

SORT_TYPE_T get_sort_type_parameter(const char* name){
  char* param_value_str = (char*)get_hash_value(parameters, name);
  SORT_TYPE_T param_value;
  BOOLEAN_T success = string_to_sort_type(param_value_str, &param_value);

  if( ! success){
    carp(CARP_FATAL, "Sort_type parameter %s has the value %s which " 
         "is not of the correct type", name, param_value_str);
  }
  return param_value;
}

ALGORITHM_TYPE_T get_algorithm_type_parameter(const char* name){
  char* param_value_str = (char*)get_hash_value(parameters, name);
  ALGORITHM_TYPE_T param_value;
  BOOLEAN_T success = string_to_algorithm_type(param_value_str, &param_value);

  if(!success){
    carp(CARP_FATAL, "Algorithm_type parameter %s has the value %s "
         "which is not of the correct type.", name, param_value_str);
  }
  return param_value;
}


SCORER_TYPE_T get_scorer_type_parameter(const char* name){
  char* param_value_str = (char*)get_hash_value(parameters, name);
  SCORER_TYPE_T param_value;
  BOOLEAN_T success = string_to_scorer_type(param_value_str, &param_value);

  if(!success){
    carp(CARP_FATAL, "Scorer_type parameter %s has the value %s " 
         "which is not of the correct type.", name, param_value_str);
  }
  return param_value;
}

ION_TYPE_T get_ion_type_parameter(const char* name){
  char* param_value_str = (char*)get_hash_value(parameters, name);
  ION_TYPE_T param_value;
  BOOLEAN_T success = string_to_ion_type(param_value_str, &param_value);

  if(!success){
    carp(CARP_FATAL, 
   "Ion_type parameter %s ahs the value %s which is not of the correct type.",
         name, param_value_str);
  }
  return param_value;
}
/**************************************************
 *   SETTERS (private)
 **************************************************
 */

BOOLEAN_T set_boolean_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 BOOLEAN_T set_value,  ///< the value to be set -in
 const char* usage,          ///< message for the usage statement
 const char* filenotes,      ///< additional informatino for the params file
 const char* foruser
  )
{
  BOOLEAN_T result;
    
  // check if parameters cah be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }

  const char* bool_str;
  if(set_value){
    bool_str = "TRUE";
  }
  else{
    bool_str = "FALSE";
  }
  result = add_or_update_hash_copy(parameters, name, bool_str);
  result = add_or_update_hash_copy(usages, name, usage);
  result = add_or_update_hash_copy(file_notes, name, filenotes);
  result = add_or_update_hash_copy(for_users, name, foruser);
  result = add_or_update_hash_copy(types, name, (void*)"BOOLEAN_T");
  return result;
}

BOOLEAN_T set_int_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 int set_value,  ///< the value to be set -in
 int min_value,  ///< the minimum accepted value -in
 int max_value,  ///< the maximum accepted value -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,  ///< additional info for param file
 const char* foruser     ///< true if should be revealed to user
  )
{
  BOOLEAN_T result;
  char buffer[PARAMETER_LENGTH];
  
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  
  //stringify default, min, and max values and set
  snprintf(buffer, PARAMETER_LENGTH, "%i", set_value);
  result = add_or_update_hash_copy(parameters, name, buffer);

  snprintf(buffer, PARAMETER_LENGTH, "%i", min_value);
  result = add_or_update_hash_copy(min_values, name, buffer);

  snprintf(buffer, PARAMETER_LENGTH, "%i", max_value);
  result = add_or_update_hash_copy(max_values, name, buffer);

  result = add_or_update_hash_copy(usages, name, usage);
  result = add_or_update_hash_copy(file_notes, name, filenotes);
  result = add_or_update_hash_copy(for_users, name, foruser);
  result = add_or_update_hash_copy(types, name, (void*)"INT_ARG");
  return result;
}

BOOLEAN_T set_double_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 double set_value,  ///< the value to be set -in
 double min_value,  ///< the value to be set -in
 double max_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,  ///< additional info for param file
 const char* foruser     ///< "true" if should be revealed to user
  )
{
  BOOLEAN_T result;
  char buffer[PARAMETER_LENGTH];
  
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  
  // convert to string
  snprintf(buffer, PARAMETER_LENGTH, "%f", set_value);
  result = add_or_update_hash_copy(parameters, name, buffer);    

  snprintf(buffer, PARAMETER_LENGTH, "%f", min_value);
  result = add_or_update_hash_copy(min_values, name, buffer);    

  snprintf(buffer, PARAMETER_LENGTH, "%f", max_value);
  result = add_or_update_hash_copy(max_values, name, buffer);    

  result = add_or_update_hash_copy(usages, name, usage);    
  result = add_or_update_hash_copy(file_notes, name, filenotes);
  result = add_or_update_hash_copy(for_users, name, foruser);
  result = add_or_update_hash_copy(types, name, (void*)"DOUBLE_ARG");    
  return result;
}

/**
 * temporary replacement for function, return name once all exe's are fixed
 * \returns TRUE if paramater value is set, else FALSE
 */ 
BOOLEAN_T set_string_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 const char* set_value,  ///< the value to be set -in
 const char* usage,
 const char* filenotes,  ///< additional comments for parameter file
 const char* foruser
  )
{
  BOOLEAN_T result;

  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  
  if( set_value == NULL ){
    set_value = "__NULL_STR";
  }

  result = add_or_update_hash_copy(parameters, name, set_value);
  result = add_or_update_hash_copy(usages, name, usage);
  result = add_or_update_hash_copy(file_notes, name, filenotes);
  result = add_or_update_hash_copy(for_users, name, foruser);
  result = add_or_update_hash_copy(types, name, (void*)"STRING_ARG");

  return result;
}

BOOLEAN_T set_mass_type_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 MASS_TYPE_T set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,  ///< additional info for param file
 const char* foruser
 )
{
  BOOLEAN_T result;
  char value_str[265] ;

  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  
  /* stringify the value */
  mass_type_to_string(set_value, value_str);
  
  result = add_or_update_hash_copy(parameters, name, value_str);
  result = add_or_update_hash_copy(usages, name, usage);
  result = add_or_update_hash_copy(file_notes, name, filenotes);
  result = add_or_update_hash_copy(for_users, name, foruser);
  result = add_or_update_hash_copy(types, name, (void*)"MASS_TYPE_T");
  return result;

}

BOOLEAN_T set_digest_type_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 DIGEST_T set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser
  )
{
  BOOLEAN_T result = TRUE;
  
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  
  /* stringify the value */
  char* value_str = digest_type_to_string(set_value);
  carp(CARP_DETAILED_DEBUG, "Setting digest param '%s' to value '%s'.", name, value_str);

  result = add_or_update_hash_copy(parameters, name, value_str);
  result = add_or_update_hash_copy(usages, name, usage);
  result = add_or_update_hash_copy(file_notes, name, filenotes);
  result = add_or_update_hash_copy(for_users, name, foruser);
  result = add_or_update_hash_copy(types, name, (void*)"DIGEST_T");
  free(value_str);
  return result;

}

BOOLEAN_T set_enzyme_type_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 ENZYME_T set_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,   ///< additional info for param file
 const char* foruser
  )
{
  BOOLEAN_T result = TRUE;
  
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  
  /* stringify the value */
  char* value_str = enzyme_type_to_string(set_value);

  result = add_or_update_hash_copy(parameters, name, value_str);
  result = add_or_update_hash_copy(usages, name, usage);
  result = add_or_update_hash_copy(file_notes, name, filenotes);
  result = add_or_update_hash_copy(for_users, name, foruser);
  result = add_or_update_hash_copy(types, name, (void*)"ENZYME_T");
  free(value_str);
  return result;

}

BOOLEAN_T set_sort_type_parameter(
  const char* name,
  SORT_TYPE_T set_value,
  const char* usage,
  const char* filenotes,
  const char* foruser)
{
  BOOLEAN_T result = TRUE;
  char value_str[SMALL_BUFFER];
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  /* stringify value */
  sort_type_to_string(set_value, value_str);
  
  result = add_or_update_hash_copy(parameters, name, value_str);
  result = add_or_update_hash_copy(usages, name, usage);
  result = add_or_update_hash_copy(file_notes, name, filenotes);
  result = add_or_update_hash_copy(for_users, name, foruser);
  result = add_or_update_hash_copy(types, name, (void*)"SORT_TYPE_T");
  return result;
}

BOOLEAN_T set_algorithm_type_parameter(
 const char* name,
 ALGORITHM_TYPE_T set_value,
 const char* usage,
 const char* filenotes,
 const char* foruser)
{
  BOOLEAN_T result = TRUE;
  char value_str[SMALL_BUFFER];
  
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  /* stringify value */
  algorithm_type_to_string(set_value, value_str);
  carp(CARP_DETAILED_DEBUG, "setting algorithm type to %s", value_str);  

  result = add_or_update_hash_copy(parameters, name, value_str);
  result = add_or_update_hash_copy(usages, name, usage);
  result = add_or_update_hash_copy(file_notes, name, filenotes);
  result = add_or_update_hash_copy(for_users, name, foruser);
  result = add_or_update_hash_copy(types, name, (void*)"ALGORITHM_TYPE_T");
  return result;
}


BOOLEAN_T set_scorer_type_parameter(
 const char* name,
 SCORER_TYPE_T set_value,
 const char* usage, 
 const char* filenotes,
 const char* foruser)
{
  BOOLEAN_T result = TRUE;
  char value_str[SMALL_BUFFER];
  
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  /* stringify value */
  scorer_type_to_string(set_value, value_str);
  carp(CARP_DETAILED_DEBUG, "setting score type to %s", value_str);  

  result = add_or_update_hash_copy(parameters, name, value_str);
  result = add_or_update_hash_copy(usages, name, usage);
  result = add_or_update_hash_copy(file_notes, name, filenotes);
  result = add_or_update_hash_copy(for_users, name, foruser);
  result = add_or_update_hash_copy(types, name, (void*)"SCORER_TYPE_T");


  return result;
}

BOOLEAN_T set_ion_type_parameter(
 const char* name,
 ION_TYPE_T set_value,
 const char* usage,
 const char* filenotes,
 const char* foruser)
{
  BOOLEAN_T result = TRUE;
  char value_str[SMALL_BUFFER];

  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  
  /* stringify value */
  ion_type_to_string(set_value, value_str);

  result = add_or_update_hash_copy(parameters, name, value_str);
  result = add_or_update_hash_copy(usages, name, usage);
  result = add_or_update_hash_copy(file_notes, name, filenotes);
  result = add_or_update_hash_copy(for_users, name, foruser);
  result = add_or_update_hash_copy(types, name, (void*)"ION_TYPE_T");
  return result;
}
/**
 * Routines that return crux enumerated types. 
 */

BOOLEAN_T string_to_param_type(const char* name, PARAMETER_TYPE_T* result ){
  BOOLEAN_T success = TRUE;
  if( name == NULL ){
    return FALSE;
  }

  int param_type = convert_enum_type_str(
                   name, -10, parameter_type_strings, NUMBER_PARAMETER_TYPES);
  (*result) = (PARAMETER_TYPE_T)param_type;

  if( param_type == -10 ){
    success = FALSE;
  }
  return success;
}

/*
 * Applies any static mods to the aa masses
 */
BOOLEAN_T update_aa_masses(){
  BOOLEAN_T success = TRUE;
  int aa;
  char aa_str[2];
  aa_str[1] = '\0';
  int alphabet_size = (int)'A' + ((int)'Z'-(int)'A'); //get_alpha_size(ALL_SIZE);
  carp(CARP_DETAILED_DEBUG, "updating masses, last is %d", alphabet_size);

  for(aa=(int)'A'; aa< alphabet_size -1; aa++){
    aa_str[0] = (char)aa;
    //aa_to_string(aa, aa_str);
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
  //carp(CARP_DEBUG, "getting aa mods, all %d of them", num_mods);
  *mods = list_of_mods;
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
  //carp(CARP_DEBUG, "getting c mods, all %d of them", num_c_mods);
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
  //carp(CARP_DEBUG, "getting n mods, all %d of them", num_n_mods);
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
  return num_mods + num_c_mods + num_n_mods;
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
char* read_mass_change(AA_MOD_T* mod, char* line, char separator){
  //carp(CARP_DEBUG, "token points to %s", line);

  aa_mod_set_mass_change(mod, atof(line));
  //mod->mass_change = atof(line);
  if( aa_mod_get_mass_change(mod) == 0){
  //if( mod->mass_change == 0){
    carp(CARP_FATAL, "The mass change is not valid for mod %s", line);
  }
  char* next = line;
  while(*next != separator){
    next++;
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

  BOOLEAN_T* aa_list = aa_mod_get_aa_list(mod);
  while( *line != '\0' && *line != ':'){
    char aa = toupper( *line );
    carp(CARP_DETAILED_DEBUG, "aa is %c", aa);

    if( aa < 'A' || aa > 'Z' ){
      carp(CARP_FATAL, "The letter '%c' in the aa list is invalid.", aa);
    }
    carp(CARP_DETAILED_DEBUG, "aa index is %d", aa - 'A');
    aa_list[aa - 'A'] = TRUE;
    //mod->aa_list[aa - 'A'] = TRUE;
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
void read_max_per_peptide(AA_MOD_T* mod, char* line){
  //carp(CARP_DETAILED_DEBUG, "token points to %s", line);
  if( *line == '\0' ){
    carp(CARP_FATAL, "Missing maximum mods per peptide for mod %s", line);
  }

  aa_mod_set_max_per_peptide(mod, atoi(line));
  //mod->max_per_peptide = atoi(line);
  if( aa_mod_get_max_per_peptide(mod) == 0 ){
  //if( mod->max_per_peptide == 0 ){
    carp(CARP_FATAL, "Maximum mods per peptide is invalid for mod %s", line);
  }

}

/**
 * \brief Set the max_distance field in the mod with the value
 * pointed to by line.
 *
 * Fails if line does not point to a valid integer.
 * \returns void
 */
void read_max_distance(AA_MOD_T* mod, char* line){
  //carp(CARP_DEBUG, "token points to %s", line);
  if( *line == '\0' ){
    carp(CARP_FATAL,
         "Missing maximum distance from protein terminus for mod %s", line);
  }

  aa_mod_set_max_distance(mod, atoi(line));
  // 0 is a valid distance, would have to check some other way
  //    carp(CARP_FATAL, "Maximum mods per peptide is invalid for mod %s", line);

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
 MOD_POSITION_T position){///< type of mod (any, c-, n-term)

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
    token = read_mass_change(cur_mod, token, ':');

    // fill in values for standard mods
    if( position == ANY_POSITION ){
      // read the aa list and set the values in mod
      token = set_aa_list(cur_mod, token, ':');

      // get max per peptide
      read_max_per_peptide(cur_mod, token);
    }// fill in values for c- or n-mod
    else{
      // get the max distance
      read_max_distance(cur_mod, token);

      // set all bools to true
      int i = 0;
      BOOLEAN_T* aa_list = aa_mod_get_aa_list(cur_mod);
      for(i=0; i<AA_LIST_LENGTH; i++){
        aa_list[i] = TRUE;
      }
      // set type to c-/n-term and max to 1
      aa_mod_set_position(cur_mod, position);
      aa_mod_set_max_per_peptide(cur_mod, 1);
    }

    //  increment counter and get next mod
    cur_index++;
    //print_mod(cur_mod);

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
  //AA_MOD_T* cur_mod = list_of_mods[num_mods]; // num_mods == 0
  int total_num_mods = 0;

  total_num_mods = read_mods(param_file, total_num_mods,
                             "mod=", ANY_POSITION);
  num_mods = total_num_mods;  // set global var

  // Read the file again to get the cmods
  rewind( param_file );

  // set cmod pointer to next in array
  list_of_c_mods = &list_of_mods[total_num_mods];

  total_num_mods = read_mods(param_file, total_num_mods, "cmod=", C_TERM);
  num_c_mods = total_num_mods - num_mods;

  // if no cmods present, don't point to the list of mods
  if( num_c_mods == 0){
    list_of_c_mods = NULL;
  }

  // Read the file again to get the nmods
  rewind( param_file );

  // set nmod pointer to next in array
  list_of_n_mods = &list_of_mods[total_num_mods];

  total_num_mods = read_mods(param_file, total_num_mods, "nmod=", N_TERM);
  num_n_mods = total_num_mods - num_mods - num_c_mods;

  // if no nmods present, don't point to the list of mods
  if( num_n_mods == 0){
    list_of_n_mods = NULL;
  }

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

