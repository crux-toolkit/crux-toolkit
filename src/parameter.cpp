/***********************************************************************//**
 * \file parameter.cpp
 * FILE: parameter.cpp
 * AUTHOR: written by Tobias Mann, CRUXified by Chris Park
 * Missed-cleavage conversion: Kha Nguyen
 * \brief General parameter handling utilities. MUST declare ALL
 * optional command parameters here inside initalialize_parameters.
 ****************************************************************************/

#include "util/crux-utils.h"
#include "io/LineFileReader.h"
#include "parameter.h"
#include "util/Params.h"
#include "util/StringUtils.h"
#include "util/WinCrux.h"
#include "model/Peptide.h"
#include <iostream>

using namespace std;

/*
 * Global variables
 */

static const char* parameter_type_strings[NUMBER_PARAMETER_TYPES] = { 
  "INT_ARG", "DOUBLE_ARG", "STRING_ARG", "MASS_TYPE_T", "DIGEST_T", 
  "ENZYME_T", 
  "bool", "SCORER_TYPE_T", "ION_TYPE_T",
  "ALGORITHM_T", "HARDKLOR_ALGORITHM_TYPE_T",
  "WINDOW_TYPE_T", "MEASURE_TYPE_T", "THRESHOLD_T", 
  "PARSIMONY_TYPE_T", "QUANT_LEVEL_TYPE_T", "DECOY_TYPE_T", "MASS_FORMAT_T"};

//one hash for parameter values, one for usage statements, one for types
// all hashes keyed on parameter/option name

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

// Global variables specifying enzymatic digestion rules.
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

bool string_to_param_type(const char*, PARAMETER_TYPE_T* );

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
  /* initialize the list of mods */                           
  for (int mod_idx = 0; mod_idx < MAX_AA_MODS; mod_idx++) {
    //initialize_aa_mod(&list_of_mods[mod_idx], mod_idx);     
    list_of_mods[mod_idx] = new AA_MOD_T(mod_idx);              
  }                                                           

  /* initialize custom enzyme variables */
  pre_list_size = 0;
  post_list_size = 0;
  pre_cleavage_list = NULL;
  post_cleavage_list = NULL;
  pre_for_inclusion = true;
  post_for_inclusion = false;

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
}

/**
 * Read the value given for custom-enzyme and enter values into global
 * params.  Correct syntax is [A-Z]|[A-Z] or {A-Z}|{A-Z}.  An X
 * indicates that any residue is legal. Sets pre/post_list size and
 * allocates memory for pre/post_cleavage_list.  Sets
 * pre/post_for_inclusion as true if [] encloses list or false if {}
 * encloses list.  For special case of [X], set p_cleavage_list as
 * empty and inclusion as false.
 *
 * Global variables: pre_cleavage_list, post_cleavage_list,
 * pre_list_size, post_list_size, pre_for_inclusion,
 * post_for_inclusion;
 *
 */
void parse_custom_enzyme(const string& rule_str){
  bool success = true;
  int len = rule_str.length();
  int idx = 0;
  int pipe_idx = 0;

  // 1. find the |
  for(idx = 0; idx < len; idx++) {
    if( rule_str[idx] == '|' ) {
      pipe_idx = idx;
      break;
    }
  }
  // Check that there isn't a second |.
  for(idx = idx+1; idx < len; idx++) {
    if( rule_str[idx] == '|' ) {
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

  // 3. Check that braces match, and set inclusion
  // pre-list.
  if (pipe_idx < 1) {
    success = false;
  } else if ( (rule_str[pre_first_idx-1] == '[') && 
              (rule_str[pre_end_idx] == ']') ) {
    pre_for_inclusion = true;
  } else if ( (rule_str[pre_first_idx-1] == '{') && 
              (rule_str[pre_end_idx] == '}') ) {
    pre_for_inclusion = false;
  } else {
    success = false;
  }

  // post list
  if (pipe_idx + 2 >= len ) {
    success = false;
  } else if ( (rule_str[post_first_idx-1] == '[') && 
              (rule_str[post_end_idx] == ']') ) {
    post_for_inclusion = true;
  } else if ( (rule_str[post_first_idx-1] == '{') && 
              (rule_str[post_end_idx] == '}') ) {
    post_for_inclusion = false;
  } else {
    success = false;
  }

  // check that braces aren't empty 
  if ( (pre_first_idx >= pre_end_idx) || (post_first_idx >= post_end_idx) ) {
    success = false;
  }

  if ( success == false ) {
    carp(CARP_FATAL, "Custom enzyme syntax '%s' is incorrect.  "
         "Must be of the form [AZ]|[AZ] or with [] replaced by {}. "
         "AZ is a list of residues (letters A-Z) required [] or prohibited {}. "
         "Use [X] to indicate that any residue is legal.",
         rule_str.c_str());
  }

  // 4. allocate lists and fill
  pre_list_size = pre_end_idx - pre_first_idx;
  pre_cleavage_list = (char*)mycalloc(pre_list_size, sizeof(char));
  for (idx = 0; idx < pre_list_size; idx++) {
    pre_cleavage_list[idx] = rule_str[pre_first_idx+idx];
  }

  post_list_size = post_end_idx - post_first_idx;
  post_cleavage_list = (char*)mycalloc(post_list_size, sizeof(char));
  for (idx = 0; idx < post_list_size; idx++) {
    post_cleavage_list[idx] = rule_str[post_first_idx+idx];
  }

  // 5. check special case of [X] (i.e., cleave before/after every
  //    residue) or {X} (i.e., do not cleave before/after any residue).
  if (strncmp( rule_str.c_str(), "[X]", pre_list_size+2) == 0) {
    free(pre_cleavage_list);
    pre_cleavage_list = NULL;
    pre_list_size = 0;
    pre_for_inclusion = false;
  }
  if (strncmp( rule_str.c_str(), "{X}", pre_list_size+2) == 0) {
    free(pre_cleavage_list);
    pre_cleavage_list = NULL;
    pre_list_size = 0;
    pre_for_inclusion = true;
  }
  if (strncmp( rule_str.c_str()+post_first_idx-1, "[X]", post_list_size+2) == 0) {
    free(post_cleavage_list);
    post_cleavage_list = NULL;
    post_list_size = 0;
    post_for_inclusion = false;
  }
  if (strncmp( rule_str.c_str()+post_first_idx-1, "{X}", post_list_size+2) == 0) {
    free(post_cleavage_list);
    post_cleavage_list = NULL;
    post_list_size = 0;
    post_for_inclusion = true;
  }

  // 6. Say what happened.
  carp(CARP_INFO, "Pre-cleavage list: %s", pre_cleavage_list);
  if (pre_for_inclusion) {
    carp(CARP_INFO, "Pre-cleavage residues are included.");
  } else {
    carp(CARP_INFO, "Pre-cleavage residues are excluded.");
  }
  carp(CARP_INFO, "Post-cleavage list: %s", post_cleavage_list);
  if (post_for_inclusion) {
    carp(CARP_INFO, "Post-cleavage residues are included.");
  } else {
    carp(CARP_INFO, "Post-cleavage residues are excluded.");
  }

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
void print_mods_parameter_file(ostream* param_file, 
                               const char* name,
                               int (*mod_getter)(AA_MOD_T***)){
  // get mod description
  char comments[PARAMETER_BUFFER] = "";
  strcat_formatted(comments, "# ", Params::GetUsage(name).c_str());
  strcat_formatted(comments, "# ", Params::GetFileNotes(name).c_str());
  int precision = Params::GetInt("mod-precision");

  // get list of mods to print
  AA_MOD_T** mod_list = NULL;
  int total_mods = (*mod_getter)(&mod_list);
  for( int mod_idx = 0 ; mod_idx < total_mods; mod_idx++){
    float mass = mod_list[mod_idx]->getMassChange();

    // standard mods have the format mass:aa_list:max
    if( strcmp(name, "mod") == 0 ){
      int max = mod_list[mod_idx]->getMaxPerPeptide();
      bool* aas_modified = mod_list[mod_idx]->getAAList();
      char aa_str[PARAMETER_BUFFER] = "";
      char* aa_str_ptr = aa_str;
      for(int aa_idx = 0; aa_idx < AA_LIST_LENGTH; aa_idx++){
        if( aas_modified[aa_idx] == true ){
          sprintf(aa_str_ptr, "%c", (aa_idx + 'A'));
          aa_str_ptr++;
        }
      }
      char buffer[1024];
      sprintf(buffer, "%s%s=%.*f:%s:%i", comments, name, precision, 
              mass, aa_str, max); 
      *param_file << buffer << endl << endl;
    } else { // nmod, cmod have the format mass:end_distance
      int distance = mod_list[mod_idx]->getMaxDistance();
      char buffer[1024];
      sprintf(buffer, "%s%s=%.*f:%i", comments, name, precision, 
              mass, distance);
      *param_file << buffer << endl << endl;
    }
  }

  // if there were no mods, print placeholder
  if( total_mods == 0 ){
    *param_file << Params::ProcessHtmlDocTags(comments) << name << "=NO MODS" << endl << endl;
  }
}

/**
 *
 * parse the parameter file given the filename
 */
void parse_parameter_file(
  const char* parameter_filename ///< the parameter file to be parsed -in
  )
{
  FILE *file;
  char *line;
  int idx;

  carp(CARP_DETAILED_DEBUG, "Parsing parameter file '%s'",parameter_filename);

  /* check if parameter file exists, if not die */
  if(access(parameter_filename, F_OK)){
    carp(CARP_FATAL, "Could not open parameter file '%s'", parameter_filename);
  }

  LineFileReader line_reader(parameter_filename);
  
  bool found_comet = false;

  while(line_reader.hasNext()) {
    string line = line_reader.next();
    bool found_equal = false;
    if (found_comet) {
      comet_enzyme_info_lines_.push_back(line);
    } else {
      if (line.find("[COMET_ENZYME_INFO]") != string::npos) {
        comet_enzyme_info_lines_.clear();
        found_comet = true;
      } else {
        size_t comment = line.find('#');
        if (comment != string::npos) {
          line = line.substr(0, comment);
        }
        size_t equals = line.find('=');
        if (equals == string::npos) {
          continue;
        }
        string option_name = StringUtils::Trim(line.substr(0, equals));
        string option_value = StringUtils::Trim(line.substr(equals + 1));
        if (!Params::Exists(option_name)) {
          carp(CARP_WARNING, "Read parameter '%s' from parameter file, but no "
                             "such parameter exists.", option_name.c_str());
          continue;
        } else if (Params::IsArgument(option_name)) {
          carp(CARP_WARNING, "Read parameter '%s' from parameter file, but is "
                             "an argument, not an option.", option_name.c_str());
          continue;
        }
        Params::Set(option_name, option_value);
      }
    }
  }
  if (comet_enzyme_info_lines_.empty()) {
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
}

/**************************************************
 *   GETTERS (public)
 **************************************************
 */

DIGEST_T get_digest_type_parameter( const char* name ){
  string param = Params::GetString(name);
  DIGEST_T digest_type = string_to_digest_type(param);
  if( digest_type == INVALID_DIGEST ){
    carp(CARP_FATAL, "Digest_type parameter %s has the value %s " 
         "which is not of the correct type.", name, param.c_str());
  }
  return digest_type;
}

ENZYME_T get_enzyme_type_parameter( const char* name ){
  string param = Params::GetString(name);
  ENZYME_T enzyme_type = string_to_enzyme_type(param);
  if( enzyme_type == INVALID_ENZYME ){
    carp(CARP_FATAL, "Enzyme_type parameter %s has the value %s " 
         "which is not of the correct type.", name, param.c_str());
  }
  return enzyme_type;
}

MASS_TYPE_T get_mass_type_parameter(
   const char* name
   ){
  string param_value_str = Params::GetString(name);
  MASS_TYPE_T param_value;
  bool success = string_to_mass_type(param_value_str, &param_value);

  if( ! success ){
    carp(CARP_FATAL, 
         "Mass_type parameter %s has the value %s which is not of "
          "the correct type", name, param_value_str.c_str());
  }
  return param_value;
}

THRESHOLD_T get_threshold_type_parameter(
  const char* name
  ){
  return string_to_threshold_type(Params::GetString(name));
}

DECOY_TYPE_T get_decoy_type_parameter(
  const char* name
  ){
  return string_to_decoy_type(Params::GetString(name));
}

DECOY_TYPE_T get_tide_decoy_type_parameter(
  const char* name
) {
  return string_to_tide_decoy_type(Params::GetString(name));
}

MASS_FORMAT_T get_mass_format_type_parameter(
  const char* name
  ){
  return string_to_mass_format(Params::GetString(name));
}

int get_max_ion_charge_parameter(
  const char* name
  ){
  string param_value_str = Params::GetString(name);
  if (param_value_str == "peptide") {
    return BILLION; //using this with min function on peptide charge.
  }
  int ans = atoi(param_value_str.c_str());
  if (ans <= 0 || ans > 6) {
    carp(CARP_FATAL,
      "Max_ion_charge parameter %s has the value %s which is not a "
      "legal value", name, param_value_str.c_str());
  }
  return ans;
}

char get_delimiter_parameter(
  const char* name
  ) {
  string param_value_str = Params::GetString(name);
  if (param_value_str == "tab") {
    return '\t'; //using this with min function on peptide charge.
  } else {
    if (param_value_str.length() != 1) {
      carp(CARP_FATAL,
        "delimiter parameter %s with value %s is not a single character or "
        "'tab'", name, param_value_str.c_str());
    }
    return param_value_str[0];
  }
}

COLTYPE_T get_column_type_parameter(const char* name){
  return string_to_column_type(Params::GetString(name));
}

COMPARISON_T get_comparison_parameter(const char* name) {
  return string_to_comparison(Params::GetString(name));
}


/**************************************************
 *   SETTERS (private)
 **************************************************
 */

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

  mod->setMassChange(atof(line));
  if( mod->getMassChange() == 0){
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

  bool* aa_list = mod->getAAList();
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

  mod->setMaxPerPeptide(atoi(line));
  if( mod->getMaxPerPeptide() == 0 ){
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
      mod->setPreventsCleavage(true);
      break;
    case 'F':
    case 'f':
      carp(CARP_DEBUG, "prevents_cleavage set to false %s", line);
      mod->setPreventsCleavage(false);
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
      mod->setPreventsXLink(true);
      break;
    case 'F':
    case 'f':
      carp(CARP_DEBUG, "prevents_xlink set to false %s",line);
      mod->setPreventsXLink(false);
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
  mod->setMaxDistance(max_distance);

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


    } else { // fill in values for c- or n-mod
      // get the max distance
      read_max_distance(cur_mod, token);

      // set all bools to true
      int i = 0;
      bool* aa_list = cur_mod->getAAList();
      for(i=0; i<AA_LIST_LENGTH; i++){
        aa_list[i] = true;
      }
      // set type to c-/n-term and max to 1
      cur_mod->setPosition(position);
      cur_mod->setMaxPerPeptide(1);
    }

    //  increment counter and get next mod
    cur_index++;

  }// repeat until end of file

  free(line);
  return cur_index;
}


void parse_modspec_line(
  const std::string& line,
  vector<AA_MOD_T*>& static_mods,
  vector<AA_MOD_T*>& variable_mods,
  MOD_POSITION_T position = ANY_POSITION
  ) {

  static_mods.clear();
  variable_mods.clear();
  if (line == "") {
    return;
  }
  int pos = 0;
  bool terminus = position != ANY_POSITION;
  const char* spec_text = line.c_str();

  while (true) {
    char c;
    int next_pos = -1;
    sscanf(spec_text + pos, " %c%n", &c, &next_pos);
    if (next_pos == -1) {
      carp(CARP_FATAL, "Expected modification specification. %s %i", spec_text, pos);
    }
    pos += next_pos - 1;
    unsigned int limit = 0;
    if (c >= '1' && c <= '9') {
      sscanf(spec_text + pos, "%u%n", &limit, &next_pos);
      if (limit == UINT_MAX) {
        carp(CARP_FATAL, "Limit too big. %s %i", spec_text, pos);
      }
      pos += next_pos;
    }
    int aa_len = -1, plus_pos = -1, delta_pos = -1, end_pos = -1;
    if (!terminus)
      sscanf(spec_text + pos, "%*[ACDEFGHIKLMNPQRSTVWY]%n%n%*[+-]%n%*[0-9.]%n",
	     &aa_len, &plus_pos, &delta_pos, &end_pos);
    else 
      sscanf(spec_text + pos, "%*[ACDEFGHIKLMNPQRSTVWYX]%n%n%*[+-]%n%*[0-9.]%n",
	     &aa_len, &plus_pos, &delta_pos, &end_pos);

    if (aa_len == -1) {
      carp(CARP_FATAL, "Expected amino acid symbol. %s %i", spec_text, pos);
    }
    assert(plus_pos != -1);
    if (delta_pos == -1) {
      carp(CARP_FATAL, "Expected '+' 'or' - and modification amount. %s %i", spec_text, pos + plus_pos);
    }
    if (end_pos == -1) {
      carp(CARP_FATAL, "Expected modification amount. %s %i", spec_text, pos+delta_pos);
    }
    if ((limit == 0) && (aa_len != 1))
      carp(CARP_FATAL, "Static modifications must be specified "
		   "for one amino acid at a time. %s %i", spec_text, pos);
    int confirm_end_pos = -1;
    double delta;
    sscanf(spec_text + pos + delta_pos, "%lg%n", &delta, &confirm_end_pos);
    if (delta_pos + confirm_end_pos != end_pos) {
      carp(CARP_FATAL, "Cannot parse modification amount. %s %i", spec_text, pos + delta_pos);
    }
    if (*(spec_text + pos + plus_pos) == '-') {
      delta *= -1;
    }
    
    AA_MOD_T *mod = new AA_MOD_T();
    mod->setMaxPerPeptide(limit);
    mod->setMassChange(delta);
    mod->setAminoAcids(spec_text + pos, aa_len);
    
    bool prevents_xlink = false;
    bool prevents_cleavage = false;
    //* Optional prevents xlink, prevent cleavage.
    //first T indicates prvents cleavage
    pos += end_pos;
    if (spec_text[pos] == ':') {
      pos++;
      if (spec_text[pos] == 'T') {
        prevents_cleavage = true;
      } else if (spec_text[pos] != 'F') {
        carp(CARP_FATAL, "Cannot parse T|F for prevents cleavage %s %i", spec_text, pos);
      }
      pos++;
    }
    if (spec_text[pos] == ':') {
      pos++;
      if (spec_text[pos] == 'T') {
        prevents_xlink = true;
      } else if (spec_text[pos] != 'F') {
        carp(CARP_FATAL, "Cannot parse T|F for prevents xlink %s %i", spec_text, pos);
      }
      pos++;
    }
    
        
    mod->setPosition(position);
    mod->setPreventsXLink(prevents_xlink);
    mod->setPreventsCleavage(prevents_cleavage);
    mod->print();

    if (limit == 0) {
      static_mods.push_back(mod);
    } else {
      variable_mods.push_back(mod);
    }
      
    if (spec_text[pos] == '\0') {
      break;
    }
    if (spec_text[pos] == ',') {
      ++pos;
    }
    
  }
}

void set_modspec() {

  vector<AA_MOD_T*> variable_mods;
  vector<AA_MOD_T*> static_mods;
  
  parse_modspec_line(
    Params::GetString("mods-spec"),
    static_mods,
    variable_mods
  );

  //Use the static mods as increases in the amino acid mass  
  for (size_t idx=0;idx<static_mods.size();idx++) {
    AA_MOD_T* static_mod = static_mods[idx];
    double delta_mass = static_mod->getMassChange();
    bool* aa_list = static_mod->getAAList();
    for (int aa_idx=0;aa_idx < AA_LIST_LENGTH;aa_idx++) {
      if (aa_list[aa_idx]) {
        carp(CARP_DEBUG, "increasing %c by %lf", (aa_idx+'A'), delta_mass);
        increase_amino_acid_mass(aa_idx+'A', delta_mass);
      }
    }
  }
  
  vector<AA_MOD_T*> monolink_mods;
  vector<AA_MOD_T*> static_mods2;
  string mono_link_str = Params::GetString("mono-link");
  parse_modspec_line(
    mono_link_str,
    static_mods2,
    monolink_mods
  );
  if (static_mods2.size() != 0) {
    carp(CARP_FATAL, "Cannot define monolink as a static mod:%s", mono_link_str.c_str());
  }

  for (size_t idx = 0;idx < monolink_mods.size();idx++) {
    monolink_mods[idx]->setMonoLink(true);
    variable_mods.push_back(monolink_mods[idx]);
  }

  carp(CARP_DEBUG, "Number of static:%d variable:%d", static_mods.size(), variable_mods.size());
  
  //Fill up expected list of modifications.
  for (size_t idx=0;idx<variable_mods.size();idx++) {
    list_of_mods[idx] = variable_mods[idx];
    list_of_mods[idx]->setModIdx(idx);
  }
  
  num_mods = variable_mods.size();
  list_of_variable_mods = &(list_of_mods[0]);
}


/**
 * \brief Read the paramter file and populate the static parameter
 * list of AA_MODS, inlcuding the list of position mods.
 *
 * Also updates the array of amino_masses.  Dies with an error if the
 * number of mods in the parameter file is greater than MAX_AA_MODS.
 * \returns void
 */
void read_mods_from_file(const char* param_filename) {
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
    carp(CARP_FATAL, "Cannot specify more than one fixed n-terminal modification.");
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
  if (num_c_mods == 0) {
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
  Params::Set("mod-precision", max_precision);

  // close file
  fclose(param_file);
  carp(CARP_DEBUG, "Finished reading mods file");
}

void incrementNumMods() {
  num_mods++;
}

void resetMods() {
  for (int i = 0; i < MAX_AA_MODS; i++) {
    delete list_of_mods[i];
    list_of_mods[i] = new AA_MOD_T(i);
  }
  num_mods = 0;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

