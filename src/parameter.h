/**
 * \file parameter.h
 * $Revision: 1.37 $
 * \brief General parameter handling utilities. All values stored here.

 * \detail MUST declare ALL optional command line parameters and
 * required command line arguments here in initalialize_parameters.
 * Parameters can only be set via parse_cmd_line_into_params_hash()
 * which takes values from the command line and from an optional
 * parameter file (provided on the command line with the --parameter
 * option). Options are checked for correct type and legal values.
 * Exits with usage statement on error.  Parameter values can be
 * retrieved with get_<type>_paramter functions.
 * 
 ****************************************************************************/
#ifndef PARAMETER_FILE_H
#define PARAMETER_FILE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <string>
#include "util/utils.h"
#include "util/crux-utils.h"
#include "io/carp.h"
#include "model/objects.h"

#include "util/modifications.h"

static const int PARAMETER_LENGTH = 1024; 
///< maximum length of parameter name and value in characters
static const int NUM_PARAMS = 512; ///< maximum number of parameters allowed
static const int MAX_LINE_LENGTH = 4096; ///< maximum line length in the parameter file
static const int MILLION = 1000000;
static const int BILLION = 1000000000;
static const int SMALL_BUFFER = 256;
static const int MAX_SET_PARAMS = 256;

// The size of the bins for discretizing the m/z axis of the
// observed spectrum.  For use with monoisotopic mass.
static const FLOAT_T BIN_WIDTH_MONO = 1.0005079;
// The size of the bins for discretizing the m/z axis of the
// observed spectrum.  For use with average mass.
static const FLOAT_T BIN_WIDTH_AVERAGE = 1.0011413;

// Global variables
// NOTE (BF mar-10-09): Could be like mod lists, but will require a
// get_parameter call for each...residue? protein?
extern char* pre_cleavage_list;
extern char* post_cleavage_list;
extern int pre_list_size;
extern int post_list_size;
extern bool pre_for_inclusion;
extern bool post_for_inclusion;

// TODO (BF 1-28-08): these should be private. move to parameter.c
/**
 * Data types of parameters.  Used for checking valid parameter input
 * from user.
 *
 * To add a new parameter type:  
 *  (for NEW types) 1. create enum, 
 *                  2. create array of strings, 
 *                  3. write string-to-type,
 *                  4. write type-to-string
 *  (for ALL types) 5. add to parameter-type enum and strings
 *                  6. write get-type-parameter
 *                  7. write set-type-parameter
 *                  8. add to the check_type_and_bounds
 *
 */
enum parameter_type {
  INT_P,             ///< parameters of type int
  DOUBLE_P,          ///< parameters of type double
  STRING_P,          ///< parameters of type char*
  MASS_TYPE_P,       ///< parameters of type MASS_TYPE_T
  DIGEST_TYPE_P,   ///< parameters of type DIGEST_T
  ENZYME_TYPE_P,     ///< parameters of type ENZYME_T
  //PEPTIDE_TYPE_P,    ///< parameters of type PEPTIDE_TYPE_T
  BOOLEAN_P,         ///< parameters of type bool
  SCORER_TYPE_P,     ///< parameters of type SCORER_TYPE_T
  ION_TYPE_P,        ///< parameters of type ION_TYPE_T
  ALGORITHM_TYPE_P,  ///< parameters of type ALGORITHM_TYPE_T
  HARDKLOR_ALGORITHM_TYPE_P, ///< parameters of type HARDKLOR_ALGORITHM_T
  WINDOW_TYPE_P,     ///< parameters of type WINDOW_TYPE_T
  MEASURE_TYPE_P,    ///< parameters of type MEASURE_TYPE_T
  THRESHOLD_P,       ///< parameters of type THRESHOLD_TYPE_T
  PARSIMONY_TYPE_P,  ///< parameters of type PARSIMONY_TYPE_T
  QUANT_LEVEL_TYPE_P,///< parameters of type QUANT_LEVEL_TYPE_T
  DECOY_TYPE_P,      ///< parameters of type DECOY_TYPE_T
  MASS_FORMAT_P,     ///< parameters of type MASS_FORMAT_T

  NUMBER_PARAMETER_TYPES  ///< leave this last, number of types
};
typedef enum parameter_type PARAMETER_TYPE_T;

/**
 * /brief Initialize parameters to default values.
 *
 * Every required argument for every executable and every option and
 * its default value must be declared here.  Allocates hash tables for
 * holding parameter values.
 */
void initialize_parameters(void);

void parse_custom_enzyme(const std::string& rule_str);

/**
 *
 * parse the parameter file given the filename
 */
void parse_parameter_file(
  const char* parameter_filename ///< the parameter file to be parsed -in
  );

void read_mods_from_file(const char* param_file);

MASS_TYPE_T get_mass_type_parameter(
  const char* name
  );

char get_delimiter_parameter(
  const char* name
  );

DIGEST_T get_digest_type_parameter(
  const char* name
  );

ENZYME_T get_enzyme_type_parameter(
  const char* name
  );

DIGEST_T get_digest_type_parameter(
  const char* name
  );

ENZYME_T get_enzyme_type_parameter(
  const char* name
  );

THRESHOLD_T get_threshold_type_parameter(
  const char* name
  );

QUANT_LEVEL_TYPE_T get_quant_level_type_parameter(
  const char* name
  );

DECOY_TYPE_T get_decoy_type_parameter(
  const char* name
  );

DECOY_TYPE_T get_tide_decoy_type_parameter(
  const char* name
  );

MASS_FORMAT_T get_mass_format_type_parameter(
  const char* name
  );

int get_max_ion_charge_parameter(
  const char* name
  );

double get_mz_bin_width();
 
double get_mz_bin_offset();

COLTYPE_T get_column_type_parameter(
  const char* name
  );

COMPARISON_T get_comparison_parameter(
  const char* name
  );

/**
 * \returns the comet enzyme info lines parsed from the file
 * or generated defaults
 */
const std::vector<std::string>& get_comet_enzyme_info_lines();

/**
 * Prints the parameters.  If lead_string is not null, preprends it to
 * each line.
 */
void print_parameters(
  char* first_line,  ///< the first line to be printed before the parameter list -in
  char* parameter_filename,  ///< the parameter file name -in
  char* lead_string,  ///< the lead string to be printed before each line -in
  FILE* outstream  ///< the output stream -out
  );


/**
 * \brief Get the pointer to the list of AA_MODs requested by the
 * user.  Does not include the c- and n-term mods.  Argument is a
 * reference to an array of pointers.  Return 0 and set mods == NULL
 * if there are no aa_mods.
 * \returns The number of items pointed to by mods
 */
int get_aa_mod_list(AA_MOD_T*** mods);

/**
 * Adds 1 to the number of variable modifications in the
 * parameters
 */
void incrementNumMods();

void resetMods();

/**
 * \brief Get the pointer to the list of AA_MODs for the peptide
 * c-terminus.  Argument is a reference to an array of
 * pointers. Return 0 and set mods==NULL if there are no c-term 
 * mods.
 *
 * \returns The number of items pointed to by mods
 */
int get_c_mod_list(AA_MOD_T*** mods);

/**
 * \brief Get the pointer to the list of AA_MODs for the peptide
 * n-terminus.  Argument is a reference to an array of
 * pointers. Return 0 and set mods==NULL if there are no n-term 
 * mods.
 *
 * \returns The number of items pointed to by mods
 */
int get_n_mod_list(AA_MOD_T*** mods);

/**
 * \brief Get the pointer to the list of AA_MODs requested by the
 * user.  Includes aa_mods, c- and n-term mods.  Argument is a
 * reference to an array of pointers.  Returns 0 and sets mods == NULL
 * if there are no aa_mods.
 * \returns The number of items pointed to by mods
 */
int get_all_aa_mod_list(AA_MOD_T*** mods);

/**
 * \returns The index of the C_TERM or N_TERM fixed modification in
 * the global list of modifications.
 */
int get_fixed_mod_index(MOD_POSITION_T p);

/**
 * \returns the number of fixed terminal modifications: 0, 1, or 2.
 */
int get_num_fixed_mods();

/**
 *  Print modifications to the given file in the format used in the
 *  parameter file.  Expected names are "mod", "cmod", "nmod".  The
 *  function to get the list of modifications should correspond to the
 *  mod name.  (i.e. for "mod", use get_aa_mod_list(), for "cmod" use
 *  get_c_mod_list(), for "nmod" use get_n_mod_list())
 */
void print_mods_parameter_file(std::ostream* param_file, 
                               const char* name,
                               int (*mod_getter)(AA_MOD_T***));


void set_modspec();

#endif
