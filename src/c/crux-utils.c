/**
 * \file crux-utils.c
 * \brief General-use functions for crux
 */


#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include "crux-utils.h"
#include "parameter.h"

/**
 * PRECISION, determines the precision of the compare float, users
 * should lower the number if need more precision
 */
#define PRECISION 0.000000005 

/**
 * the maximum error in terms of Units in the Last Place. 
 * This specifies how big an error we are willing to accept in terms of the value of the least significant 
 * digit of the floating point numbers representation. 
 * MAX_ULPS can also be interpreted in terms of how many representable floats 
 * we are willing to accept between A and B. This function will allow MAX_ULPS-1 floats between A and B.
 */
#define MAX_ULPS 2

/* Functions for converting custom types to and from strings */

/**
 * The string version of isotopic mass type (average, mono)
 */
static char* mass_type_strings[NUMBER_MASS_TYPES] = {"average", "mono"};

BOOLEAN_T string_to_mass_type(char* name, MASS_TYPE_T* result){
  BOOLEAN_T success = TRUE;
  //this is copied from parameter.c::get_peptide_mass_type
  int mass_type = convert_enum_type_str(
                          name, -10, mass_type_strings, NUMBER_MASS_TYPES);

  (*result) = (MASS_TYPE_T)mass_type;

  if( mass_type < 0 ){
    success = FALSE;
  }
  return success;
}

BOOLEAN_T mass_type_to_string(MASS_TYPE_T type, char* type_str){
  BOOLEAN_T success = TRUE;
  if( (int)type > NUMBER_MASS_TYPES ){
    success = FALSE;
    type_str = NULL;
  }

  //  type_str = mass_type_strings[type];
  strcpy(type_str, mass_type_strings[type]);

  return success;
}

// replace peptide type
/**
 * The string versions of digest types
 */
static char* digest_type_strings[NUMBER_DIGEST_TYPES] =
  {"invalid", "full-digest", "partial-digest", "non-specific-digest"};

DIGEST_T string_to_digest_type(char* name){
  int clev_int = convert_enum_type_str(name, -10, 
                                       digest_type_strings, 
                                       NUMBER_DIGEST_TYPES);
  if( clev_int < 0 ){
    clev_int = 0;
  }

  return (DIGEST_T)clev_int;
}

char* digest_type_to_string(DIGEST_T type){
  if( (int)type > NUMBER_DIGEST_TYPES){
    return NULL;
  }

  char* type_str = my_copy_string(digest_type_strings[type]);

  return type_str;
}

/**
 * The string version of enzyme types
 */
static char* enzyme_type_strings[NUMBER_ENZYME_TYPES] = 
  {"invalid", "no-enzyme", "trypsin", "chymotrypsin", "elastase",
   "clostripain", "cyanogen-bromide", "iodosobenzoate", 
   "proline-endopeptidase", "staph-protease", "aspn", 
   "modified-chymotrypsin", "elastase-trypsin-chymotrypsin",
   "custom-enzyme"};

ENZYME_T string_to_enzyme_type(char* name){
  int enz_int = convert_enum_type_str(name, -10, 
                                      enzyme_type_strings, 
                                      NUMBER_ENZYME_TYPES);
  if( enz_int < 0 ){
    enz_int = 0;
  }

  return (ENZYME_T)enz_int;
}

char* enzyme_type_to_string(ENZYME_T type){
  if( (int)type > NUMBER_ENZYME_TYPES){
    return NULL;
  }

  char* type_str = my_copy_string(enzyme_type_strings[type]);

  return type_str;
}


/**
 * The string version of peptide cleavage type
 */
/*
static char* peptide_type_strings[NUMBER_PEPTIDE_TYPES] = 
{"tryptic", "partial", "N_TRYPTIC", "C_TRYPTIC", "NOT_TRYPTIC", "all"};

BOOLEAN_T string_to_peptide_type(char* name, PEPTIDE_TYPE_T* result){

  BOOLEAN_T success = TRUE;
  //this is copied from parameter.c::get_peptide_mass_type
  int pep_type = convert_enum_type_str(
                     name, -10, peptide_type_strings, NUMBER_PEPTIDE_TYPES);
  (*result) = (PEPTIDE_TYPE_T)pep_type;

  if( pep_type < 0 ){
    success = FALSE;
  }
  return success;
}

BOOLEAN_T peptide_type_to_string(PEPTIDE_TYPE_T type, char* type_str){
  BOOLEAN_T success = TRUE;
  if( (int)type > NUMBER_PEPTIDE_TYPES){
    success = FALSE;
    type_str = NULL;
  }

  //type_str = peptide_type_strings[type];
  strcpy(type_str, peptide_type_strings[type]);

  return success;
}
*/
/**
 * The string version of peptide sort types
 */
static char* sort_type_strings[NUMBER_SORT_TYPES] =
  { "none", "mass", "length", "lexical" };

BOOLEAN_T string_to_sort_type(char* name, SORT_TYPE_T* result){
  BOOLEAN_T success = TRUE;

  int sort_type = convert_enum_type_str(
                        name, -10, sort_type_strings, NUMBER_SORT_TYPES);
  (*result) = (SORT_TYPE_T)sort_type;

  if( sort_type < 0){
    success = FALSE;
  }
  return success;
}

BOOLEAN_T sort_type_to_string(SORT_TYPE_T type, 
                              char* type_str){
  BOOLEAN_T success = TRUE;
  if( (int)type > NUMBER_SORT_TYPES ){
    success = FALSE;
    type_str = NULL;
  }
  strcpy(type_str, sort_type_strings[type]);
  return success;
}

/*
 * The string version of ion types
 */
static char* ion_type_strings[NUMBER_ION_TYPES] = {
  "a", "b", "c", "x", "y", "z", "p", "by", "bya", "all" };

BOOLEAN_T string_to_ion_type(char* name, ION_TYPE_T* result){
  BOOLEAN_T success = TRUE;

  int ion_type = convert_enum_type_str(
                       name, -10, ion_type_strings, NUMBER_ION_TYPES);
  (*result) = (ION_TYPE_T)ion_type;

  if( ion_type < 0){
    success = FALSE;
  }
  return success;
}

BOOLEAN_T ion_type_to_string(ION_TYPE_T type,
                             char* type_str){
  BOOLEAN_T success = TRUE;
  if( (int)type > NUMBER_ION_TYPES ){
    success = FALSE;
    type_str = NULL;
  }
  strcpy(type_str, ion_type_strings[type]);
  return success;
}

/*
 * The string version of ALGORITHM_TYPE_T
 */
static char* algorithm_type_strings[NUMBER_ALGORITHM_TYPES] = 
  {"percolator", "rczar", "curve-fit",
   //"qvalue",
   "none", "all"};

BOOLEAN_T string_to_algorithm_type(char* name, ALGORITHM_TYPE_T* result){
  BOOLEAN_T success = TRUE;

  int algorithm_type = convert_enum_type_str(name, -10,
                                             algorithm_type_strings,
                                             NUMBER_ALGORITHM_TYPES);
  (*result) = (ALGORITHM_TYPE_T)algorithm_type;

  if(algorithm_type < 0){
    success = FALSE;
  }
  return success;
}

BOOLEAN_T algorithm_type_to_string(ALGORITHM_TYPE_T type, char* type_str){
  BOOLEAN_T success = TRUE;
  if( (int)type > NUMBER_ALGORITHM_TYPES){
    success = FALSE;
    type_str = NULL;
  }
  strcpy(type_str, algorithm_type_strings[type]);
  return success;
}

/*
 * The string version of SCORER_TYPE_T
 */
static char* scorer_type_strings[NUMBER_SCORER_TYPES] = 
  {"sp", "xcorr", "dotp", "logp_exp_sp", "logp_bonf_exp_sp", 
   "logp_evd_xcorr", "logp_bonf_evd_xcorr", "logp_weibull_sp", 
   "sp-pvalue",  //"sp-logp", 
   "logp_weibull_xcorr", 
   "xcorr-pvalue", //"xcorr-logp", 
   "q_value", "percolator_score", 
   "qvalue"};//"logp_qvalue_weibull_xcorr" };
//TODO: this should probably be changed, these strings are the option args
//Instead could have an if block in string_to_type

BOOLEAN_T string_to_scorer_type(char* name, SCORER_TYPE_T* result){
  BOOLEAN_T success = TRUE;

  int scorer_type = convert_enum_type_str(name, -10, scorer_type_strings,
                                          NUMBER_SCORER_TYPES);
  (*result) = (SCORER_TYPE_T)scorer_type;

  if( scorer_type < 0){
    success = FALSE;
  }
  return success;
}

BOOLEAN_T scorer_type_to_string(SCORER_TYPE_T type, char* type_str){
  BOOLEAN_T success = TRUE;
  if( (int)type > NUMBER_SCORER_TYPES){
    success = FALSE;
    type_str = NULL;
  }
  strcpy(type_str, scorer_type_strings[type]);
  return success;
}


/**
 * returns a heap allocated copy of the src string
 */
char* my_copy_string(char* src){
  if( src == NULL ){
    return NULL;
  }
  int length = strlen(src) +1; // +\0
  char* copy = 
    (char *)mymalloc(sizeof(char)*length);
  return strncpy(copy, src, length);  
}

/**
 * returns copy of the src string upto the specified length
 * includes a null terminating \\0 character
 * the string is heap allocated thus, user must free
 */
char* copy_string_part(char* src, int length){
  char* copy = (char*)mycalloc(length+1, sizeof(char));
  strncpy(copy, src, length);
  copy[length] = '\0';
  return copy;
}

/**
 * \returns the 0 if equal, 1 if float_a is larger, -1 if float_b is larger
 * compare the absolute value of the difference of two numbers with an 
 * appropriate epsilon to get relations.
 * Multiplying the epsilon by sum of the comparands adjusts the comparison 
 * to the range of the numbers, allowing a single epsilon to be used for many, 
 * or perhaps all compares.
 */
inline int compare_float(FLOAT_T float_a, FLOAT_T float_b){
  FLOAT_T EPSILON = PRECISION;
  FLOAT_T sum = float_a + float_b;
  // a == b
  if( fabsf(float_a - float_b) <= fabsf(sum)* EPSILON ){
    return 0;
  }
  // a > b
  else if((float_a - float_b) > fabsf(sum)* EPSILON){
    return 1;
  }
  // a < b
  else{
    return -1;
  }
}

/**
 *\returns TRUE if float_a is between the interaval of min and max, else FALSE
 */
inline BOOLEAN_T compare_float_three(FLOAT_T float_a, FLOAT_T min, FLOAT_T max){
  if(compare_float(float_a, min) == -1 ||
     compare_float(float_a, max) ==  1){
    return FALSE;
  }
  return TRUE;
}

/**
 * parses the filename and path  
 * returns an array A, with A[0] the filename and A[1] the path to the filename
 * returns A[1] NULL if only a filename was passed in
 * ex) ../../file_name => returns filename , ../../
 *     file_name => returns filename, NULL
 *\returns A heap allocated array of both filename and path
 */
char** parse_filename_path(char* file){
  int len = strlen(file);
  int end_idx = len;
  int end_path = -1;  // index of where the last "/" is located
  char* path = NULL;
  char* filename = NULL;
  char** result = (char**)mycalloc(2, sizeof(char*));

  for(; end_idx > 0; --end_idx){
    if(strncmp(&file[end_idx - 1], "/", 1) == 0){
      end_path = end_idx;
      break;
    }
  }
  // copy path, if there is a "/" in the file
  if(end_path != -1){
    path = copy_string_part(file, end_path);
  }
  // copy filename
  filename = copy_string_part(&file[end_idx], len); 
  
  // set result with filename and path
  result[0] = filename;
  result[1] = path;
  
  return result;
}

/**
 * \brief Parses the filename, path, and file extension of given string.
 *  
 * The array returned, A, contains the filename (A[0]) striped of the
 * given file extension and the path (A[1]).  Path is NULL if no / in
 * given name.  Filename is unchanged if it does not include the
 * extension. e.g. Given "../../filname.ext" and ".ext" returns
 * A[0]="filename" A[1]="../../".  Given "filename" returns A[0] =
 * NULL and A[1] = "filename". 
 *
 * \returns A heap allocated array of filename striped of extension and
 * path.
 */
char** parse_filename_path_extension(
       char* file, ///< filename and path to parse -in
       char* extension ///< extension to look for (includes leading .) --in
){

  carp(CARP_DETAILED_DEBUG, "Given path/file %s and ext %s", file, extension);
  char** file_path_array = parse_filename_path(file);
  char* trimmed_filename = file_path_array[0];

  // look for extension
  if( extension != NULL ){

    carp(CARP_DETAILED_DEBUG, "File trimmed of path is %s", trimmed_filename);
    if( ! suffix_compare(trimmed_filename, extension) ){
        return file_path_array;
    }
    int file_len = strlen(trimmed_filename);
    int ext_len = strlen(extension);

    /* compare_suffix replaces the following
    int file_idx = file_len;
    int ext_idx = ext_len;

    if( ext_len > file_len ){
      carp(CARP_ERROR, 
           "Cannot parse file extension.  The file extension %s is longer "
           "than the filename '%s'", extension, file);
      return file_path_array;
    }
    carp(CARP_DETAILED_DEBUG, "Name len %d ext len %d", file_idx, ext_idx);

    //compare name and ext from end of strings backwards
    for(ext_idx = ext_idx; ext_idx > -1; ext_idx--){
      carp(CARP_DETAILED_DEBUG, "Name[%d]='%d', ext[%d]='%d'", 
           file_idx, trimmed_filename[file_idx], ext_idx, extension[ext_idx]);
      // if they stop matching, don't change filename
      if( extension[ext_idx] != trimmed_filename[file_idx--]){
        return file_path_array;
      }
    }
    */
    // after comparing the whole extension, it matched
    trimmed_filename[file_len - ext_len] = '\0';
    carp(CARP_DETAILED_DEBUG, "Final trimmed filename %s", trimmed_filename);

  }

  return file_path_array;
}

/**
 * parses the filename
 * ex) ../../file_name => returns filename
 *\returns A heap allocated array of filename
 */
char* parse_filename(char* file){
  int len = strlen(file);
  int end_idx = len;
  int end_path = -1;  // index of where the last "/" is located
  char* filename = NULL;
  
  for(; end_idx > 0; --end_idx){
    if(strncmp(&file[end_idx - 1], "/", 1) == 0){
      end_path = end_idx;
      break;
    }
  }
  
  // copy filename
  filename = copy_string_part(&file[end_idx], len); 
  
  return filename;
}



/**
 * convert the integer into a string
 * \returns a heap allocated string
 */
char* int_to_char(unsigned int i){
  unsigned int digit = i / 10;
  char* int_string = (char*)mycalloc(digit+2, sizeof(char));
  sprintf(int_string, "%d", i);
  return int_string;
}
 
/**
 * convert the integer into a string
 * \returns a heap allocated string
 */
char* signed_int_to_char(int i){
  int digit = abs(i)/ 10;
  char* int_string = (char*)mycalloc(digit+2, sizeof(char));
  sprintf(int_string, "%d", i);
  return int_string;
}
/**
 * Gives the peptide type as defined by the string
 * Returns FALSE if the string is not a valid type
 */
//BOOLEAN_T string_to_peptide_type

/**
 *prints the peptide type given it's enum value
 */
/*
void print_peptide_type(PEPTIDE_TYPE_T peptide_type, FILE* file){
  if(peptide_type == TRYPTIC){
    fprintf(file, "%s", "TRYPTIC");
  }
  else if(peptide_type == PARTIALLY_TRYPTIC){
    fprintf(file, "%s", "PARTIALLY_TRYPTIC");
  }
  else if(peptide_type == N_TRYPTIC){
    fprintf(file, "%s", "N_TRYPTIC");
  }
  else if(peptide_type == C_TRYPTIC){
    fprintf(file, "%s", "C_TRYPTIC");
  }
  else if(peptide_type == NOT_TRYPTIC){
    fprintf(file, "%s", "NOT_TRYPTIC");
  }
  else if(peptide_type == ANY_TRYPTIC){
    fprintf(file, "%s", "ANY_TRYPTIC");
  }
}
*/
/**
 * given two strings return a concatenated third string
 * \returns a heap allocated string that concatenates the two inputs
 */
char* cat_string(char* string_one, char* string_two){
  int len_one = strlen(string_one);
  int len_two = strlen(string_two);
  
  char* result = (char*)mycalloc(len_one + len_two + 1, sizeof(char));
  strncpy(result, string_one, len_one);
  strncpy(&result[len_one], string_two, len_two);
  return result;
}

/**
 * Adds the fileroot parameter to a string as a prefix.
 * Given a pointer to pointer to a string, if the fileroot parameter is set
 * the memory for the string is reallocated, and the fileroot string
 * is added as a prefix.
 */
void prefix_fileroot_to_name(char** name) {
  char* fileroot = get_string_parameter("fileroot");
  if (fileroot != NULL) {
    int len_name = strlen(*name);
    int len_root = strlen(fileroot);
    *name = myrealloc(*name, len_root + len_name + 2);
    memmove(*name + len_root + 1, *name, len_name + 1);
    strcpy(*name, fileroot);
    (*name)[len_root] = '.';
    free(fileroot);
  };
}

/**
 * \brief Check if the string has the correct prefix
 * \returns TRUE if the string starts with the given prefix, else FALSE
 */
BOOLEAN_T prefix_compare(
  char* string, ///< The string to check
  char* prefix  ///< The prefix to find in the string
  )
{
  int len = strlen(string);
  int len_prefix = strlen(prefix);

  if(len_prefix > len){
    return FALSE;
  }
  
  if(strncmp(string, prefix, len_prefix) == 0){
    return TRUE;
  }
  
  return FALSE;
}

/**
 * \brief Check if the string has the correct suffix
 * \returns TRUE if the end of the string matches the given suffix, else FALSE
 */
BOOLEAN_T suffix_compare(
  char* string, ///< The string to check
  char* suffix  ///< The suffix to find in the string
  )
{
    int string_len = strlen(string);
    int suffix_len = strlen(suffix);
    int string_idx = string_len;
    int suffix_idx = suffix_len;

    if( suffix_len > string_len ){
      return FALSE;
    }

    //compare name and ext from end of strings backwards
    for(suffix_idx = suffix_idx; suffix_idx > -1; suffix_idx--){
      //carp(CARP_DETAILED_DEBUG, "Name[%d]='%d', ext[%d]='%d'", 
      //   string_idx, string[string_idx], suffix_idx, suffix[suffix_idx]);
      // if they stop matching, don't change filename
      if( suffix[suffix_idx] != string[string_idx--]){
        return FALSE;
      }
    }

  return TRUE;
}

/**
 * Given the path and the filename return a file with path
 * "path/filename".  Returns filename unchanged if path = NULL.
 * \returns a heap allocated string, "path/filename"
 */
char* get_full_filename(char* path, char* filename){

  char* result = NULL;
  if( path == NULL ){
    result = my_copy_string(filename);
  }else{
    // TODO (BF 26-Feb-08) don't add second / if path already ends in /
    char* ready_path = cat_string(path, "/");
    result = cat_string(ready_path, filename);
    free(ready_path);
  }

  return result;
}

/**
 * \brief Decide if a file name is a decoy csm file
 * \returns TRUE if name ends in -decoy-#.csm 
 */
BOOLEAN_T name_is_decoy(char* name){
  char* name_end = strrchr(name, '\0');
  char* last_d = strrchr(name, 'd');
  int decoy_len = strlen("decoy-1.csm");

  carp(CARP_DEBUG, "Name is %s", name);
  // where is the last d in the name
  if( last_d == NULL ){
    carp(CARP_DEBUG, "Found no 'd'");
    return FALSE;
  }
  if( (name_end - last_d) == decoy_len 
      && *(last_d-1) == '.'
      && strncmp(last_d, "decoy-", 6)==0 ){
    carp(CARP_DEBUG, "Name is a decoy file");
    return TRUE;
  } 

  carp(CARP_DEBUG, "Found a 'd' but name is not a decoy file");
  return FALSE;
}

/**
 * returns the file size of the given filename
 */
long get_filesize(char *FileName){
    struct stat file;
    // return file size
    if(!stat(FileName,&file)){
      return file.st_size;
    }
    return 0;
}

/**
 * \brief A function for creating a directory to hold output files from crux.
 * 
 * Tries to create a the named directory for use as the output directory for crux.
 * If the overwrite option is true, an existing directory wtih that
 * name will not cause an error. 
 * 
 * \returns 0 if successful, -1 if an error occured.
*/
int create_output_directory(
  char *output_folder, // Name of output folder.
  BOOLEAN_T overwrite  // Whether or not to overwrite an existing dir 
) 
{

  int result = -1;
  BOOLEAN_T path_is_directory = FALSE;
  BOOLEAN_T path_exists = FALSE;
  struct stat stat_buffer;

  // Does the output directory alredy exist?
  if (stat(output_folder, &stat_buffer)) {
    if (errno == ENOENT) {
      // stat failed because the path doesn't exist.
      path_exists = FALSE;
      path_is_directory = FALSE;
    }
    else {
      // stat failed for some other reason
      carp(
        CARP_ERROR,
        "Unable to check for status of output directory '%s': %s.\n",
        output_folder,
        strerror(errno)
      );
      result = -1;
    }
  }
  else {
    path_exists = TRUE;
    path_is_directory = S_ISDIR(stat_buffer.st_mode);
  }

  if (path_exists) {
    if (!path_is_directory) {
      carp(
        CARP_ERROR,
        "A non-directory file named '%s' already exists,\n"
        "so that name can't be used for an output directory.\n",
        output_folder
      );
      result = -1;
    }
    else {
      if (!overwrite) {
        carp(
          CARP_WARNING,
          "The output directory '%s' already exists.\nExisting files will not"
          " be overwritten.\n",
          output_folder
        );
        result = -1;
      }
      else {
        carp(
          CARP_WARNING,
          "The output directory '%s' already exists.\nExisting files will"
          " be overwritten.\n",
          output_folder
        );
        result = 0;
      }
    }
  }
  else {
    // The directory doesn't exist, so we can create it.
    // Does this accomodate the case where one or more of the
    // parent directories doesn't exit?
    int dir_access = S_IRWXU + S_IRWXG + S_IRWXO;
    if (mkdir(output_folder, dir_access)) {
      // mkdir failed
      carp(
        CARP_ERROR,
        "Unable to create output directory '%s': %s.\n",
        output_folder,
        strerror(errno)
      );
      result = -1;
    }
    else {
      result = 0;
      carp(
        CARP_INFO,
        "Writing results to output directory '%s'.",
        output_folder
      );
    }
  }
  return result;
} 

/**
 * returns whether the given filename is a directory.
 * Returns TRUE if a directory, FALSE otherwise.
 * Terminates program if unable to determine status of file.
 */
BOOLEAN_T is_directory(char *FileName){
    struct stat file;
    if(stat(FileName,&file) == 0){
      // return directory status
      return S_ISDIR(file.st_mode);
    }
    else {
      char *error = strerror(errno);
      carp(
        CARP_FATAL, 
        "stat failed. Unable to determine status of %s. Error: %s.", 
        FileName,
        error
      );
      return FALSE; // Avoid compiler warning
    }
}

/**
 * deletes a given directory and it's files inside.
 * assumes that there's no sub directories, only files
 * \returns TRUE if successfully deleted directory
 */
BOOLEAN_T delete_dir(char* dir) {
  struct dirent **namelist =NULL;
  int num_file =0;
  int result;
  char* cwd = getcwd(NULL, 0); //gnu lib

  // does the directory to remove exist?, if so move into it..
  if(chdir(dir) == -1){
    carp(CARP_DETAILED_DEBUG, "Could not find directory '%s' to remove", dir);
    return FALSE;
  }

  // collect all files in dir
  num_file = scandir(".", &namelist, 0, alphasort);

  // delete all files in temp dir
  while(num_file--){
    remove(namelist[num_file]->d_name);
    free(namelist[num_file]);
  }
  free(namelist);

  //chdir(".."); // assumes the directory to delete is in cwd
  chdir(cwd);
  result = rmdir(dir);
  if(result == FALSE){
    return FALSE;
  }
  
  return TRUE;
}

/**
 * \brief Take a filename, strip its leading path information (if
 * any) and file extension (if any).  Add a new path (if given) and a
 * new suffix (exension).
 *
 * If given ../dir/filename.ext, .new-ext, .ext, otherdir would return
 * otherdir/filename.new-ext 
 * \returns A heap allocated filename
 */
char* generate_name_path(
  char* filename,
  char* old_suffix,
  char* new_suffix,
  char* new_path
  ){

  carp(CARP_DEBUG, "Generate name given filename '%s', old suffix '%s', " \
       "new suffix '%s', new path '%s'", 
       filename, old_suffix, new_suffix, new_path);

  // parse path, filename, extension
  char** name_path = parse_filename_path_extension(filename, old_suffix);

  // add the extension
  char* new_name = cat_string(name_path[0], new_suffix);

  // add path to new filename
  char* full_new_name = get_full_filename(new_path, new_name);

  // cleanup
  free(name_path[0]);
  free(name_path[1]);
  free(new_name);

  carp(CARP_DEBUG, "Final name is '%s'", full_new_name);
  return full_new_name;
}
/**
 * given a fasta_file name it returns a name with the name_tag add to the end
 * Suffix may be NULL
 * format: suffix_myfasta_nameTag
 * \returns A heap allocated file name of the given fasta file
 */
char* generate_name(
  char* fasta_filename,
  char* name_tag,
  char* file_extension,
  char* suffix
  )
{
  int len = strlen(fasta_filename);
  int end_path = len;  // index of where the "." is located in the file
  char* name = NULL;
  char* after_suffix = NULL;
  int suffix_length = 0;
  char** file_n_path = NULL;
  int length = 0;

  // cut off the file extension if needed
  int end_idx;
  for(end_idx = len; end_idx > 0; --end_idx){
    if(strcmp(&fasta_filename[end_idx - 1], file_extension) == 0){
      end_path = end_idx - 1;
      break;
    }
  }
  
  // check suffix
  if(suffix != NULL){
    suffix_length = strlen(suffix);
    file_n_path = parse_filename_path(fasta_filename);
  }

  name = (char*)mycalloc(
      suffix_length + end_path + strlen(name_tag) + 1, sizeof(char));
  after_suffix = name;
  
  // if suffix exit add to top
  if(suffix_length != 0){
    length = strlen(file_n_path[1]);
    if(file_n_path[1] != NULL){
      strncpy(name, file_n_path[1], length);
      after_suffix = &name[length];
    }    
    strncpy(after_suffix, suffix, suffix_length);
    after_suffix = &after_suffix[suffix_length];
    
    length = strlen(file_n_path[0]);
    
    strncpy(after_suffix, file_n_path[0], (length - (len-end_path)));
    after_suffix = &after_suffix[(length - (len-end_path))];
    
    free(file_n_path[0]);
    free(file_n_path[1]);
    free(file_n_path);
  }
  else{
    strncpy(after_suffix, fasta_filename, end_path);
  }
  
  strcat(after_suffix, name_tag);
  return name;
}

/**
 * \brief Create the correct filename for a binary psm file, 
 * search.target.csm for target search and search.decoy-#.csm for 
 * decoy searches.
 *
 * Adds the appropriate
 * extension depending on the file index (0=target, 1=first decoy,
 * 2=second decoy, etc).
 * \returns A heap allocated char* with the new filename.
 */
char* generate_psm_filename(int file_index) {///< target/decoy index -in
  carp(CARP_DEBUG, "Given index %d", file_index);

  char* fullname = mymalloc(sizeof(char) * 30);
  if( file_index == 0 ){
    sprintf(fullname, "search.target.csm");
  }else{
    sprintf(fullname, "search.decoy-%i.csm", file_index);
  }
  prefix_fileroot_to_name(&fullname);

  return fullname;

}

/**
 * checks if each AA is an AA
 *\returns TRUE if sequence is valid else, FALSE
 */
BOOLEAN_T valid_peptide_sequence( char* sequence){
  // iterate over all AA and check if with in range
  while(sequence[0] != '\0'){
    if(sequence[0] < 65 || sequence[0] > 90 ){
      return FALSE;
    }
    ++sequence;
  }
  return TRUE;
}

/**
 * \brief Open and create a file of the given name in the given
 * directory.
 *
 * Assumes the directory exists.  Fails if file can't be opened or if
 * file exists and overwrite is false.
 *\returns A file handle to the newly created file.
 */
FILE* create_file_in_path(
  char* filename,  ///< the filename to create & open -in
  char* directory,  ///< the directory to open the file in -in
  BOOLEAN_T overwrite  ///< replace file (T) or die if exists (F)
  )
{
  char* file_full_path = get_full_filename(directory, filename);
  // FIXME CEG consider using stat instead
  FILE* file = fopen(file_full_path, "r"); //to test if file exists
  if( file != NULL ){  
    //The file exists, are we allowed to overwrite it?
    fclose(file);
    file = NULL;
    if( ! overwrite ){
        // Not allowed to overwrite, we must die.
        carp(
          CARP_FATAL, 
          "The file '%s' already exists and cannot be overwritten. " \
            "Use --overwrite T to replace or choose a different output file name",
          file_full_path
        );
    }
    else {
      // Allowed to overwrite, send warning message.
      carp(
        CARP_WARNING, 
        "The file '%s' already exists and will be overwriten.",
        file_full_path
      );
    }
  }
  
  file = fopen(file_full_path, "w+"); //read and write, replace existing

  if(file == NULL){
    carp(CARP_FATAL, "Failed to create and open file: %s", file_full_path);
  }
  
  free(file_full_path);

  return file;
}

/**
 *\returns a heap allocated feature name array for the algorithm type
 */
char** generate_feature_name_array(
  ALGORITHM_TYPE_T algorithm ///< the algorithm's feature name to produce -in
)
{
  char** name_array = NULL;

  switch(algorithm){
    case PERCOLATOR_ALGORITHM:
    case RCZAR_ALGORITHM:
    case QVALUE_ALGORITHM:
    case ALL_ALGORITHM:
    case NO_ALGORITHM:
      name_array = (char**)mycalloc(20, sizeof(char *));
      name_array[0] =  my_copy_string("XCorr");
      name_array[1] =  my_copy_string("DeltCN");
      name_array[2] =  my_copy_string("DeltLCN");
      name_array[3] =  my_copy_string("Sp");
      name_array[4] =  my_copy_string("lnrSp");
      name_array[5] =  my_copy_string("dM");
      name_array[6] =  my_copy_string("absdM");
      name_array[7] =  my_copy_string("Mass");
      name_array[8] =  my_copy_string("ionFrac");
      name_array[9] =  my_copy_string("lnSM");
      name_array[10] =  my_copy_string("enzN");
      name_array[11] =  my_copy_string("enzC");
      name_array[12] =  my_copy_string("enzInt");
      name_array[13] =  my_copy_string("pepLen");
      name_array[14] =  my_copy_string("charge1");
      name_array[15] =  my_copy_string("charge2");
      name_array[16] =  my_copy_string("charge3");
      name_array[17] =  my_copy_string("numPep");
      name_array[18] =  my_copy_string("numProt");
      name_array[19] =  my_copy_string("pepSite");
  }
  
  return name_array;
}

/**
 * User define our upper and our lower bounds.
 * The random number will always be 
 * between low and high, inclusive.
 * There is no seeding in this function, user must do it for themselves
 *\returns a random number between the interval user provides
 */
int get_random_number_interval(
  int low, ///< the number for lower bound -in
  int high ///< the number for higher bound -in
  )
{  
  return (rand() % (high - low + 1) + low);
}

/**
 *\returns the number of digits in the number
 */
int get_number_digits(
  int number ///< the number to count digits
  )
{
  int idx = 0;
  for(; number >= 10; ++idx){
    number = number/10;    
  }

  return ++idx;
}

void swap_quick(
  FLOAT_T* a,
  int idx,
  int jdx
  )
{
  FLOAT_T temp = 0;
  temp = a[idx];
  a[idx] = a[jdx];
  a[jdx] = temp;
}
 
int Random(int i, int j) {
  return i + rand() % (j-i+1);
}

void quick_sort(FLOAT_T a[], int left, int right) {
  int last = left, i;

  if (left >= right) return;
  
  swap_quick(a,left,Random(left,right));
  for (i = left + 1; i <= right; i++)
    if (a[i] > a[left]) /// CHECK THIS!!
      swap_quick(a,++last,i);
  swap_quick(a,left,last);
  quick_sort(a,left,last-1);
  quick_sort(a,last+1,right);
}

void quicksort(FLOAT_T a[], int array_size){
  quick_sort(a, 0, array_size-1);
}

/**
 * \brief Shuffle an array of floats.  Uses the Knuth algorithm.  Uses
 * get_random_number_interval() to generate random numbers. 
 */
void shuffle_floats(FLOAT_T* array, int size){
  if( array == NULL ){
    carp(CARP_ERROR, "Cannot shuffle NULL array.");
    return;
  }

  int idx, switch_idx;
  int last_element_idx = size - 1;
  FLOAT_T temp_value;
  for(idx=0; idx < size; idx++){
    switch_idx = get_random_number_interval(idx, last_element_idx);
    temp_value = array[idx];
    array[idx] = array[switch_idx];
    array[switch_idx] = temp_value;
  }
}

/**
 * \brief Comparison function for reverse sorting floats.
 * \returns -1,0,1 if a is <,=,> b
 */
int compare_floats_descending(const void* a, const void* b){

  FLOAT_T diff = ( *(float*)b - *(float*)a);
  if( diff < 0 ){
    return -1;
  }else if( diff > 0 ){
    return 1;
  }else{
    return 0;
  }
  
}


/**
 * Fits a three-parameter Weibull distribution to the input data. 
 * Implementation of Weibull distribution parameter estimation from 
 * http:// www.chinarel.com/onlincebook/LifeDataWeb/rank_regression_on_y.htm
 * \returns eta, beta, c (which in this case is the amount the data should
 * be shifted by) and the best correlation coefficient
 */
void fit_three_parameter_weibull(
    FLOAT_T* data, ///< the data to be fit -in
    int fit_data_points, ///< the number of data points to fit -in
    int total_data_points, ///< the total number of data points to fit -in
    FLOAT_T min_shift, ///< the minimum shift to allow -in
    FLOAT_T max_shift, ///< the maximum shift to allow -in
    FLOAT_T step,      ///< step for shift -in
    FLOAT_T* eta,      ///< the eta parameter of the Weibull dist -out
    FLOAT_T* beta,      ///< the beta parameter of the Weibull dist -out
    FLOAT_T* shift,     ///< the best shift -out
    FLOAT_T* correlation   ///< the best correlation -out
    ){
  
  FLOAT_T correlation_tolerance = 0.1;
  
  FLOAT_T best_eta = 0.0;
  FLOAT_T best_beta = 0.0;
  FLOAT_T best_shift = 0.0;
  FLOAT_T best_correlation = 0.0;

  FLOAT_T cur_eta = 0.0;
  FLOAT_T cur_beta = 0.0;
  FLOAT_T cur_correlation = 0.0;
  FLOAT_T cur_shift;

  // FIXME remove below
  int idx;
  for (idx=0; idx<fit_data_points; idx++){
    carp(CARP_DETAILED_DEBUG, "X[%i]=%.6f", idx, data[idx]);
  }

  for (cur_shift = max_shift; cur_shift > min_shift ; cur_shift -= step){

    fit_two_parameter_weibull(data, fit_data_points, total_data_points, 
        cur_shift, &cur_eta, &cur_beta, &cur_correlation);

    if (cur_correlation > best_correlation){
      *eta = best_eta = cur_eta;
      *beta = best_beta = cur_beta;
      *shift = best_shift = cur_shift;
      *correlation = best_correlation = cur_correlation;
    } else if (cur_correlation < best_correlation - correlation_tolerance){
      *eta = best_eta;
      *beta = best_beta;
      *shift = best_shift;
      *correlation = best_correlation;
      carp(CARP_DETAILED_DEBUG, "Stat: Mu, Corr = %.6f, %.6f\n", cur_shift, cur_correlation);
      carp(CARP_DETAILED_DEBUG, "Stat: Eta, Beta, Shift = %.6f, %.6f, %.6f", 
          best_eta, best_beta, best_shift);
      return;
    }
  }
}

/**
 * Fits a two-parameter Weibull distribution to the input data. 
 *
 * Called by the three parameter weibull fitting function to see if
 * the proposed shift gives the best correlation.  If there are too
 * few data points, sets correlation to 0 (minimum value).
 * http:// www.chinarel.com/onlincebook/LifeDataWeb/rank_regression_on_y.htm
 * \returns eta, beta and the correlation coefficient.
 */
void fit_two_parameter_weibull(
    FLOAT_T* data, ///< the data to be fit. should be in descending order -in
    int fit_data_points, ///< the number of data points to fit -in
    int total_data_points, ///< the total number of data points -in
    FLOAT_T shift, ///< the amount by which to shift our data -in
    FLOAT_T* eta,      ///< the eta parameter of the Weibull dist -out
    FLOAT_T* beta,      ///< the beta parameter of the Weibull dist -out
    FLOAT_T* correlation ///< the best correlation -out
    ){

  FLOAT_T* X = calloc(sizeof(float) , total_data_points); //hold data here

  // transform data into an array of values for fitting
  // shift (including only non-neg data values) and take log
  int idx;
  for(idx=0; idx < fit_data_points; idx++){
    FLOAT_T score = data[idx] + shift; // move right by shift
    if (score <= 0.0){
      carp(CARP_DEBUG, "Reached negative score at idx %i", idx);
      fit_data_points = idx;
      break;
    } 
    X[idx] = log(score);
    // carp(CARP_DEBUG, "X[%i]=%.6f=ln(%.6f)", idx, X[idx], score);
  }

  FLOAT_T* F_T = mymalloc(sizeof(float) * total_data_points);
  for(idx=0; idx < fit_data_points; idx++){
    int reverse_idx = total_data_points - idx;
    // magic numbers 0.3 and 0.4 are never changed
    F_T[idx] = (reverse_idx - 0.3) / (total_data_points + 0.4);
    //carp(CARP_DEBUG, "F[%i]=%.6f", idx, F_T[idx]);
  }

  FLOAT_T* Y   = mymalloc(sizeof(float) * total_data_points);
  for(idx=0; idx < fit_data_points; idx++){
    Y[idx] = log( -log(1.0 - F_T[idx]) );
    //carp(CARP_DEBUG, "Y[%i]=%.6f", idx, Y[idx]);
  }

  int N = fit_data_points; // rename for formula's sake
  FLOAT_T sum_Y  = 0.0;
  FLOAT_T sum_X  = 0.0;
  FLOAT_T sum_XY = 0.0;
  FLOAT_T sum_XX = 0.0;
  for(idx=0; idx < fit_data_points; idx++){
    sum_Y  += Y[idx];
    sum_X  += X[idx];
    sum_XX += X[idx] * X[idx];
    sum_XY += X[idx] * Y[idx];
  }
  carp(CARP_DETAILED_DEBUG, "sum_X=%.6f", sum_X);
  carp(CARP_DETAILED_DEBUG, "sum_Y=%.6f", sum_Y);
  carp(CARP_DETAILED_DEBUG, "sum_XX=%.6f", sum_XX);
  carp(CARP_DETAILED_DEBUG, "sum_XY=%.6f", sum_XY);

  FLOAT_T b_num    = sum_XY - (sum_X * sum_Y / N);
  carp(CARP_DETAILED_DEBUG, "b_num=%.6f", b_num);
  FLOAT_T b_denom  = sum_XX - sum_X * sum_X / N;
  carp(CARP_DETAILED_DEBUG, "b_denom=%.6f", b_denom);
  FLOAT_T b_hat    = b_num / b_denom;

  FLOAT_T a_hat    = (sum_Y - b_hat * sum_X) / N;
  *beta = b_hat;
  *eta  = exp( - a_hat / *beta );

  FLOAT_T c_num   = 0.0;
  FLOAT_T c_denom_X = 0.0;
  FLOAT_T c_denom_Y = 0.0;
  FLOAT_T mean_X = sum_X / N;
  FLOAT_T mean_Y = sum_Y / N;
  for (idx=0; idx < N; idx++){
    FLOAT_T X_delta = X[idx] - mean_X; 
    FLOAT_T Y_delta = Y[idx] - mean_Y;
    c_num += X_delta * Y_delta;
    c_denom_X += X_delta * X_delta;
    c_denom_Y += Y_delta * Y_delta;
  }
  FLOAT_T c_denom = sqrt(c_denom_X * c_denom_Y);
  if (c_denom == 0.0){
    //carp(CARP_FATAL, "Zero denominator in correlation calculation!");
    carp(CARP_DETAILED_DEBUG, "Zero denominator in correlation calculation!");
    *correlation = 0.0; // min value
    *eta = 0;
    *beta = 0;
  }
  *correlation = c_num / c_denom;

  carp(CARP_DETAILED_DEBUG, "eta=%.6f", *eta);
  carp(CARP_DETAILED_DEBUG, "beta=%.6f", *beta);
  carp(CARP_DETAILED_DEBUG, "correlation=%.6f", *correlation);

  free(F_T);
  free(Y);
  free(X);
}




