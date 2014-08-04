/****************************************************************//**
 * \file utils.cpp
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 9-8-97
 * PROJECT: shared
 * COPYRIGHT: 1997-2001 Columbia University
 * DESCRIPTION: Various useful generic utilities.
 ********************************************************************/
#include <stdarg.h>
#include <string.h>
#ifdef _MSC_VER
#include <Winsock2.h>
#else
#include <sys/time.h>
#include <unistd.h> 
#endif
#include <assert.h>
#include <errno.h>
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/uniform_int_distribution.hpp"
#include "utils.h"
#include "carp.h"
#include "WinCrux.h"

#ifdef INCLUDE_GETLINE

/*********************************************************
 This function replaces the GNU extension of the same name.
 Reads a line from the given stream.
 *********************************************************/
int getline(char **lineptr, size_t *n, FILE *stream) {

  const size_t BUFFSIZE = 100;

  // Check the input values.
  if (lineptr == NULL || stream == NULL) {
    errno = EINVAL;
    return -1;
  }

  // Read the first character from the stream.
  size_t index = 0;
  int c = fgetc(stream);
  int e = ferror(stream);
  if (c == EOF || e) {
    return -1;
  }

  // Allocate an buffer if needed.
  if (*lineptr == NULL) {
    *lineptr = (char *) mymalloc((*n + BUFFSIZE) * sizeof(char));
    *n += BUFFSIZE;
  }

  // Copy from the stream to the buffer until we find a line end.
  while(c != '\n' && c != EOF && !e) {
    (*lineptr)[index++] = c;
    if (index > *n - 1) {
      // Out of space, expand the buffer. 
      *lineptr = (char *) myrealloc(*lineptr, *n + BUFFSIZE);
      *n += BUFFSIZE;
    }
    c = fgetc(stream);
    e = ferror(stream);
  }

  // Reached end of line, end of file, or read error.
  if (!e) {
    
    if (c != EOF) {
      (*lineptr)[index++] = c;
      if (index > (*n - 1)) {
        *lineptr = (char *) myrealloc(*lineptr, *n + 1);
        (*n)++;
      }
    }

    // Terminate the string.
    (*lineptr)[index] = 0;
    
    // Return the length of the string
    // without the terminating null.
    return index;
  }
  else {
    // Some sort of read error
    return -1;
  }
}

#endif

/***********************************************************************
 * Return the value to replace a missing value -- NaN.
 ***********************************************************************/
double NaN
  (void)
{
  return atof("NaN");
}

/***********************************************************************
 * Return elapsed time in microseconds since the last call.
 ***********************************************************************/
double wall_clock(){
  static int first_call = 1;
  static double first_time;
  double t;
  struct timeval tp;

  gettimeofday(&tp, NULL);
  if(first_call == 1){
    first_time = (1E6*((double)tp.tv_sec)+
                 ((double)tp.tv_usec));
    t = first_time;
    first_call = 0;
  } else {
    t = (1E6*((double)tp.tv_sec)+
        ((double)tp.tv_usec));
    t = t - first_time;
  }
  return (double) t;
}

/************************************************************************
 * See .h file for description.
 ************************************************************************/
bool open_file
  (const char *    filename,      /* Name of the file to be opened. */
   const char *    file_mode,     /* Mode to be passed to fopen. */
   bool allow_stdin,         /* If true, filename "-" is stdin. */
   const char *    file_description,   
   const char *    content_description,
   FILE **         afile)         /* Pointer to the open file. */
{
  if (filename == NULL) {
    carp(CARP_ERROR, "No %s filename specified.\n", file_description);
    return(false);
  } else if ((allow_stdin) && (strcmp(filename, "-") == 0)) {
    if (strchr(file_mode, 'r') != NULL) {
        carp(CARP_INFO, "Reading %s from stdin.\n", content_description);
      *afile = stdin;
    } else if (strchr(file_mode, 'w') != NULL) {
        carp(CARP_INFO, "Writing %s to stdout.\n", content_description);
      *afile = stdout;
    } else {
      carp(CARP_INFO, "Sorry, I can't figure out whether to use stdin ");
      carp(CARP_INFO, "or stdout for %s.\n", content_description);
      return(false);
    }
  } else if ((*afile = fopen(filename, file_mode)) == NULL) {
    carp(CARP_INFO, "Error opening file %s.\n", filename);
    return(false);
  }
  return(true);
}

/********************************************************************
 * void mymalloc, mycalloc, myrealloc
 * 
 * See .h file for descriptions.
 ********************************************************************/
void *mymalloc
  (size_t size)
{
  void * temp_ptr;

  if (size == 0)
    size++;

  temp_ptr = malloc(size);

  if (temp_ptr == NULL) {
    carp(CARP_FATAL, "Memory exhausted.  Cannot allocate %d bytes.", (int)size);
  }
    

  return(temp_ptr);
}

void *mycalloc
  (size_t nelem,
   size_t size)
{
  void * temp_ptr;

  /* Make sure we allocate something. */
  if (size == 0) {
    size = 1;
  }
  if (nelem == 0) {
    nelem = 1;
  }

  temp_ptr = calloc(nelem, size);

  if (temp_ptr == NULL) {
    carp(CARP_FATAL, "Memory exhausted.  Cannot allocate %d bytes.", (int)size);
  }

  return(temp_ptr);
}

void * myrealloc
  (void * ptr,
   size_t  size)
{
  void * temp_ptr;

  /* Make sure we allocate something. */
  if (size == 0)
    size = 1;
  assert(size > 0);

  /* Some non-ANSI systems complain about reallocating NULL pointers. */
  if (ptr == NULL) {
    temp_ptr = malloc(size);
  } else {
    temp_ptr = realloc(ptr, size);
  }

  if (temp_ptr == NULL) {
    carp(CARP_FATAL, "Memory exhausted.  Cannot allocate %d bytes.", (int)size);
  }

  return(temp_ptr);
}

/********************************************************************
 * fwrite with a check to make sure it was successful (useful for NFS problems)
 ********************************************************************/
bool myfwrite
  (const void *ptr, 
   size_t size, 
   size_t nitems, 
   FILE *stream){

  size_t ret = fwrite(ptr, size, nitems, stream);
  if (nitems != ret){
    carp(CARP_ERROR, "Problem writing %i items", nitems);
    return false;
  }
  return true;
}

/**********************************************************************
 * Compute a logarithm.
 **********************************************************************/
PROB_T my_log
  (PROB_T x)
{
  if (x > 0.0) {
    return(LOG_VALUE(log(x)));
  } else if (x < 0.0) {
    carp(CARP_FATAL, "Tried to take the log of a negative value (%g).", x);
  } /* else if (x == 0.0) */
  return(LOG_ZERO);
}

/* The lookup table. */
static const int LOG_PRECISION = 100000;
static PROB_T log_table[(int) LOG_PRECISION + 2];

/**********************************************************************
 * Set up lookup table for log(x), 0 <= x <= 1.
 **********************************************************************/
void init_log_prob
  (void)
{
  int    i_table;
  PROB_T table_value;

  log_table[0] = LOG_ZERO;
  for (i_table = 1; i_table <= LOG_PRECISION; i_table++) {
    table_value = (double)(i_table / LOG_PRECISION);
    log_table[i_table] = log(table_value);
  }
  log_table[i_table] = 0;  /* For use in iterpolation when x=1 */
}

/**********************************************************************
 * Efficiently find log(x), when 0 < x <= 1.  Doesn't check bounds.
 **********************************************************************/
PROB_T log_prob
  (PROB_T value)
{
  const PROB_T scaled_value = value * LOG_PRECISION;
  const int    log_index = (int)scaled_value;
  const PROB_T decimal_part = scaled_value - log_index;
  const PROB_T lower_value = log_table[log_index];
  const PROB_T upper_value = log_table[log_index+1];
  const PROB_T interpolation = decimal_part * (lower_value - upper_value);

  if (value == 0.0) {
    return(LOG_ZERO);
  }
  return(lower_value + interpolation);
}


/**************************************************************************
 * See .h file for description.
 **************************************************************************/
bool is_zero
  (double    value,
   bool log_form)
{
  if ((log_form) && (value < LOG_SMALL)) {
    return(true);
  } else if ((!log_form) && (value == 0.0)) {
    return(true);
  } else {
    return(false);
  }
}

/**************************************************************************
 * See .h file for description.
 **************************************************************************/
bool almost_equal
  (double value1,
   double value2,
   double slop)
{
  if ((value1 - slop > value2) || (value1 + slop < value2)) {
    return(false);
  } else {
    return(true);
  }
}

/*************************************************************************
 * Convert a boolean to and from a "true" or "false" string.
 *************************************************************************/
char* boolean_to_string
 (bool the_boolean)
{
  static char * true_or_false;
  static bool first_time = true;

  if (first_time) {
    true_or_false = (char *)mymalloc(sizeof(char) * 6);
    first_time = false;
  }

  if (the_boolean) {
    strcpy(true_or_false, "true");
  } else {
    strcpy(true_or_false, "false");
  }
  return(true_or_false);
}

bool boolean_from_string
  (char* true_or_false)
{
  if (strcmp(true_or_false, "true") == 0) {
    return(true);
  } else if (strcmp(true_or_false, "false") == 0) {
    return(false);
  } else {
    carp(CARP_FATAL, "Invalid input to boolean_from_string (%s)\n", true_or_false);
  }
  return(false); /* Unreachable. */
}


/**************************************************************************
 * Does a given character appear in a given string?
 **************************************************************************/
bool char_in_string
  (const char* a_string,
   char        a_char)
{
  int  i_string;    /* Index into the string. */
  char string_char; /* Character appearing at that index. */

  i_string = 0;
  string_char = a_string[i_string];
  while (string_char != '\0') {
    if (string_char == a_char) {
      return(true);
    }
    i_string++;
    string_char = a_string[i_string];
  }
  return(false);
}

/**************************************************************************
 * Generic functions for converting between integer and string
 * representations of an enumerated type.
 *
 * Assumes that the longest string representation of the enumerated
 * type does not exceed 100 characters.
 *
 * Assumes that the zeroth enumerated type element is invalid.
 **************************************************************************/
const char * convert_enum_type
  (int     enum_type, /* The enumerated type object to be converted. */
   const char *  enum_strs[], /* String values associated with this type. */
   int     num_enums) /* Number of values of the type. */
{
  if ((enum_type <= 0) || (enum_type >= num_enums)) {
    carp(CARP_FATAL, "Illegal enumerated type value (%d).", enum_type);
  }

  return(enum_strs[enum_type]);
}
    
int convert_enum_type_str
  (const char *  enum_type_str, /* String to be converted. */
   int     default_value, /* Value to return if string not found. */
   const char ** enum_strs,     /* String values associated with this type. */
   int     num_enums)     /* Number of values of the type. */
{
  int i_enum;

  /* If no string was given, return the default. */
  if (enum_type_str == NULL) {
    return(default_value);
  }

  /* Search for the value corresponding to the given string. */
  for (i_enum = 0; i_enum < num_enums; i_enum++) {
    if (strcmp(enum_type_str, enum_strs[i_enum]) == 0) {
      return(i_enum);
    }
  }
  //  carp(CARP_FATAL, "Illegal value (%s).", enum_type_str);
  //  exit(1);
  return( default_value ); 
}

/************************************************************************//**
 * \brief Get the name of the CPU.
 ****************************************************************************/
static const int MAX_HOST_NAME = 100;
const char* hostname
 ()
{
  static char the_hostname[MAX_HOST_NAME];
#ifdef _MSC_VER
  WSADATA wsaData;
  WSAStartup(MAKEWORD(2, 2), &wsaData);
#endif
  int result = gethostname(the_hostname, MAX_HOST_NAME);
  if (result != 0) {
    snprintf(the_hostname, MAX_HOST_NAME, "unknown");
  }

  return(the_hostname);
}

/************************************************************************//**
 * \brief Get the current date and time.
 ****************************************************************************/
const char* date_and_time
  ()
{
  FILE *           date_stream;
  static char      the_date[MAX_HOST_NAME];
  static bool first_time = true;

  if (first_time) {
    // FIXME: we should be able to get the date
    // from a function call rather than having
    // to call a shell command.
#ifdef _MSC_VER
    date_stream = (FILE *)popen("date /T", "r");
#else
    date_stream = (FILE *)popen("date", "r");
#endif
    if( fgets(the_date, MAX_HOST_NAME, date_stream) == NULL ){ return NULL; }
    pclose(date_stream);
  }

  /* Remove the EOL. */
  assert(the_date[strlen(the_date)-1] == '\n');
  the_date[strlen(the_date)-1] = '\0';

  first_time = false;
  return(the_date);
}

/****************************************************************************
 * Copy a string, with allocation.
 ****************************************************************************/
void copy_string
 (char** target,
  char*  source)
{
  if (source == NULL) {
    *target = NULL;
  } else {
    *target = (char *)mycalloc(strlen(source) + 1, sizeof(char));
    strcpy(*target, source);
  }
}

/************************************************************************
 * Copy an array of integers.
 ************************************************************************/
void copy_int_array
 (int  nelems,
  int* source,
  int* target)
{
  int i;

  for (i = 0; i < nelems; i++)
    target[i] = source[i];
}

/**
 * parses a file of length max_lines and returns an array of strings
 */
char** parse_file(
  char* file_name,
  int max_lines,
  int* num_lines
  ){

  FILE *infile;
  if (open_file(file_name, "r", 1, "input", "", &infile) == 0)
    exit(1);

  size_t buf_length = 1024;
  char** lines = (char**) mycalloc(max_lines, sizeof(char*));
  int line_idx = 0;
  int length = 0;
  while ((length = getline(&lines[line_idx], &buf_length, infile)) != -1){
    char* line = lines[line_idx];
    if (line[length-2] == '\n' || line[length-2] == '\r'){
      line[length-2] = '\0';
    } else if (line[length-1] == '\n' || line[length-1] == '\r'){
      line[length-1]='\0';
    }
    line_idx++;
    if (line_idx >= max_lines){
      carp(CARP_FATAL, "Number of lines in %s exceeds maximum of %i!", 
          file_name, max_lines);
    }
  }
  free(lines[line_idx]);
  fclose(infile);
  *num_lines = line_idx;

  return lines;
}

boost::mt19937& get_mt19937() {
  static boost::mt19937 mt19937_;
  return mt19937_;
}

/**
 * Returns an integer in the range between 0 and UNIFORM_INT_DISTRIBUTION_MAX
 */
int myrandom() {
  static boost::random::uniform_int_distribution<>
    uniform_int_distribution_(0, UNIFORM_INT_DISTRIBUTION_MAX);
  return uniform_int_distribution_(get_mt19937());
}

/**
 * Returns an integer in the range [0, max)
 */
int myrandom_limit(int max) {
  return myrandom() % max;
}

void mysrandom(unsigned seed) {
  get_mt19937().seed(seed);
}

#ifdef MAIN


int main (int argc, char *argv[])
{
  FILE *outfile;

  if (argc != 2) {
    carp(CARP_FATAL, "USAGE: utils <filename>");
  }

  if (open_file(argv[1], "w", 1, "output", "", &outfile) == 0)
    exit(1);

  // double double_array[8] = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0};
  // fwrite(double_array, sizeof(double), 8, outfile);
  int int_array[8] = {1,2,3,4,5,6,7,8};
  bool a = myfwrite(int_array, sizeof(int), 0, outfile);
  if (a == true){ 
    printf("myfwrite succeeded\n");
  } else {
    printf("myfwrite failed\n");
  }
  fclose(outfile);

   return(0);
}

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
