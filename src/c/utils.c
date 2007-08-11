/********************************************************************
 * FILE: utils.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 9-8-97
 * PROJECT: shared
 * COPYRIGHT: 1997-2001 Columbia University
 * DESCRIPTION: Various useful generic utilities.
 ********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>
#include <assert.h>
#include "utils.h"


/***********************************************************************
 * Return the value to replace a missing value -- NaN.
 ***********************************************************************/
double NaN
  (void)
{
  return atof("NaN");
}

/**********************************************************************
 * See .h file for description.
 **********************************************************************/
#ifdef NOCLOCK
double myclock() {return(0);}
#else

#ifdef crayc90
/* No problem on the CRAY. */
#include <time.h>
double myclock() {return((double)clock());}

#else
int getrusage(int who, struct rusage *rusage);

double myclock()
{
  static BOOLEAN_T first_time = TRUE;
  static double    start_time;
  double           elapsed;
  struct rusage    ru;

  if (first_time) {
    getrusage(RUSAGE_SELF, &ru);
    start_time = (ru.ru_utime.tv_sec * 1.0E6) + ru.ru_utime.tv_usec;
    first_time = FALSE;
    return 0;

  } else {
    getrusage(RUSAGE_SELF, &ru);
    elapsed = (ru.ru_utime.tv_sec * 1.0E6) + ru.ru_utime.tv_usec -
      start_time;
    return elapsed;
  }
}
#endif /* crayc90 */
#endif /* NOCLOCK */

/***********************************************************************
 * Return elapsed time in microseconds since the last call.
 ***********************************************************************/
double wall_clock(){
  struct timeval tp;
  static double first_time;
  double t;
  static int first_call = 1;

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
BOOLEAN_T open_file
  (char *    filename,            /* Name of the file to be opened. */
   char *    file_mode,           /* Mode to be passed to fopen. */
   BOOLEAN_T allow_stdin,         /* If true, filename "-" is stdin. */
   char *    file_description,   
   char *    content_description,
   FILE **         afile)               /* Pointer to the open file. */
{
  if (filename == NULL) {
    fprintf(stderr, "Error: No %s filename specified.\n", file_description);
    return(FALSE);
  } else if ((allow_stdin) && (strcmp(filename, "-") == 0)) {
    if (strchr(file_mode, 'r') != NULL) {
      fprintf(stderr, "Reading %s from stdin.\n", content_description);
      *afile = stdin;
    } else if (strchr(file_mode, 'w') != NULL) {
      fprintf(stderr, "Writing %s to stdout.\n", content_description);
      *afile = stdout;
    } else {
      fprintf(stderr, "Sorry, I can't figure out whether to use stdin ");
      fprintf(stderr, "or stdout for %s.\n", content_description);
      return(FALSE);
    }
  } else if ((*afile = fopen(filename, file_mode)) == NULL) {
    fprintf(stderr, "Error opening file %s.\n", filename);
    return(FALSE);
  }
  return(TRUE);
}

/*************************************************************************
 * Run a program from a given directory with given arguments.  Return
 * the resulting pipe.
 *************************************************************************/
static FILE* run_program
  (char*      program,     // The program to run in the pipe.
   char*      directory,   // Directory where program resides.
   char*      arguments,   // Program arguments.
   char*      type)        // Read ("r") or write ("w").
{
  char* command;
  FILE* return_value;

  // Allocate space for the command.
  command = (char*)mymalloc(sizeof(char) * (strlen(directory)
					    + strlen(program) 
					    + strlen(arguments) + 2));

  // Formulate the command.  Deal with directories possibly ending
  // with a slash or not.
  if (strlen(directory) == 0) {
    sprintf(command, "%s %s", program, arguments);
  } else if (directory[strlen(directory) - 1] == '/') {
    sprintf(command, "%s%s %s", directory, program, arguments);
  } else {
    sprintf(command, "%s/%s %s", directory, program, arguments);
  }

  // Run the program.
  return_value = popen(command, type);
  myfree(command);
  return(return_value);
}


/*************************************************************************
 * Attempt to run a given program in a given directory with the given
 * arguments, and check that it gives the expected one-line reply.
 *************************************************************************/
static BOOLEAN_T try_to_run
  (char*      program,          // The program to run in the pipe.
   char*      directory,        // Directory to look in.
   char*      test_arguments,   // Arguments used when searching for program.
   char*      expected_reply)   // Expected reply from search.
{
  char* reply;
  FILE* pipe;
  BOOLEAN_T return_value;

  // Allocate space for the reply.
  reply = (char*)mymalloc(sizeof(char) * (strlen(expected_reply) + 1));

  // Run the command.
  pipe = run_program(program, directory, test_arguments, "r");

  // Check the pipe.
  if (pipe == NULL) {
    return_value = FALSE;
  } else {

    // Read from the pipe.
    if (fgets(reply, strlen(expected_reply) + 1, pipe) == NULL) {
      return_value = FALSE;
    } else {
      return_value = (strcmp(reply, expected_reply) == 0);
    }
  }

  // Close the pipe.
  if (pclose(pipe) == -1) {
    return_value = FALSE;
  }

  myfree(reply);
  return(return_value);
}


/*************************************************************************
 * Open a read-only pipe using a given command line.
 *************************************************************************/
FILE* open_command_pipe
  (char*     program,          // The program to run in the pipe.
   char*     directory,        // Directory to look in.
   char*     test_arguments,   // Arguments used when searching for program.
   char*     expected_reply,   // Expected reply from search.
   char*     real_arguments,   // Arguments used when running the program.
   BOOLEAN_T stdout_on_error,  // If command fails, return STDOUT?
   char*     error_message)    // Error or warning if command fails.
{
  FILE* return_value;

  // Try to run the command with no directory specified.
  if (try_to_run(program, "", test_arguments, expected_reply)) {
    return_value = run_program(program, "", real_arguments, "w");
  }

  // Try to run the program in the specified directory.
  else if (try_to_run(program, directory, test_arguments, expected_reply)) {
    return_value = run_program(program, directory, real_arguments, "w");

  } else {

    // If we failed, print the error message.
    fprintf(stderr, error_message);
    if (stdout_on_error) {
      return_value = stdout;
    } else {
      exit(1);
    }
  }

  return(return_value);
}


/********************************************************************
 * See .h file for description.
 ********************************************************************/
void die
  (char *format, 
   ...)
{
  va_list  argp;

  fprintf(stderr, "FATAL: ");
  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  fprintf(stderr, "\n");
  fflush(stderr);

#ifdef DEBUG
  abort();
#else
  exit(1);
#endif
}


/**************************************************************************
 * See .h file for description.
 **************************************************************************/
void myassert
  (BOOLEAN_T die_on_error,
   BOOLEAN_T test,
   char * const    format,
   ...)
{
  va_list  argp;

  if (!test) {

    if (die_on_error) {
      fprintf(stderr, "FATAL: ");
    } else {
      fprintf(stderr, "WARNING: ");
    }

    /* Issue the error message. */
    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
    fprintf(stderr, "\n");
    fflush(stderr);
    
    if (die_on_error) {
#ifdef DEBUG
      abort();
#else
      exit(1);
#endif
    }
  }      
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
    die("Memory exhausted.  Cannot allocate %d bytes.", (int)size);
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

  if (temp_ptr == NULL)
    die("Memory exhausted.  Cannot allocate %d bytes.", (int)size);

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

  if (temp_ptr == NULL)
    die("Memory exhausted.  Cannot allocate %d bytes.", (int)size);

  return(temp_ptr);
}

#ifdef MYRAND
#define MY_RAND_MAX 4096

/********************************************************************
 * Primary function for the built-in random number generator. 
 ********************************************************************/
static double my_rand
  (long seed)
{
  static long stored_seed = 0;

  /* If this is the first call, just set the seed. */
  if (stored_seed == 0) {
    stored_seed = seed;
  }

  /* Otherwise, create a new pseudorandom number. */
  else {
    stored_seed = abs((stored_seed / 3) * stored_seed + 7718);
  }

  /* Make sure the pseudorandom number is in the right range. */
  return((double)(stored_seed % MY_RAND_MAX) / (double)MY_RAND_MAX);
}
#else
/* The stupid include file doesn't have these prototypes. */
void srand48();
double drand48();

#endif

/********************************************************************
 * See .h file for description.
 ********************************************************************/
void my_srand
  (long seed)
{
#ifdef MYRAND
  my_rand(seed);
#else
  srand48(seed);
#endif
}

/********************************************************************
 * See .h file for description.
 ********************************************************************/
double my_drand
  (void)
{
#ifdef MYRAND
  return(my_rand(0));
#else
  return(drand48());
#endif
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
    die("Tried to take the log of a negative value (%g).", x);
  } /* else if (x == 0.0) */
  return(LOG_ZERO);
}

/* The lookup table. */
#define LOG_PRECISION 1.0e5
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
BOOLEAN_T is_zero
  (double    value,
   BOOLEAN_T log_form)
{
  if ((log_form) && (value < LOG_SMALL)) {
    return(TRUE);
  } else if ((!log_form) && (value == 0.0)) {
    return(TRUE);
  } else {
    return(FALSE);
  }
}

/**************************************************************************
 * See .h file for description.
 **************************************************************************/
BOOLEAN_T almost_equal
  (double value1,
   double value2,
   double slop)
{
  if ((value1 - slop > value2) || (value1 + slop < value2)) {
    return(FALSE);
  } else {
    return(TRUE);
  }
}

/*************************************************************************
 * Convert a boolean to and from a "true" or "false" string.
 *************************************************************************/
char* boolean_to_string
 (BOOLEAN_T the_boolean)
{
  static char * true_or_false;
  static BOOLEAN_T first_time = TRUE;

  if (first_time) {
    true_or_false = (char *)mymalloc(sizeof(char) * 6);
    first_time = FALSE;
  }

  if (the_boolean) {
    strcpy(true_or_false, "true");
  } else {
    strcpy(true_or_false, "false");
  }
  return(true_or_false);
}

BOOLEAN_T boolean_from_string
  (char* true_or_false)
{
  if (strcmp(true_or_false, "true") == 0) {
    return(TRUE);
  } else if (strcmp(true_or_false, "false") == 0) {
    return(FALSE);
  } else {
    die("Invalid input to boolean_from_string (%s)\n", true_or_false);
  }
  return(FALSE); /* Unreachable. */
}


/**************************************************************************
 * Does a given character appear in a given string?
 **************************************************************************/
BOOLEAN_T char_in_string
  (const char* a_string,
   char        a_char)
{
  int  i_string;    /* Index into the string. */
  char string_char; /* Character appearing at that index. */

  i_string = 0;
  string_char = a_string[i_string];
  while (string_char != '\0') {
    if (string_char == a_char) {
      return(TRUE);
    }
    i_string++;
    string_char = a_string[i_string];
  }
  return(FALSE);
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
char * convert_enum_type
  (int     enum_type, /* The enumerated type object to be converted. */
   char *  enum_strs[], /* String values associated with this type. */
   int     num_enums) /* Number of values of the type. */
{
  if ((enum_type <= 0) || (enum_type >= num_enums)) {
    die("Illegal enumerated type value (%d).", enum_type);
  }

  return(enum_strs[enum_type]);
}
    
int convert_enum_type_str
  (char *  enum_type_str, /* String to be converted. */
   int     default_value, /* Value to return if first arg is null. */
   char ** enum_strs,     /* String values associated with this type. */
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
  die("Illegal value (%s).", enum_type_str);
  return(0); /* Unreachable. */
}

/****************************************************************************
 * Get the name of the CPU.
 ****************************************************************************/
#define HOST_LENGTH 100
const char* hostname
  ()
{
  FILE *           hostname_stream;
  static char      the_hostname[HOST_LENGTH];
  static BOOLEAN_T first_time = TRUE;
  int              num_scanned;

  if (first_time) {
    hostname_stream = (FILE *)popen("hostname", "r"); /* SGI needs cast. */
    num_scanned = fscanf(hostname_stream, "%s", the_hostname);
    assert(num_scanned == 1);
    num_scanned = pclose(hostname_stream);
    assert(num_scanned == 0);
  }
  return(the_hostname);
}

/****************************************************************************
 * Get the current date and time.
 ****************************************************************************/
const char* date_and_time
  ()
{
  FILE *           date_stream;
  static char      the_date[HOST_LENGTH];
  static BOOLEAN_T first_time = TRUE;

  if (first_time) {
    date_stream = (FILE *)popen("date", "r"); /* SGI needs cast. */
    fgets(the_date, HOST_LENGTH, date_stream);
    pclose(date_stream);
  }

  /* Remove the EOL. */
  assert(the_date[strlen(the_date)-1] == '\n');
  the_date[strlen(the_date)-1] = '\0';

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

#ifdef MAIN


int main (int argc, char *argv[])
{
  FILE *infile;
  char word[1000];
  long seed;
  int i, j;

  if (argc != 2) {
    die("USAGE: utils <filename>");
  }

  if (open_file(argv[1], "r", 1, "input", "", &infile) == 0)
    exit(1);

  while (fscanf(infile, "%s", word) == 1)
    printf("%s ", word);

  fclose(infile);

  /* Test the random number generator. */
  seed = time(0);
  my_srand(seed);
  printf("\nSome random numbers (seed=%ld): \n", seed);
  for (i = 0; i < 10; i++) {
    for (j = 0; j < 10; j++) {
      printf("%6.4f ", my_drand());
    }
    printf("\n");
  }
  return(0);
}

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
