/*************************************************************************//**
 * \file carp.cpp
 * \brief Object for representing a single protein.
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "carp.h"
#include "crux-utils.h"
#include "parameter.h"
#include "utils.h"

using namespace std;

/**
 * Constants
 */
static int G_verbosity; 
static FILE* log_file = NULL;

HASH_T * messages_;
unsigned int hash_size_ = 1000;

void set_verbosity_level(int verbosity){
  G_verbosity = verbosity;
}

int get_verbosity_level(void){
  return G_verbosity;
}

/**
 * Open log file for carp messages.
 *
 * Parameters must have been processed before calling this function.
 */
void open_log_file(char **log_file_name) {
  char* output_dir = get_string_parameter("output-dir");
  bool overwrite = get_boolean_parameter("overwrite");
  prefix_fileroot_to_name(log_file_name);
  log_file = create_file_in_path(*log_file_name, output_dir, overwrite);
  free(output_dir);
}

/**
 * Print command line to log file.
 *
 * Parameters must have been processed before calling this function.
 */
void log_command_line(int argc, char *argv[]) {
  // Command line arguments were shifted, shift back.
  ++argc;
  --argv;
  if (log_file != NULL) {
    fprintf(log_file, "COMMAND LINE: ");
    int i = 0;
    for (i = 0; i < argc; ++i) {
      fprintf(log_file, "%s%c", argv[i], i < (argc - 1) ? ' ' : '\n');
    }
  }
}

static void carp_print(const char *string) {
  fprintf(stderr, "%s", string);
  if (log_file != NULL) {
    fprintf(log_file, "%s", string);
  }
}

/**
 * Print message to log file.
 *
 * Print severity level and message to log file.
 * The term 'carp' is used because 'log' is already used 
 * by the math library. 
 *
 * Verbosity of CARP_FATAL will cause the 
 * program to exit with status code 1.
 *
 */
void carp( int verbosity, const char* format, ...) {
  if (verbosity <= G_verbosity){
    va_list  argp;

    if (verbosity == CARP_WARNING){
      carp_print("WARNING: ");
    }
    else if (verbosity == CARP_ERROR){
      carp_print("ERROR: ");
    }
    else if (verbosity == CARP_FATAL){
      carp_print("FATAL: ");
    }
    else if (verbosity == CARP_INFO){
      carp_print("INFO: ");
    }
    else if (verbosity == CARP_DETAILED_INFO){
      carp_print("DETAILED INFO: ");
    }
    else if (verbosity == CARP_DEBUG){
      carp_print("DEBUG: ");
    }
    else if (verbosity == CARP_DETAILED_DEBUG){
      carp_print("DETAILED DEBUG: ");
    } 
    else {
      carp_print("UNKNOWN: ");
    }

    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
    if (log_file != NULL) { 
      va_start(argp, format); //BF: added to fix segfault
      vfprintf(log_file, format, argp);
      va_end(argp);
    }
    carp_print("\n");
    fflush(stderr);
    if (log_file != NULL) {
      fflush(log_file);
    }
  } 
  if (verbosity == CARP_FATAL) {
    // Fatal carps cause the program to exit
#ifdef DEBUG
    abort(); // Dump core in DEBUG mode.  Use 'make CXXFLAGS=-DDEBUG"'
#else
    exit(1);
#endif
  }
}

void carp ( int verbosity, string& msg) {

  string carp_msg = msg;
  carp(verbosity, "%s", carp_msg.c_str());
}

/**
 * Similar to carp(), but multiple similar messages will
 * be printed only once.
 */
void warn_once(const char * msg1, const char * msg2_format, ...) {

  // Create hash table if not exist
  if (messages_ == NULL) {
    messages_ = new_hash(hash_size_);
  }

  // Look up msg1 in the hash table
  if (get_hash_value(messages_, msg1) == NULL) {
    // if msg1 does not exist in the hash table, then we will save it to the hash table
    // we also have to make sure that the hash table still has enough space
    if (hash_size(messages_) == hash_size_) {
      hash_size_ *= 2;
      HASH_T * temp = new_hash(hash_size_);
      HASH_ITERATOR_T * it = new_hash_iterator(messages_);
      while (hash_iterator_has_next(it)) {
        char * key = hash_iterator_next(it);
        add_hash(temp, key, key);
      }
      free_hash_iterator(it);
      free_hash(messages_);
      messages_ = temp;
    }
    
    // now we add msg1 to our hash table
    add_hash(messages_, msg1, msg1);
    // and we print msg1 to stderr and log msg1 to log_file
    carp(CARP_WARNING, msg1);
  }
  
  va_list argp;
  // msg2 is only printed to the log file
  if (log_file != NULL) { 
      va_start(argp, msg2_format); //BF: added to fix segfault
      vfprintf(log_file, msg2_format, argp);
      va_end(argp);
      fflush(log_file);
  }
  
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

