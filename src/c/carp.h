/**
 * \file carp.h 
 * $Revision: 1.6 $
 * \brief Provides methods for logging error messages, and setting verbosity level.
 *****************************************************************************/
#ifndef CARP_H 
#define CARP_H

#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/** 
 * Verbosity level for a fatal error (e.g., could not open an input file)
 */
#define CARP_FATAL 0    

/** 
 * Verbosity level for a serious, not fatal, error
 * (e.g., could not close a file handle)
 */
#define CARP_ERROR 10   

/**
 * Verbosity level for a warning (e.g., a spectrum has no peaks)
 */
#define CARP_WARNING 20 

/**
 * Verbosity level for informational message (e.g., processed X lines of file)
 */
#define CARP_INFO 30    

/**
 * Verbosity level for detailed informational message (e.g, on spectrum 1000 )
 */
#define CARP_DETAILED_INFO 40    

/** 
 * Verbosity level for a debugging message
 */
#define CARP_DEBUG 50   

/** 
 * Verbosity level for very detailed debugging message
 */
#define CARP_DETAILED_DEBUG 60  

/** 
 * The maximum verbosity level
 */
#define CARP_MAX 100 

//SJM: Code to help with optimization of verbosity code.
//Allow for compilation to remove carp commands
//in the preprocessor.

#define IF_CARP(x,y) if (get_verbosity_level() >= x) {y;}
#define CRUX_DEBUG 

#ifdef CRUX_DEBUG
#define IF_CARP_DEBUG(y) IF_CARP(CARP_DEBUG,y)
#define IF_CARP_DETAILED_DEBUG(y) IF_CARP(CARP_DETAILED_DEBUG,y)
#else
#define IF_CARP_DEBUG(y) y;
#define IF_CARP_DETAILED_DEBUG(y) y;
#endif


#include <stdio.h>
#include "utils.h"

/**
 * Sets the global verbosity level. Expects an integer greater than or equal
 * to CARP_FATAL and less than CARP_MAX
 */
void set_verbosity_level(int verbosity);

/**
 * \returns the current verbosity level.
 */
int get_verbosity_level(void);

/**
 * Open log file for carp messages.
 *
 * Parameters must have been processed before calling this function.
 */
void open_log_file(char **log_file_name);

/**
 * Print command line to log file.
 *
 * Parameters must have been processed before calling this function.
 */
void log_command_line(int argc, char *argv[]);

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
void carp(
  int verbosity, 
  const char* format,
  ...
);

/**
 * \def carp_once( verbosity, msg, ...)
 *
 * \brief Print message to log file, just once.
 *
 * Similar to the carp function, below, but will only print the message once.
 */
#define carp_once( verbosity, msg, ... ) \
{ \
  static BOOLEAN_T _carp_; \
  if (!_carp_) carp(verbosity, msg, ## __VA_ARGS__); \
  _carp_ = 1;\
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

#ifdef __cplusplus
}
#endif

#endif
