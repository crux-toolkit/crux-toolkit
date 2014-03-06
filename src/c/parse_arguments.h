/**
 * \file parse_arguments.h
 * \brief  Support for parsing of arguments from command line
 */

/* CREATE DATE: 5/22/2004
 
 * AUTHOR: Charles E. Grant
 * MODIFIED: Chris Park

 * PROJECT: utilities
 
 * COPYRIGHT: 2004, University of Washington
 
 * DESCRIPTION: Support for parsing of arguments from command line
 
 * CAVEATS:  It is expected that command line arguments will be give in
 * the form "cmd [option1] [option2] ... [optionn] required1 required2 ..."
 * Each option will start with '-' and may be followed by a single
 * option value. Any argument that does not start with '-' and is
 * not an option value will be taken to be a required argument
 **************************************************************************/
#ifndef PARSE_ARGUMENTS_H
#define PARSE_ARGUMENTS_H

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define DOXYGEN_SHOULD_SKIP_THIS

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include "carp.h"
#include "objects.h"
#include "hash.h"
#include "WinCrux.h"

enum argument_type { FLAG_ARG, INT_ARG, LONG_ARG, DOUBLE_ARG, 
                     STRING_ARG, BOOLEAN_ARG };
enum argument_error { NO_ERROR, UNKNOWN_OPTION, MISSING_VALUE,
                      INVALID_VALUE, UNKNOWN_REQ_ARG, MISSING_REQ_ARG,
                      TOO_MANY_REQ_ARGS, TOO_MANY_OPT_ARGS};



int parse_arguments_set_opt(const char * name, const char * usage, 
                                void * container, enum argument_type type,
                                bool print=true);
int parse_arguments_set_req(const char * name, const char * usage, 
                              void * container, enum argument_type type,
                              bool print=true);
int parse_arguments(int argc, char * argv[], int die_on_error);
int parse_arguments_into_hash(int argc, char * argv[], HASH_T* h, 
                              int die_on_error);
int parse_arguments_get_error(/*const*/ char ** s);
char * parse_arguments_get_usage(const char * name);
enum argument_type string_to_argument_type(char* arg_type_str);

/**
 * updates all the parameters in the parameter file with the 
 * higher precedence command line parameters
 * returns TRUE is sucessful, else FALSE
 */
//BOOLEAN_T update_parameter(void);

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
#endif /* PARSE_ARGUMENTS_H */
