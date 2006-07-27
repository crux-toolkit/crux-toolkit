/**************************************************************************
 * FILE parse_arguments.h
 *
 * CREATE DATE: 5/22/2004
 
 * AUTHOR: Charles E. Grant
 
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

enum argument_type { FLAG_ARG, INT_ARG, LONG_ARG, DOUBLE_ARG, STRING_ARG };
enum argument_error { NO_ERROR, UNKNOWN_OPTION, MISSING_VALUE,
                       INVALID_VALUE, UNKNOWN_REQ_ARG, MISSING_REQ_ARG,
                        TOO_MANY_REQ_ARGS, TOO_MANY_OPT_ARGS};

/**
 * The argument struct holds information about a command line argument.
 */
typedef struct {
  const char *name;
  const char *usage;
  void *container;
  enum argument_type type;
} argument;

int parse_arguments_set_opt(const char * name, const char * usage, 
                                void * container, enum argument_type type);
int parse_arguments_set_req(const char * name, const char * usage, 
                              void * container, enum argument_type type);
int parse_arguments(int argc, char * argv[], int die_on_error);
int parse_arguments_get_error(const char ** s);
char * parse_arguments_get_usage(const char * name);

#endif /*PARSE_ARGUMENTS*/
