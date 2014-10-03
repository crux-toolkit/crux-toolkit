/**
 * \file parse_arguments.cpp
 * \brief A central storage location for parameter values, reading in
 * from command line and parameter file.
 */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include "carp.h"
#include "parameter.h"
#include "parse_arguments.h"
#include "objects.h"


/* Limits on the number of arguments */
static const int MAX_OPT_ARGS = 500;
static const int MAX_REQ_ARGS = 25;

/* Make the error message bufer large enough
 * to accomodate the longest allowed argument
 * and the longest allowed error message, and 
 * a little slack for nulls and conjucntion text */
static const int MAX_VALUE_STRLEN = 50; 
static const int MAX_ARG_LENGTH = 250;
static const int MAX_MESSAGE_BUFFER = 510;



/**
 * The argument struct holds information about a command line argument.
 */
typedef struct {
  bool command_line; ///<  the value come from the command line
  bool print; ///< print this option?
  const char *name;  ///< the name of arguemt
  const char *usage; ///< the type of argument
  void *container;  ///< A pointer to storage for the parsed value of the option. 
  enum argument_type type; ///< arguemnt type, (int, char ...?)
} argument;


/* These variables are for error handling */
enum argument_error error = NO_ERROR;
char message[MAX_MESSAGE_BUFFER];

int argument_count = 0;
char ** arguments = NULL;
int optional_count = 0;
argument optional[MAX_OPT_ARGS];
int required_count = 0;
int required_index = 0;
argument required[MAX_REQ_ARGS];
char * usage = NULL;

/* Forward declarations */

int assign_value_from_required(/*const*/ argument * req,  /*const*/ char * value);
int assign_value_from_option(/*const*/ argument * option,  int *index);
int assign_value_from_required_to_hash(/*const*/ argument * req,  
                                       /*const*/ char * value, HASH_T* h);
int assign_value_from_option_to_hash(/*const*/ argument * option,  
                                     int *index, 
                                     HASH_T* h);
/*const*/ argument * get_next_req_argument();
/*const*/ argument * find_option(char * name);
/*const*/ char * base_name(/*const*/ char *s);
void build_message(const char * arg);
size_t get_usage_size(const char * name);
int sprintf_option_default_value(argument * o);
const char * get_option_value_type(argument * o);
int is_numeric(/*const*/ char * s);

/***********************************************************************
 * Function:    parse_arguments_set_opt
 * 
 * Description: define the properties of a command line option
 * 
 * Parameters:
 *   name       A pointer to a string defining the name of the option.
 *
 *   usage      A pointer to a string defining a usage message for the option
 *
 *   container  A pointer to storage for the parsed value of the option. 
 *  
 *   type       The type of the option value (flag, int, long, double, char*)
 *
 *  Caveats:    At most MAX_OPT_ARGS options can be registered.
 *
 *  Returns     1 if successful, 0 if the allowed number of options has
 *              been exceeded.
 * 
 ***********************************************************************/
int parse_arguments_set_opt(const char * name, const char * usage, void * container, 
                enum argument_type type, bool print) {
        
  int result = 0;

  if (optional_count < MAX_OPT_ARGS) {
    optional[optional_count].name = name;
    optional[optional_count].usage = usage;
    optional[optional_count].container = container;
    optional[optional_count].type = type;
    optional[optional_count].print = print;
    optional[optional_count].command_line = false;
    optional_count++;
    result = 1;
  } else {
    error = TOO_MANY_OPT_ARGS;
    build_message(NULL);
  }

  return result;

}

/***********************************************************************
 * Function:    parse_arguments_set_req
 * 
 * Description: define the properties of a required command line argument
 * 
 * Parameters:
 *   name       A pointer to a string defining the name of the option.
 *
 *   usage      A pointer to a string defining a usage message for the option
 *
 *   container  A pointer to storage for the parsed value of the option. 
 *  
 *   type       The type of the option value (flag, int, long, double, char*)
 *
 *  Caveats:    At most MAX_REQ_ARGS options can be registered.
 *
 *  Returns     1 if successful, 0 if the allowed number of required args has
 *              been exceeded.
 * 
 ***********************************************************************/
int parse_arguments_set_req(const char * name, const char * usage, void * container, 
                enum argument_type type, bool print) {

  int result = 0;
  
  if (required_count < MAX_REQ_ARGS) {
    required[required_count].name = name;
    required[required_count].usage = usage;
    required[required_count].container = container;
    required[required_count].type = type;
    required[required_count].command_line = false;
    required[required_count].print = true;
    required_count++;
    result = 1;
  } else {
    error = TOO_MANY_REQ_ARGS;
    build_message(NULL);
  }

  return result;
}

/***********************************************************************
 * Function:      parse_arguments
 * 
 * Description:   Run through the command line arguments interpreting them
 *                using the definitions established by parse_arguments_set_opt
 *                and parse_arguments_set_req.
 * 
 * Parameters:
 *   argc         The number of command line arguments
 *  
 *   argv         The array of command line arguments
 *
 *   die_on_error If non-zero, an error while parsing will cause a usage
 *                message to be printed and the programe to exit. Otherwise
 *                an error message will be recorded and the function will
 *                return 0.
 *
 * 
 * Returns      1 if parse succeded
 *              0 if parse encountered an error
 ***********************************************************************/
int parse_arguments(int argc, char * argv[], int die_on_error) {

  int i;
  int n;
  int result = 0;
  /*const*/ argument * option;
  /*const*/ argument * req;

  argument_count = argc;
  arguments = argv;

  option = NULL;
  req = NULL;
  for (i = 1; i < argc; i++) {
    n = strlen(argv[i]);
    if (argv[i][0] == '-' && argv[i][1] == '-' && n > 1) {
      if ((option = find_option(&(argv[i][1]))) != NULL) {
        error = (argument_error)assign_value_from_option(option, &i);
        if (error != NO_ERROR) {
          /* Missing or incorrect value */
          build_message(argv[i]);
          break;
        }
      } else {
        /* Invalid option */
        error = UNKNOWN_OPTION;
        build_message(argv[i]);
        break;
      }
    } else {
      if ((req = get_next_req_argument()) != NULL) {
        error = (argument_error)assign_value_from_required(req, argv[i]);
        if (error != NO_ERROR) {
          /* Should never reach here */
          break;
        }
      } else {
        /* Too many command line arguments */
        error = UNKNOWN_REQ_ARG;
        build_message(argv[i]);
        break;
      }
    }
  }
  if (error == NO_ERROR && required_index != required_count) {
    error = MISSING_REQ_ARG;
    build_message(required[required_index].name);
  } else if (error == NO_ERROR) {
    result = 1;
  }
  
  if (result != 1 && die_on_error) {
    carp(
      CARP_FATAL, 
      "%s: %s\n%s", 
      base_name(argv[0]), 
      message,
      parse_arguments_get_usage(base_name(argv[0]))
    );
  }
  return result;
}


/***********************************************************************
 * Function:      parse_arguments_into_hash
 * 
 * Description:   Run through the command line arguments interpreting them
 *                using the definitions established by parse_arguments_set_opt
 *                and parse_arguments_set_req.  Ignore the container 
 *                variable and instead store parameter values by inserting
 *                them into the given hash table, using the agument name as
 *                the key.
 * 
 * Parameters:
 *   argc         The number of command line arguments
 *  
 *   argv         The array of command line arguments
 *
 *   HASH_T       The hash table where values are stored
 *
 *   die_on_error If non-zero, an error while parsing will cause a usage
 *                message to be printed and the programe to exit. Otherwise
 *                an error message will be recorded and the function will
 *                return 0.
 *
 * 
 * Returns      1 if parse succeded
 *              0 if parse encountered an error
 ***********************************************************************/
int parse_arguments_into_hash(int argc, char * argv[], 
                              HASH_T* hash, int die_on_error) {

  carp(CARP_DETAILED_DEBUG, "Parsing arguments, inserting values into hash");
  int i;
  int n;
  int result = 0;
  /*const*/ argument * option;
  /*const*/ argument * req;

  argument_count = argc;
  arguments = argv;

  option = NULL;
  req = NULL;
  for (i = 1; i < argc; i++) {
    n = strlen(argv[i]);
    if (argv[i][0] == '-' && argv[i][1] == '-' && n > 1) {
      if ((option = find_option(&(argv[i][1]))) != NULL) {
        // error = assign_value_from_option(option, &i);
        error = (argument_error)assign_value_from_option_to_hash(option, &i, hash);
        if (error != NO_ERROR) {
          /* Missing or incorrect value */
          build_message(argv[i]);
          break;
        }
      } else {
        /* Invalid option */
        error = UNKNOWN_OPTION;
        build_message(argv[i]);
        break;
      }
    } else {
      if ((req = get_next_req_argument()) != NULL) {
        //error = assign_value_from_required(req, argv[i]);
        error = (argument_error)assign_value_from_required_to_hash(req, argv[i], hash);
        if (error != NO_ERROR) {
          /* Should never reach here */
          break;
        }
      } else {
        /* Too many command line arguments */
        error = UNKNOWN_REQ_ARG;
        build_message(argv[i]);
        break;
      }
    }
  }
  if (error == NO_ERROR && required_index != required_count) {
    error = MISSING_REQ_ARG;
    build_message(required[required_index].name);
  } else if (error == NO_ERROR) {
    result = 1;
  }
  
  if (result != 1 && die_on_error) {
    carp(
      CARP_FATAL, 
      "%s: %s\n%s", 
      base_name(argv[0]), 
      message,
      parse_arguments_get_usage(base_name(argv[0]))
    );
  }
  
  // For testing so that cmd line can be parsed more than once
  required_index = 0;
  return result;
}

/***********************************************************************
 * Function:    get_next_req_argument
 * 
 * Description: Provides the next required argument to be assigned a value 
 *              (required arguments are assigned sequentially)
 * 
 * Parameters:  None
 * 
 * Returns:      A pointer to the next required argument, NULL if no required
 *               arguments remain.
 ***********************************************************************/
/*const*/ argument * get_next_req_argument() {

  /*const*/ argument *result = NULL;

  if (required_index < required_count) {
    result = &(required[required_index]);
    required_index++;
  }

  return result;
}

/***********************************************************************
 * Function:    find_option
 * 
 * Description:  Find a valid command line option whose name matches the
 *               the given string.
 * 
 * Parameters:
 *   name        A pointer to a null terminated string which should be the
 *               name of a valid option i.e. matching the name of one of the
 *               options.
 * 
 * Returns:      A pointer to the matching optional argument, NULL if no match found.
 ***********************************************************************/
/*const*/ argument * find_option(char * name) {

  /*const*/ argument * option_found = NULL;
  int i;

  /* Skip over extra '-' for long style options */
  if (name[0] == '-') {
     name++;
  }

  /* Linear search for matching option */
  for (i = 0; i < optional_count; i++) {
    if (0 == strncmp(name, optional[i].name, MAX_ARG_LENGTH)) {
      option_found = &(optional[i]);
    }
  }
  return option_found;
}

/***********************************************************************
 * Function:    assign_value_from_required
 * 
 * Description:  Assign the value from the provided required argument to the 
 *               container assigned to it.
 * 
 * Parameters:
 *   req         A pointer to the argument struct describing this optional
 *               argument.
 *
 *   index       An integer giving the current index into the array of
 *               arguments.
 * 
 * Returns:      An integer should always be NO_ERROR
 ***********************************************************************/
int assign_value_from_required(/*const*/ argument * req,  /*const*/ char * value) {

  switch (req->type) {
    case FLAG_ARG:
      *((int *) req->container) = 1;
      break;
    case INT_ARG:
      *((int *) req->container) = atoi(value);
      break;
    case LONG_ARG:
      *((long *) req->container) = atol(value);
      break;
    case DOUBLE_ARG:
      *((double *) req->container) = atof(value);
      break;
    case STRING_ARG:
      //if container free container?
      *((/*const*/ char **) req->container) = value;
      break;
    case BOOLEAN_ARG:
      *((bool*) req->container) = atoi(value);
      break;
  }
  
  // yes this required value came from the command line
  req->command_line = true;

  return NO_ERROR;
}

/***********************************************************************
 * Function:    assign_value_from_required_to_hash
 * 
 * Description:  Assign the value from the provided required argument to the 
 *               container assigned to it.
 * 
 * Parameters:
 *   req         A pointer to the argument struct describing this optional
 *               argument.
 *
 *   index       An integer giving the current index into the array of
 *               arguments.
 *
 *   hash        A pointer to a hash table into which values should be entered
 * 
 * Returns:      An integer should always be NO_ERROR
 ***********************************************************************/
int assign_value_from_required_to_hash(/*const*/ argument * req,  
                                       /*const*/ char * value,
                                       HASH_T* hash) {

  carp(CARP_DETAILED_DEBUG, 
       "Assigning required '%s' of type '%i' to value '%s'", 
       req->name, (int)req->type, value);
  switch (req->type) {
    case FLAG_ARG:
      //      *((int *) req->container) = 1;
      update_hash_value(hash, req->name, (void*)"1");
      break;
    case INT_ARG:
    case LONG_ARG:
    case DOUBLE_ARG:
    case STRING_ARG:
    case BOOLEAN_ARG:
      //      *((int *) req->container) = atoi(value);
      add_or_update_hash(hash, req->name, value);
      break;
  }
  /* BF: there was no type checking unlike option.  why? */
  
  // yes this required value came from the command line
  req->command_line = true;

  return NO_ERROR;
}

/***********************************************************************
 * Function:    assign_value_from_option
 * 
 * Description:  Assign the value from the provided option to the 
 *               container assigned to it.
 * 
 * Parameters:
 *   option      A pointer to the argument struct describing this optional
 *               argument.
 *
 *   index       A pointer to an integer giving the current index into the 
 *               array of arguments.
 * 
 * Returns:      An integer from the argument_error enumeration
 ***********************************************************************/
int assign_value_from_option(/*const*/ argument * option,  int *index) {

  int more_args = 0;
  int next_arg_is_not_option = 1;
  enum argument_error result = NO_ERROR;
  /*const*/ char * value = NULL;

  if (*index < argument_count -1) {
     more_args = 1;
     if (arguments[(*index) + 1][1] == '-') {
       next_arg_is_not_option = 0;
     }
  }

  switch (option->type) {
    case FLAG_ARG:
      /* No value for this argument */
      *((int *) option->container) = 1;
      result = NO_ERROR;
      break;
    case INT_ARG:
      if (more_args && next_arg_is_not_option) {
        /* Next argument should be value */
        (*index)++;
        value = arguments[*index];
        if (is_numeric(value)) {
          *((int *) option->container) = atoi(value);
          result = NO_ERROR;
        } else {
          /* Value not numeric */
          (*index)--;
          result = INVALID_VALUE;
        }
      } else {
        /* Missing value */
        result = MISSING_VALUE;
      }
      break;
    case LONG_ARG:
      if (more_args && next_arg_is_not_option) {
        /* Next argument should be value */
        (*index)++;
        value = arguments[*index];
        if (is_numeric(value)) {
          *((long *) option->container) = atol(value);
          result = NO_ERROR;
        } else {
          /* Value not numeric */
          (*index)--;
          result = INVALID_VALUE;
        }
      } else {
        /* Missing value */
        result = MISSING_VALUE;
      }
      break;
    case DOUBLE_ARG:
      if (more_args && next_arg_is_not_option) {
        /* Next argument should be value */
        (*index)++;
        value = arguments[*index];
        if (is_numeric(value)) {
          *((double *) option->container) = atof(value);
          result = NO_ERROR;
        } else {
          /* Value not numeric */
          (*index)--;
          result = INVALID_VALUE;
        }
      } else {
        /* Missing value */
        result = MISSING_VALUE;
      }
      break;
    case STRING_ARG:
      if (more_args && next_arg_is_not_option) {
        /* Next argument should be value */
        (*index)++;
        value = arguments[*index];
        *((/*const*/ char **) option->container) = value;
        result = NO_ERROR;
      } else {
        /* Missing value */
        result = MISSING_VALUE;
      }
      break;
    case BOOLEAN_ARG:
      if (more_args && next_arg_is_not_option) {
        /* Next argument should be value (T or F) */
        (*index)++;
        value = arguments[*index];
        *((/*const*/ char **) option->container) = value;
        if (value[0] != 'T' || value[0] != 'F'){
          result = INVALID_VALUE;
        }else{
          result = NO_ERROR;
        }
      } else {
        /* Missing value */
        result = MISSING_VALUE;
      }
      break;
  }
  
  // yes this optional value came from the command line
  option->command_line = true;

  return result;
}

/***********************************************************************
 * Function:    assign_value_from_option_to_hash
 * 
 * Description:  Assign the value from the provided option to the 
 *               container assigned to it.
 * 
 * Parameters:
 *   option      A pointer to the argument struct describing this optional
 *               argument.
 *
 *   index       A pointer to an integer giving the current index into the 
 *               array of arguments.
 *
 *   hash        A pointer to a hash into which values will be inserted.
 * 
 * Returns:      An integer from the argument_error enumeration
 ***********************************************************************/
int assign_value_from_option_to_hash(/*const*/ argument * option,  
                                     int *index, 
                                     HASH_T* hash) {
  carp(CARP_DETAILED_DEBUG, "Assigning option '%s' to value '%s'", 
       option->name, arguments[(*index)+1]);
  int more_args = 0;
  int next_arg_is_not_option = 1;
  enum argument_error result = NO_ERROR;
  /*const*/ char * value = NULL;

  if (*index < argument_count -1) {
     more_args = 1;
     std::string arg_value = arguments[(*index) + 1];
     if (arg_value.compare(0, 2, "--") == 0) {
       next_arg_is_not_option = 0;
     }
  }

  switch (option->type) {
    case FLAG_ARG:
      /* No value for this argument */
      //      *((int *) option->container) = 1;
      update_hash_value(hash, option->name , (void*)"1");
      result = NO_ERROR;
      break;
    case INT_ARG:
      if (more_args && next_arg_is_not_option) {
        /* Next argument should be value */
        (*index)++;
        value = arguments[*index];
        if (is_numeric(value)) {
          //          *((int *) option->container) = atoi(value);
          update_hash_value(hash, option->name, value);
          result = NO_ERROR;
        } else {
          /* Value not numeric */
          (*index)--;
          result = INVALID_VALUE;
        }
      } else {
        /* Missing value */
        result = MISSING_VALUE;
      }
      break;
    case LONG_ARG:
      if (more_args && next_arg_is_not_option) {
        /* Next argument should be value */
        (*index)++;
        value = arguments[*index];
        if (is_numeric(value)) {
          //          *((long *) option->container) = atol(value);
          update_hash_value(hash, option->name, value);
          result = NO_ERROR;
        } else {
          /* Value not numeric */
          (*index)--;
          result = INVALID_VALUE;
        }
      } else {
        /* Missing value */
        result = MISSING_VALUE;
      }
      break;
    case DOUBLE_ARG:
      if (more_args && next_arg_is_not_option) {
        /* Next argument should be value */
        (*index)++;
        value = arguments[*index];
        if (is_numeric(value)) {
          //          *((double *) option->container) = atof(value);
          update_hash_value(hash, option->name, value);
          result = NO_ERROR;
        } else {
          /* Value not numeric */
          (*index)--;
          result = INVALID_VALUE;
        }
      } else {
        /* Missing value */
        result = MISSING_VALUE;
      }
      break;
    case STRING_ARG:
      if (more_args && next_arg_is_not_option) {
        /* Next argument should be value */
        (*index)++;
        value = arguments[*index];
        //        *((/*const*/ char **) option->container) = value;
        update_hash_value(hash, option->name, value);
        result = NO_ERROR;
      } else {
        /* Missing value */
        result = MISSING_VALUE;
      }
      break;
      /* BF added, should remove and use string for these */
    case BOOLEAN_ARG:
      if (more_args && next_arg_is_not_option) {
        /* Next argument should be value (T or F) */
        (*index)++;
        value = arguments[*index];
        *((/*const*/ char **) option->container) = value;
        if (value[0] != 'T' || value[0] != 'F'){
          result = INVALID_VALUE;
        }else{
          result = NO_ERROR;
        }
      } else {
        /* Missing value */
        result = MISSING_VALUE;
      }
      break;
  }
  
  // yes this optional value came from the command line
  option->command_line = true;

  return result;
}

/***********************************************************************
 * Function:     parse_arguments_get_error
 * 
 * Description:  Provides information about an error in a call to 
 *               parse_arguments
 * 
 * Parameters:
 *   message     A pointer to a char pointer, which will be set to
 *               point at a text description of the error.
 *
 * Returns:      An member of the argument_error enumeration
 *               describing the error that occured
 ***********************************************************************/
int parse_arguments_get_error(/*const*/ char ** s) {
  *s = message;
  return error;
}

/***********************************************************************
 * Function:     get_usage_size()
 * 
 * Description:  Calculates the memory needed for the usage string
 * 
 * Parameters:  
 *   name        A pointer to a null terminated string containing the
 *               name of the program
 *
 * Returns:      the amount of memory needed for the usage string
 ***********************************************************************/
size_t get_usage_size(const char * name) {
  size_t memory_used = 0;
  int i = 0;
  memory_used += 100; /* Fixed text and slack */
  memory_used += strlen(name);
  for (i = 0; i < required_count; i++) {
    /* The argument name appears twice in the usage string */
    memory_used += 2 * strlen(required[i].name);
    memory_used += strlen(required[i].usage);
    /* Spaces, newlines, brackets and misc for each arg */
    memory_used += 12;
  }
  for (i = 0; i < optional_count; i++) {
    memory_used += strlen(optional[i].name);
    memory_used += strlen(optional[i].usage);
    /* default value, type, spaces, newlines, brackets and misc for each arg */
    memory_used += 100;
  }
  return memory_used;
}

/***********************************************************************
 * Function:     get_option_value_type
 * 
 * Description:  Returns a string describing the type
 *               of an option value
 * 
 * Parameters:   
 *   o           A pointer to an argument struct
 *
 * Returns:      A pointer to a null terminated string
 *               describing the type of the option value
 *               "", "int", "long", "double", "string"
 ***********************************************************************/
const char * get_option_value_type(argument * o) {

  const char *v = "";
  if (o) {
    switch (o->type){
      case FLAG_ARG:
        break;
      case INT_ARG:
        v = "int";
        break;
      case LONG_ARG:
        v = "long";
        break;
      case DOUBLE_ARG:
        v = "double";
        break;
      case STRING_ARG:
        v = "string";
        break;
      case BOOLEAN_ARG:
        v = "boolean";
        break;
    }
  }

  return v;
}

/***********************************************************************
 * Function:     sprintf_option_default_value
 * 
 * Description:  Print the default value of an option at the 
 *               '\0' marker of the usage string.
 * 
 * Parameters:  
 *   o           A pointer to an argument struct
 *
 * Caveats:      The displayed default value will be truncated
 *               at MAX_VALUE_STRLEN characters.
 *
 * Returns:      An int indicating the number of characters
 *               written to the usage string.
 ***********************************************************************/
int sprintf_option_default_value(argument * o) {

  size_t offset = strlen(usage);
  char * destination = usage + offset;
  int result = 0;

  if (o) {
    switch (o->type){
      case FLAG_ARG:
        result = snprintf(destination, MAX_VALUE_STRLEN, "(default = %d)", 
                        *((int *) o->container));
        break;
      case INT_ARG:
        result = snprintf(destination, MAX_VALUE_STRLEN, "(default = %d)", 
                        *((int *) o->container));
        break;
      case LONG_ARG:
        result = snprintf(destination, MAX_VALUE_STRLEN, "(default = %ld)", 
                        *((long *) o->container));
        break;
      case DOUBLE_ARG:
        result = snprintf(destination, MAX_VALUE_STRLEN, "(default = %g)",
                        *((double *) o->container));
        break;
      case STRING_ARG:
        if (*((char **)o->container) != NULL) {
          result = snprintf(destination, MAX_VALUE_STRLEN, "(default = %s)", 
                          *((char **) o->container));
        } else {
          result = snprintf(destination, MAX_VALUE_STRLEN, "(default = none)");
        }
        break;
      case BOOLEAN_ARG:
        if (*((char **)o->container) != NULL) {
          result = snprintf(destination, MAX_VALUE_STRLEN, "(default = %s)", 
                          *((char **) o->container));
        } else {
          result = snprintf(destination, MAX_VALUE_STRLEN, "(default = none)");
        }
        break;
    }
  }

  return result;
}

/***********************************************************************
 * Function:     parse_arguments_get_usage
 * 
 * Description:  Builds a usage string from the information
 *               passed in calls to parse_arguments_set_opt
 *               and parse_arguments_set_req.
 * 
 * Parameters:  
 *   name        A pointer to a string containing the name of the program.
 *               We use this rather then the 0th command-line argument to
 *               avoid having to parse the command line before building
 *               the usage string and to avoid issuess with links and aliases
 *               disguising the real name of the program.
 *
 *
 * Note:         Assumes default value is in the usage string, does not
 *               look in container for value.
 *
 * Caveats:      Caller is responsible for freeing the memory associated with
 *               the usage string.
 *
 * Returns:      A pointer to a null terminated string containing the
 *               usage message. If this is null it indicates malloc
 *               failed. We must be out of memory and printing a usage
 *               statment is the least of our worries.
 ***********************************************************************/
char * parse_arguments_get_usage(const char * name) {

  int i = 0;
  /* Calculate the size of the buffer we'll need.                     */
  /* This should protect strcat and sprintf against buffer overflows. */
  size_t usage_size = get_usage_size(name);
  usage = (char *) malloc(usage_size);
  if (usage) {

    /* Add the command line summary */
    *usage = 0;
    strcat(usage, "\nUSAGE:\n\n");
    strcat(usage, "   crux ");
    strcat(usage, name);
    if (optional_count > 0) {
      strcat(usage, " [options]");
    }
    strcat(usage, " ");
    for (i = 0; i < required_count; i++) {
      strcat(usage, "<");
      strcat(usage, required[i].name);
      strcat(usage, "> ");
    }
    strcat(usage, "\n\n");

    /* Add the required argument usage comments */
    strcat(usage, "REQUIRED ARGUMENTS:\n\n");
    for (i = 0; i < required_count; i++) {
      std::stringstream ss;
      ss << "<" << required[i].name << "> " << required[i].usage << '\n';
      strcat_formatted(usage, "   ", ss.str().c_str());
    }
    strcat(usage, "\n");

    /* Add the optional argument usage comments */
    strcat(usage, "OPTIONAL ARGUMENTS:\n\n");
    for (i = 0; i < optional_count; i++) {
      if (optional[i].print) {
        strcat(usage, "  [--");
        strcat(usage, optional[i].name);
        if (optional[i].type != FLAG_ARG) {
          strcat(usage, " <");
          strcat(usage, get_option_value_type(&optional[i]));
          strcat(usage, ">");
        }
        strcat(usage, "]\n");
        strcat_formatted(usage, "     ", optional[i].usage);
      }
    }
    strcat(usage, "\nAdditional parameters are documented in "
                 "the online documentation.");
    strcat(usage, "\n");

  }

  return usage;
}

/***********************************************************************
 * Function:    build_message()
 * 
 * Description: Builds a null terminated string describing the most recent 
 *              error encountered by simple_getopt.
 * 
 * Parameter:
 *   arg        A pointer to a null terminated string to be included
 *              in the message.
 * 
 * Caveats:     Message will be truncated to fit in fixed buffer
 *              If you change this function make sure you adjust
 *              the size of the message buffer to accomodate your
 *              your changes without overflowing.
 ***********************************************************************/
void build_message(const char * arg) {
  
  /* Keep the error messages well below MAX_MESSAGE_BUFFER */
  /* in size, to make sure we don't truncate them          */
  /* The error message are in 1-1 correspondence with the  */
  /* argument_error enumeration                            */

  const char * error_messages[] = {
    "No error.",
    "Invalid option (%s).",
    "The option %s is missing its value.",
    "The value for the option %s should be numeric.",
    "The argument %s was not expected.",
    "The required argument <%s> is missing.",
    "Too many required arguments defined. Only %d allowed.",
    "Too many optional arguments defined. Only %d allowed."
  };
      
  switch (error) {
    case TOO_MANY_REQ_ARGS:
      snprintf(message, MAX_MESSAGE_BUFFER - 1, error_messages[error],
               MAX_REQ_ARGS);
      break;
    case TOO_MANY_OPT_ARGS:
      snprintf(message, MAX_MESSAGE_BUFFER - 1, error_messages[error],
               MAX_OPT_ARGS);
      break;
    default:
      snprintf(message, MAX_MESSAGE_BUFFER - 1, error_messages[error], arg);
      break;
  };
}

/***********************************************************************
 * Function:    is_numeric(const char * s)
 * 
 * Description: Tests whether a string represent a base 10 number.
 * 
 * Parameter:
 *   s          Pointer to null terminated string
 * 
 * Caveats:     A numeric is a string of the form ^[+,-]?d+\.?d*$
 *              or ^[+,-]?\.d+$
 *
 * Returns      1 if the string was numeric, 0 otherwise
 ***********************************************************************/
int is_numeric(/*const*/ char * s) {

  int result = 0;
  int found_decimal = 0;

  if (*s == '+' || *s == '-' || *s == '.' || isdigit(*s)) {
    if (*s == '.') { 
      found_decimal = 1; 
    }
    s++;
    
    while (*s) {
      if (*s == '.' || isdigit(*s)) {
        if (*s == '.') {
          if(found_decimal){
            break;
          }
          found_decimal = 1;
        } 
      } 
      else {
        break;
      }
      s++;
    }
    
    if (*s) {
      /* if *s != 0 must have broken out of the while */
      result = 0;
    } else {
      result = 1;
    }
  }

  return result;
}

/***********************************************************************
 * Function:    base_name(const char * s)
 * 
 * Description: Tests whether a string represent a base 10 number.
 * 
 * Parameter:
 *   s          Pointer to null terminated string
 * 
 * Returns      Pointer to a null terminated string containing
 *              the base name of a path
 ***********************************************************************/
/*const*/ char * base_name(/*const*/ char *s) {

  size_t l = 0;
  /*const*/ char * p = NULL;

  l = strlen(s);
  p = s + l;
  while (p != s) {
    if (*p != '/') {
      p--;
    } else {
      p++;
      break;
    }
  }
  return p;
}

/**
 * updates all the parameters in the parameter file with the 
 * higher precedence command line parameters
 * returns true is sucessful, else false
 */
/*bool update_parameter(){
  int array_idx = 0;
  argument* optional_arg = NULL;
  argument* required_arg = NULL;
  char dest[PARAMETER_LENGTH];
  int result = 0;
  char* bool = "false";
  
  // look at every optional arg and update parameter file
  // if the valuse were set from command line
  for(array_idx = 0; array_idx < optional_count; ++array_idx){
    optional_arg = &optional[array_idx];
    
    // skip the arguments that are set from default
    if(!optional_arg->command_line){
      continue;
    }
    
    // convert to string
    switch (optional_arg->type){
      case FLAG_ARG:
        if(*((int *) optional_arg->container) == 1){
          bool = "true";
        }
        result = snprintf(dest, PARAMETER_LENGTH, "%s", 
                          bool);
        break;
      case INT_ARG:
        result = snprintf(dest, PARAMETER_LENGTH, "%d", 
                        *((int *) optional_arg->container));
        break;
      case LONG_ARG:
        result = snprintf(dest, PARAMETER_LENGTH, "%ld", 
                        *((long *) optional_arg->container));
        break;
      case DOUBLE_ARG:
        result = snprintf(dest, PARAMETER_LENGTH, "%g",
                        *((double *) optional_arg->container));
        break;
      case STRING_ARG:
        result = snprintf(dest, PARAMETER_LENGTH, "%s", 
                          *((char **) optional_arg->container));
        break;
      case BOOLEAN_ARG:
        result = snprintf(dest, PARAMETER_LENGTH, "%s", 
                          *((char **) optional_arg->container));
        break;
    }

    // update value in paramters
    if(!set_options_command_line(optional_arg->name, dest, false) || result < 1){
      fprintf(stderr,"failed to update parameters from command line\n");
      return false;
    } 
  }

  // look at every required arg and update parameter file
  for( array_idx = 0; array_idx < required_count; ++array_idx){
    required_arg = &required[array_idx];
        
    // convert to string
    switch (required_arg->type){
      case FLAG_ARG:
        if(*((int *)required_arg->container) == 1){
          bool = "true";
        }
        else{
          bool = "false";
        }
        result = snprintf(dest, PARAMETER_LENGTH, "%s", 
                          bool);
        break;
      case INT_ARG:
        result = snprintf(dest, PARAMETER_LENGTH, "%d", 
                        *((int *) required_arg->container));
        break;
      case LONG_ARG:
        result = snprintf(dest, PARAMETER_LENGTH, "%ld", 
                        *((long *) required_arg->container));
        break;
      case DOUBLE_ARG:
        result = snprintf(dest, PARAMETER_LENGTH, "%g",
                        *((double *) required_arg->container));
        break;
      case STRING_ARG:
        result = snprintf(dest, PARAMETER_LENGTH, "%s", 
                          *((char **) required_arg->container));
        break;
          case BOOLEAN_ARG:
        result = snprintf(dest, PARAMETER_LENGTH, "%s", 
                          *((char **) required_arg->container));
        break;
}

    // update value in paramters
    if(!set_options_command_line(required_arg->name, dest, true) || result < 1){
      fprintf(stderr,"failed to update parameters from command line\n");
      return false;
    } 
  }
  return true;
}
*/

enum argument_type string_to_argument_type(char* arg_type_string){

  enum argument_type type = STRING_ARG;
    
  if( strcmp(arg_type_string, "FLAG_ARG") == 0){
    type = FLAG_ARG;
  }
  else if( strcmp(arg_type_string, "INT_ARG") == 0){
    type = INT_ARG;
  }
  else if( strcmp(arg_type_string, "LONG_ARG") == 0){
    type = LONG_ARG;
  }
  else if( strcmp(arg_type_string, "DOUBLE_ARG") == 0){
    type = DOUBLE_ARG;
  }
  else if( strcmp(arg_type_string, "STRING_ARG") == 0){
    type = STRING_ARG;
  }
  else if( strcmp(arg_type_string, "BOOLEAN_ARG") == 0){
    type = BOOLEAN_ARG;
  }

  return type;
}
