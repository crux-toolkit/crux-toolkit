/******************************************************************************
 * FILE: parameter-file.h
 * AUTHOR: ??
 * CREATE DATE: ??
 * DESCRIPTION: General parameter handling utilities.
 *****************************************************************************/
#ifndef PARAMETER_FILE_H
#define PARAMETER_FILE_H
#include "utils.h"

#define PARAMETER_LENGTH 1024
#define NUM_PARAMS 512
#define MAX_LINE_LENGTH 4096

void parse_parameter_file
(char* filename);

BOOLEAN_T get_boolean_parameter
(char* name, 
 BOOLEAN_T default_value);

int get_int_parameter
(char* name,
 int default_value);

double get_double_parameter
(char* name, 
 double default_value);

char* get_string_parameter
(char* name);

void print_parameters
(char* first_line,
 char* parameter_filename,
 char* lead_string,
 FILE* outstream);

void check_unused_parameters ();

#endif
