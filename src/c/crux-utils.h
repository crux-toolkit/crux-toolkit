/**
 * \file crux-utils.h
 * $Revision: 1.4 $
 * $Author: aklammer $
 * \brief Utilities for the crux project
 */
#ifndef CRUX_UTILS_H
#define CRUX_UTILS_H
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "utils.h"

char* my_copy_string(char* src);


/**
 * returns copy of the src string upto the specified length
 * includes a null terminating `\0' character
 * the string is heap allocated thus, user must free
 */
char* copy_string_part(char*src, int length);

#endif
