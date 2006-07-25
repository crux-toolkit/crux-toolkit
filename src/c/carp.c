/*****************************************************************************
 * \file carp.c
 * $Revision: 1.1 $
 * \brief: Object for representing a single protein.
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "utils.h"

/**
 * Constants
 */
static int G_verbosity; 

void set_verbosity_level(int verbosity){
  G_verbosity = verbosity;
}

int get_verbosity_level(void){
  return G_verbosity;
}

void carp(
    int verbosity, 
    char* format,
    ...){
  if (verbosity <= G_verbosity){
    va_list  argp;
    fprintf(stderr, "WARNING: ");
    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
    fprintf(stderr, "\n");
    fflush(stderr);
  } 
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

