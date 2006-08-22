/*****************************************************************************
 * \file carp.c
 * $Revision: 1.3 $
 * \brief: Object for representing a single protein.
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "utils.h"
#include "carp.h"

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

BOOLEAN_T carp(
    int verbosity, 
    char* format,
    ...){
  if (verbosity <= G_verbosity){
    va_list  argp;

    if(verbosity == CARP_WARNING){
      fprintf(stderr, "WARNING: ");
    }
    else if(verbosity == CARP_ERROR){
      fprintf(stderr, "ERROR: ");
    }
    else if(verbosity == CARP_FATAL){
      fprintf(stderr, "FATAL: ");
    }
    else if(verbosity == CARP_INFO){
      fprintf(stderr, "INFO: ");
    }
    else{
      fprintf(stderr, "DEBUG: ");
    }

    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
    fprintf(stderr, "\n");
    fflush(stderr);
  } 
  return TRUE;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

