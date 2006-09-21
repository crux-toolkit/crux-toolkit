/**
 * \file ion.h
 * $Revision: 1.2 $
 * \brief Object for representing one ion in an ion_series.
 *
 */
#ifndef ION_H
#define ION_H
#include <stdio.h>
#include <stdlib.h>
#include "objects.h"

/**
 * \returns an ION_T object
 */
ION_T* new_ion (
  ION_TYPE_T type,   ///< intensity for the new ion -in 
  int cleavage_idx, ///< index into the peptide amide bonds of this ion
  int charge, ///< charge of the ion
  char* peptide ///< location for the new ion -in
  );

/**
 * frees A ION_T object
 */
void free_ion (
  ION_T* garbage_ion ///< the ion to free -in
  );

/**
 * \returns the location of ION_T object
 */
float get_ion_mass(
  ION_T* working_ion///< return the location of this ion -in 
  );


/**
 * sets the mass of the ION_T object
 * while this can be calculated from the char*, cleavage_idx and
 * modifications, it allows some optimizations if we allow it to be set
 * instead
 */
void set_ion_mass(
  ION_T* working_ion, ///<set the m/z location of this ion -out
  float mass ///< the m/z location -in
  );

/**
 * prints the location and fields of ION_T object to the file, in the
 * following format:
 *
 * m/z \\t mass \\t charge \\t ion-series \\t  ...
 *  peptide-bond-index \\t modifications \n
 *
 * Where:
 *
 * m/z - is the ion's mass-to-charge
 * mass - is the ion's (charged) mass
 * charge - is the ion's charge e.g. 1,2,3
 * ion-series - is one of (b,y,p)
 * bond-index - is in [1...n), where n is peptide length
 * modifications - is one of (none|nh3|h2o)
 *
 * if the ion has more than one modification, each will be printed on a
 * separate line, with the necessary number of tabs to right justify
 */
void print_ion(
  ION_T* ion, ///< print this ion -in
  FILE* file ///< to this file -in
  );

/**
 * Adds the given ION_MODIFICATION to this ion
 */
void add_modification(
  ION_T* ion,///< ion to which to add the modification -mod
  ION_MODIFICATION_T modification ///< add this modification to the ion -in
  );

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
