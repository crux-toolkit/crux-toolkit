/**
 * \file spectrum.h 
 * $Revision: 1.20 $
 * \brief Object for representing one spectrum.
 *****************************************************************************/
#ifndef SPECTRUM_H
#define SPECTRUM_H

#include "utils.h"
#include "peak.h"
#include <stdio.h>

/**
 * \typedef SPECTRUM_T 
 * \brief A spectrum
 */
typedef struct spectrum SPECTRUM_T;

/**
 * \typedef PEAK_ITERATOR_T 
 * \brief An object to iterate over the peaks in a spectrum
 */
typedef struct peak_iterator PEAK_ITERATOR_T;

/**
 * The enum for spectrum type (MS1, MS2, MS3)
 */
enum _spectrum_type { MS1, MS2, MS3 };

/**
 * \typedef SPECTRUM_TYPE_T 
 * \brief The typedef for spectrum type (MS1, MS2, MS3)
 */
typedef enum _spectrum_type SPECTRUM_TYPE_T;

/**
 * \returns An (empty) spectrum object.
 */
SPECTRUM_T* allocate_spectrum(void);

/**
 * \returns A new spectrum object, populated with the user specified parameters.
 */
SPECTRUM_T* new_spectrum(
  int               first_scan,         ///< The number of the first scan
  int               last_scan,          ///< The number of the last scan
  SPECTRUM_TYPE_T   spectrum_type,      ///< The type of spectrum.
  float             precursor_mz,       ///< The m/z of the precursor (for MS-MS spectra)
  int*              possible_z,         ///< The possible charge states of this spectrum
  int               num_possible_z,     ///< The number of possible charge states of this spectrum  
  char*             filename);          ///< Optional filename

/**
 * Frees an allocated spectrum object.
 */
void free_spectrum (SPECTRUM_T* spectrum);

/**
 * Prints a spectrum object to file.
 */
void print_spectrum(SPECTRUM_T* spectrum, FILE* file);

/**
 * Prints a spectrum object to STDOUT.
 */
void print_spectrum_stdout(SPECTRUM_T* spectrum);

/**
 * Copies spectrum object src to dest.
 */
void copy_spectrum(
  SPECTRUM_T* src,
  SPECTRUM_T* dest);

/**
 * Parses a spectrum from file.
 * \returns TRUE if success. FALSE is failure.
 */
BOOLEAN_T parse_spectrum_file(
  SPECTRUM_T* spectrum,
  FILE* file);

/**
 * \returns TRUE if success. FALSE if failure.
 */
BOOLEAN_T parse_spectrum(
  SPECTRUM_T* spectrum,
  char*      filename);

/** 
 * Access routines of the form get_<object>_<field> and set_<object>_<field>. 
 * FIXME Chris, could you create the get and set methods for the object fields?
 */

/**
 * Additional get and set methods
 */

/**
 * \returns the number of the first scan
 */
int get_spectrum_first_scan(SPECTRUM_T* spectrum);

/**
 * \sets the number of the first scan
 */
void set_spectrum_first_scan(SPECTRUM_T* spectrum, int first_scan);

/**
 * \returns the number of the last scan
 */
int get_spectrum_last_scan(SPECTRUM_T* spectrum);

/**
 * \sets the number of the last scan
 */
void set_spectrum_last_scan(SPECTRUM_T* spectrum, int last_scan);

/**
 * \returns the spectrum_id
 */
int get_spectrum_id(SPECTRUM_T* spectrum);

/**
 * \sets the spectrum_id
 */
void set_spectrum_id(SPECTRUM_T* spectrum, int id);

/**
 * \returns the spectrum_type
 */
SPECTRUM_TYPE_T get_spectrum_spectrum_type(SPECTRUM_T* spectrum);

/**
 * \sets the spectrum_type
 */
void set_spectrum_spectrum_type(SPECTRUM_T* spectrum, SPECTRUM_TYPE_T spectrum_type);

/**
 * \returns the m/z of the precursor
 */
float get_spectrum_precursor_mz(SPECTRUM_T* spectrum);

/**
 * \sets the m/z of the precursor
 */
void set_spectrum_precursor_mz(SPECTRUM_T* spectrum, float precursor_mz);

/**
 * \returns the possible charge states of this spectrum
 * returns an int* to a heap allocated copy of the src spectrum
 * thus, the user must free the memory
 * number of possible charge states can be gained by 
 * the get_num_possible_z function.
 */
int* get_spectrum_possible_z(SPECTRUM_T* spectrum);
 
/**
 * \sets the possible charge states of this spectrum
 * the function copies the possible_z into a heap allocated memory
 * num_possible_z must match the array size of possible_z 
 * frees the memory of the possible_z that is replaced
 * updates the number of possible charge states field
 */
void set_spectrum_possible_z(SPECTRUM_T* spectrum, 
                             int* possible_z, 
                             int num_possible_z);


/**
 * \sets the possible charge states of this spectrum
 * this function should only be used when possible_z is set to NULL
 * to change existing possible_z use set_spectrum_possible_z()
 * the function copies the possible_z into a heap allocated memory
 * num_possible_z must match the array size of possible_z 
 * updates the number of possible charge states field
 */
void set_spectrum_new_possible_z(SPECTRUM_T* spectrum, 
                                 int* possible_z, 
                                 int num_possible_z);

/**
 * \returns the number of possible charge states of this spectrum
 */
int get_spectrum_num_possible_z(SPECTRUM_T* spectrum);

/**
 * \returns the minimum m/z of all peaks
 */
float get_spectrum_min_peak_mz(SPECTRUM_T* spectrum);

/**
 * \returns the maximum m/z of all peaks
 */
float get_spectrum_max_peak_mz(SPECTRUM_T* spectrum);

/**
 * \returns the number of peaks
 */
int get_spectrum_num_peaks(SPECTRUM_T* spectrum);

/**
 * \returns the sum of intensities in all peaks
 */
double get_spectrum_total_energy(SPECTRUM_T* spectrum);

/**
 * \returns the filename of the ms2 file the spectrum was parsed
 * returns a char* to a heap allocated copy of the filename
 * user must free the memory
 */
char* get_spectrum_filename(SPECTRUM_T* spectrum);

/**
 * \sets the filename of the ms2 file the spectrum was parsed
 * copies the value from arguement char* filename into a heap allocated memory
 * frees memory for the filename that is replaced
 */
void set_spectrum_filename(SPECTRUM_T* spectrum, char* filename);


/**
 * \sets the filename of the ms2 file the spectrum was parsed
 * this function should be used only the first time the filename is set(set to NULL)
 * to change existing filename use set_spectrum_filename
 * copies the value from arguement char* filename into a heap allocated memory
 */
void set_spectrum_new_filename(SPECTRUM_T* spectrum, char* filename);


/**
 * \returns The intensity of the peak with the maximum intensity.
 */
float get_spectrum_max_peak_intensity(SPECTRUM_T* spectrum);

/**
 * \returns The mass of the charged precursor ion, according to the formula 
 * mass = m/z * charge
 */
float get_spectrum_mass(SPECTRUM_T* spectrum, int charge);

/**
 * \returns The mass of the neutral precursor ion, according to the formula 
 * mass = m/z * charge - mass_H * charge
 */
float get_spectrum_neutral_mass(SPECTRUM_T* spectrum, int charge);

/**
 * \returns The mass of the singly charged precursor ion, according to the formula 
 * mass = m/z * charge - (mass_H * (charge - 1))
 */
float get_spectrum_singly_charged_mass(SPECTRUM_T* spectrum, int charge);

/**
 * Adds a peak to the spectrum.
 */
BOOLEAN_T add_peak(SPECTRUM_T* spectrum, PEAK_T* peak);

/**
 * Adds a peak to the spectrum given an intensity and location.
 * Calls add_peak function
 */
BOOLEAN_T add_peak_to_spectrum(
  SPECTRUM_T* spectrum, 
  float intensity,
  float location_mz );



/******************************************************************************/

/**
 * Instantiates a new peak_iterator from a spectrum.
 * \returns a PEAK_ITERATOR_T object.
 */
PEAK_ITERATOR_T* new_peak_iterator(SPECTRUM_T* spectrum);        

/**
 * Frees an allocated peak_iterator object.
 */
void free_peak_iterator(PEAK_ITERATOR_T* peak_iterator);

/**
 * The basic iterator functions. 
 * \returns TRUE if there are additional peaks to iterate over, FALSE if not.
 */
BOOLEAN_T peak_iterator_has_next(PEAK_ITERATOR_T* peak_iterator);

/**
 * \returns The next peak object in the spectrum, in order of m/z.
 */
PEAK_T* peak_iterator_next(PEAK_ITERATOR_T* peak_iterator);

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif

/** \mainpage The crux API documentation page.
 * \section Introduction
 * Welcome to crux, a C software package for analysis of tandem mass
 * spectrometry data. Click on the links above to see documentation for
 * crux objects and their user interfaces.
 */
