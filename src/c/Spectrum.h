/**
 * \file spectrum.h 
 * $Revision: 1.43 $
 * \brief Object for representing one spectrum.
 *****************************************************************************/
#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <stdio.h>
#include <vector>
#include <string>
#include "utils.h"
#include "objects.h"
#include "peak.h"

/**
 * \class Spectrum 
 * \brief A mass spectrum

 * A mass spectrum consists mainly of a list of peak objects along with
 * some identifying information. A single spectrum is generated from one 
 * or more "scans" of the mass spectrometer; each scan is identified by 
 * a unique increasing positive integer. The range of scans that
 * generated a particular spectrum are indicated by the member variables 
 * "first_scan" and "last_scan". In addition to scan information, 
 * a tandem fragmentation mass spectrum has information 
 * about the m/z of the intact ion that generated the spectrum, which is
 * indicated by "precursor_mz" member variable.
 * Also, while the m/z of particular spectrum is known, the charge state of
 * the originating ion is unknown; the possible charge states of the 
 * precursor ion is stored "possible_z" and "num_possible_z". 
 * Finally, some summary information that can be derived from the spectrum
 * peaks but is convenient to have is stored as "min_peak_mz",
 * "max_peak_mz", and "total_energy".
 */
class Spectrum{
  friend class PeakIterator; // allow the iterator access to the peaks
 protected:
  // member variables
  int              first_scan_;    ///< The number of the first scan
  int              last_scan_;     ///< The number of the last scan
  FLOAT_T          precursor_mz_;  ///< The m/z of precursor (MS-MS spectra)
  std::vector<int> possible_z_; ///< The possible charge states of this spectrum
  std::vector<PEAK_T*>  peaks_;         ///< The spectrum peaks
  FLOAT_T          min_peak_mz_;   ///< The minimum m/z of all peaks
  FLOAT_T          max_peak_mz_;   ///< The maximum m/z of all peaks
  double           total_energy_;  ///< The sum of intensities in all peaks
  std::string           filename_;      ///< Optional filename
  std::vector<std::string> i_lines_v_;  ///< store i lines
  std::vector<std::string> d_lines_v_;  ///< store d lines
  bool             has_peaks_;  ///< Does the spectrum contain peak information
  bool             sorted_by_mz_; ///< Are the spectrum peaks sorted by m/z...
  bool             sorted_by_intensity_; ///< ... or by intensity?
  bool             has_mz_peak_array_; ///< Is the mz_peak_array populated.
  PEAK_T**         mz_peak_array_;  ///< Allows rapid peak retrieval by mz.

  // constants
  /**
   * m/z resolution.  I.e., 5 == 0.2 m/z units
   */
  static const int MZ_TO_PEAK_ARRAY_RESOLUTION = 5;
  static const int MAX_PEAK_MZ = 5000; ///< Maximum possible m/z value.
  static const int MAX_CHARGE = 6;     ///< Maximum allowed charge.
  
  // private methods
  /**
   * Parses the 'S' line of a spectrum
   * \returns TRUE if success. FALSE is failure.
   */
  bool parse_S_line
    (char* line, ///< 'S' line to parse -in
     int buf_length ///< line length -in
     );

  /**
   * Parses the 'Z' line of the a spectrum
   * \returns TRUE if success. FALSE is failure.
   */
  bool parse_Z_line(char* line);  ///< 'Z' line to parse -in

  /**
   * Parses the 'D' line of the a spectrum
   * \returns TRUE if success. FALSE is failure.
   */
  bool parse_D_line(char* line);  ///< 'D' line to parse -in

  /**
   * Parses the 'I' line of the a spectrum
   * \returns TRUE if success. FALSE is failure.
   */
  bool parse_I_line(char* line);  ///< 'I' line to parse -in

  /**
   * Updates num_peaks, min_peak_mz, max_peak_mz, total_energy fields.
   */
  void update_fields
    (FLOAT_T intensity,///< the intensity of the peak that has been added -in
     FLOAT_T location  ///< the location of the peak that has been added -in
     );

 public:
  /**
   * Default constructor.
   */
  Spectrum();

  /**
   * Constructor initializes spectrum with given values.
   */
  Spectrum
    (int               first_scan,         ///< number of the first scan -in
     int               last_scan,          ///< number of the last scan -in
     FLOAT_T           precursor_mz,       ///< m/z of the precursor
     const std::vector<int>& possible_z,   ///< possible charge states
     const char*       filename
     );
  
  /**
   * Copy constructor.  Deep copy--allocates new peaks peak array. 
   */
  Spectrum(const Spectrum& old_spec);

  /**
   * Default destructor.
   */
  ~Spectrum();

  /**
   * Prints a spectrum object to file.
   */
  void print(FILE* file); ///< output file to print at -out

  /**
   * Prints a spectrum with the given intensities instead of the
   * observed peaks.  Assumes intensities are in m/z bins from 0 to
   * max_mz_bin.  Only prints non-zero intensities.
   */
  void print_processed_peaks
    (int cur_charge,       ///< print at this charge state
     FLOAT_T* intensities, ///< intensities of new peaks
     int max_mz_bin,       ///< num_bins in intensities
     FILE* file);          ///< print to this file

  /**
   * Prints a spectrum object to file in xml format.
   */
  void print_xml
    (FILE* file,           ///< output file to print at -out
     int charge,            ///< charge used for the search -in
     int index);            ///< used to output index to file

  /**
   * Prints a spectrum object to file in sqt format.
   */
  void print_sqt
    (FILE* file,           ///< output file to print at -out
     int num_matches,      ///< number of peptides compared to this spec -in
     int charge            ///< charge used for the search -in
     );

  /**
   * Parses a spectrum from a file, either mgf or ms2.
   */
  static Spectrum* parse_file(FILE* file, const char* filename);

  /**
   * Parses a spectrum from file, deleting any existing values in the
   * object.
   * Skips Header line "H"
   * \returns The newly allocated spectrum or NULL on error or EOF.
   */
  static Spectrum* parse_ms2(FILE* file, const char* filename); 

  /**
   * Parse a spectrum from a file in mgf format.
   * \returns A newly allocated spectrum or NULL on error or EOF.
   */
  static Spectrum* parse_mgf(FILE* file, const char* filename); 

  /**
   * Transfer values from an MSToolkit spectrum to the crux Spectrum.
   */
  // FIXME: pass an MSToolkit spectrum instead of void*
  bool parse_mstoolkit_spectrum(void* mst_spectrum, const char* filename);

  /**
   * Parse the spectrum from the tab-delimited result file
   *\returns the parsed spectrum , else returns NULL for failed parse
   */
  static Spectrum* parse_tab_delimited(MatchFileReader& file); 
  
  /**
   * Normalize peak intensities so that they sum to unity.
   */
  void sum_normalize();

  /**
   * Populate peaks with rank information.
   */
  void rank_peaks();

  /**
   * \returns The number of the first scan.
   */
  int get_first_scan();

  /**
   * \returns The number of the last scan.
   */
  int get_last_scan();

  /**
   * \returns The m/z of the precursor.
   */
  FLOAT_T get_precursor_mz();

  /**
   * \returns The a const reference to a vector of the possible charge
   * states of this spectrum.
   */
  const std::vector<int>& get_possible_z();

  /**
   * Considers the spectrum-charge parameter and returns the
   * appropriate charge states that should be searched for this
   * spectrum: all of them or the one selected by the parameter.
   * /returns A vector of charge states to consider for this spectrum.
   */ 
  std::vector<int> get_charges_to_search();
  
  /**
   * \returns The number of possible charge states of this spectrum.
   */
  int get_num_possible_z();

  /**
   * \returns The minimum m/z of all peaks.
   */
  FLOAT_T get_min_peak_mz();
  
  /**
   * \returns The maximum m/z of all peaks.
   */
  FLOAT_T get_max_peak_mz();
  
  /**
   * \returns The number of peaks.
   */
  int get_num_peaks();

  /**
   * \returns The closest PEAK_T within 'max' of 'mz' in 'spectrum'
   * NULL if no peak within 'max'
   * This should lazily create the data structures within the
   * spectrum object that it needs.
   */
  PEAK_T* get_nearest_peak
    (FLOAT_T mz, ///< the mz of the peak around which to sum intensities -in
     FLOAT_T max ///< the maximum distance to get intensity -in
     );
  
  /**
   * \returns The sum of intensities in all peaks.
   */
  double get_total_energy();

  /**
   * \returns The intensity of the peak with the maximum intensity.
   */
  FLOAT_T get_max_peak_intensity();

  /**
   * \returns The mass of the charged precursor ion, according to the formula 
   * mass = m/z * charge
   */
  FLOAT_T get_mass(int charge); ///< the charge of precursor ion -in

  /**
   * \returns The mass of the neutral precursor ion, according to the formula 
   * mass = m/z * charge - mass_H * charge
   */
  FLOAT_T get_neutral_mass(int charge); ///< the charge of precursor ion -in

  /**
   * \returns The mass of the singly charged precursor ion, according
   * to the formula mass = m/z * charge - (mass_H * (charge - 1))
   */
  FLOAT_T get_singly_charged_mass(int charge); ///< the charge of the precursor ion -in

  /**
   * Adds a possible charge(z) to the spectrum.
   */
  bool add_possible_z(int charge);  ///< charge to add

  /**
   * Adds a peak to the spectrum given a intensity and location.
   * Calls update_fields.
   * \returns TRUE if successfully added.
   */
  bool add_peak
    (FLOAT_T intensity,  ///< the intensity of peak to add -in
     FLOAT_T location_mz ///< the location of peak to add -in
     );

  /**
   * Creates and fills mz_peak_array_, the array of pointers to peaks
   * in the Spectrum's vector of peaks.  Peaks in the array are
   * indexed by ???
   */
  void populate_mz_peak_array();

};    


/**
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
