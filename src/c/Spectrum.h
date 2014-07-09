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
#include "Peak.h"

#include "Spectrum.h"
#include "MSToolkit/Spectrum.h"
#include "pwiz/data/msdata/SpectrumInfo.hpp"
#include "SpectrumZState.h"

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

namespace Crux {

class Spectrum{
 protected:
  // member variables
  int              first_scan_;    ///< The number of the first scan
  int              last_scan_;     ///< The number of the last scan
  FLOAT_T          precursor_mz_;  ///< The m/z of precursor (MS-MS spectra)
  std::vector<SpectrumZState> zstates_;
  std::vector<SpectrumZState> ezstates_;
  std::vector<Peak*>  peaks_;         ///< The spectrum peaks
  FLOAT_T          min_peak_mz_;   ///< The minimum m/z of all peaks
  FLOAT_T          max_peak_mz_;   ///< The maximum m/z of all peaks
  double           total_energy_;  ///< The sum of intensities in all peaks
  FLOAT_T          lowest_sp_;  ///< The lowest Sp score (for writing SQT)
  std::string           filename_;      ///< Optional filename
  std::string           stripped_filename_; ///< filename, no path or extension
  std::vector<std::string> i_lines_v_;  ///< store i lines
  std::vector<std::string> d_lines_v_;  ///< store d lines
  bool             has_peaks_;  ///< Does the spectrum contain peak information
  bool             sorted_by_mz_; ///< Are the spectrum peaks sorted by m/z...
  bool             sorted_by_intensity_; ///< ... or by intensity?
  bool             has_mz_peak_array_; ///< Is the mz_peak_array populated.
  Peak         **mz_peak_array_;  ///< Allows rapid peak retrieval by mz.

  // constants
  /**
   * m/z resolution.  I.e., 5 == 0.2 m/z units
   */
  static const int MZ_TO_PEAK_ARRAY_RESOLUTION = 5;
  static const int MAX_PEAK_MZ = 5000; ///< Maximum possible m/z value.
  static const int MAX_CHARGE = 6;     ///< Maximum allowed charge.
  
  // private methods
  /**
   * Updates num_peaks, min_peak_mz, max_peak_mz, total_energy fields.
   */
  void updateFields
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

  void copyFrom(Spectrum* src);


  /**
   * Default destructor.
   */
  ~Spectrum();

  /**
   * \returns the peak iterator that signifies the start of the peaks 
   * in the spectrum
   */
  PeakIterator begin() const;

  /**
   * \returns the peak iterator that signifies the end of the peaks 
   * in the spectrum
   */
  PeakIterator end() const;

  /**
   * Prints a spectrum object to file.
   */
  void print(FILE* file); ///< output file to print at -out

  /**
   * Prints a spectrum with the given intensities instead of the
   * observed peaks.  Assumes intensities are in m/z bins from 0 to
   * max_mz_bin.  Only prints non-zero intensities.
   */
  void printProcessedPeaks
    (SpectrumZState& zstate,       ///< print at this charge state
     FLOAT_T* intensities, ///< intensities of new peaks
     int max_mz_bin,       ///< num_bins in intensities
     FILE* file);          ///< print to this file

  /**
   * Prints a spectrum object to file in sqt format.
   */
  void printSqt
    (FILE* file,           ///< output file to print at -out
     int num_matches,      ///< number of peptides compared to this spec -in
     SpectrumZState& zstate            ///< charge used for the search -in
     );

  /**
   * Transfer values from an MSToolkit spectrum to the crux Spectrum.
   */
  bool parseMstoolkitSpectrum(MSToolkit::Spectrum* mst_spectrum, 
                                const char* filename = NULL);

  /**
   * Transfer values from a proteowizard Spectrum object to the
   * crux spectrum.
   */
  bool parsePwizSpecInfo(
    const pwiz::msdata::SpectrumPtr& pwiz_spectrum,
    int firstScan,
    int lastScan
  );

  /**
   * Parse the spectrum from the tab-delimited result file
   *\returns the parsed spectrum , else returns NULL for failed parse
   */
  static Spectrum* parseTabDelimited(MatchFileReader& file); 
  
  /**
   * Normalize peak intensities so that they sum to unity.
   */
  void sumNormalize();

  /**
   * Sort peaks
   */
  void sortPeaks(PEAK_SORT_TYPE_T type);

  /**
   * Populate peaks with rank information.
   */
  void rankPeaks();

  /**
   * \returns The number of the first scan.
   */
  int getFirstScan() const;

  /**
   * \returns The number of the last scan.
   */
  int getLastScan() const;

  /**
   * \returns The m/z of the precursor.
   */
  FLOAT_T getPrecursorMz() const;

  /**
   * \returns The a const reference to a vector of the possible charge
   * states of this spectrum. If EZ states are available, return those.
   */
  const std::vector<SpectrumZState>& getZStates() const;

  /**
   * \returns the ZState at the requested index
   */
  const SpectrumZState& getZState(int idx);

  /**
   * Considers the spectrum-charge parameter and returns the
   * appropriate charge states that should be searched for this
   * spectrum: all of them or the one selected by the parameter.
   * /returns A vector of charge states to consider for this spectrum.
   */ 
  //std::vector<int> getChargesToSearch();
 
  std::vector<SpectrumZState> getZStatesToSearch();

  
  /**
   * \returns The number of possible charge states of this spectrum.
   */
  unsigned int getNumZStates() const;

  /**
   * \returns The minimum m/z of all peaks.
   */
  FLOAT_T getMinPeakMz();
  
  /**
   * \returns The maximum m/z of all peaks.
   */
  FLOAT_T getMaxPeakMz();
  
  /**
   * \returns The number of peaks.
   */
  int getNumPeaks() const;

  /**
   * \returns The closest PEAK_T within 'max' of 'mz' in 'spectrum'
   * NULL if no peak within 'max'
   * This should lazily create the data structures within the
   * spectrum object that it needs.
   */
  Peak * getNearestPeak
    (FLOAT_T mz, ///< the mz of the peak around which to sum intensities -in
     FLOAT_T max ///< the maximum distance to get intensity -in
     );
  
  /**
   * \returns The PEAK_T within 'max' of 'mz' in 'spectrum'
   * that is the maximum intensity.
   * NULL if no peak within 'max'
   * This should lazily create the data structures within the
   * spectrum object that it needs.
   */
  Peak* getMaxIntensityPeak(
    FLOAT_T mz, ///< the mz of the peak to find
    FLOAT_T max ///< the maximum distance to get intensity -in
  );


  /**
   * \returns The sum of intensities in all peaks.
   */
  double getTotalEnergy();

  /**
   * Sets the total ion current.
   */
  void setTotalEnergy(
    FLOAT_T tic ///< the total ion current
  );

  /**
   * Sets the lowest Sp score.
   */
  void setLowestSp(
    FLOAT_T sp ///< the lowest Sp score for this spectrum
  );

  /**
   * \returns The intensity of the peak with the maximum intensity.
   */
  FLOAT_T getMaxPeakIntensity();

  /**
   * \returns The mass of the singly charged precursor ion, according
   * to the formula mass = m/z * charge - (mass_H * (charge - 1))
   */
  //FLOAT_T getSinglyChargedMass(int charge); ///< the charge of the precursor ion -in

  /**
   * Adds a possible charge(z) to the spectrum.
   */
  //bool addPossibleZ(int charge);  ///< charge to add

  /**
   * Adds a peak to the spectrum given a intensity and location.
   * Calls update_fields.
   * \returns TRUE if successfully added.
   */
  bool addPeak
    (FLOAT_T intensity,  ///< the intensity of peak to add -in
     FLOAT_T location_mz ///< the location of peak to add -in
     );

  /**
   * Creates and fills mz_peak_array_, the array of pointers to peaks
   * in the Spectrum's vector of peaks.  Peaks in the array are
   * indexed by ???
   */
  void populateMzPeakArray();

  /**
   *if ms2 file dose not have any Z line then assignZState will create it  
   */
  bool assignZState();

  /**
   * \returns The name of the file this spectrum came from or an empty
   * string, if unavailable.
   */
  const char* getFilename();
};    

}

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
