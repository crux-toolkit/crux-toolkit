/*************************************************************************//**
 * \file Spectrum.cpp
 * AUTHOR: Chris Park, cpp-ified by Barbara Frewen
 * CREATE DATE:  June 22 2006, turned into a class Sept 21, 2010
 * \brief code to support working with spectra
 ****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "objects.h"
#include "Spectrum.h"
#include "utils.h"
#include "mass.h"
#include "parameter.h"
#include "Scorer.h"
#include "carp.h"
#include <vector>
#include <string>
#include "DelimitedFile.h"
#include "MatchFileReader.h"
#include "MSToolkit/Spectrum.h"

using namespace std;
using namespace Crux;
namespace pzd = pwiz::msdata;

/**
 * Default constructor.
 */
Spectrum::Spectrum() :
   first_scan_(0),
   last_scan_(0),
   precursor_mz_(0),
   min_peak_mz_(0),
   max_peak_mz_(0),
   total_energy_(0),
   lowest_sp_(0),
   has_peaks_(false),
   sorted_by_mz_(false),
   sorted_by_intensity_(false),
   has_mz_peak_array_(false)
{
  mz_peak_array_ = NULL;
}

/**
 * Constructor initializes spectrum with given values.
 */ 
Spectrum::Spectrum
(int               first_scan,   ///< The number of the first scan -in
 int               last_scan,    ///< The number of the last scan -in
 FLOAT_T           precursor_mz, ///< The m/z of the precursor 
 const vector<int>& possible_z,  ///< The possible charge states 
 const char*       filename      ///< Optional filename
 ) : 
   first_scan_(first_scan),
   last_scan_(last_scan),
   precursor_mz_(precursor_mz), 
   min_peak_mz_(0),
   max_peak_mz_(0),
   total_energy_(0),
   filename_(filename),
   has_peaks_(false),
   sorted_by_mz_(false),
   sorted_by_intensity_(false),
   has_mz_peak_array_(false)
 {
  mz_peak_array_ = NULL;

  for (unsigned int idx=0;idx<possible_z.size();idx++) {
    SpectrumZState zstate;
    zstate.setMZ(precursor_mz, possible_z.at(idx));
    zstates_.push_back(zstate);
  }


}

/**
 * Default destructor.
 */
Spectrum::~Spectrum()
{
  free_peak_vector(peaks_);
  
  if(has_mz_peak_array_){
    delete [] mz_peak_array_;
  }
}

/**
 * \returns the peak iterator that signifies the start of the peaks 
 * in the spectrum
 */
PeakIterator Spectrum::begin() const {

  return peaks_.begin();
}

/**
 * \returns the peak iterator that signifies the end of the peaks 
 * in the spectrum
 */
PeakIterator Spectrum::end() const {
  return peaks_.end();
}

/**
 * Prints a spectrum object to file in ms2 format.
 */
void Spectrum::print(FILE* file) ///< output file to print at -out
{
  int mass_precision = get_int_parameter("mass-precision");
  fprintf(file, "S\t%06d\t%06d\t%.*f\n", 
         first_scan_,
         last_scan_,
         mass_precision,
         precursor_mz_);

  // print 'I' line
  for(size_t line_idx = 0; line_idx < i_lines_v_.size(); line_idx++){
    fprintf(file, "%s\n", (i_lines_v_[line_idx]).c_str());
  }
  
  // print 'Z', 'D' line
  for(size_t z_idx = 0; z_idx < zstates_.size(); z_idx++){
    fprintf(file, "Z\t%d\t%.*f\n", zstates_[z_idx].getCharge(), mass_precision,
            zstates_[z_idx].getSinglyChargedMass());
    // are there any 'D' lines to print?
    if(z_idx < d_lines_v_.size() ){
      fprintf(file, "%s", d_lines_v_[z_idx].c_str());
    }
  }
  
  // print 'EZ' line
  for (size_t ez_idx = 0; ez_idx < ezstates_.size(); ez_idx++) {
    fprintf(file, "I\tEZ\t%d\t%.4f\t%.4f\t%.4f\n", ezstates_[ez_idx].getCharge(),
      ezstates_[ez_idx].getSinglyChargedMass(),
      ezstates_[ez_idx].getRTime(),
      ezstates_[ez_idx].getArea());
  }

  if (zstates_.size() == 0 && ezstates_.size() != 0) {
    for (size_t ez_idx = 0; ez_idx < ezstates_.size(); ez_idx++) {
      fprintf(file, "Z\t%d\t%.*f\n", ezstates_[ez_idx].getCharge(), mass_precision,
      ezstates_[ez_idx].getSinglyChargedMass());
    }
  }

  // print peaks
  for(int peak_idx = 0; peak_idx < (int)peaks_.size(); ++peak_idx){
    fprintf(file, "%.*f %.4f\n",
            mass_precision,
            peaks_[peak_idx]->getLocation(),
            peaks_[peak_idx]->getIntensity());
  }
}

/**
 * Prints a spectrum in ms2 format with the given intensities instead of the
 * observed peaks.  Assumes intensities are in m/z bins from 0 to
 * max_mz_bin.  Only prints non-zero intensities.
 */
void Spectrum::printProcessedPeaks(
  SpectrumZState& zstate,           ///< print at this charge state
  FLOAT_T* intensities, ///< intensities of new peaks
  int max_mz_bin,       ///< num_bins in intensities
  FILE* file){          ///< print to this file

  int mass_precision = get_int_parameter("mass-precision");

  // print S line
  fprintf(file, "S\t%06d\t%06d\t%.*f\n", 
          first_scan_,
          last_scan_,
          mass_precision,
          precursor_mz_);

  // print I line(s)
  for(size_t line_idx = 0; line_idx < i_lines_v_.size(); line_idx++){
    fprintf(file, "%s\n", (i_lines_v_[line_idx]).c_str());
  }

  // print 'Z', 'D' line
  if( zstate.getCharge() != 0 ){  // print only one charge state
    fprintf(file, "Z\t%d\t%.*f\n", zstate.getCharge(), mass_precision,
            zstate.getSinglyChargedMass());
    // TODO find associated Z line and print
  } else {  // print all charge states

    for(size_t z_idx = 0; z_idx < zstates_.size(); z_idx++){
      fprintf(file, "Z\t%d\t%.*f\n", zstates_[z_idx].getCharge(), mass_precision,
              zstates_[z_idx].getSinglyChargedMass());
      // are there any 'D' lines to print?
      if(z_idx < d_lines_v_.size()){
        fprintf(file, "%s", d_lines_v_[z_idx].c_str());
      }
    }
  }

  // print peaks
  for(int bin_idx = 0; bin_idx < max_mz_bin; bin_idx++){
    if( intensities[bin_idx] != 0 ){
      fprintf(file, "%d %.*f\n", bin_idx, mass_precision, 
              intensities[bin_idx]); 
    }
  }
  return;
}


/**
 * Prints a spectrum object to file in sqt format.
 */
void Spectrum::printSqt(
  FILE* file,           ///< output file to print to -out
  int num_matches,      ///< number of peptides compared to this spec -in
  SpectrumZState& zstate            ///< charge used for the search -in
  ){

  fprintf(file,
          "S\t%d\t%d\t%d\t%.1f\t%s\t%.*f\t%.2f\t%.*f\t%d\n", 
          first_scan_, 
          last_scan_,
          zstate.getCharge(), 
          0.0, // FIXME dummy <process time>
          "server", // FIXME dummy <server>
          get_int_parameter("mass-precision"),
          zstate.getSinglyChargedMass(), //this is used in search
          total_energy_,
          get_int_parameter("precision"),
          lowest_sp_,
          num_matches);
}

/**
 * Copy constructor.  Deep copy--allocates new peaks for peak vector.
 */
 Spectrum::Spectrum(
  const Spectrum& old_spectrum ///< the spectrum to take values from
  ) :
 first_scan_(old_spectrum.first_scan_),
 last_scan_(old_spectrum.last_scan_),
 precursor_mz_(old_spectrum.precursor_mz_),
 zstates_(old_spectrum.zstates_),
 min_peak_mz_(old_spectrum.min_peak_mz_),
 max_peak_mz_(old_spectrum.max_peak_mz_),
 total_energy_(old_spectrum.total_energy_),
 filename_(old_spectrum.filename_),
 i_lines_v_(old_spectrum.i_lines_v_),
 d_lines_v_(old_spectrum.d_lines_v_),
 has_peaks_(old_spectrum.has_peaks_),
 sorted_by_mz_(old_spectrum.sorted_by_mz_),
 sorted_by_intensity_(old_spectrum.sorted_by_intensity_),
 has_mz_peak_array_(old_spectrum.has_mz_peak_array_)
{

  // copy each peak
  for(int peak_idx=0; peak_idx < (int)old_spectrum.peaks_.size(); ++peak_idx){
    this->addPeak(old_spectrum.peaks_[peak_idx]->getIntensity(),
                  old_spectrum.peaks_[peak_idx]->getLocation());
  }

  /*  Should we do this??
  if( old_spectrum.mz_peak_array ){
    populateMzPeakArray();
  }
  */
}

void Spectrum::copyFrom(Spectrum *src) {

 first_scan_ = src->first_scan_;
 last_scan_ = src->last_scan_;
 precursor_mz_ = src->precursor_mz_;
 zstates_ = src->zstates_;
 min_peak_mz_ = src->min_peak_mz_;
 max_peak_mz_ = src->max_peak_mz_;
 filename_ = src->filename_;
 i_lines_v_  = src->i_lines_v_;
 d_lines_v_ = src-> d_lines_v_;
 has_peaks_ = src-> has_peaks_;
 sorted_by_mz_ = src->sorted_by_mz_;
 sorted_by_intensity_ = src->sorted_by_intensity_;
 has_mz_peak_array_ = src->has_mz_peak_array_;
 // copy each peak
 for(int peak_idx=0; peak_idx < (int)src->peaks_.size(); ++peak_idx){
   this->addPeak(src->peaks_[peak_idx]->getIntensity(),
                  src->peaks_[peak_idx]->getLocation());
  }

  /*  Should we do this??
  if( old_spectrum.mz_peak_array ){
    populateMzPeakArray();
  }
  */

}

/**
 * Transfer values from an MSToolkit spectrum to the crux Spectrum.
 * \returns TRUE if success. FALSE is failure.
 */
bool Spectrum::parseMstoolkitSpectrum
  (MSToolkit::Spectrum* mst_spectrum, ///< the input MSToolkit spectrum -in
  const char* filename ///< filename of the spectrum
  ) {

  // clear any existing values
  zstates_.clear();

  free_peak_vector(peaks_);
  i_lines_v_.clear();
  d_lines_v_.clear();
  if( mz_peak_array_ ){ delete [] mz_peak_array_; }

  MSToolkit::Spectrum* mst_real_spectrum = (MSToolkit::Spectrum*)mst_spectrum;

  //set first_scan, last_scan, and precursor_mz.
  first_scan_ = mst_real_spectrum->getScanNumber();
  last_scan_ = mst_real_spectrum->getScanNumber();
  precursor_mz_ = mst_real_spectrum->getMZ();

  // setfilename of empty spectrum
  filename_ = filename;

  //add all peaks.
  for(int peak_idx = 0; peak_idx < (int)mst_real_spectrum->size(); peak_idx++){
    this->addPeak(mst_real_spectrum->at(peak_idx).intensity,
                   mst_real_spectrum->at(peak_idx).mz);
  }
  
  //add possible charge states.
  if(  mst_real_spectrum->sizeZ() > 0 ){
    for (int z_idx = 0; z_idx < mst_real_spectrum -> sizeZ(); z_idx++) {
      SpectrumZState zstate;
      zstate.setSinglyChargedMass(
        mst_real_spectrum->atZ(z_idx).mz,
        mst_real_spectrum->atZ(z_idx).z);
      zstates_.push_back(zstate);
    }
  } else { // if no charge states detected, decide based on spectrum
    assignZState(); 
  }

  return true;
}

/**
 * Transfer values from a proteowizard SpectrumInfo object to the
 * crux spectrum.
 */
bool Spectrum::parsePwizSpecInfo(
  const pzd::SpectrumPtr& pwiz_spectrum,
  int firstScan,
  int lastScan
){
  // clear any existing values
  zstates_.clear();
  ezstates_.clear();
  free_peak_vector(peaks_);
  i_lines_v_.clear();
  d_lines_v_.clear();
  if( mz_peak_array_ ){ free(mz_peak_array_); }

  // assign new values
  first_scan_ = firstScan;
  last_scan_ = lastScan;

  // get peaks
  int num_peaks = pwiz_spectrum->defaultArrayLength;
  vector<double>& mzs = pwiz_spectrum->getMZArray()->data;
  vector<double>& intensities = pwiz_spectrum->getIntensityArray()->data;
  for(int peak_idx = 0; peak_idx < num_peaks; peak_idx++){
    addPeak(intensities[peak_idx], mzs[peak_idx]);
  }
  has_peaks_ = true;

  // get precursor m/z and charge
  // is there exactly one precursor?
  if( pwiz_spectrum->precursors.size() != 1 ){  
    carp(CARP_FATAL, 
         "Spectrum %d has more than one precursor.", first_scan_);
  }
  // get the isolation window as the precursor m/z
  pzd::IsolationWindow iso_window = 
                       pwiz_spectrum->precursors[0].isolationWindow;
  bool have_precursor_mz = false;
  if (iso_window.hasCVParam(pzd::MS_isolation_window_target_m_z)) {
    precursor_mz_ =
      iso_window.cvParam(pzd::MS_isolation_window_target_m_z).valueAs<double>();
    have_precursor_mz = true;
  }

  // each charge state(s) stored in selectedIon(s)
  // is there at least one selected ion?
  vector<pzd::SelectedIon> ions = pwiz_spectrum->precursors[0].selectedIons;
  if( ions.empty() ){
    carp(CARP_FATAL, "No selected ions in spectrum %d.", first_scan_);
  }

  bool knownCharge = ions[0].hasCVParam(pzd::MS_charge_state);
  bool hasZLine = false;
  const vector<pwiz::data::UserParam>& specUserParams = pwiz_spectrum->userParams;
  if (!knownCharge) {
    for (vector<pwiz::data::UserParam>::const_iterator i = specUserParams.begin();
         i != specUserParams.end();
         ++i) {
      if (i->name == "ms2 file charge state") {
        hasZLine = true;
        break;
      }
    }
  }

  // determined charge states will be stored
  // one per selected ion
  if (knownCharge) {
    carp(CARP_DEBUG, "MS_charge_state");
    // get each charge state and possibly the associated mass
    for(size_t ion_idx = 0; ion_idx < ions.size(); ion_idx++){
      int charge = ions[ion_idx].cvParam(pzd::MS_charge_state).valueAs<int>();
      carp(CARP_DEBUG, "Charge:%d", charge);
      if (ions[ion_idx].hasCVParam(pzd::MS_accurate_mass)) {
        //bullseye-determined charge states
        FLOAT_T accurate_mass = 
          ions[ion_idx].cvParam(pzd::MS_accurate_mass).valueAs<FLOAT_T>();
        carp(CARP_DEBUG, "accurate mass:%f charge:%i",accurate_mass, charge);
        SpectrumZState zstate;
        zstate.setSinglyChargedMass(accurate_mass, charge);
        ezstates_.push_back(zstate);
      } else if (ions[ion_idx].hasCVParam(pzd::MS_selected_ion_m_z)) {
        FLOAT_T mz =
          ions[ion_idx].cvParam(pzd::MS_selected_ion_m_z).valueAs<FLOAT_T>();
        carp(CARP_DEBUG, "mz:%g", mz);
        //if we don't have a precursor set yet, set it now.
        if (!have_precursor_mz) {
          precursor_mz_ = mz;
        }
        SpectrumZState zstate;
        zstate.setMZ(mz, charge);
        zstates_.push_back(zstate);
      } else {
        ostringstream oss;
        oss << "Cannot find precursor mass! CVParams:"<<endl;

        for (vector<pwiz::data::CVParam>::iterator iter = ions[0].cvParams.begin();
             iter != ions[0].cvParams.end();
             ++iter) {
          oss<<"id:"<<(int)iter->cvid<<" value:"<<iter->value;
        }
        string err_string = oss.str();
        carp(CARP_FATAL, err_string);
      }
    }
  } else if (hasZLine) {
    for (vector<pwiz::data::UserParam>::const_iterator i = specUserParams.begin();
         i != specUserParams.end();
         ++i) {
      if (i->name == "ms2 file charge state") {
        const string& zLine = i->value; // "<charge> <m/z>"
        size_t separator = zLine.find(" ");
        SpectrumZState zstate;
        zstate.setMZ(boost::lexical_cast<double>(zLine.substr(separator + 1)),
                     boost::lexical_cast<int>(zLine.substr(0, separator)));
        zstates_.push_back(zstate);
      }
    }
  } else if (ions[0].hasCVParam(pzd::MS_possible_charge_state)) {
    // possible charge states will all be stored in the first selected ion
    carp(CARP_DEBUG, "charges stored ion");
    vector<pzd::CVParam> charges = 
      ions[0].cvParamChildren(pzd::MS_possible_charge_state);
    for (size_t charge_idx = 0; charge_idx < charges.size(); charge_idx++) {
      SpectrumZState zstate;
      zstate.setMZ(precursor_mz_, charges[charge_idx].valueAs<int>());
      zstates_.push_back(zstate);
    }
  } else { // we have no charge information
    assignZState(); //do choose charge and add +1 or +2,+3
  }

  return true;
}

/**
 * Adds a peak to the spectrum given a intensity and location
 * calls update_spectrum_fields to update num_peaks, min_peak ...
 */
bool Spectrum::addPeak
( FLOAT_T intensity, ///< the intensity of peak to add -in
  FLOAT_T location_mz ///< the location of peak to add -in
  )
{

  /* TODO : why do we return a bool here? Do we need to test
   * for success?
   */ 
  if (intensity > 0) {

    Peak *peak = new Peak(intensity, location_mz);
    peaks_.push_back(peak);

    updateFields(intensity, location_mz);
    has_peaks_ = true;
  }
  return true;

}

/**
 * Creates and fills mz_peak_array_, the array of pointers to peaks
 * in the Spectrum's vector of peaks.  Peaks in the array are
 * indexed by ???
 */
void Spectrum::populateMzPeakArray()
{
  if (has_mz_peak_array_ == true){
    return;
  }
  
  int array_length = MZ_TO_PEAK_ARRAY_RESOLUTION * MAX_PEAK_MZ;
  mz_peak_array_ = new Peak * [array_length];
  for (int peak_idx = 0; peak_idx < array_length; peak_idx++){
    mz_peak_array_[peak_idx] = NULL;
  }
  for(int peak_idx = 0; peak_idx < (int)peaks_.size(); peak_idx++){
    Peak * peak = peaks_[peak_idx];
    FLOAT_T peak_mz = peak->getLocation();
    int mz_idx = (int) (peak_mz * MZ_TO_PEAK_ARRAY_RESOLUTION);
    if (mz_peak_array_[mz_idx] != NULL){
      carp(CARP_INFO, "Peak collision at mz %.3f = %i", peak_mz, mz_idx);
      if (mz_peak_array_[mz_idx]->getIntensity() < peak->getIntensity()) {
        mz_peak_array_[mz_idx] = peak;
      }
    } else {
      mz_peak_array_[mz_idx] = peak; 
    }
  }
  has_mz_peak_array_ = true;
}

/**
 * \returns The closest intensity within 'max' of 'mz' in 'spectrum'
 * NULL if no peak.
 * This should lazily create the data structures within the
 * spectrum object that it needs.
 * TODO: reimplement with faster peak lookup
 */
Peak * Spectrum::getNearestPeak(
  FLOAT_T mz, ///< the mz of the peak around which to sum intensities -in
  FLOAT_T max ///< the maximum distance to get intensity -in
  )
{
  this->populateMzPeakArray(); // for rapid peak lookup by mz

  FLOAT_T min_distance = BILLION;
  int min_mz_idx = (int)((mz - max) * MZ_TO_PEAK_ARRAY_RESOLUTION + 0.5);
  min_mz_idx = min_mz_idx < 0 ? 0 : min_mz_idx;
  int max_mz_idx = (int)((mz + max) * MZ_TO_PEAK_ARRAY_RESOLUTION + 0.5);
  int absolute_max_mz_idx = MAX_PEAK_MZ * MZ_TO_PEAK_ARRAY_RESOLUTION - 1;
  max_mz_idx = max_mz_idx > absolute_max_mz_idx 
    ? absolute_max_mz_idx : max_mz_idx;
  Peak * peak = NULL;
  Peak * nearest_peak = NULL;
  int peak_idx;
  for (peak_idx=min_mz_idx; peak_idx < max_mz_idx + 1; peak_idx++){
    if ((peak = mz_peak_array_[peak_idx]) == NULL){
      continue;
    }
    FLOAT_T peak_mz = peak->getLocation();
    FLOAT_T distance = fabs(mz - peak_mz);
    if (distance > max){
      continue;
    }
    if (distance < min_distance){
      nearest_peak = peak;
      min_distance = distance;
    }
  }
  return nearest_peak;
}

/**
 * \returns The PEAK_T within 'max' of 'mz' in 'spectrum'
 * that is the maximum intensity.
 * NULL if no peak within 'max'
 * This should lazily create the data structures within the
 * spectrum object that it needs.
 */
Peak* Spectrum::getMaxIntensityPeak(
  FLOAT_T mz, ///< the mz of the peak to find
  FLOAT_T max ///< the maximum distance to get intensity -in
  ) {

  FLOAT_T max_intensity = -BILLION;
  Peak* max_intensity_peak = NULL;

  for (PeakIterator peak_iter = begin();
    peak_iter != end();
    ++peak_iter) {

    Peak* peak = *peak_iter;
    FLOAT_T peak_mz = peak->getLocation();
    FLOAT_T distance = fabs(mz - peak_mz);
    FLOAT_T intensity = peak->getIntensity();
    if ((distance <= max) && (intensity > max_intensity)){
      max_intensity_peak = peak;
      max_intensity = intensity;
    }
  }
  return max_intensity_peak;



}




/**
 * Updates num_peaks, min_peak_mz, max_peak_mz, total_energy.
 */
void Spectrum::updateFields(
  FLOAT_T intensity, ///< the intensity of the peak that has been added -in
  FLOAT_T location ///< the location of the peak that has been added -in
  )
{
 
  // is new peak the smallest peak
  if(peaks_.size() == 1 || min_peak_mz_ > location){
    min_peak_mz_ = location;
  }
  // is new peak the largest peak
  if(peaks_.size() == 1 || max_peak_mz_ < location){
    max_peak_mz_ = location;
  }
  // update total_energy
  total_energy_ += intensity;
}

/**
 * \returns The number of the first scan.
 */
int Spectrum::getFirstScan() const
{
  return first_scan_;
}

/**
 * \returns The number of the last scan.
 */
int Spectrum::getLastScan() const
{
  return last_scan_;
}

/**
 * \returns The m/z of the precursor.
 */
FLOAT_T Spectrum::getPrecursorMz() const
{
  return precursor_mz_;
}

/**
 * \returns The minimum m/z of all peaks.
 */
FLOAT_T Spectrum::getMinPeakMz()
{
  return min_peak_mz_;
}

/**
 * \returns The maximum m/z of all peaks.
 */
FLOAT_T Spectrum::getMaxPeakMz()
{
  return max_peak_mz_;
}

/**
 * \returns The number of peaks.
 */
int Spectrum::getNumPeaks() const
{
  return (int)peaks_.size();
}


/**
 * \returns The sum of intensities in all peaks.
 */
double Spectrum::getTotalEnergy()
{
  return total_energy_;
}

/**
 * Sets the total ion current.
 */
void Spectrum::setTotalEnergy(FLOAT_T tic)
{
  total_energy_ = tic;
}

/**
 * Sets the lowest Sp score.
 */
void Spectrum::setLowestSp(FLOAT_T sp)
{
  lowest_sp_ = sp;
}

/**
 * \returns A read-only reference to the vector of possible chare
 * states for this spectrum.  If EZ states are available, return those.
 */
const vector<SpectrumZState>& Spectrum::getZStates() const {
  if (ezstates_.size() != 0) {
    return ezstates_;
  } else {
    return zstates_;
  }
}


/**
 *  Considers the spectrum-charge parameter and returns the
 *  appropriate charge states that should be searched for this
 *  spectrum: all of them or the one selected by the parameter.
 * /returns A vector of charge states to consider for this spectrum.
 */ 
vector<SpectrumZState> Spectrum::getZStatesToSearch() {

  vector<SpectrumZState> select_zstates;
  const char* charge_str = get_string_parameter_pointer("spectrum-charge");

  
  if( strcmp( charge_str, "all") == 0){ // return full array of charges
    select_zstates = getZStates();
  } else { // return a single charge state.

    int param_charge = atoi(charge_str);
    
    if( (param_charge < 1) || (param_charge > MAX_CHARGE) ){
      carp(CARP_FATAL, "spectrum-charge option must be 1,2,3,.. %d or 'all'.  "
           "'%s' is not valid", MAX_CHARGE, charge_str);
    }

    for (unsigned int zstate_idx=0;zstate_idx < getNumZStates();zstate_idx++) {
      if (getZState(zstate_idx).getCharge() == param_charge) {
        select_zstates.push_back(getZState(zstate_idx));
      }
    }
  }

  return select_zstates;

}

/**
 * \returns the ZState at the requested index
 */
const SpectrumZState& Spectrum::getZState(
  int idx ///< the zstate index
) {
  return getZStates().at(idx);
}


/**
 * \returns The number of possible charge states of this spectrum.
 */
unsigned int Spectrum::getNumZStates() const {
  return getZStates().size();
}

/**
 * \returns The intensity of the peak with the maximum intensity.
 */
FLOAT_T Spectrum::getMaxPeakIntensity()
{
  FLOAT_T max_intensity = -1;

  for(int peak_idx = 0; peak_idx < (int)peaks_.size(); ++peak_idx){
    if (max_intensity <= peaks_[peak_idx]->getIntensity()) {
      max_intensity = peaks_[peak_idx]->getIntensity();
    }
  }
  return max_intensity; 
}


/**
 * Parse the spectrum from the tab-delimited result file.
 *\returns The parsed spectrum, else returns NULL for failed parse.
 */
Spectrum* Spectrum::parseTabDelimited(
  MatchFileReader& file ///< output stream -out
  ) {

  Spectrum* spectrum = new Spectrum();

  spectrum->first_scan_ = file.getInteger(SCAN_COL);
  spectrum->last_scan_ = spectrum->first_scan_;

  spectrum->precursor_mz_ = file.getFloat(SPECTRUM_PRECURSOR_MZ_COL);
  //Is it okay to assign an individual spectrum object for each charge?

  int charge = file.getInteger(CHARGE_COL);

  FLOAT_T neutral_mass = file.getFloat(SPECTRUM_NEUTRAL_MASS_COL);
  
  SpectrumZState zstate;

  zstate.setNeutralMass(neutral_mass, charge);

  spectrum->zstates_.push_back(zstate);


  /*
  TODO : Implement these in the tab delimited file?
  spectrum -> min_peak_mz = file.getFloat("spectrum min peak mz");
  spectrum -> max_peak_mz = file.getFloat("spectrum max peak mz");
  spectrum -> num_peaks = file.getInteger("spectrum num peaks");
  spectrum -> total_energy = file.getInteger("spectrum total energy");
  */

  spectrum->has_peaks_ = false;
  return spectrum;

}

/**
 * Normalize peak intensities so that they sum to unity.
 */
void Spectrum::sumNormalize()
{
  for(int peak_idx = 0; peak_idx < (int)peaks_.size(); peak_idx++){
    Peak * peak = peaks_[peak_idx];
    FLOAT_T new_intensity = peak->getIntensity() / total_energy_;
    peak->setIntensity(new_intensity);
  }
}

/**
 * Sort peaks
 */
void Spectrum::sortPeaks(PEAK_SORT_TYPE_T type)
{
  if ((type == _PEAK_LOCATION && sorted_by_mz_) ||
      (type == _PEAK_INTENSITY && sorted_by_intensity_)) {
    return;
  }
  sort_peaks(peaks_, type);
  sorted_by_mz_ = (type == _PEAK_LOCATION);
  sorted_by_intensity_ = (type == _PEAK_INTENSITY);
}

/**
 * Populate peaks with rank information.
 */
void Spectrum::rankPeaks()
{
  sort_peaks(peaks_, _PEAK_INTENSITY);
  sorted_by_intensity_ = true;
  sorted_by_mz_ = false;
  int rank = (int)peaks_.size();
  for(int peak_idx = 0; peak_idx < (int) peaks_.size(); peak_idx++){
    Peak * peak = peaks_[peak_idx];
    FLOAT_T new_rank = rank/(float)peaks_.size();
    rank--;
    peak->setIntensityRank(new_rank);
  }

}

/**
 * \returns The name of the file (no path or extension) this spectrum
 * came from or an empty string, if filename is unavailable.
 */
const char* Spectrum::getFilename(){
  
  if( filename_.empty() ){
    return "";
  }

  // store the filename stripped of path and extension
  if( stripped_filename_.empty() ){
    char** name_ext_array = parse_filename_path_extension(filename_.c_str(), 
                                                          NULL);
    stripped_filename_ = name_ext_array[0];
    free(name_ext_array[0]);
    if(name_ext_array[1] != NULL ){
      free(name_ext_array[1]);
    }
    free(name_ext_array);
  }

  return stripped_filename_.c_str();
}

/**
 * \Determine charge state for a spectrum without Z line 
 * /return true if it can determine cahrge state and return false if it can't create z line 
 */
bool Spectrum::assignZState(){
  carp_once(CARP_WARNING,
       "Spectrum %i has no charge state.\nCalculating charge",
       first_scan_);
  
  CHARGE_STATE_T charge = choose_charge(precursor_mz_, peaks_);
  SpectrumZState zstate;
  switch(charge){
  case SINGLE_CHARGE_STATE:
    zstate.setMZ(precursor_mz_, 1);
    zstates_.push_back(zstate);
    return true; 
  case MULTIPLE_CHARGE_STATE:
    zstate.setMZ(precursor_mz_, 2);
    zstates_.push_back(zstate);
    zstate.setMZ(precursor_mz_, 3);
    zstates_.push_back(zstate);
    return true;    
  case INVALID_CHARGE_STATE: 
  case NUMBER_CHARGE_STATE:
    carp(CARP_ERROR, "Could not determine charge state for spectrum %d.", 
         first_scan_);
    return false; 
 }
  return true; 
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

