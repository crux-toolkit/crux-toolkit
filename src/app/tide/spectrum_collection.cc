// Benjamin Diament
//
// See .h file.

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <functional>
#include "spectrum.pb.h"
#include "spectrum_collection.h"
#include "mass_constants.h"
#include "max_mz.h"
#include "records.h"
#include "records_to_vector-inl.h"
#include "util/mass.h"
#include "util/Params.h"

using namespace std;
using google::protobuf::uint64;

#define CHECK(x) GOOGLE_CHECK((x))

// Integerization constant for the XCorr p-value calculation.
#define EVIDENCE_INT_SCALE 500.0

Spectrum::Spectrum(const pb::Spectrum& spec) {
  spectrum_number_ = spec.spectrum_number();
  precursor_m_z_ = spec.precursor_m_z();
  rtime_ = spec.rtime();
  for (int i = 0; i < spec.charge_state_size(); ++i)
    charge_states_.push_back(spec.charge_state(i));
  int size = spec.peak_m_z_size();
  CHECK(size == spec.peak_intensity_size());
  ReservePeaks(size);
  uint64 total = 0;
  double m_z_denom = spec.peak_m_z_denominator();
  double intensity_denom = spec.peak_intensity_denominator();
  for (int i = 0; i < size; ++i) {
    CHECK(spec.peak_m_z(i) > 0);
    total += spec.peak_m_z(i); // deltas of m/z are stored
    peak_m_z_.push_back(total / m_z_denom);
    peak_intensity_.push_back(spec.peak_intensity(i) / intensity_denom);
  }
}

// A spectrum can have multiple precursor charges assigned.  This
// reports the maximum such charge state.
int Spectrum::MaxCharge() const {
  vector<int>::const_iterator i = max_element(charge_states_.begin(), charge_states_.end());
  return i != charge_states_.end() ? *i : 1;
}

// Report maximum intensity peak in the given m/z range.
double Spectrum::MaxPeakInRange( double min_range, double max_range ) const {
  double return_value = 0.0;

  for (int i = 0; i < this->Size(); ++i) {
    double mz = peak_m_z_[i];
    if ( (min_range <= mz) && (mz <= max_range) ) {
      double intensity = peak_intensity_[i];
      if (intensity > return_value) {
        return_value = intensity;
      }
    }
  }
  return(return_value);
}

static inline bool IsInt(double x) {
  // See whether x is quite close to an integer (within 0.001).
  return fabs(x - uint64(x+0.5)) < 0.001;
}

static inline bool CheckDenom(const vector<double>& vals, int denom) {
  // See whether all vals can be accommodated by denom when rendered as a
  // fraction.
  double d_denom = denom;
  for (int i = 0; i < vals.size(); ++ i) {
    if (!IsInt(vals[i] * d_denom))
      return false;
  }
  return true;
}

static inline int GetDenom(const vector<double>& vals) {
  // See how much precision is given in the vals array. Not especially fast,
  // but fast enough. Used only for converting spectrum input format.
  const int kMaxPrecision = 10000; // store at most 3 digits of precision
  for (int precision = 1; precision < kMaxPrecision; precision *= 10)
    if (CheckDenom(vals, precision))
      return precision;
  return kMaxPrecision;
}

void Spectrum::FillPB(pb::Spectrum* spec) {
  spec->Clear();
  spec->set_spectrum_number(spectrum_number_);
  if (precursor_m_z_ > 0)
    spec->set_precursor_m_z(precursor_m_z_);
  spec->set_rtime(rtime_);
  for (int i = 0; i < NumChargeStates(); ++i)
    spec->add_charge_state(ChargeState(i));
  int size = peak_m_z_.size();
  CHECK(size == peak_intensity_.size());
  int m_z_denom = GetDenom(peak_m_z_);
  int intensity_denom = GetDenom(peak_intensity_);
  spec->set_peak_m_z_denominator(m_z_denom);
  spec->set_peak_intensity_denominator(intensity_denom);
  uint64 last = 0;
  for (int i = 0; i < size; ++i) {
    uint64 val = uint64(peak_m_z_[i]*m_z_denom + 0.5);
    CHECK(val > last);
    spec->add_peak_m_z(val - last);
    last = val;
    spec->add_peak_intensity(uint64(peak_intensity_[i]*intensity_denom + 0.5));
  }
}

void Spectrum::SortIfNecessary() {
  if (adjacent_find(peak_m_z_.begin(), peak_m_z_.end(), greater<double>())
      == peak_m_z_.end())
    return;

  // TODO: eliminate copy operations
  int size = Size();
  vector< pair<double, double> > pairs;
  for (int i = 0; i < size; ++i)
    pairs[i] = make_pair(peak_m_z_[i], peak_intensity_[i]);
  sort(pairs.begin(), pairs.begin() + size);
  for (int i = 0; i < size; ++i) {
    peak_m_z_[i] = pairs[i].first;
    peak_intensity_[i] = pairs[i].second;
  }
}

// Do Morpheus-style simple(-istic?) deisotoping.  "For each
// peak, lower m/z peaks are considered. If the reference peak
// lies where an expected peak would lie for a charge state from
// one to the charge state of the precursor, within mass
// tolerance, and is of lower abundance, the reference peak is
// considered to be an isotopic peak and removed."
bool Spectrum::Deisotope(int index, double deisotope_threshold) const {
  if (deisotope_threshold == 0.0) {
    return false;
  }
  double location = M_Z(index);
  double intensity = Intensity(index);
  int maxCharge = MaxCharge();
  for (int fragCharge = 1; fragCharge < maxCharge; fragCharge++) {
    double isotopic_peak = location - (ISOTOPE_SPACING / fragCharge);
    double ppm_difference = (location * deisotope_threshold) / 1e6;
    double isotopic_intensity = MaxPeakInRange(isotopic_peak - ppm_difference,
                                               isotopic_peak + ppm_difference);

    if (intensity < isotopic_intensity) {
      carp(CARP_DETAILED_DEBUG,
           "Removing isotopic peak (%g, %g) because of peak in [%g, %g] with intensity %g.",
           location, intensity, isotopic_peak - ppm_difference,
           isotopic_peak + ppm_difference, isotopic_intensity);
      return true;
    }
  }
  return false;
}

/* Calculates vector of cleavage evidence for an observed spectrum, using XCorr
 * b/y/neutral peak sets and heights.
 *
 * Written by Jeff Howbert, May, 2013 (as function createEvidenceArrayObserved).
 * Extended and modified by Jeff Howbert, October, 2013.
 * Ported to and integrated with Tide by Jeff Howbert, November, 2013.
 */
vector<double> Spectrum::CreateEvidenceVector(
  double binWidth,
  double binOffset,
  int charge,
  double pepMassMonoMean,
  int maxPrecurMass,
  long int* num_range_skipped,
  long int* num_precursors_skipped,
  long int* num_isotopes_skipped,
  long int* num_retained
) const {
  // TODO need to review these constants, decide which can be moved to parameter file
  const double maxIntensPerRegion = 50.0;
  const double BYHeight = 50.0;
  const double NH3LossHeight = 10.0;
  const double COLossHeight = 10.0;    // for creating a ions on the fly from b ions
  const double H2OLossHeight = 10.0;
  const double FlankingHeight = BYHeight / 2;;
  // TODO end need to review
  int numPeaks = Size();
  double experimentalMassCutoff = PrecursorMZ() * charge + 50.0;
  double maxIonMass = 0.0;
  double maxIonIntens = 0.0;

  // Find max ion mass and max ion intensity
  bool skipPreprocess = Params::GetBool("skip-preprocessing");
  bool remove_precursor = !skipPreprocess && Params::GetBool("remove-precursor-peak");
  double precursorMZExclude = Params::GetDouble("remove-precursor-tolerance");
  double deisotope_threshold = Params::GetDouble("deisotope");
  set<int> peakSkip;
  for (int ion = 0; ion < numPeaks; ion++) {
    double ionMass = M_Z(ion);
    double ionIntens = Intensity(ion);
    if (ionMass >= experimentalMassCutoff) {
      peakSkip.insert(ion);
      if (num_range_skipped) {
        (*num_range_skipped)++;
      }
      continue;
    } else if (remove_precursor && ionMass > PrecursorMZ() - precursorMZExclude &&
               ionMass < PrecursorMZ() + precursorMZExclude) {
      peakSkip.insert(ion);
      if (num_precursors_skipped) {
        (*num_precursors_skipped)++;
      }
      continue;
    } else if (deisotope_threshold != 0.0 && Deisotope(ion, deisotope_threshold)) {
      peakSkip.insert(ion);
      if (num_isotopes_skipped) {
        (*num_isotopes_skipped)++;
      }
      continue;
    }
    if (num_retained) {
      (*num_retained)++;
    }
    if (maxIonMass < ionMass) {
      maxIonMass = ionMass;
    }
    if (maxIonIntens < ionIntens) {
      maxIonIntens = ionIntens;
    }
  }

  // 10 bin intensity normalization 
  int regionSelector = (int)floor(MassConstants::mass2bin(maxIonMass) / (double)NUM_SPECTRUM_REGIONS);
  vector<double> intensObs(maxPrecurMass, 0);
  vector<int> intensRegion(maxPrecurMass, -1);
  for (int ion = 0; ion < numPeaks; ion++) {
    if (peakSkip.find(ion) != peakSkip.end()) {
      continue;
    }
    double ionMass = M_Z(ion);
    double ionIntens = Intensity(ion);
    int ionBin = MassConstants::mass2bin(ionMass);
    int region = (int)floor((double)(ionBin) / (double)regionSelector);
    if (region >= NUM_SPECTRUM_REGIONS) {
      region = NUM_SPECTRUM_REGIONS - 1;
    }
    intensRegion[ionBin] = region;
    if (intensObs[ionBin] < ionIntens) {
      intensObs[ionBin] = ionIntens;
    }
  }

  maxIonIntens = sqrt(maxIonIntens);
  for (vector<double>::iterator i = intensObs.begin(); i != intensObs.end(); i++) {
    *i = sqrt(*i);
    if (*i <= 0.05 * maxIonIntens) {
      *i = 0.0;
    }
  }

  vector<double> maxRegion(NUM_SPECTRUM_REGIONS, 0);
  for (int i = 0; i < maxPrecurMass; i++) {
    int reg = intensRegion[i];
    if (reg >= 0 && maxRegion[reg] < intensObs[i]) {
      maxRegion[reg] = intensObs[i];
    }
  }
  for (int i = 0; i < maxPrecurMass; i++) {
    int reg = intensRegion[i];
    if (reg >= 0 && maxRegion[reg] > 0.0) {
      intensObs[i] *= (maxIntensPerRegion / maxRegion[reg]);
    }
  }

  // ***** Adapted from tide/spectrum_preprocess2.cc.
  // TODO replace, if possible, with call to
  // static void SubtractBackground(double* observed, int end).
  // Note numerous small changes from Tide code.
  vector<double> partial_sums;
  partial_sums.reserve(maxPrecurMass);
  double total = 0.0;
  for (vector<double>::const_iterator i = intensObs.begin(); i != intensObs.end(); i++) {
    partial_sums.push_back(total += *i);
  }
  const double multiplier = 1.0 / (MAX_XCORR_OFFSET * 2.0 + 1.0);
  for (int i = 0; i < maxPrecurMass; ++i) {
    int right = std::min(maxPrecurMass - 1, i + MAX_XCORR_OFFSET);
    int left = std::max(0, i - MAX_XCORR_OFFSET - 1);
    intensObs[i] -= multiplier * (partial_sums[right] - partial_sums[left]);
  }

  bool flankingPeaks = Params::GetBool("use-flanking-peaks");
  bool nlPeaks = Params::GetBool("use-neutral-loss-peaks");
  int binFirst = MassConstants::mass2bin(30);
  int binLast = MassConstants::mass2bin(pepMassMonoMean - 47);
  vector<double> evidence(maxPrecurMass, 0);
  for (int i = binFirst; i <= binLast; i++) {
    // b ion
    double bIonMass = (i - 0.5 + binOffset) * binWidth;
    int ionBin = MassConstants::mass2bin(bIonMass);
    evidence[i] += intensObs[ionBin] * BYHeight;
    for (int j = 2; j < charge; j++) {
      evidence[i] += intensObs[MassConstants::mass2bin(bIonMass, j)] * BYHeight;
    }
    // y ion
    double yIonMass = pepMassMonoMean + 2 * MASS_H_MONO - bIonMass;
    ionBin = MassConstants::mass2bin(yIonMass);
    evidence[i] += intensObs[ionBin] * BYHeight;
    for (int j = 2; j < charge; j++) {
      evidence[i] += intensObs[MassConstants::mass2bin(yIonMass, j)] * BYHeight;
    }
    if (flankingPeaks) {
      // flanking peaks for b ions
      ionBin = MassConstants::mass2bin(bIonMass, 1);
      evidence[i] += intensObs[ionBin + 1] * FlankingHeight;
      evidence[i] += intensObs[ionBin - 1] * FlankingHeight;
      for (int j = 2; j < charge; j++) {
        evidence[i] += intensObs[MassConstants::mass2bin(bIonMass, j) + 1] * FlankingHeight;
        evidence[i] += intensObs[MassConstants::mass2bin(bIonMass, j) - 1] * FlankingHeight;
      }
      // flanking peaks for y ions
      ionBin = MassConstants::mass2bin(yIonMass, charge);
      evidence[i] += intensObs[ionBin + 1] * FlankingHeight;
      evidence[i] += intensObs[ionBin - 1] * FlankingHeight;
      for (int j = 2; j < charge; j++) {
        evidence[i] += intensObs[MassConstants::mass2bin(yIonMass, j) + 1] * FlankingHeight;
        evidence[i] += intensObs[MassConstants::mass2bin(yIonMass, j) - 1] * FlankingHeight;
      }
    }
    if (nlPeaks) {
      // NH3 loss from b ion
      double ionMassNH3Loss = bIonMass - MASS_NH3_MONO;
      ionBin = MassConstants::mass2bin(ionMassNH3Loss);
      evidence[i] += intensObs[ionBin] * NH3LossHeight;
      for (int j = 2; j < charge; j++) {
        evidence[i] += intensObs[MassConstants::mass2bin(ionMassNH3Loss, j)] * NH3LossHeight;
      }
      // NH3 loss from y ion
      ionMassNH3Loss = yIonMass - MASS_NH3_MONO;
      ionBin = MassConstants::mass2bin(ionMassNH3Loss);
      evidence[i] += intensObs[ionBin] * NH3LossHeight;
      for (int j = 2; j < charge; j++) {
        evidence[i] += intensObs[MassConstants::mass2bin(ionMassNH3Loss, j)] * NH3LossHeight;
      }
      // CO and H2O loss from b ion
      double ionMassCOLoss = bIonMass - MASS_CO_MONO;
      double ionMassH2OLoss = bIonMass - MASS_H2O_MONO;
      evidence[i] += intensObs[MassConstants::mass2bin(ionMassCOLoss)] * COLossHeight;
      evidence[i] += intensObs[MassConstants::mass2bin(ionMassH2OLoss)] * H2OLossHeight;
      for (int j = 2; j < charge; j++) {
        evidence[i] += intensObs[MassConstants::mass2bin(ionMassCOLoss, j)] * COLossHeight;
        evidence[i] += intensObs[MassConstants::mass2bin(ionMassH2OLoss, j)] * H2OLossHeight;
      }
      // H2O loss from y ion
      ionMassH2OLoss = yIonMass - MASS_H2O_MONO;
      evidence[i] += intensObs[MassConstants::mass2bin(ionMassH2OLoss)] * H2OLossHeight;
      for (int j = 2; j < charge; j++) {
        evidence[i] += intensObs[MassConstants::mass2bin(ionMassH2OLoss, j)] * H2OLossHeight;
      }
    }
  }
  return evidence;
}

vector<int> Spectrum::CreateEvidenceVectorDiscretized(
  double binWidth,
  double binOffset,
  int charge,
  double pepMassMonoMean,
  int maxPrecurMass,
  long int* num_range_skipped,
  long int* num_precursors_skipped,
  long int* num_isotopes_skipped,
  long int* num_retained
) const {
  vector<double> evidence =
    CreateEvidenceVector(binWidth, binOffset, charge, pepMassMonoMean, maxPrecurMass,
                         num_range_skipped, num_precursors_skipped, num_isotopes_skipped, num_retained);
  vector<int> discretized;
  discretized.reserve(evidence.size());
  for (vector<double>::const_iterator i = evidence.begin(); i != evidence.end(); i++) {
    discretized.push_back((int)floor(*i / EVIDENCE_INT_SCALE + 0.5));
  }
  return discretized;
}

void SpectrumCollection::ReadMS(istream& in, bool ms1) {
  // Parse MS2 file format.
  // Not very fast: uses scanf. CONSIDER speed-up.
  static const int kMaxLine = 1000;
  char line[kMaxLine];
  Spectrum* spectrum = NULL;
  while (in.getline(line, kMaxLine)) {
    switch(line[0]) {
    case 'S': {
        if (spectrum)
	  spectra_.push_back(spectrum);
	int specnum1, specnum2;
	double precursor_m_z = 0;
	int ok1 = 0;
	int ok2 = 0;
	sscanf(line, "S %d %d%n%lf%n", &specnum1, &specnum2, &ok1,
	       &precursor_m_z, &ok2);
	CHECK((ok1 > 0) && (ms1 != (ok2 > 0)));
	CHECK(specnum1 == specnum2);
	spectrum = new Spectrum(specnum1, precursor_m_z);
      }
      break;
    case 'I': {
        int pos = 0;
	double rtime = -1;
	sscanf(line, "I RTime %n%lf", &pos, &rtime);
	if (pos > 0) {
	  CHECK(rtime >= 0);
	  spectrum->SetRTime(rtime);
	}
      }
      break;
    case 'Z':
      int charge;
      CHECK(1 == sscanf(line, "Z %d", &charge));
      spectrum->AddChargeState(charge);
      break;
    case '0': case '1': case '2': case '3': case '4':
    case '5': case '6': case '7': case '8': case '9':
      double location, intensity;
      CHECK(2 == sscanf(line, "%lf %lf", &location, &intensity));
      if (spectrum->Size() > 0) // check that m_z values are increasing.
	CHECK(location > spectrum->M_Z(spectrum->Size() - 1));
      if (intensity > 0)
	spectrum->AddPeak(location, intensity);
      break;
    default:
      break;
    }
  }
  if (spectrum)
    spectra_.push_back(spectrum);
}

bool SpectrumCollection::ReadSpectrumRecords(const string& filename,
					     pb::Header* header) {
  pb::Header tmp_header;
  if (header == NULL)
    header = &tmp_header;
  HeadedRecordReader reader(filename, header);
  if (header->file_type() != pb::Header::SPECTRA)
    return false;
  pb::Spectrum pb_spectrum;
  while (!reader.Done()) {
    reader.Read(&pb_spectrum);
    spectra_.push_back(new Spectrum(pb_spectrum));
  }
  if (!reader.OK()) {
    for (int i = 0; i < spectra_.size(); ++i)
      delete spectra_[i];
    spectra_.clear();
    return false;
  }
  return true;
}

void SpectrumCollection::MakeSpecCharges() {
  // Create one entry in the spec_charges_ array for each
  // (spectrum, charge) pair.
  int spectrum_index = 0;
  vector<Spectrum*>::iterator i = spectra_.begin();
  for (; i != spectra_.end(); ++i) {
    for (int j = 0; j < (*i)->NumChargeStates(); ++j) {
      int charge = (*i)->ChargeState(j);
      double neutral_mass = (((*i)->PrecursorMZ() - MASS_PROTON)
			     * charge);
      spec_charges_.push_back(SpecCharge(neutral_mass, charge, *i,
                                         spectrum_index));
    }
    spectrum_index++;
  }
}

double SpectrumCollection::FindHighestMZ() const {
  // Return the maximum MZ seen across all input spectra.
  double highest = 0;
  vector<Spectrum*>::const_iterator i = spectra_.begin();
  for (; i != spectra_.end(); ++i) {
    CHECK((*i)->Size() > 0) << "ERROR: spectrum " << (*i)->SpectrumNumber()
			    << " has no peaks.\n";
    double last_peak = (*i)->M_Z((*i)->Size() - 1);
    if (last_peak > highest)
      highest = last_peak;
  }
  return highest;
}

void SpectrumCollection::Sort() {
  MakeSpecCharges();
  sort(spec_charges_.begin(), spec_charges_.end());
}
