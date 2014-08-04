// Benjamin Diament
//
// See .h file.

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include "spectrum.pb.h"
#include "spectrum_collection.h"
#include "mass_constants.h"
#include "records.h"
#include "records_to_vector-inl.h"

using namespace std;
using google::protobuf::uint64;

#define CHECK(x) GOOGLE_CHECK((x))

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
      double neutral_mass = (((*i)->PrecursorMZ() - MassConstants::proton)
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
