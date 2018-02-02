// Benjamin Diament
//
// The Spectrum class represents a single spectrum, unprocessed. Instances can
// be read from and to the pb::Spectrum protocol buffer type (q.v.), which is
// relatively compact.
//
// The SpectrumCollection class represents all the spectra in the input.
// Initialize with ReadMS2() or ReadSpectrumRecords(). ReadMS2() takes an MS2
// format, ReadSpectrumRecords() takes a file of records of spectrum.proto.
//
// SpectrumCollection::Sort() creates one entry in the spec_charges_ array for
// each (spectrum, charge) pair, e.g. for the case where a spectrum has
// multiple possible charge states. The spectra have to be sorted by neutral
// mass in order to perform the rolling-window join. (See
// active_peptide_queue.{h,cc}.)
//
// SpectrumCollection::FindHighestMZ() returns the maximum MZ seen across all
// input spectra. This is cached by the MaxMZ class.

#ifndef SPECTRUM_COLLECTION_H
#define SPECTRUM_COLLECTION_H

#include <iostream>
#include <vector>
#include "header.pb.h"
#include "spectrum.pb.h"

using namespace std;

// Number of m/z regions in XCorr normalization.
#define NUM_SPECTRUM_REGIONS 10

class Spectrum {
 public:
  // Manual instantiation and specification
  Spectrum(int spectrum_number, double precursor_m_z)
    : spectrum_number_(spectrum_number), precursor_m_z_(precursor_m_z) {
  }
  void ReservePeaks(int num) {
    peak_m_z_.reserve(num);
    peak_intensity_.reserve(num);
  }
  void SetRTime(double rtime) { rtime_ = rtime; }
  void AddChargeState(int charge_state) {
    charge_states_.push_back(charge_state);
  }
  void AddPeak(double m_z, double intensity) {
    peak_m_z_.push_back(m_z);
    peak_intensity_.push_back(intensity);
  }
  
  explicit Spectrum(const pb::Spectrum& spec); // Instantiation from PB
  void FillPB(pb::Spectrum* spec);

  int SpectrumNumber() const { return spectrum_number_; }
  double PrecursorMZ() const { return precursor_m_z_; }
  double RTime() const { return rtime_; }

  int NumChargeStates() const { return charge_states_.size(); }
  int ChargeState(int index) const { return charge_states_[index]; }

  int Size() const { return peak_m_z_.size(); } // number of peaks
  double M_Z(int index) const { return peak_m_z_[index]; }
  double Intensity(int index) const { return peak_intensity_[index]; }

  void SortIfNecessary();

  bool Deisotope(int index, double deisotope_threshold) const;

  std::vector<double> CreateEvidenceVector(
    double binWidth,
    double binOffset,
    int charge,
    double pepMassMonoMean,
    int maxPrecurMass,
    long int* num_range_skipped = NULL,
    long int* num_precursors_skipped = NULL,
    long int* num_isotopes_skipped = NULL,
    long int* num_retained = NULL) const;
  std::vector<int> CreateEvidenceVectorDiscretized(
    double binWidth,
    double binOffset,
    int charge,
    double pepMassMonoMean,
    int maxPrecurMass,
    long int* num_range_skipped = NULL,
    long int* num_precursors_skipped = NULL,
    long int* num_isotopes_skipped = NULL,
    long int* num_retained = NULL) const;

  int MaxCharge() const;
  double MaxPeakInRange( double min_range, double max_range ) const;
  
 private:
  int spectrum_number_;
  double rtime_;
  double precursor_m_z_;
  vector<int> charge_states_;

  vector<double> peak_m_z_;
  vector<double> peak_intensity_;
};

class SpectrumCollection {
 public:
  ~SpectrumCollection() {
    for (int i = 0; i < spectra_.size(); ++i)
      delete spectra_[i];
  }

  void ReadMS(istream& in, bool ms1);
  bool ReadSpectrumRecords(const string& filename, pb::Header* header = NULL);
  void Sort();
  int Size() const { return(spectra_.size()); } // number of spectra

  template<typename BinaryPredicate>
  void Sort(BinaryPredicate Predicate) {
    MakeSpecCharges();
    sort(spec_charges_.begin(), spec_charges_.end(), Predicate);
  }

  double FindHighestMZ() const;

  struct SpecCharge {
    double neutral_mass;
    int charge;
    Spectrum* spectrum;
    int spectrum_index;

    SpecCharge(double neutral_mass_param, int charge_param,
               Spectrum* spectrum_param, int spectrum_index_param)
    : neutral_mass(neutral_mass_param), charge(charge_param),
      spectrum(spectrum_param), spectrum_index(spectrum_index_param) {
    }

    bool operator<(const SpecCharge& other) const {
      return (neutral_mass < other.neutral_mass);
    }
  };

  const vector<SpecCharge>* SpecCharges() const { return &spec_charges_; }
  vector<Spectrum*>* Spectra() { return &spectra_; }

 private:
  void MakeSpecCharges();

  vector<Spectrum*> spectra_;
  vector<SpecCharge> spec_charges_;
};

#endif // SPECTRUM_COLLECTION_H
