#ifndef SPECFEATURES_H
#define SPECFEATURES_H

#include<iostream>
#include<fstream>
#include <map>
#include <sstream>
#include <vector>
#include <assert.h>
#include <cstring>
#include <math.h>
#include "carp.h"
#include "SpectrumCollectionFactory.h"

using namespace std;

#define NUM_IONS 2
//typedef enum {B_ION, Y_ION} ion_type;
//number of flanking peaks
#define NUM_FL 2
#define NUM_NL  3

//typedef enum {H2O,NH3,CO} nl_type;
#define NUM_AA  26


class SpecFeaturesGenerator
{
 public:
  SpecFeaturesGenerator();
  ~SpecFeaturesGenerator();
  void clear();
  void read_ms2_file(const string& filename);
  void initialize_aa_tables();

  /*
   * model m3 consists of 3 features:
   * by-ions
   * all flanking peaks
   * all neutral losses
   */
  void get_spec_features_m3(int scan, int ch, string &peptide, double *features);
  void get_spec_features_m7(int scan, int ch, string &peptide, double *features);

 protected:
  void read_spectrum();
  void process_observed_spectrum();
  void normalize_each_region(
    double max_intensity_overall,
    vector<double> &max_intensity_per_region,
		int region_selector
  );
  void shift_peaks();
  void get_observed_spectrum(int scan);
  void clear_tspec(double **tspec,int num_features);
  void allocate_tspec(double ***tspec, int num_features);
  void zero_out_tspec(double **tspec, int num_features);
  void add_intensity(double *tspec, int bin, double intensity);
  double sproduct(double *tspec);

  Crux::SpectrumCollection* spectra_;
  Crux::Spectrum* spectrum_;

  double precursor_mz_;
  int charge_;

  //for the processed spectrum
  int max_mz_;
  vector<double> peaks_;
  //for the raw spectrum
  vector<double> mz_values_;
  vector<double> intens_values_;

  double **ts_m3_;
  double **ts_m7_;
 
  //mass info for neutral losses
  // While GNU g++ implements an extension that
  // allows intializing static const doubles, this
  // is not allowed by the strict C++ standard, so we
  // have to initialze the doubles in the implementation
  // to support compilers other than g++.
  static const double mass_h2o_mono;
  static const double mass_nh3_mono;
  static const double mass_co_mono;
  static const double proton_mass;
  static const double bin_width_mono;
  static const int num_regions = 10;
  static const int max_per_region = 50;
  static const int max_xcorr_offset = 75;

  vector<double> aa_masses_mono;
  vector<double> nl_masses_mono;

};

#endif //SPECFEATURES_H
