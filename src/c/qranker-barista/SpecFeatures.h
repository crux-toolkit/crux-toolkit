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
  int open_ms2_file_for_reading(string &ms2_filename);

  void read_ms2_file();
  void read_spectrum(string &tempstr);
  void process_observed_spectrum();
  void normalize_each_region(double max_intensity_overall, vector<double> &max_intensity_per_region,
			     int region_selector);
  void shift_peaks();

  void read_processed_ms2_file();
  void read_processed_spectrum(string &tempstr);
  void save_spec_positions(string &out_dir);
  void save_retention_times(string &out_dir);
  void load_spec_positions(string &in_dir);
  void initialize_aa_tables();

  void get_observed_spectrum(string &spec);
  void get_processed_observed_spectrum(string &spec);
  void clear_tspec(double **tspec,int num_features);
  void allocate_tspec(double ***tspec, int num_features, int max_mz);
  void zero_out_tspec(double **tspec, int num_features, int max_mz);
  void add_intensity(double *tspec, int bin, double intensity);
  double sproduct(double *tspec, vector<double> &ospec);
  
  /*
   * model m3 consists of 3 features:
   * by-ions
   * all flanking peaks
   * all neutral losses
   */
  //void get_spec_features_m3(string &spec, string &peptide, double *features);
  void get_spec_features_m3(int scan, int ch, string &peptide, double *features);
  void get_spec_features_m6(int scan, int ch, string &peptide, double *features);
  void get_spec_features_m7(int scan, int ch, string &peptide, double *features);

  void get_spec_features_aa_end(int scan, int ch, string &peptide, double *features);
  void get_spec_features_aa_mid(int scan, int ch, string &peptide, double *features);

 protected:
  ifstream f_ms2;
  streamoff pos_in_file;

  int first_scan;
  int last_scan;
  double precursor_mz;
  int charge;
  vector<int> all_charges_of_spec;
  double rtime;

  //map<string,unsigned long> spec_to_pos_in_file;
  map<string,streamoff> spec_to_pos_in_file;
  map <int, double> scan_to_rtime;

  //for the processed spectrum
  int max_mz;
  vector<double>peaks;
  //for the raw spectrum
  vector<double> mz_values;
  vector<double> intens_values;

  double **ts_main_ion;
  double **ts_m3;
  double **ts_m6;
  double **ts_m7;
 
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
  vector <double> nl_masses_mono;

};

#endif //SPECFEATURES_H
