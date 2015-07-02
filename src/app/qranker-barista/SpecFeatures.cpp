#include "SpecFeatures.h"

const double SpecFeaturesGenerator::mass_h2o_mono = 18.01056;
const double SpecFeaturesGenerator::mass_nh3_mono = 17.02655;
const double SpecFeaturesGenerator::mass_co_mono = 27.9949;
const double SpecFeaturesGenerator::proton_mass = 1.00727646688;
const double SpecFeaturesGenerator::bin_width_mono = 1.0005079;

SpecFeaturesGenerator :: SpecFeaturesGenerator():
  max_mz_(1025),
  ts_m3_(NULL), ts_m7_(NULL),
  spectra_(NULL), spectrum_(NULL)
{
  mz_values_.reserve(1000);
  intens_values_.reserve(1000);
  peaks_.resize(max_mz_, 0.0);
}

SpecFeaturesGenerator :: ~SpecFeaturesGenerator()
{
  if (spectra_) {
    delete spectra_;
  }
  if (spectrum_) {
    delete spectrum_;
  }
  clear();
}

void SpecFeaturesGenerator :: clear()
{
  max_mz_ = 0;
  clear_tspec(ts_m3_, 3); ts_m3_ = NULL;
  clear_tspec(ts_m7_, 7); ts_m7_ = NULL;
}

void SpecFeaturesGenerator :: read_spectrum()
{
  precursor_mz_ = spectrum_->getPrecursorMz();
  
  mz_values_.clear();
  intens_values_.clear();

  for (PeakIterator i = spectrum_->begin(); i != spectrum_->end(); ++i) {
    Peak* peak = *i;
    FLOAT_T mz = peak->getLocation();
    mz_values_.push_back(mz);
    intens_values_.push_back(peak->getIntensity());
    /*int mz_bin = (int)(mz / bin_width_mono + 0.5);
    if (mz_bin >= max_mz_) {
      max_mz_ = mz_bin + 1;
    }*/
  }
}


void SpecFeaturesGenerator :: shift_peaks()
{
#if 1
  int n = peaks_.size();
  vector<double> sums(n);
  vector<double> integral(n+1);
  integral[0] = 0;
  for (int i=0; i<n; i++)
    integral[i+1] = integral[i] + peaks_[i];
  for (int i=0; i<n; i++)
    {
      int ilo = i - max_xcorr_offset;
      int ihi = i + max_xcorr_offset+1;
      sums[i] = integral[ihi<=n ? ihi : n] - integral[ilo>=0 ? ilo : 0]; 
    }
#else
  vector<double> sums;
  sums.resize(peaks_.size(),0.0);
  
  for(unsigned int idx = 0; idx < peaks_.size(); idx++)
    {
      for(int sub_idx = idx-max_xcorr_offset; sub_idx <= idx+max_xcorr_offset; sub_idx++)
	if(sub_idx > -1 && sub_idx < peaks_.size())
	  sums[idx] += peaks_[sub_idx];
    }
#endif
  for(unsigned int idx = 0; idx < peaks_.size(); idx++)
    peaks_[idx] -= (sums[idx]/(max_xcorr_offset*2.0+1));

}

void SpecFeaturesGenerator :: normalize_each_region(double max_intensity_overall, vector<double> &max_intensity_per_region,
						    int region_selector)
{
  int region_idx = 0;
  double max_intensity = max_intensity_per_region[region_idx];
  for (int i = 0; i < (int)peaks_.size(); i++)
    {
    
      if(i >= (region_idx+1)*region_selector && (region_idx+1)<num_regions)
	{
	  region_idx++;
	  max_intensity = max_intensity_per_region[region_idx];
	}
      // Don't normalize if no peaks in region, and for compatibility
      // with SEQUEST drop peaks with intensity less than 1/20 of
      // the overall max intensity.
      if((max_intensity != 0)
	&& (peaks_[i] > 0.05 * max_intensity_overall))
	{
	  // normalize intensity to max 50
	  peaks_[i] = (peaks_[i] /max_intensity) * max_per_region;
	}
      // no more peaks beyond the 10 regions mark, exit
      if(i > num_regions * region_selector){
	return;
      }
    }
}


void SpecFeaturesGenerator :: process_observed_spectrum()
{
  double experimental_mass_cut_off = precursor_mz_ * charge_ + 50.0;

  if ((int)peaks_.size() != max_mz_)
    peaks_.resize(max_mz_);
  peaks_.assign(peaks_.size(), 0.0);

  assert(mz_values_.size() == intens_values_.size());
  double max_peak_location = 0.0;
  double max_peak_intensity = 0.0;
  for (unsigned int i = 0; i < mz_values_.size(); i++) {
    double peak_location = mz_values_[i];
    if (peak_location < experimental_mass_cut_off &&
        peak_location > max_peak_location)
      max_peak_location = peak_location;
  }
  int region_selector = (int)max_peak_location / num_regions;
  vector<double> max_peak_intensity_per_region;
  max_peak_intensity_per_region.resize(num_regions,0.0);
  
  for (unsigned int i = 0; i < mz_values_.size(); i++)
  {
    double peak_location = mz_values_[i];
    double peak_intensity = intens_values_[i];
    //if above experimental mass, skip
    if (peak_location > experimental_mass_cut_off)
      continue;
    // skip all peaks within precursor ion mz +/- 15
    if (peak_location < precursor_mz_ + 15 &&
        peak_location > precursor_mz_ - 15)
      continue;

    //get the bin and the region
    int mz = (int)(peak_location / bin_width_mono + 0.5);

    int region = mz / region_selector;
    // don't let index beyond array
    if (region >= num_regions) {
      continue;
    }

    //update max mz over all spectra
    /*if (mz >= max_mz_) {
      max_mz_ = mz + 1;
      peaks_.resize(max_mz_,0.0);
    }*/
    // sqrt the original intensity
    peak_intensity = sqrt(peak_intensity);
    if (peak_intensity > max_peak_intensity)
      max_peak_intensity = peak_intensity;
    if (peaks_[mz] < peak_intensity) {
      peaks_[mz] = peak_intensity;
      // check if this peak is max intensity in the region(one out of 10)
      if (max_peak_intensity_per_region[region] < peak_intensity) {
        max_peak_intensity_per_region[region] = peak_intensity;
      }
    }
  }
  normalize_each_region(max_peak_intensity, max_peak_intensity_per_region,
                        region_selector);
  shift_peaks();
}

void SpecFeaturesGenerator :: read_ms2_file(const string& filename)
{
  spectra_ = SpectrumCollectionFactory::create(filename.c_str());
  spectra_->parse();

  // TODO hack for finding max m/z out of all spectra
  for (SpectrumIterator i = spectra_->begin(); i != spectra_->end(); ++i) {
    int mz_bin = (int)((*i)->getMaxPeakMz() / bin_width_mono + 0.5);
    if (mz_bin >= max_mz_) {
      max_mz_ = mz_bin + 1;
    }
  }
}

/*************************************************************************************/

void SpecFeaturesGenerator :: get_observed_spectrum(int scan)
{
  if (spectrum_) {
    delete spectrum_;
    spectrum_ = NULL;
  }

  for (SpectrumIterator i = spectra_->begin(); i != spectra_->end(); ++i) {
    if (scan == (*i)->getFirstScan()) {
      spectrum_ = new Crux::Spectrum();
      spectrum_->copyFrom(*i);
      break;
    }
  }
  // TODO this doesn't work because subclasses of spectrumcollection override
  //spectrum_ = spectra_->getSpectrum(scan);

  // not found in MS2 file
  if (!spectrum_)
    carp(CARP_FATAL, "Spectrum \"%d\" not found in MS2 file", scan);

  read_spectrum();
  process_observed_spectrum();
}


void SpecFeaturesGenerator :: clear_tspec(double **tspec, int num_features)
{
  if (tspec) {
    for (int i = 0; i < num_features; i++)
	    if (tspec[i]) {
	      delete [] tspec[i];
        tspec[i] = NULL;
      }
    delete [] tspec;
    tspec = NULL;
  }
}

void SpecFeaturesGenerator :: allocate_tspec(double ***tspec, int num_features)
{
  (*tspec) = new double*[num_features];
  for(int i = 0; i < num_features; i++)
    (*tspec)[i] = new double[max_mz_];
}

void SpecFeaturesGenerator :: zero_out_tspec(double **tspec, int num_features)
{
  for(int i = 0; i < num_features; i++)
    memset(tspec[i], 0, sizeof(double) * max_mz_);
}

void SpecFeaturesGenerator :: add_intensity(double *tspec, int idx, double intensity)
{
  if(idx > -1 && idx < max_mz_ &&
     tspec[idx] < intensity) {
	  tspec[idx] = intensity;
  }
}

double SpecFeaturesGenerator :: sproduct(double *tspec)
{
  double sm = 0;
  for (unsigned int i = 0; i < peaks_.size(); i++)
    sm += tspec[i] * peaks_[i];
  return sm;
}

void SpecFeaturesGenerator :: get_spec_features_m3(int scan, int ch, string &peptide, double *features)
{
  charge_ = ch;
  get_observed_spectrum(scan);
  
  //allocate the theoretical spectrum if necessary
  if (!ts_m3_)
    allocate_tspec(&ts_m3_, 3);
  zero_out_tspec(ts_m3_, 3);
  
  double mz = 0;
  int idx = 0;
  double fragment_mass;
  
  //start from the beginning of the peptide and do the B-ION
  fragment_mass = 0.0;
  for(unsigned int i = 0; i < peptide.size()-1;i++)
    {
      char aa = peptide.at(i);
      if(aa != '.')
	{
	  //do the B-ION
	  fragment_mass += aa_masses_mono[aa-'A'];
	  //add the fragment
	  for (int charg = 1; charg < ch; charg++)
	    {
	      //add the ion itself
	      mz = (fragment_mass+charg*proton_mass)/charg;
	      idx = (int)(mz/bin_width_mono+0.5);
	      add_intensity(ts_m3_[0], idx, 1.0);
	      
	      //add its flanking peaks
	      add_intensity(ts_m3_[1], idx-1, 1.0);
	      add_intensity(ts_m3_[1], idx+1, 1.0);
	      
	      //add its neutral losses nh3 and h2o
	      double mz_nh3 = mz-(mass_nh3_mono/charg);
	      idx = (int)(mz_nh3/bin_width_mono+0.5);
	      add_intensity(ts_m3_[2], idx, 1.0);
	      double mz_h2o = mz-(mass_h2o_mono/charg);
	      idx = (int)(mz_h2o/bin_width_mono+0.5);
	      add_intensity(ts_m3_[2], idx, 1.0);

	      //add neutral loss co
	      double mz_co = mz-(mass_co_mono/charg);
	      idx = (int)(mz_co/bin_width_mono+0.5);
	      add_intensity(ts_m3_[2], idx, 1.0);
	    }
	}
    }

    //do the Y-ION
    fragment_mass = mass_h2o_mono;
    for(unsigned int i = peptide.size(); i > 1;i--)
      {
	char aa = peptide.at(i-1);
	if(aa != '.')
	  {
	    fragment_mass += aa_masses_mono[aa-'A'];
	    //add the fragment
	    for (int charg = 1; charg < ch; charg++)
	      {
		//add the ion itself
		mz = (fragment_mass+charg*proton_mass)/charg;
		idx = (int)(mz/bin_width_mono+0.5);
		add_intensity(ts_m3_[0], idx, 1.0);
		
		//add its flanking peaks
		add_intensity(ts_m3_[1], idx-1, 1.0);
		add_intensity(ts_m3_[1], idx+1, 1.0);
		
		//add its neutral losses
		double mz_nh3 = mz-(mass_nh3_mono/charg);
		idx = (int)(mz_nh3/bin_width_mono+0.5);
		add_intensity(ts_m3_[2], idx, 1.0);
		//double mz_h2o = mz-(mass_h2o_mono/charg);
		//idx = (int)(mz_h2o/bin_width_mono+0.5);
		//add_intensity(ts_m3[2], idx, 1.0);
	      }
	  }
      }

  features[0] = sproduct(ts_m3_[0]);
  features[1] = sproduct(ts_m3_[1]);
  features[2] = sproduct(ts_m3_[2]);

  //cout << features[0] << " " << features[1] << " " << features[2] << endl;

}

void SpecFeaturesGenerator :: get_spec_features_m7(int scan, int ch, string &peptide, double *features)
{
  charge_ = ch;
  get_observed_spectrum(scan);
  
  //allocate the theoretical spectrum if necessary
  if (!ts_m7_)
    allocate_tspec(&ts_m7_, 7);
  zero_out_tspec(ts_m7_, 7);

  /*
   * 0. b-ion
   * 1. y_ion
   * 2. flanking
   * 3. h2o_nl
   * 4. co2_nl
   * 5. nh3_nl_b
   * 6. nh3_nl_y
   */
  
  double mz = 0;
  int idx = 0;
  double fragment_mass;
  
  //start from the beginning of the peptide and do the B-ION
  fragment_mass = 0.0;
  for(unsigned int i = 0; i < peptide.size()-1;i++)
    {
      char aa = peptide.at(i);
      if(aa != '.')
	{
	  //do the B-ION
	  fragment_mass += aa_masses_mono[aa-'A'];
	  //add the fragment
	  for (int charg = 1; charg < ch; charg++)
	    {
	      //add the ion itself
	      mz = (fragment_mass+charg*proton_mass)/charg;
	      idx = (int)(mz/bin_width_mono+0.5);
	      //0 in array
	      add_intensity(ts_m7_[0], idx, 1.0);
	      
	      //add its flanking peaks (2 in array)
	      add_intensity(ts_m7_[2], idx-1, 1.0);
	      add_intensity(ts_m7_[2], idx+1, 1.0);
	      
	      //add its neutral losses nh3(5 in array) 
	      double mz_nh3 = mz-(mass_nh3_mono/charg);
	      idx = (int)(mz_nh3/bin_width_mono+0.5);
	      add_intensity(ts_m7_[5], idx, 1.0);

	      //add h2o (3 in array)
	      double mz_h2o = mz-(mass_h2o_mono/charg);
	      idx = (int)(mz_h2o/bin_width_mono+0.5);
	      add_intensity(ts_m7_[3], idx, 1.0);

	      //add neutral loss co (4 in array)
	      double mz_co = mz-(mass_co_mono/charg);
	      idx = (int)(mz_co/bin_width_mono+0.5);
	      add_intensity(ts_m7_[4], idx, 1.0);
	    }
	}
    }

    //do the Y-ION
    fragment_mass = mass_h2o_mono;
    for(unsigned int i = peptide.size(); i > 1;i--)
      {
	char aa = peptide.at(i-1);
	if(aa != '.')
	  {
	    fragment_mass += aa_masses_mono[aa-'A'];
	    //add the fragment
	    for (int charg = 1; charg < ch; charg++)
	      {
		//add the ion itself ( 1 in array)
		mz = (fragment_mass+charg*proton_mass)/charg;
		idx = (int)(mz/bin_width_mono+0.5);
		add_intensity(ts_m7_[1], idx, 1.0);
		
		//add its flanking peaks ( 2 in array)
		add_intensity(ts_m7_[2], idx-1, 1.0);
		add_intensity(ts_m7_[2], idx+1, 1.0);
		
		//add its neutral losses (6 in array)
		double mz_nh3 = mz-(mass_nh3_mono/charg);
		idx = (int)(mz_nh3/bin_width_mono+0.5);
		add_intensity(ts_m7_[6], idx, 1.0);
	      }
	  }
      }

  features[0] = sproduct(ts_m7_[0]);
  features[1] = sproduct(ts_m7_[1]);
  features[2] = sproduct(ts_m7_[2]);
  features[3] = sproduct(ts_m7_[3]);
  features[4] = sproduct(ts_m7_[4]);
  features[5] = sproduct(ts_m7_[5]);
  features[6] = sproduct(ts_m7_[6]);

  //cout << features[0] << " " << features[1] << " " << features[2] << " " << features[3] << " " << features[4] << " " << features[5] << " " << features[6] << endl;

}





void SpecFeaturesGenerator :: initialize_aa_tables()
{
  aa_masses_mono.resize(NUM_AA,0.0);
  aa_masses_mono['A' - 'A'] = 71.03711;
  aa_masses_mono['B' - 'A'] = 114.53494;
  aa_masses_mono['C' - 'A'] = 103.00919;
  aa_masses_mono['D' - 'A'] = 115.02694;
  aa_masses_mono['E' - 'A'] = 129.04259;
  aa_masses_mono['F' - 'A'] = 147.06841;
  aa_masses_mono['G' - 'A'] = 57.02146;
  aa_masses_mono['H' - 'A'] = 137.05891;
  aa_masses_mono['I' - 'A'] = 113.08406;
  aa_masses_mono['K' - 'A'] = 128.09496;
  aa_masses_mono['L' - 'A'] = 113.08406;
  aa_masses_mono['M' - 'A'] = 131.04049;
  aa_masses_mono['N' - 'A'] = 114.04293;
  aa_masses_mono['O' - 'A'] = 114.07931;
  aa_masses_mono['P' - 'A'] = 97.05276;
  aa_masses_mono['Q' - 'A'] = 128.05858;
  aa_masses_mono['R' - 'A'] = 156.10111;
  aa_masses_mono['S' - 'A'] = 87.03203;
  aa_masses_mono['T' - 'A'] = 101.04768;
  aa_masses_mono['U' - 'A'] = 150.04344;
  aa_masses_mono['V' - 'A'] = 99.06841;
  aa_masses_mono['W' - 'A'] = 186.07931;
  aa_masses_mono['X' - 'A'] = 113.08406;
  aa_masses_mono['Y' - 'A'] = 163.06333;
  aa_masses_mono['Z' - 'A'] = 128.55059;
  
  //nl_masses_mono.resize(NUM_NL,0.0);
  //nl_masses_mono[H2O] = mass_h2o_mono;
  //nl_masses_mono[NH3] = mass_nh3_mono;
  //nl_masses_mono[CO] = mass_co_mono;

}

