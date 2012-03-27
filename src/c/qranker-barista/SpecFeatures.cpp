#include "SpecFeatures.h"

const double SpecFeaturesGenerator::mass_h2o_mono = 18.01056;
const double SpecFeaturesGenerator::mass_nh3_mono = 17.02655;
const double SpecFeaturesGenerator::mass_co_mono = 27.9949;
const double SpecFeaturesGenerator::proton_mass = 1.00727646688;
const double SpecFeaturesGenerator::bin_width_mono = 1.0005079;

SpecFeaturesGenerator :: SpecFeaturesGenerator(): max_mz(1025),ts_main_ion((double**)0),
						  ts_m3((double**)0),ts_m6((double**)0),ts_m7((double**)0)
{
  peaks.resize(max_mz,0.0);
}

SpecFeaturesGenerator :: ~SpecFeaturesGenerator()
{
  //if(f_ms2.is_open())
  //f_ms2.close();
  clear();
}

void SpecFeaturesGenerator :: clear()
{
  spec_to_pos_in_file.clear();
  max_mz = 0;
  clear_tspec(ts_m3,3); ts_m3 = 0;
  clear_tspec(ts_m6,6); ts_m6 = 0;
  clear_tspec(ts_m7,7); ts_m7 = 0;
  clear_tspec(ts_main_ion, NUM_AA*2); ts_main_ion = 0;
  
  if(f_ms2.is_open())
    f_ms2.close();
}

int SpecFeaturesGenerator :: open_ms2_file_for_reading(string &ms2_filename)
{
  
  //f_ms2 = fopen(ms2_filename.c_str(),"r");
  
  f_ms2.open(ms2_filename.c_str());
  if(!f_ms2.is_open())
    {
      cout << "could not open file " << ms2_filename << " for reading" << endl;
      return 0;
    }
  
  return 1;
}

void SpecFeaturesGenerator :: save_spec_positions(string &out_dir)
{
  ostringstream fname;
  fname << out_dir << "/spec_to_pos_in_file.txt";
  ofstream f_spec_to_pos_in_file(fname.str().c_str());
  f_spec_to_pos_in_file << max_mz << endl;
  for(map<string,streamoff> :: iterator it = spec_to_pos_in_file.begin(); it != spec_to_pos_in_file.end(); it++)
    f_spec_to_pos_in_file << it->first << " " << it->second << endl;
  f_spec_to_pos_in_file.close();
}


void SpecFeaturesGenerator :: save_retention_times(string &out_dir)
{
  ostringstream fname;
  fname << out_dir << "/scan_to_rtime.txt";
  ofstream f_scan_to_rtimes(fname.str().c_str());
  for(map<int,double> :: iterator it = scan_to_rtime.begin(); it != scan_to_rtime.end(); it++)
    f_scan_to_rtimes << it->first << " " << it->second << endl;
  f_scan_to_rtimes.close();
}


void SpecFeaturesGenerator :: load_spec_positions(string &in_dir)
{
  ostringstream fname;
  fname << in_dir << "/spec_to_pos_in_file.txt";
  cout << fname.str() << endl;
  ifstream f_spec_to_pos_in_file(fname.str().c_str());
  f_spec_to_pos_in_file >> max_mz;
  string spec;
  unsigned long line_num;
  while(!f_spec_to_pos_in_file.eof())
    {
      f_spec_to_pos_in_file >> spec;
      f_spec_to_pos_in_file >> line_num;
      spec_to_pos_in_file[spec] = line_num;
    }
  f_spec_to_pos_in_file.close();
}

void SpecFeaturesGenerator :: read_processed_spectrum(string &tempstr)
{
  //read the info on the S line
  pos_in_file = f_ms2.tellg();
  f_ms2 >> first_scan;
  f_ms2 >> last_scan;
  f_ms2 >> precursor_mz;
  
  f_ms2 >> tempstr;
  while(!f_ms2.eof())
    {
      if((tempstr.compare("D") == 0) || (tempstr.compare("I") == 0))
	getline(f_ms2,tempstr);
      else if (tempstr.compare("S") == 0)
	break;
      else if (tempstr.compare("Z") == 0)
	{
	  f_ms2 >> charge;
	  getline(f_ms2,tempstr);
	}
      else
	{
	  istringstream mz(tempstr);
	  int x;
	  mz >> x;
	  f_ms2 >> tempstr;
	  istringstream intens(tempstr);
	  double y;
	  intens >> y;
	  if(x >= max_mz)
	    {
	      max_mz = x+1;
	      peaks.resize(max_mz,0.0);
	    }
	  peaks[x] = y;
	}
      f_ms2 >> tempstr;
      
    }
}

void SpecFeaturesGenerator :: read_processed_ms2_file()
{
  string line;
  string tempstr;
  int num_spec_read = -1;
  f_ms2 >> tempstr;
  while(!f_ms2.eof())
    {
      if(tempstr.compare("H") == 0)
	{
	  getline(f_ms2,line);
	  f_ms2 >> tempstr;
	}
      if (tempstr.compare("S") == 0)
	{
	  num_spec_read++;
	  //if((num_spec_read%1000) == 0)
	  //cout << num_spec_read << endl;
	  if(num_spec_read > 0)
	    {
	      ostringstream spec;
	      spec << first_scan << "." << charge;
	      spec_to_pos_in_file[spec.str()] = pos_in_file;
	    }
	  read_processed_spectrum(tempstr);
	}
    }
  //get the last spectrum 
  if(num_spec_read > 0)
    {
      ostringstream spec;
      spec << first_scan << "." << charge;
      spec_to_pos_in_file[spec.str()] = pos_in_file;
    } 
  
  //cout << "num_spec_read " << num_spec_read << endl;
  //cout << max_mz << "\n";
}



void SpecFeaturesGenerator :: read_spectrum(string &tempstr)
{
  //read the info on the S line
  pos_in_file = f_ms2.tellg();
  f_ms2 >> first_scan;
  f_ms2 >> last_scan;
  f_ms2 >> precursor_mz;
  
  f_ms2 >> tempstr;
  all_charges_of_spec.erase(all_charges_of_spec.begin(),all_charges_of_spec.end());
  mz_values.erase(mz_values.begin(),mz_values.end());
  intens_values.erase(intens_values.begin(),intens_values.end());
  while(!f_ms2.eof())
    {
      if(tempstr.compare("D") == 0)
	getline(f_ms2,tempstr);
      else if(tempstr.compare("I") == 0)
	{
	  f_ms2 >> tempstr;
	  if(tempstr.compare("RTime") == 0)
	    {
	      f_ms2 >> rtime;
	      //cout << rtime << endl;
	    }
	  getline(f_ms2,tempstr);
	}
      else if (tempstr.compare("S") == 0)
	break;
      else if (tempstr.compare("Z") == 0)
	{
	  
	  f_ms2 >> charge;
	  //cout << "scan " << first_scan << " " << charge << " " << endl;
	  all_charges_of_spec.push_back(charge);
	  getline(f_ms2,tempstr);
	}
      else
	{
	  istringstream mz(tempstr);
	  double x;
	  mz >> x;
	  f_ms2 >> tempstr;
	  istringstream intens(tempstr);
	  double y;
	  intens >> y;
	  mz_values.push_back(x);
	  intens_values.push_back(y);
	  //update max_mz
	  int mz_bin = (int)(x/bin_width_mono+0.5);
	  if(mz_bin >= max_mz)
	    max_mz = mz_bin+1;
	}
      f_ms2 >> tempstr;
    }
}


void SpecFeaturesGenerator :: shift_peaks()
{
#if 1
  int n = peaks.size();
  vector<double> sums(n);
  vector<double> integral(n+1);
  integral[0] = 0;
  for (int i=0; i<n; i++)
    integral[i+1] = integral[i] + peaks[i];
  for (int i=0; i<n; i++)
    {
      int ilo = i - max_xcorr_offset;
      int ihi = i + max_xcorr_offset+1;
      sums[i] = integral[ihi<=n ? ihi : n] - integral[ilo>=0 ? ilo : 0]; 
    }
#else
  vector<double> sums;
  sums.resize(peaks.size(),0.0);
  
  for(unsigned int idx = 0; idx < peaks.size(); idx++)
    {
      for(int sub_idx = idx-max_xcorr_offset; sub_idx <= idx+max_xcorr_offset; sub_idx++)
	if(sub_idx > -1 && sub_idx < peaks.size())
	  sums[idx] += peaks[sub_idx];
    }
#endif
  for(unsigned int idx = 0; idx < peaks.size(); idx++)
    peaks[idx] -= (sums[idx]/(max_xcorr_offset*2.0+1));

}

void SpecFeaturesGenerator :: normalize_each_region(double max_intensity_overall, vector<double> &max_intensity_per_region,
						    int region_selector)
{
  int region_idx = 0;
  double max_intensity = max_intensity_per_region[region_idx];
  for (int i = 0; i < (int)peaks.size(); i++)
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
	&& (peaks[i] > 0.05 * max_intensity_overall))
	{
	  // normalize intensity to max 50
	  peaks[i] = (peaks[i] /max_intensity) * max_per_region;
	}
      // no more peaks beyond the 10 regions mark, exit
      if(i > num_regions * region_selector){
	return;
      }
    }
}


void SpecFeaturesGenerator :: process_observed_spectrum()
{
  double experimental_mass_cut_off = precursor_mz*charge+50.0;

  /*
  if(experimental_mass_cut_off > 512){
    int x = (int)experimental_mass_cut_off / 1024;
    double y = experimental_mass_cut_off - (1024 * x);
    max_mz = x * 1024;

    if(y > 0){
      max_mz += 1024;
    }
  }
  */

  if((int)peaks.size() != max_mz)
    peaks.resize(max_mz);
  peaks.assign(peaks.size(),0.0);


  assert(mz_values.size() == intens_values.size());
  double max_peak_location = 0.0;
  double max_peak_intensity = 0.0;
  for(unsigned int i = 0; i < mz_values.size(); i++)
    {
      double peak_location = mz_values[i];
      if(peak_location < experimental_mass_cut_off)
	{
	  if(peak_location > max_peak_location)
	    max_peak_location = peak_location;
	  
	}
    }
  int region_selector = (int) max_peak_location/num_regions;
  vector<double> max_peak_intensity_per_region;
  max_peak_intensity_per_region.resize(num_regions,0.0);
  
  for(unsigned int i = 0; i < mz_values.size(); i++)
    {
      double peak_location = mz_values[i];
      double peak_intensity = intens_values[i];
      //if above experimental mass, skip
       if(peak_location > experimental_mass_cut_off)
	 continue;
       // skip all peaks within precursor ion mz +/- 15
       if(peak_location < precursor_mz + 15 &&  peak_location > precursor_mz - 15)
	 continue;
       
       //get the bin and the region
       int mz = (int)(peak_location/bin_width_mono+0.5);
              	 
       int region = mz/region_selector;
       // don't let index beyond array
       if(region >= num_regions){
	 continue;
       }

       //update max mz over all spectra
       if(mz >= max_mz)
	 {
	   max_mz = mz+1;
	   peaks.resize(max_mz,0.0);
	 }
       // sqrt the original intensity
       peak_intensity = sqrt(peak_intensity);
       if(peak_intensity > max_peak_intensity)
	 max_peak_intensity = peak_intensity;
       if(peaks[mz] < peak_intensity)
	 {
	   peaks[mz] = peak_intensity;
	   // check if this peak is max intensity in the region(one out of 10)
	   if(max_peak_intensity_per_region[region] < peak_intensity){
	     max_peak_intensity_per_region[region] = peak_intensity;
	   }
	 }
    }
  normalize_each_region(max_peak_intensity, max_peak_intensity_per_region,region_selector);
  shift_peaks();
}



void SpecFeaturesGenerator :: read_ms2_file()
{
  string line;
  string tempstr;
  mz_values.reserve(1000);
  intens_values.reserve(1000);

  int num_spec_read = -1;
  f_ms2 >> tempstr;
  while(!f_ms2.eof())
    {
      if(tempstr.compare("H") == 0)
	{
	  getline(f_ms2,line);
	  f_ms2 >> tempstr;
	}
      if (tempstr.compare("S") == 0)
	{
	  num_spec_read++;
	  //if((num_spec_read%5000) == 0)
	  //cout << "spectrum number " << num_spec_read << " " << pos_in_file  << endl;
	  if(num_spec_read > 0)
	    {
	      scan_to_rtime[first_scan] = rtime;
	      ostringstream spec;
	      for(unsigned int k = 0; k < all_charges_of_spec.size();k++)
		{
		  int ch = all_charges_of_spec[k];
		  spec << first_scan << "." << ch;
		  spec_to_pos_in_file[spec.str()] = pos_in_file;
		  spec.str("");
		}
	    }
	  read_spectrum(tempstr);
	}
    }

  //get the last spectrum 
  if(num_spec_read > 0)
    {
      scan_to_rtime[first_scan] = rtime;
      ostringstream spec;
      for(unsigned int k = 0; k < all_charges_of_spec.size();k++)
	{
	  int ch = all_charges_of_spec[k];
	  spec << first_scan << "." << ch;
	  spec_to_pos_in_file[spec.str()] = pos_in_file;
	  spec.str("");
	}
      all_charges_of_spec.clear();
    }
      
  //cout << "num_spec_read " << num_spec_read << endl;
  //cout << max_mz << "\n";
}



/*************************************************************************************/

void SpecFeaturesGenerator :: get_processed_observed_spectrum(string &spec)
{
  if(f_ms2.fail())
    f_ms2.clear();
  int pos_in_file = spec_to_pos_in_file[spec];
  //cout << pos_in_file << endl;
  f_ms2.seekg(pos_in_file,ios::beg);
  
  if((int)peaks.size() < max_mz)
    peaks.resize(max_mz);
  peaks.assign(peaks.size(),0.0);
  string tempstr;
  read_processed_spectrum(tempstr);
  ostringstream spec_read;
  spec_read << first_scan << "." << charge;
  assert(spec.compare(spec_read.str()) == 0);
}


void SpecFeaturesGenerator :: get_observed_spectrum(string &spec)
{
  if(f_ms2.fail())
    f_ms2.clear();
  pos_in_file = spec_to_pos_in_file[spec];
  f_ms2.seekg(pos_in_file,ios::beg);
  pos_in_file = f_ms2.tellg();
    
  string tempstr;
  read_spectrum(tempstr);
  ostringstream spec_read;
  int flag = 0;
  for(unsigned int i = 0; i < all_charges_of_spec.size();i++)
    {
      int ch = all_charges_of_spec[i];
      spec_read << first_scan << "." << ch;
      if(spec.compare(spec_read.str()) == 0)
	{
	  flag = 1;
	  charge = ch;
	}
      spec_read.str("");
    }
  assert(flag==1);
  process_observed_spectrum();
}


void SpecFeaturesGenerator :: clear_tspec(double **tspec,int num_features)
{
  if(tspec)
    {
      for (int i = 0; i < num_features;i++)
	if(tspec[i])
	  delete [] tspec[i];
      delete [] tspec;
    }
}

void SpecFeaturesGenerator :: allocate_tspec(double ***tspec, int num_features, int max_mz)
{
  (*tspec) = new double*[num_features];
  for(int i = 0; i < num_features;i++)
    (*tspec)[i] = new double[max_mz];
}

void SpecFeaturesGenerator :: zero_out_tspec(double **tspec, int num_features, int max_mz)
{
  for(int i = 0; i < num_features;i++)
    memset(tspec[i],0,sizeof(double)*max_mz);
}

void SpecFeaturesGenerator :: add_intensity(double *tspec, int idx, double intensity)
{
  if(idx > -1 && idx < max_mz)
    {
      if(tspec[idx] < intensity)
	tspec[idx] = intensity;
    }
}


double SpecFeaturesGenerator :: sproduct(double *tspec, vector<double> &ospec)
{
  double sm = 0;
  for(unsigned int i = 0; i < ospec.size();i++)
    sm+=tspec[i]*ospec[i];
  return sm;
}


void SpecFeaturesGenerator :: get_spec_features_m3(int scan, int ch , string &peptide, double *features)
{         
  ostringstream scan_stream;
  scan_stream << scan << "." << ch;
  string spec = scan_stream.str();

  get_observed_spectrum(spec);
  
  //allocate the theoretical spectrum if necessary
  if(!ts_m3)
    allocate_tspec(&ts_m3,3,max_mz);
  zero_out_tspec(ts_m3,3,max_mz);
  
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
	      add_intensity(ts_m3[0], idx, 1.0);
	      
	      //add its flanking peaks
	      add_intensity(ts_m3[1], idx-1, 1.0);
	      add_intensity(ts_m3[1], idx+1, 1.0);
	      
	      //add its neutral losses nh3 and h2o
	      double mz_nh3 = mz-(mass_nh3_mono/charg);
	      idx = (int)(mz_nh3/bin_width_mono+0.5);
	      add_intensity(ts_m3[2], idx, 1.0);
	      double mz_h2o = mz-(mass_h2o_mono/charg);
	      idx = (int)(mz_h2o/bin_width_mono+0.5);
	      add_intensity(ts_m3[2], idx, 1.0);

	      //add neutral loss co
	      double mz_co = mz-(mass_co_mono/charg);
	      idx = (int)(mz_co/bin_width_mono+0.5);
	      add_intensity(ts_m3[2], idx, 1.0);
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
		add_intensity(ts_m3[0], idx, 1.0);
		
		//add its flanking peaks
		add_intensity(ts_m3[1], idx-1, 1.0);
		add_intensity(ts_m3[1], idx+1, 1.0);
		
		//add its neutral losses
		double mz_nh3 = mz-(mass_nh3_mono/charg);
		idx = (int)(mz_nh3/bin_width_mono+0.5);
		add_intensity(ts_m3[2], idx, 1.0);
		//double mz_h2o = mz-(mass_h2o_mono/charg);
		//idx = (int)(mz_h2o/bin_width_mono+0.5);
		//add_intensity(ts_m3[2], idx, 1.0);
	      }
	  }
      }

  features[0] = sproduct(ts_m3[0],peaks);
  features[1] = sproduct(ts_m3[1],peaks);
  features[2] = sproduct(ts_m3[2],peaks);

  //cout << features[0] << " " << features[1] << " " << features[2] << endl;

}



void SpecFeaturesGenerator :: get_spec_features_m6(int scan, int ch, string &peptide, double *features)
{
  ostringstream scan_stream;
  scan_stream << scan << "." << ch;
  string spec = scan_stream.str();
  
  get_observed_spectrum(spec);
  
  //allocate the theoretical spectrum if necessary
  if(!ts_m6)
    allocate_tspec(&ts_m6,6,max_mz);
  zero_out_tspec(ts_m6,6,max_mz);

  /*
   * 0. b-ion
   * 1. y_ion
   * 2. flanking
   * 3. h2o_nl
   * 4. co2_nl
   * 5. nh3_nl
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
	      add_intensity(ts_m6[0], idx, 1.0);
	      
	      //add its flanking peaks (2 in array)
	      add_intensity(ts_m6[2], idx-1, 1.0);
	      add_intensity(ts_m6[2], idx+1, 1.0);
	      
	      //add its neutral losses nh3(5 in array) 
	      double mz_nh3 = mz-(mass_nh3_mono/charg);
	      idx = (int)(mz_nh3/bin_width_mono+0.5);
	      add_intensity(ts_m6[5], idx, 1.0);

	      //add h2o (3 in array)
	      double mz_h2o = mz-(mass_h2o_mono/charg);
	      idx = (int)(mz_h2o/bin_width_mono+0.5);
	      add_intensity(ts_m6[3], idx, 1.0);

	      //add neutral loss co (4 in array)
	      double mz_co = mz-(mass_co_mono/charg);
	      idx = (int)(mz_co/bin_width_mono+0.5);
	      add_intensity(ts_m6[4], idx, 1.0);
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
		add_intensity(ts_m6[1], idx, 1.0);
		
		//add its flanking peaks ( 2 in array)
		add_intensity(ts_m6[2], idx-1, 1.0);
		add_intensity(ts_m6[2], idx+1, 1.0);
		
		//add its neutral losses (5 in array)
		double mz_nh3 = mz-(mass_nh3_mono/charg);
		idx = (int)(mz_nh3/bin_width_mono+0.5);
		add_intensity(ts_m6[5], idx, 1.0);
	      }
	  }
      }

  features[0] = sproduct(ts_m6[0],peaks);
  features[1] = sproduct(ts_m6[1],peaks);
  features[2] = sproduct(ts_m6[2],peaks);
  features[3] = sproduct(ts_m6[3],peaks);
  features[4] = sproduct(ts_m6[4],peaks);
  features[5] = sproduct(ts_m6[5],peaks);

  //cout << features[0] << " " << features[1] << " " << features[2] << " " << features[3] << " " << features[4] << " " << features[5] << endl;

}



void SpecFeaturesGenerator :: get_spec_features_m7(int scan, int ch, string &peptide, double *features)
{
  
  ostringstream scan_stream;
  scan_stream << scan << "." << ch;
  string spec = scan_stream.str();

  get_observed_spectrum(spec);
  
  //allocate the theoretical spectrum if necessary
  if(!ts_m7)
    allocate_tspec(&ts_m7,7,max_mz);
  zero_out_tspec(ts_m7,7,max_mz);

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
	      add_intensity(ts_m7[0], idx, 1.0);
	      
	      //add its flanking peaks (2 in array)
	      add_intensity(ts_m7[2], idx-1, 1.0);
	      add_intensity(ts_m7[2], idx+1, 1.0);
	      
	      //add its neutral losses nh3(5 in array) 
	      double mz_nh3 = mz-(mass_nh3_mono/charg);
	      idx = (int)(mz_nh3/bin_width_mono+0.5);
	      add_intensity(ts_m7[5], idx, 1.0);

	      //add h2o (3 in array)
	      double mz_h2o = mz-(mass_h2o_mono/charg);
	      idx = (int)(mz_h2o/bin_width_mono+0.5);
	      add_intensity(ts_m7[3], idx, 1.0);

	      //add neutral loss co (4 in array)
	      double mz_co = mz-(mass_co_mono/charg);
	      idx = (int)(mz_co/bin_width_mono+0.5);
	      add_intensity(ts_m7[4], idx, 1.0);
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
		add_intensity(ts_m7[1], idx, 1.0);
		
		//add its flanking peaks ( 2 in array)
		add_intensity(ts_m7[2], idx-1, 1.0);
		add_intensity(ts_m7[2], idx+1, 1.0);
		
		//add its neutral losses (6 in array)
		double mz_nh3 = mz-(mass_nh3_mono/charg);
		idx = (int)(mz_nh3/bin_width_mono+0.5);
		add_intensity(ts_m7[6], idx, 1.0);
	      }
	  }
      }

  features[0] = sproduct(ts_m7[0],peaks);
  features[1] = sproduct(ts_m7[1],peaks);
  features[2] = sproduct(ts_m7[2],peaks);
  features[3] = sproduct(ts_m7[3],peaks);
  features[4] = sproduct(ts_m7[4],peaks);
  features[5] = sproduct(ts_m7[5],peaks);
  features[6] = sproduct(ts_m7[6],peaks);

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


void SpecFeaturesGenerator :: get_spec_features_aa_end(int scan, int ch, string &peptide, double *features)
{
  
  ostringstream scan_stream;
  scan_stream << scan << "." << ch;
  string spec = scan_stream.str();
  
  get_observed_spectrum(spec);
  
  //allocate the theoretical spectrum if necessary
  if(!ts_main_ion)
    allocate_tspec(&ts_main_ion,NUM_AA*2,max_mz);
  zero_out_tspec(ts_main_ion,NUM_AA*2,max_mz);
  
  double mz = 0;
  int idx = 0;
  int aa_ind = 0;
  int offset;
  double fragment_mass;
  
  //start from the beginning of the peptide and do the B-ION
  offset = 0;
  fragment_mass = 0.0;
  for(unsigned int i = 0; i < peptide.size()-1;i++)
    {
      char aa = peptide.at(i);
      aa_ind = aa-'A';
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
	      add_intensity(ts_main_ion[offset+aa_ind], idx, 1.0);
	    }
	}
    }

  //do the Y-ION
  offset = NUM_AA;
  fragment_mass = mass_h2o_mono;
  for(unsigned int i = peptide.size(); i > 1;i--)
    {
      char aa = peptide.at(i-1);
      aa_ind = aa-'A';
      if(aa != '.')
	  {
	    fragment_mass += aa_masses_mono[aa-'A'];
	    //add the fragment
	    for (int charg = 1; charg < ch; charg++)
	      {
		//add the ion itself
		mz = (fragment_mass+charg*proton_mass)/charg;
		idx = (int)(mz/bin_width_mono+0.5);
		add_intensity(ts_main_ion[offset+aa_ind], idx, 1.0);
	      }
	  }
      }

  for(int i = 0; i < 2*NUM_AA;i++)
    features[i] = sproduct(ts_main_ion[i],peaks);

}


void SpecFeaturesGenerator :: get_spec_features_aa_mid(int scan, int ch, string &peptide, double *features)
{
  ostringstream scan_stream;
  scan_stream << scan << "." << ch;
  string spec = scan_stream.str();
  
  get_observed_spectrum(spec);
  
  //allocate the theoretical spectrum if necessary
  if(!ts_main_ion)
    allocate_tspec(&ts_main_ion,NUM_AA*2,max_mz);
  zero_out_tspec(ts_main_ion,NUM_AA*2,max_mz);
  
  double mz = 0;
  int idx = 0;
  int aa_ind = 0;
  int offset;
  double fragment_mass;
  
  //start from the beginning of the peptide and do the B-ION
  offset = 0;
  fragment_mass = 0.0;
  for(unsigned int i = 0; i < peptide.size()-1;i++)
    {
      char aa = peptide.at(i);
      aa_ind = aa-'A';
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
	      add_intensity(ts_main_ion[offset+aa_ind], idx, 1.0);
	      for (unsigned int j = 0; j < i ; j++)
		{
		  char bb = peptide.at(j);
		  aa_ind = bb-'A';
		  add_intensity(ts_main_ion[offset+aa_ind], idx, 1.0);
		}
	    }
	}
    }

  //do the Y-ION
  offset = NUM_AA;
  fragment_mass = mass_h2o_mono;
  for(unsigned int i = peptide.size(); i > 1;i--)
    {
      char aa = peptide.at(i-1);
      aa_ind = aa-'A';
      if(aa != '.')
	  {
	    fragment_mass += aa_masses_mono[aa-'A'];
	    //add the fragment
	    for (int charg = 1; charg < ch; charg++)
	      {
		//add the ion itself
		mz = (fragment_mass+charg*proton_mass)/charg;
		idx = (int)(mz/bin_width_mono+0.5);
		add_intensity(ts_main_ion[offset+aa_ind], idx, 1.0);
		for (unsigned int j = peptide.size(); j >i ; j--)
		  {
		    char bb = peptide.at(j-1);
		    aa_ind = bb-'A';
		    add_intensity(ts_main_ion[offset+aa_ind], idx, 1.0);
		  }
	      }
	  }
      }

  for(int i = 0; i < 2*NUM_AA;i++)
    features[i] = sproduct(ts_main_ion[i],peaks);

}






//int main()
//{
  /*
  SpecFeaturesGenerator sp;
  string ms2_fn = "yeast/yeast-01_processed.ms2";
  sp.open_ms2_file_for_reading(ms2_fn);

  //sp.read_processed_ms2_file();
  //string out_dir = "yeast";
  //sp.save_spec_positions(out_dir);
  //sp.clear();
  string in_dir = "yeast";
  sp.load_spec_positions(in_dir);
  
  
  sp.initialize_aa_tables();
  double *feat = new double[3];
  string peptide = "ELPIVTR";
  string spec = "591.2";
  //sp.get_spec_features_aa_end(spec, peptide, feat);
  sp.get_spec_features_m7(591, 2, peptide, feat);

  //string peptide = "R.HEGEHGLYFR.S";
  peptide = "HEGEHGLYFR";
  spec = "591.3";
  sp.get_spec_features_m7(591,3, peptide, feat);
  delete [] feat;
  */

/*
  SpecFeaturesGenerator sp;
  string ms2_fn = "yeast/yeast-01.ms2";
  sp.open_ms2_file_for_reading(ms2_fn);
  //sp.read_ms2_file();
  //string out_dir = "yeast";
  //sp.save_spec_positions(out_dir);
  //sp.clear();
  string in_dir = "yeast";
  sp.load_spec_positions(in_dir);
  
  
  sp.initialize_aa_tables();
  double *feat = new double[7];
  
  string peptide = "ELPIVTR";
  string spec = "591.2";
  sp.get_spec_features_m7(591,2, peptide, feat);


  //string peptide = "R.HEGEHGLYFR.S";
  //peptide = "HEGEHGLYFR";
  //spec = "591.3";
  //sp.get_spec_features_m7(591,3, peptide, feat);

  delete [] feat;



}
*/

