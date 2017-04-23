//TODO - change cerr/couts to carp.

#include "xlink_assign_ions.h"

#include "XLinkMatch.h"
#include "SelfLoopPeptide.h"
#include "XLinkPeptide.h"

#include "model/objects.h"
#include "model/Scorer.h"
#include "model/IonSeries.h"
#include "model/Ion.h"
#include "io/SpectrumCollectionFactory.h"
#include "io/DelimitedFile.h"
#include "util/Params.h"

#include <math.h>
#include <assert.h>
#include <ctype.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include <iostream>
#include <fstream>


using namespace Crux;

XLinkAssignIons::XLinkAssignIons() {
}

XLinkAssignIons::~XLinkAssignIons() {
}

int XLinkAssignIons::main(int argc, char** argv) {
  /* Get Arguments */
  string peptideAStr = Params::GetString("peptide A");
  string peptideBStr = Params::GetString("peptide B");
  char* peptideA = my_copy_string(peptideAStr.c_str());
  char* peptideB = my_copy_string(peptideBStr.c_str());
  
  int posA     = Params::GetInt("pos A");
  int posB     = Params::GetInt("pos B");
  int charge   = Params::GetInt("charge state"); 
  int scan_num = Params::GetInt("scan number"); 

  string ms2_file = Params::GetString("ms2 file");

  XLinkPeptide::setLinkerMass(Params::GetDouble("link mass"));
 
  // create new ion series
  
  // a single peptide linked to itself
  if (strcmp(peptideB, "NULL") == 0) {
    cout << "B is null" << endl; 
    peptideB = NULL;
  }

  // read ms2 file
  Crux::SpectrumCollection* collection = SpectrumCollectionFactory::create(ms2_file);
  collection->parse();

  
  //cout << "lp " << lp << endl; 
  Spectrum* spectrum = NULL;
  Spectrum* current_spectrum = NULL;

  //TODO allow a binary search on both mgf and ms2 files.

  for (SpectrumIterator spectrum_iterator = collection->begin();
    spectrum_iterator != collection->end();
    ++spectrum_iterator) {

    current_spectrum = *spectrum_iterator;
    if (current_spectrum->getFirstScan() == scan_num) {
      spectrum = current_spectrum;
      break;
    }
  }

  // search for spectrum with correct scan number
  if (spectrum == NULL) {
    carp(CARP_ERROR, "Failed to find spectrum with scan number %d.", scan_num);
    delete collection;
    exit(1);
  }

  //created linked peptide.
  XLinkMatch* xlp = NULL;
  if (peptideB != NULL) {
    xlp = new XLinkPeptide(peptideA, peptideB, posA-1, posB-1);
  } else {
    xlp = new SelfLoopPeptide(peptideA, posA-1, posB-1); 
  }
  
  IonConstraint* ion_constraint = IonConstraint::newIonConstraintSequestXcorr(charge);
  ion_constraint->setUseIonType(BY_ION, false);
  ion_constraint->setUseIonType(BYA_ION, false);
  ion_constraint->setUseIonType(ALL_ION, false);
  ion_constraint->setUseIonType(A_ION, Params::GetBool("use-a-ions"));
  ion_constraint->setUseIonType(B_ION, Params::GetBool("use-b-ions"));
  ion_constraint->setUseIonType(C_ION, Params::GetBool("use-c-ions"));
  ion_constraint->setUseIonType(X_ION, Params::GetBool("use-x-ions"));
  ion_constraint->setUseIonType(Y_ION, Params::GetBool("use-y-ions"));
  ion_constraint->setUseIonType(Z_ION, Params::GetBool("use-z-ions"));
  IonSeries* ion_series = new IonSeries(ion_constraint, charge);
  xlp->predictIons(ion_series, charge);
  print_spectrum(spectrum, xlp, ion_series);

  // free heap
  delete collection;
  delete ion_series;
  delete ion_constraint;
  delete xlp;
  //free_spectrum(spectrum);
  return 0;
}

void XLinkAssignIons::print_spectrum(
  Spectrum* spectrum,
  XLinkMatch* xlp,
  IonSeries* ion_series
  ) {
 
  carp(CARP_INFO, "Number of peaks = %d", spectrum->getNumPeaks());

  FLOAT_T bin_width = Params::GetDouble("mz-bin-width");
      
  map<Peak*, Ion*> matched;
  double matched_intensity = 0;
      
  for (IonIterator ion_iterator = ion_series->begin();
       ion_iterator != ion_series->end();
        ++ion_iterator) {
    
    Ion* ion = *ion_iterator;
    Peak* peak = spectrum->getMaxIntensityPeak(ion->getMassZ(), bin_width);  
    if (peak != NULL) {
      if (matched.find(peak) == matched.end()) {
        matched_intensity += peak->getIntensity();
        matched[peak] = ion;
      } else {
        carp(CARP_WARNING, "Peak at %g matches more than one theoretical peak.", peak->getLocation());
      }
    }
  }

  double total_intensity = spectrum->getTotalEnergy();
  double frac_intensity = matched_intensity / total_intensity;
  carp(CARP_INFO, "Fraction matched intensity = %lf / %lf = %lf",
       matched_intensity, total_intensity, frac_intensity);
      

  //now print out the spectrum

  DelimitedFile result_file;

  unsigned int mz_obs_col = result_file.addColumn("m/z obs");
  unsigned int int_col = result_file.addColumn("intensity");
  unsigned int int_matched_col = result_file.addColumn("matched intensity");
  unsigned int mz_calc_col = result_file.addColumn("m/z calc");
  unsigned int mass_obs_col = result_file.addColumn("mass obs");
  unsigned int mass_calc_col = result_file.addColumn("mass calc");
  unsigned int ppm_col = result_file.addColumn("ppm(Mass)");
  unsigned int ion_col = result_file.addColumn("ion");
  unsigned int seq_col = result_file.addColumn("sequence");
      
  for (PeakIterator peak_iter = spectrum->begin();
    peak_iter != spectrum->end();
    ++peak_iter) {

    Peak * peak = *peak_iter;
    if (peak->getIntensity() > 0) { 
      unsigned int row_idx = result_file.addRow();
      result_file.setValue(mz_obs_col, row_idx, peak->getLocation());
      result_file.setValue(int_col, row_idx, peak->getIntensity());
      if (matched.find(peak) != matched.end()) {
        Ion* ion = matched[peak];
        double mz_calc = ion->getMassZ();
        double mz_obs = peak->getLocation();
        int charge = ion->getCharge();
        double mass_calc = ion->getMassFromMassZ();
        double mass_obs = (mz_obs - MASS_PROTON) * (double)charge;
        double ppm =  fabs(mass_calc - mass_obs) / mass_calc * 1e6;
         
        result_file.setValue(int_matched_col, row_idx, peak->getIntensity()); 
        result_file.setValue(mz_calc_col, row_idx, mz_calc);
        result_file.setValue(mass_obs_col, row_idx, mass_obs);
        result_file.setValue(mass_calc_col, row_idx, mass_calc);
        result_file.setValue(ppm_col, row_idx, ppm);
              
        ostringstream ion_string_oss;
        switch(ion->getType()) {
          case A_ION:
            ion_string_oss << "a";
            break;
          case B_ION:
            ion_string_oss << "b";
            break;
          case C_ION:
            ion_string_oss << "c";
            break;
          case X_ION:
            ion_string_oss << "x";
            break;
          case Y_ION:
            ion_string_oss << "y";
            break;
          case Z_ION:
            ion_string_oss << "z";
            break;
          default:
            ion_string_oss << "u";
        }
        ion_string_oss << "(+" << charge << ")";

        result_file.setValue(ion_col, row_idx, ion_string_oss.str());
        result_file.setValue(seq_col, row_idx, xlp->getIonSequence(ion));
      }
    }
  }
  cout << result_file;

}

void XLinkAssignIons::processParams() {
  set_modspec();
}

string XLinkAssignIons::getName() const {
  return "xlink-assign-ions";
}

string XLinkAssignIons::getDescription() const {
  return
    "Given a spectrum and a pair of cross-linked peptides, assign theoretical "
    "ion type labels to peaks in the observed spectrum.";
}

vector<string> XLinkAssignIons::getArgs() const {
  string arr[] = {
    "peptide A",
    "peptide B",
    "pos A",
    "pos B",
    "link mass",
    "charge state",
    "scan number",
    "ms2 file"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector<string> XLinkAssignIons::getOptions() const {
  string arr[] = {
    "verbosity",
    "spectrum-parser",
    "fragment-mass",
    "max-ion-charge",
    "mz-bin-width",
    "precision"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector< pair<string, string> > XLinkAssignIons::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("stdout",
    "tab-delimited text in which each row is a peak in the observed spectrum, "
    "and the columns are <ol>"
    "<li>The m/z value.</li>"
    "<li>The observed intensity.</li>"
    "<li>The matched intensity.</li>"
    "<li>The calculated m/z value of the theoretical peak.</li>"
    "<li>The mass associated with the observed peak.</li>"
    "<li>The mass difference (in ppm) between the observed and theoretical peaks.</li>"
    "<li>The ion type, specified as b or y, followed by the charge state in parentheses.</li>"
    "<li>The amino acid sequence of the fragment.</li>"));
  return outputs;
}

