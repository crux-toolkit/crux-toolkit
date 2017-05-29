#include "xhhc_score_peptide_spectrum.h"
#include "xhhc.h"
#include "LinkedIonSeries.h"
#include "xhhc_scorer.h"
#include "LinkedPeptide.h"
#include "XHHC_Peptide.h"

#include "XLinkPeptide.h"
#include "SelfLoopPeptide.h"
#include "XLinkScorer.h"

#include "model/IonConstraint.h"
#include "model/Scorer.h"
#include "io/SpectrumCollectionFactory.h"
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

XLinkScoreSpectrum::XLinkScoreSpectrum() {
}

XLinkScoreSpectrum::~XLinkScoreSpectrum() {
}

int XLinkScoreSpectrum::main(int argc, char** argv) {
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

  FLOAT_T link_mass = Params::GetDouble("link mass");
  XLinkPeptide::setLinkerMass(link_mass);
  string scoremethod(Params::GetString("xlink-score-method"));


  XLinkMatch* candidate = NULL;
  
  // a single peptide linked to itself
  if (strcmp(peptideB, "NULL") == 0) {
    cout << "B is null" << endl;
    if (scoremethod != "composite") {
      carp(CARP_FATAL, "For composite scoring, provide a linked peptide, "
           "not a self loop");
    }
    candidate = new SelfLoopPeptide(peptideA, posA, posB);
  } else {
    candidate = new XLinkPeptide(peptideA, peptideB, posA, posB);
  }

  // read ms2 file
  Crux::SpectrumCollection* collection = SpectrumCollectionFactory::create(ms2_file);
  
  // search for spectrum with correct scan number
  Spectrum* spectrum = collection->getSpectrum(scan_num);
  if( spectrum == NULL ) {
    carp(CARP_ERROR, "Failed to find spectrum with scan_num: %d", scan_num);
    delete collection;
    exit(1);
  }
  
  XLinkScorer xlink_scorer(spectrum, charge);
  
  IonConstraint* ion_constraint = xlink_scorer.getIonConstraintXCorr();
  ion_constraint->setUseIonType(BY_ION, false);
  ion_constraint->setUseIonType(BYA_ION, false);
  ion_constraint->setUseIonType(ALL_ION, false);
  ion_constraint->setUseIonType(A_ION, Params::GetBool("use-a-ions"));
  ion_constraint->setUseIonType(B_ION, Params::GetBool("use-b-ions"));
  ion_constraint->setUseIonType(C_ION, Params::GetBool("use-c-ions"));
  ion_constraint->setUseIonType(X_ION, Params::GetBool("use-x-ions"));
  ion_constraint->setUseIonType(Y_ION, Params::GetBool("use-y-ions"));
  ion_constraint->setUseIonType(Z_ION, Params::GetBool("use-z-ions"));

  if (scoremethod == "composite") {
    FLOAT_T score = xlink_scorer.scoreCandidate(candidate);   

    cout << score << endl;

  } else if (scoremethod == "modification") {
    
    XLinkPeptide* xlp = (XLinkPeptide*)candidate;
    FLOAT_T deltaB = xlp->getXLinkablePeptide(1).getMass(MONO) + link_mass;
    FLOAT_T scoreA = xlink_scorer.scoreXLinkablePeptide(xlp->getXLinkablePeptide(0), 0, deltaB);
    
    FLOAT_T deltaA = xlp->getXLinkablePeptide(0).getMass(MONO) + link_mass;
    FLOAT_T scoreB = xlink_scorer.scoreXLinkablePeptide(xlp->getXLinkablePeptide(1), 0, deltaA);

    if (scoreA > scoreB)
      cout << scoreA << "\t" << scoreB << endl;
    else
      cout << scoreB << "\t" << scoreA << endl;

  } else if (scoremethod == "concatenation") {

    vector<double> scores;
    double score1 = get_concat_score(peptideA, peptideB, posA, charge, spectrum);
    scores.push_back(score1);

    double score2 = get_concat_score(peptideB, peptideA, posB, charge, spectrum);
    scores.push_back(score2);

    int lengthA = string(peptideA).length();
    int lengthB = string(peptideB).length();

    double score3 = get_concat_score(peptideA, peptideB, lengthA + posB, charge, spectrum);
    scores.push_back(score3);

    double score4 = get_concat_score(peptideB, peptideA, lengthB + posA, charge, spectrum);
    scores.push_back(score4);

    sort(scores.begin(), scores.end(), less<double>());
    cout <<scores[0];
    for (int i=1; i < 4; i++) {
      cout << "\t" << scores[i];
    }

    cout << endl;
  } else {
    carp(CARP_ERROR, "Unknown score method (%s).", scoremethod.c_str());
  }
  // free heap
  delete collection;
  delete spectrum;

  return 0;
}


double XLinkScoreSpectrum::get_concat_score(char* peptideA, char* peptideB, int link_site, int charge, Spectrum* spectrum) {
  string lpeptide = string(peptideA) + string(peptideB); 
  
  IonConstraint* ion_constraint = IonConstraint::newIonConstraintSmart(XCORR, charge);
  
  IonSeries* ion_series = new IonSeries(lpeptide.c_str(), charge, ion_constraint);
  
  ion_series->predictIons();
  
  //modify ions.
  
  //int pepA_begin = 0;
  int pepB_begin = string(peptideA).length();
  int llength = lpeptide.length();
  
  for (IonIterator ion_iterator = ion_series->begin();
    ion_iterator != ion_series->end();
    ++ion_iterator) {

    Ion* ion = *ion_iterator;
    //check to see if if is the cterm of 1st peptide.
    int ion_charge = ion->getCharge();
    int cleavage_idx = ion->getCleavageIdx();
    ION_TYPE_T ion_type = ion->getType();

    //if contains cterm of 1st peptide, modify by -OH 
    carp(CARP_DEBUG, "====================");
    if (ion_type == B_ION) {
      carp(CARP_DEBUG, "B-ion");
      carp(CARP_DEBUG, "%s", lpeptide.substr(0, cleavage_idx).c_str());
    } else if (ion_type == Y_ION) {
      carp(CARP_DEBUG, "Y-ion");
      carp(CARP_DEBUG, "%s", lpeptide.substr(llength-cleavage_idx, llength).c_str());
    } else {
      continue;
    }

    carp(CARP_DEBUG, "cleavage idx:%d", cleavage_idx);

    bool cterm_1st = false;
    if (ion_type == B_ION) {
      carp(CARP_DEBUG, "B-Ion");
      if (cleavage_idx >= pepB_begin) {
        cterm_1st = true;
      }
    } else if (ion_type == Y_ION) {
      carp(CARP_DEBUG, "Y-ion");
      if (cleavage_idx > (llength - pepB_begin)) {
        cterm_1st = true;
      }
    }

    bool nterm_2nd = false;
    if (ion_type == B_ION) {
      if (cleavage_idx > pepB_begin) {
        nterm_2nd = true;
      }
    } else if (ion_type == Y_ION) {
      if (cleavage_idx >= (llength-pepB_begin)) {
        nterm_2nd = true;
      }
    }

    bool has_link_site = false;
    if (ion_type == B_ION) {
      if (cleavage_idx > link_site) {
        has_link_site = true;
      }
    } else if (ion_type == Y_ION) {
      if (cleavage_idx >= (llength- link_site)) {
        has_link_site = true;
      }
    }
      
    carp(CARP_DEBUG, "cterm:%d", cterm_1st);
    carp(CARP_DEBUG, "nterm:%d", nterm_2nd);
    carp(CARP_DEBUG, "has site:%d", has_link_site);
      

    //if it contains the cterm of the 1st peptide, modify by -OH
    if (cterm_1st) {
      FLOAT_T old_mass = (ion->getMassZ() - MASS_H_MONO) * (FLOAT_T)ion_charge;
      FLOAT_T new_mass = old_mass + MASS_H2O_MONO - MASS_H_MONO;
      FLOAT_T new_mz = (new_mass + (FLOAT_T)ion_charge) / (FLOAT_T)ion_charge;
      ion->setMassZ(new_mz);
    }
    //if contains the nterm of 2nd peptide, modify by -H
    if (nterm_2nd) {
      FLOAT_T old_mass = (ion->getMassZ() - MASS_H_MONO) * (FLOAT_T)ion_charge;
      FLOAT_T new_mass = old_mass + MASS_H_MONO;
      FLOAT_T new_mz = (new_mass + (FLOAT_T)ion_charge) / (FLOAT_T)ion_charge;
      ion->setMassZ(new_mz);
    }
    //if contains the link site, modify by link mass.
    if (has_link_site) {
      FLOAT_T old_mass = (ion->getMassZ() - MASS_H_MONO) * (FLOAT_T)ion_charge;
      FLOAT_T new_mass = old_mass + LinkedPeptide::getLinkerMass();
      FLOAT_T new_mz = (new_mass + (FLOAT_T)ion_charge) / (FLOAT_T)ion_charge;
      ion->setMassZ(new_mz);
    }
  }

  Scorer* scorer = new Scorer(XCORR); 

  // calculate the score
  FLOAT_T score = scorer->scoreSpectrumVIonSeries(spectrum, ion_series);
  return score;
}

FLOAT_T* XLinkScoreSpectrum::get_observed_raw(Spectrum* spectrum, int charge) {
  FLOAT_T peak_location = 0;
  int mz = 0;
  FLOAT_T intensity = 0;
  FLOAT_T bin_width = BIN_WIDTH_MONO;
  FLOAT_T precursor_mz = spectrum->getPrecursorMz();
  FLOAT_T experimental_mass_cut_off = precursor_mz*charge + 50;

  // set max_mz and malloc space for the observed intensity array
  FLOAT_T sp_max_mz = 512;

  if(experimental_mass_cut_off > 512) {
    int x = (int)experimental_mass_cut_off / 1024;
    FLOAT_T y = experimental_mass_cut_off - (1024 * x);
    sp_max_mz = x * 1024;

    if(y > 0) {
      sp_max_mz += 1024;
    }
  }

  // DEBUG
  // carp(CARP_INFO, "experimental_mass_cut_off: %.2f sp_max_mz: %.3f", experimental_mass_cut_off, scorer->sp_max_mz);
  FLOAT_T* observed = (FLOAT_T*)mycalloc((int)sp_max_mz, sizeof(FLOAT_T));
  
  // DEBUG
  // carp(CARP_INFO, "max_peak_mz: %.2f, region size: %d",get_spectrum_max_peak_mz(spectrum), region_selector);
  
  for (PeakIterator peak_iterator = spectrum->begin();
    peak_iterator != spectrum->end();
    ++peak_iterator) {

    Peak *peak = *peak_iterator;
    peak_location = peak->getLocation();
    
    // skip all peaks larger than experimental mass
    if(peak_location > experimental_mass_cut_off) {
      continue;
    }
    
    // skip all peaks within precursor ion mz +/- 15
    if(peak_location < precursor_mz + 15 &&  peak_location > precursor_mz - 15) {
      continue;
    }
    
    // map peak location to bin
    mz = (int)(peak_location / bin_width + 0.5);

    // get intensity
    // sqrt the original intensity
    intensity = peak->getIntensity();

    // set intensity in array with correct mz, only if max peak in the bin
    if(observed[mz] < intensity) {
      observed[mz] = intensity;
      }
    }    
  

  return observed;

}



void XLinkScoreSpectrum::print_spectrum(Spectrum* spectrum, LinkedIonSeries& ion_series) {


      Scorer* scorer = new Scorer(XCORR);
      scorer->createIntensityArrayObserved(spectrum, ion_series.getCharge());

      FLOAT_T* observed_raw = get_observed_raw(spectrum, ion_series.getCharge());
      FLOAT_T* observed_processed = scorer->getIntensityArrayObserved();


      FLOAT_T max_mz = scorer->getSpMaxMz();

        


      XHHC_Scorer xhhc_scorer(max_mz);

      FLOAT_T* theoretical = (FLOAT_T*)mycalloc((size_t)max_mz, sizeof(FLOAT_T));
      xhhc_scorer.hhcCreateIntensityArrayTheoretical(ion_series, theoretical);


      
      ofstream fout("ion_match.out");
      for (int i = 0; i < max_mz; i++) {
        fout << i << "\t" 
             << observed_raw[i] << "\t"
             << observed_processed[i] << "\t" 
             << theoretical[i] << "\t" 
             << (theoretical[i] != 0) << "\t" 
             << (theoretical[i] > 25) << endl;
      }

      fout.close();
}

void XLinkScoreSpectrum::processParams() {
  set_modspec();
}

string XLinkScoreSpectrum::getName() const {
  return "xlink-score-spectrum";
}

string XLinkScoreSpectrum::getDescription() const {
  return
    "[[nohtml:Given a cross-linked peptide and a spectrum "
    "calculate the corresponding XCorr score a number of different ways.]]"
    "[[html:Given a cross-linked peptide and a spectrum "
    "calculate the corresponding XCorr score a number of different ways, "
    "depending upon the xlink-score-method parameter:"
    "<ul><li>composite &ndash; compute a combined XCorr score</li>"
    "<li>modification - score the two peptides "
    "separately, treating the second peptide as a variable modification "
    "on the first peptide</li><li>concatenated - score as the concatenation "
    "of the two peptides.  Note that this mode gives four scores, "
    "corresponding to the two relative orders of the peptides (AB and BA) "
    "and the modification appearing on the first or second peptide.</li></ul>]]";
}

vector<string> XLinkScoreSpectrum::getArgs() const {
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

vector<string> XLinkScoreSpectrum::getOptions() const {
  string arr[] = {
    "verbosity",
    "use-flanking-peaks",
    "xlink-score-method",
    "use-a-ions",
    "use-b-ions",
    "use-c-ions",
    "use-x-ions",
    "use-y-ions",
    "use-z-ions"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector< pair<string, string> > XLinkScoreSpectrum::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("stdout",
    "XCorr score(s) in descending order"));
  return outputs;
}

