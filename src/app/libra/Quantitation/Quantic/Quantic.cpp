#include "Quantic.h"
/* ******************************************************************************

Program       : Quantic
Author        : David Shteynberg <dshteynb  AT systemsbiology.org>
Date          : 04.20.2018
SVN Info      : $Id$

Copyright (C) 2018 David Shteynberg

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

David Shteynberg
Institute for Systems Biology
401 Terry Avenue North
Seattle, WA  98109  USA

******************************************************************************** */
#define DEFAULT_PPM 1
#define MAX_CH 99
#define PROTONMASS 1.007276466
// int nTermMod(vector<int>& ntermod) {  
//   for (int t=0; t< ntermod.size(); t++) {
//     if (ntermod[t]) 
//       return 1;
//   }
//   return 0; 
// }


bool compareIntens (const double& a, const double& b)  {
  return a > b;

}
bool compareNeutQuants (const NeutQuant* a, const NeutQuant* b)  {
  return a->loss_ < b->loss_;

}

Quantic::Quantic(string& pep, int charge, double calc_neut_mass, 
		 cRamp* cramp, long scan, 
		  vector<string>& modaas, vector<double>& shift,  vector<vector<double>*>& neutlosses,
			     TPP_HASHMAP_T<char, double>* stat_mods, TPP_HASHMAP_T<char, vector<double>*>* var_mods, 
			     TPP_HASHMAP_T<char, double>* stat_prot_termods, TPP_HASHMAP_T<char, vector<double>*>* var_prot_termods, 
#ifdef MSVC
		HANDLE* mutex,
#else
		pthread_mutex_t* mutex,
#endif
		 bool is_nterm_pep, bool is_cterm_pep, TPP_HASHMAP_T<int, double>& pos_mod_hash ) { 
                  // bool massdiff_mode, bool labile_mode, bool direct_mode) {

  diaMode_ = false;
  
  mutex_ = mutex;
  
  stat_mods_ = stat_mods;

  pep_coupled_ = false;

  ppm_tol_ = DEFAULT_PPM;
  
  mass_diff_ = 0.;

  top_peaks_ = 0;
  
  stat_prot_termods_ = stat_prot_termods;


  var_mods_ = var_mods;

  var_prot_termods_ = var_prot_termods;

  has_nterm_mod_ = false;
  has_cterm_mod_ = false;
  
  is_nterm_pep_ = is_nterm_pep;
  is_cterm_pep_ = is_cterm_pep;
  
  calc_neut_mass_ = calc_neut_mass;

  shift_ = shift;
  neutlosses_ = neutlosses;
  
  
  verbose_ = 0;

  pep_ = new PeptideUser();
  
  //Peptide::defaultTables();

  // if (scan == 11842) {
  ////  cerr << "DDS: DEBUG THIS ... " << endl;
  //}

  bool nterm = false,  cterm = false;
  for (int j = 0; j < shift_.size(); j++) {

    nmods_.push_back(0);

    neutlossquants_.push_back(new vector<NeutQuant*>());
    
    neutquants_.push_back(new vector<double>());
    neutquantmeans_.push_back(new vector<double>());
    neutquantstdvs_.push_back(new vector<double>());
    neutquantvect_.push_back(new vector<vector<double>*>());
    mods_.push_back(new vector<int>());
    //  nomods_.push_back(new vector<vector<int>*>());
    
    label_.push_back( new strp_hash());

     if (modaas[j].find('c') != string::npos) {
      ctermod_.push_back(1);
      cterm = true;
      has_cterm_mod_ = true;
     }
     else if (stat_mods_->find('c') != stat_mods_->end()) {
       ctermod_.push_back(0);
       cterm = true;
       has_cterm_mod_ = true;
     }
     else {
       ctermod_.push_back(0);
     }
     
     if (modaas[j].find('n') != string::npos) {
       ntermod_.push_back(1);
       nterm = true;
       has_nterm_mod_ = true;
     }
     else if (stat_mods_->find('n') != stat_mods_->end() || pep.substr(0,1) == "n") {
       ntermod_.push_back(0);
       nterm = true;
       has_nterm_mod_ = true;
     }
     else {
       ntermod_.push_back(0);
     }

     for (int a = 0; a <  modaas[j].length(); a++) {
       string userToken("");
       string userType("");
       stringstream newModTypess;
     
       newModTypess.precision(6);
       newModTypess << "Qtc_" << modaas[j][a] << "_" << shift[j];

       userType = newModTypess.str();
#ifdef MSVC
	      WaitForSingleObject( *mutex_,    // handle to mutex
				   INFINITE);      //#include "windows.h"
#else	
	      pthread_mutex_lock( &*mutex_ );//#include <pthread.h>
#endif 
	      pep_->processUserMod(modaas[j][a], shift[j], userType, userToken, true);
#ifdef MSVC
		ReleaseMutex( *mutex_);      //#include "windows.h"
#else
		pthread_mutex_unlock( &*mutex_ );//#include <pthread.h>
#endif

     }
  }

  if (pep.substr(0,1) != "n" && nterm) {
    pep = string("n") + pep;
    has_nterm_mod_ = true;
  } 

  if (pep.substr(pep.find_last_of(string("ABCDEFGHIJKLMNOPQRSTUVWXYZnc")),1) != "c" && cterm) {
  //  if (pep.substr(pep.length()-1,1) != "c" && cterm) {
    pep = pep + string("c");
    has_cterm_mod_ = true;
  } 
  string newpep = "";
  for (int p = 0 ; p < pep.length(); p++) {
    char aa = pep.substr(p,1).c_str()[0];
    char next = p < pep.length()-1 ? pep.substr(p+1,1).c_str()[0] : '\0';

    if (stat_mods_->find(aa) != stat_mods_->end() &&  next != '[') {
      stringstream ss;
#ifdef MSVC
	      WaitForSingleObject( *mutex_,    // handle to mutex
				   INFINITE);      //#include "windows.h"
#else	
	      pthread_mutex_lock( &*mutex_ );//#include <pthread.h>
#endif 
      ss << aa << '[' << (int)(PeptideUser::getAAMonoisotopicMass(aa) + (*stat_mods_)[aa] + 0.5) << ']';
#ifdef MSVC
		ReleaseMutex( *mutex_);      //#include "windows.h"
#else
		pthread_mutex_unlock( &*mutex_ );//#include <pthread.h>
#endif
      newpep += ss.str();
    }
    else {
      newpep += aa;
    }


  }
  char aa; //TEST = 'T';
  double sh;//TEST =  79.966331;	
  
  //pep = newpep;
 
  //  int ret = Peptide::processNewMod(aa, sh, userType, userToken);  

  int ret;

  for (TPP_HASHMAP_T<char, double>::iterator itr = stat_mods_->begin(); itr != stat_mods_->end(); itr++ ) {
    aa = itr->first;
    sh = itr->second;
    string userToken("");
    string userType("");
    stringstream newModTypess;
    
    newModTypess.precision(6);
    newModTypess << "Qtc_" << aa << "_" << sh;
    
    userType = newModTypess.str();

    ostringstream ss;
    ss.str("");
    ss << aa << '[' << (int)(sh + PeptideUser::getAAMonoisotopicMass(aa) + 0.5) << ']';
     
    userToken = ss.str();

#ifdef MSVC
    WaitForSingleObject( *mutex_,    // handle to mutex
			 INFINITE);      //#include "windows.h"
#else	
    pthread_mutex_lock( &*mutex_ );//#include <pthread.h>
#endif 
    ret = pep_->processUserMod(aa, sh, userType, userToken, true);	
#ifdef MSVC
		ReleaseMutex( *mutex_);      //#include "windows.h"
#else
		pthread_mutex_unlock( &*mutex_ );//#include <pthread.h>
#endif
  }

  for (TPP_HASHMAP_T<char, vector<double>*>::iterator itr = var_mods_->begin(); itr != var_mods_->end(); itr++ ) {
    aa = itr->first;

    for (int x = 0; itr->second && x < (*itr->second).size(); x++) {
      sh = (*itr->second)[x];
      string userToken("");
      string userType("");
      stringstream newModTypess;
      
      newModTypess.precision(6);
      newModTypess << "Qtc_" << aa << "_" << sh;
      
      userType = newModTypess.str();



      ostringstream ss;
      ss.str("");
      //ss.precision(0);
      ss << aa << '[' << (int)(sh + PeptideUser::getAAMonoisotopicMass(aa) + 0.5) << ']';
      
      userToken = ss.str();

      if ((*stat_mods_).find(aa)!= (*stat_mods_).end()) {
	sh += (*stat_mods_)[aa];
	
      }
      else if ((*stat_prot_termods_).find(aa)!= (*stat_mods_).end()) {
	sh += (*stat_prot_termods_)[aa];
	
      }
#ifdef MSVC
	      WaitForSingleObject( *mutex_,    // handle to mutex
				   INFINITE);      //#include "windows.h"
#else	
	      pthread_mutex_lock( &*mutex_ );//#include <pthread.h>
#endif 
      ret = pep_->processUserMod(aa, sh, userType, userToken, true);	
#ifdef MSVC
		ReleaseMutex( *mutex_);      //#include "windows.h"
#else
		pthread_mutex_unlock( &*mutex_ );//#include <pthread.h>
#endif
    }
  }

  
  for (TPP_HASHMAP_T<char, double>::iterator itr = stat_prot_termods_->begin(); itr != stat_prot_termods_->end(); itr++ ) {
    aa = itr->first;
    sh = itr->second;
    string userToken("");
    string userType("");

    stringstream newModTypess;
      
      newModTypess.precision(6);
      newModTypess << "Qtc_" << aa << "_" << sh;
      
      userType = newModTypess.str();



      ostringstream ss;
      ss.str("");
      //      ss.precision(0);
      ss << aa << '[' << (int)(sh + PeptideUser::getAAMonoisotopicMass(aa) + 0.5) << ']';
      
      userToken = ss.str();

    //stringstream newModTypess;
	
    //newModTypess.precision(6);
    //newModTypess << "PTMPro_" << aa << "_" << sh;
    
    //userType = newModTypess.str();
#ifdef MSVC
    WaitForSingleObject( *mutex_,    // handle to mutex
			 INFINITE);      //#include "windows.h"
#else	
    pthread_mutex_lock( &*mutex_ );//#include <pthread.h>
#endif 
    ret = pep_->processUserMod(aa, sh, userType, userToken, true);	
#ifdef MSVC
		ReleaseMutex( *mutex_);      //#include "windows.h"
#else
		pthread_mutex_unlock( &*mutex_ );//#include <pthread.h>
#endif
  }

  for (TPP_HASHMAP_T<char, vector<double>*>::iterator itr = var_prot_termods_->begin(); itr != var_prot_termods_->end(); itr++ ) {
    aa = itr->first;

    for (int x = 0; x < itr->second->size(); x++) {
      sh = (*itr->second)[x];
      string userToken("");
      string userType("");
      stringstream newModTypess;
      
      newModTypess.precision(6);
      newModTypess << "Qtc_" << aa << "_" << sh;
      
      userType = newModTypess.str();



      ostringstream ss;
      ss.str("");
      //ss.precision(0);
      ss << aa << '[' << (int)(sh + PeptideUser::getAAMonoisotopicMass(aa) + 0.5) << ']';
      
      userToken = ss.str();

      //stringstream newModTypess;
	
      //newModTypess.precision(6);
      //newModTypess << "PTMPro_" << aa << "_" << sh;
      
      //userType = newModTypess.str();
	
      if ((*stat_mods_).find(aa)!= (*stat_mods_).end()) {
	sh += (*stat_mods_)[aa];
	
      }
      else if ((*stat_prot_termods_).find(aa)!= (*stat_mods_).end()) {
	sh += (*stat_prot_termods_)[aa];
	
      }
#ifdef MSVC
	      WaitForSingleObject( *mutex_,    // handle to mutex
				   INFINITE);      //#include "windows.h"
#else	
	      pthread_mutex_lock( &*mutex_ );//#include <pthread.h>
#endif 
      ret = pep_->processUserMod(aa, sh, userType, userToken, true);	
#ifdef MSVC
		ReleaseMutex( *mutex_);      //#include "windows.h"
#else
		pthread_mutex_unlock( &*mutex_ );//#include <pthread.h>
#endif
    }
  }
  

 
  pep_unmod_ = NULL;
  pep_str_ = new string(pep);
  pep_unmod_str_ = "";

  
  string::size_type pos = 0;
  int pep_pos = 0;

  NAA_ = 1;
  while (pos != string::npos) {
    string token = pep_->nextAAToken(pep, pos, pos);

    if (token.find('n') == std::string::npos) {
      pep_pos++;
    }

    if (pos_mod_hash.find(pep_pos) != pos_mod_hash.end()) {
      double tot_mass = pos_mod_hash[pep_pos];
      string userToken("");
      string userType("");
      


      double sh = tot_mass - PeptideUser::getAAMonoisotopicMass(token[0]);

      for (int s = 0; s < shift_.size(); s++) {
	if (fabs(shift_[s]-sh) < 0.005) {
	  mods_[s]->push_back(pep_pos);
	  nmods_[s]++;	 
	}
      }
      
      stringstream newModTypess;
      
      newModTypess.precision(6);
      newModTypess << "Qtc_" << token[0] << "_" << sh;
      userType = newModTypess.str();

      //ostringstream ss;
      //ss.str("");
      //ss << token[0] << '[' << (int)(sh + Peptide::getAAMonoisotopicMass(token[0]) + 0.5) << ']';
      
      userToken = token; //ss.str();

#ifdef MSVC
      WaitForSingleObject( *mutex_,    // handle to mutex
			   INFINITE);      //#include "windows.h"
#else	
      pthread_mutex_lock( &*mutex_ );//#include <pthread.h>
#endif 
      pep_->processUserMod(token[0], sh, userType, userToken, true);
#ifdef MSVC
      ReleaseMutex( *mutex_);      //#include "windows.h"
#else
      pthread_mutex_unlock( &*mutex_ );//#include <pthread.h>
#endif
    }
    
    
    
    pep_unmod_str_ += token.substr(0,1);
    NAA_++;
  }

  pep_->init(newpep, charge);

 
  cramp_ = cramp;
  scan_ = scan;
  entry_ = NULL;
  charge_ = charge;
  etd_ = false;
  unknown_mod_ = false;
  bool unspec_mod = false;
  bool stat_tok = false;
  string unspec_tok = "";
  mz_tol_ = 0.1;
 
  //COUNT modsites of all types
  int aa_pos = -1;
  int has_nterm_token = 0;
  int has_cterm_token = 0;



  

}



Quantic::~Quantic() {
    if (entry_) delete entry_;
    if (pep_str_) delete pep_str_;   
    entry_ = NULL;
    pep_str_ = NULL;
    
    //delete pep_;
    delete pep_unmod_;
    neutlossquants_.clear();
    neutlosses_.clear();
    neutquants_.clear();
    neutquantmeans_.clear();
    neutquantstdvs_.clear();

    for (int i=0; i<neutquantvect_.size(); i++) {
      for (int ii=0; ii<neutquantvect_[i]->size(); ii++)  {
	(*neutquantvect_[i])[ii]->clear();
      }
      neutquantvect_[i]->clear();
    }
    neutquantvect_.clear();

}

bool Quantic::init() {
  if (unknown_mod_) {
    return false;
  }
    rampScanInfo* scanInfo = NULL;
    rampPeakList* peaks = NULL;
    PeptideUser* pep = pep_;
  
#ifdef MSVC
	      WaitForSingleObject( *mutex_,    // handle to mutex
				   INFINITE);      //#include "windows.h"
#else	
	      pthread_mutex_lock( &*mutex_ );//#include <pthread.h>
#endif 
    scanInfo = cramp_->getScanHeaderInfo(scan_);
    //read the peaks
    peaks = cramp_->getPeakList(scan_);	    
#ifdef MSVC
		ReleaseMutex( *mutex_);      //#include "windows.h"
#else
		pthread_mutex_unlock( &*mutex_ );//#include <pthread.h>
#endif

    double precursorMz = scanInfo->m_data.precursorMZ;
    int peakCount = peaks->getPeakCount();  
    int precursorCharge = scanInfo->m_data.precursorCharge;
    double precursorIntensity = scanInfo->m_data.precursorIntensity;
    double totIonCurrent = scanInfo->m_data.totIonCurrent;
    double retentionTime = scanInfo->m_data.retentionTime;
    string fragType(scanInfo->m_data.activationMethod);

    if (fragType.find("ETD",0)==0) {
      etd_ = true;
    }
    double collisionEnergy = scanInfo->m_data.collisionEnergy;

    entry_ = new SpectraSTLibEntry(pep, "", "Normal", NULL, fragType);

    if (!(fragType.empty()) && entry_->getFragType().empty()) {
      entry_->setFragType(fragType);
    }
    
    // will overwrite the retention time from pepXML file with that from the mzXML file  
    stringstream rtss;
    rtss.precision(1);
    rtss << fixed << retentionTime << ',' << retentionTime << ',' << retentionTime; 
    entry_->setOneComment("RetentionTime", rtss.str());
    
    stringstream ticss;
    ticss.precision(2);
    ticss << totIonCurrent;
    entry_->setOneComment("TotalIonCurrent", ticss.str());
    
    stringstream pintss;
    pintss.precision(2);
    pintss << precursorIntensity;
    entry_->setOneComment("PrecursorIntensity", pintss.str());
    
    //MH: Change from >0 to >=0 because sometimes CE is 0.
    if (collisionEnergy >= 0) {
      stringstream cess;
      cess.precision(1);
      cess << fixed << collisionEnergy;
      entry_->setOneComment("CollisionEnergy", cess.str());
    } else {
      stringstream errss;
      errss << "Scan #" << scan_ << " has collision energy below zero: " << collisionEnergy;
    }
    
    
    if (precursorCharge < 1) precursorCharge = 0;
    
    //  cout << "inserting peaks " << peakCount << endl;
    // create the peak list and read the peaks one-by-one
    //if (entry_->getPeakList()->capacity() < peakCount)
    //  entry_->getPeakList()->reserve(peakCount);


    for (int j = 0; j < peakCount; j++) {
      double mz = peaks->getPeak(j)->mz;
      float intensity = (float)(peaks->getPeak(j)->intensity);
      
      if (intensity > 0.0) {
	entry_->getPeakList()->insert(mz, intensity, "", "");
      }
    }
    
    if (precursorIntensity > 0.0) {
      entry_->getPeakList()->setPrecursorIntensity(precursorIntensity);
    }
    
    stringstream maxss;
    maxss.precision(2);
    maxss << entry_->getPeakList()->getOrigMaxIntensity();
    entry_->setOneComment("OrigMaxIntensity", maxss.str());
    
    delete scanInfo;
    //delete peaks;
    
    entry_->annotatePeaks(true);
    
    
    return pep_->isGood();
}

void Quantic::computeMassDiffs(vector<double>* massdiffs) {
  setMZTolerance(0.5);
  processPeakList() ;
  if (!pep_ ||
      !entry_->getPeakList() ||
      !entry_->getPeptidePtr()->isModsSet) {
    return;
  }
  PeptideUser* mpep = pep_;
   
  if (nions_.empty()) {
    if (etd_) {
      nions_.push_back('c');
    }
    else {
      nions_.push_back('b');
      nions_.push_back('a');
    }
  }

  if (cions_.empty()) {
    if (etd_) {
      cions_.push_back('z');
    }
    else {
      cions_.push_back('y');
    }
  }
  
  char ion1;// = 'b';
  char ion2;// = 'y';

  double mzTOL = mz_tol_;

  //if (etd_) {
  //  ion1 = 'c';
  //  ion2 = 'z';
  //}
  
  double mz, foundInt, foundMZ;
  
  Peak* foundPK = NULL;
  
  map<double, double>* foundMZs = new map<double,double>();
  for (int iso = 0; 0 && iso < 7; iso++) {
    mz = iso*(PROTONMASS/ charge_) + (mpep->monoisotopicNeutralM() + charge_ * PROTONMASS) / charge_ ;
    foundPK = peakList_->findNearestPeakPtr<double>(mz, mzTOL, foundMZs);
    if (foundPK) {
      double wt = 1;
      
      foundMZ = foundPK->mz;
      
      
      foundInt = foundPK->intensity;
      
      massdiffs->push_back(mz-foundMZ);
      
      foundMZs->insert(make_pair(foundMZ, wt*foundInt));
      
    }
  }
  
  
  unsigned int max_ch = (unsigned int)mpep->charge < MAX_CH ?  (unsigned int)mpep->charge : MAX_CH;
  
  string::size_type pos = 0;
  while (pos++ < pep_unmod_str_.length()) {
    string token = pep_->getModToken(*pep_unmod_str_.substr(pos,1).c_str(),
				     pep_->mods.find(pos-nterm_)->second);
    
    foundPK = NULL;
    mz = -1;
    map<string, vector<double>>::iterator im_itr;
   
    im_itr = mpep->AAMonoisotopicImmoniumTable->find(token);
    
    int SZ = 0;

    if (im_itr != mpep->AAMonoisotopicImmoniumTable->end()) {
      SZ = im_itr->second.size();
    }

    for ( int im=0; im < SZ; im++) {
      
      
    	mz = (*mpep->AAMonoisotopicImmoniumTable)[token][im];

      if (mz>0) {
    	foundPK = peakList_->findPeakPtr<double>(mz, mzTOL, foundMZs);
	
    	if (foundPK) {
    	  foundMZ = foundPK->mz;
	
	  
	  foundInt = foundPK->intensity;
	  if ((foundMZs->find(foundMZ)) == foundMZs->end()) {
	    double wt = 1;
	    massdiffs->push_back(mz-foundMZ);
    	    foundMZs->insert(make_pair(foundMZ, wt*foundInt));
	    
	  }
	}
      }
    }
  }

  for (int i=1; i<mpep->NAA() ; i++) {
    for (unsigned int ch = 1; ch <= max_ch ; ch++) {
      for (int ni=0; ni<nions_.size(); ni++) {
	ion1 = nions_[ni];
	foundPK = NULL;
	if (ion1 == 'a' && i == 1) {
	  continue;
	}
	
	foundPK = NULL;
	
	mz = mpep->monoisotopicMZFragment(ion1, i, ch);
	if (mz>0) {
	  foundPK = peakList_->findPeakPtr<double>(mz, mzTOL, foundMZs);
	  
	  foundInt = 0;
	  if (foundPK) {
	    foundMZ = foundPK->mz;
	    
	    
	    foundInt = foundPK->intensity;
	    
	    if (( foundMZs->find(foundMZ)) == foundMZs->end()) {
	      //if (1) {
	      double wt = 1;// - sqrt(fabs(foundMZ-mz)/mzTOL);
	      
	      massdiffs->push_back(mz-foundMZ);
	      foundMZs->insert(make_pair(foundMZ, wt*foundInt));
	    }
	    
	  }
	  
	}
      }
      for (int ci=0; ci<cions_.size(); ci++) {
	ion2 = cions_[ci];	
	foundPK = NULL;
	mz = mpep->monoisotopicMZFragment(ion2, mpep->NAA() - i, ch);
	if (mz>0) {
	  foundPK = peakList_->findPeakPtr<double>(mz, mzTOL, foundMZs);
	  
	  foundInt = 0;
	  if (foundPK) {
	    foundMZ = foundPK->mz;
	    
	    
	    foundInt = foundPK->intensity;
	    
	    if ((foundMZs->find(foundMZ)) == foundMZs->end()) {
	      //if (1) {
	      double wt = 1;// - sqrt(fabs(foundMZ-mz)/mzTOL);
	      massdiffs->push_back(mz-foundMZ);
	      foundMZs->insert(make_pair(foundMZ, wt*foundInt));
	    }
	  }
	}
      }
    }
  }
}

void Quantic::processPeakList() {
  peakList_ = entry_->getPeakList();
  //peakList_->subtractNoise();
  //double retainedFraction = peakList_->clean(true, true, 0.0); // de-isotope, remove near-precursor ions, do not remove light ions

  //  for (int ix=0; ix<peakList->getNumPeaks(); ix++) 
  Peak* pk = new Peak();
  int n =  peakList_->getNumPeaks() ;

  if (n>0) {
	  peakList_->getNthLargestPeak(peakList_->getNumPeaks(), (*pk));

	  minInt_ = (*pk).intensity;

	  peakList_->getNthLargestPeak(1, (*pk));



	  maxInt_ = (*pk).intensity;
	  //peakList_->clean(true,true, 0);
  }

  n =  peakList_->getNumPeaks() ;

  if (n>0) {
	  peakList_->getNthLargestPeak(peakList_->getNumPeaks(), (*pk));

	  minInt_ = (*pk).intensity;

	  peakList_->getNthLargestPeak(1, (*pk));


	  maxInt_ = (*pk).intensity;
  }
  else {
	  minInt_ = -1;
	  maxInt_ = -1;
  }

  TIC_ =  peakList_->getTotalIonCurrent();
  delete pk;
}






void Quantic::setPrecision(ostringstream & stream , double& value) {
  stream.precision(0);
  stream.setf(ios::fixed); 
  return;
}



void Quantic::evaluatePeptideAnnotatedTIC(map<double, double>* avoidMZ) {
 
  if (!pep_ ||
      !entry_->getPeakList() ||
      !entry_->getPeptidePtr()->isModsSet) {
    return;
  }

  PeptideUser* mpep = pep_;

  
  

  int nDECOY = 0;
  double mzTOL = mz_tol_;
  
  if (mz_tol_ <= 0) {
    mzTOL = 0.1;
  }

  int idx;


  double tmpFor = 0;
  double tmpAgainst = 0;
  double minFor = 0;
  double minAgainst = 0;
  
  annot_tic_ = 0;
  annot_str_ = "";
  theor_str_ = "";

  vector<double> forEvidenceByType;
  vector<int> forCountByType;

  //vector<double> quantEvid;

  vector<NeutQuant*> quantEvid;



  int nsites =0;
  int nmods = 0;
  int ntypes = 0;

   // check b ions from left+1 to pos, y ions from len-pos to len-left-1
  
  
  int left = -1;
  int right = -1;

  int N = 0;


  int forCount = 0;
  int againstCount = 0;

  double forEvidence = minInt_;
  double againstEvidence = minInt_;
  double ratio = maxInt_/minInt_; 
  double altRatio = maxInt_/minInt_; //for the other position

  double decoyForEvidence = 0;//minInt_;
  double decoyAgainstEvidence = 0;//minInt_;

  if (nions_.empty()) {
    if (etd_) {
      nions_.push_back('c');
    }
    else {
      nions_.push_back('b');
      nions_.push_back('a');
    }
  }

  if (cions_.empty()) {
    if (etd_) {
      cions_.push_back('z');
    }
    else {
      cions_.push_back('y');
    }
  }
  
  char ion1;// = 'b';
  char ion2;// = 'y';


  //if (etd_) {
  //  ion1 = 'c';
  //  ion2 = 'z';
  //}

  double siteProbSum = 0;
  double Xsq = 0;
  left = -1;
  right = -1;
 
  forEvidence = 0; 
  double directEvidence = 0;
  decoyForEvidence = 0;
  
  double wtSum = 1;

  double mz, foundInt, foundMZ;

  Peak* foundPK = NULL;

  
      
  map<double, double>* foundMZs = new map<double,double>();
  for (int type = 0; type < shift_.size(); type++) {
    for (int n = 0; n < neutlosses_[type]->size(); n++) {
      neutlossquants_[type]->push_back(new NeutQuant((*neutlosses_[type])[n], 0, 0, 0));
      neutquants_[type]->push_back(0);
      neutquantmeans_[type]->push_back(0);
      neutquantstdvs_[type]->push_back(0);
      neutquantvect_[type]->push_back(new vector<double>());
    }
  }
  for (int iso = 0; iso < 7; iso++) {
    mz = iso*(PROTONMASS/ charge_) + (mpep->monoisotopicNeutralM() + charge_ * PROTONMASS) / charge_ ;
    foundPK = peakList_->findPeakPtr<double>(mz, mzTOL, foundMZs);

    theor_str_ += "M,"+numberToString(iso)+","+
      numberToString(mz,3) +";";

    if (foundPK) {
      double wt = 1;

      foundMZ = foundPK->mz;
      
      if (!avoidMZ || avoidMZ->find(foundMZ)==avoidMZ->end()) {
	 directEvidence += foundPK->intensity;
	 foundInt = foundPK->intensity;
	 
	 forEvidence += foundInt;
	 foundMZs->insert(make_pair(foundMZ, wt*foundInt));
	 annot_tic_ += wt*foundInt;
	 annot_str_ += "M,"+numberToString(iso)+","+
	   numberToString(foundMZ,3) +","+
	   numberToString(foundMZ-mz,3) +","+
	   numberToString(foundInt,1) +";";
	 forCount++;
      }
      
    }
 
    //BEGIN UNFRAG PRECURSOR NEUTRAL LOSSES
    double sav_mz = mz;
    int ch = charge_;
    double foundEvid, prevEvid;
    double wt = 1;
    for (int type = 0; type < shift_.size(); type++) {
      if (nmods_[type] > 0 && type < neutlosses_.size()) {
	for (int n = 0; n < neutlosses_[type]->size(); n++) {

	  if (n+1>quantEvid.size()) {
	    // quantEvid.push_back(0);
	    quantEvid.push_back(new NeutQuant((*neutlosses_[type])[n], 0, 0, 0));
	  }
	  else {
	    quantEvid[n]->quant_ = 0;
	  }

	  foundPK = NULL;
	  
	  mz = sav_mz+(*neutlosses_[type])[n]/ch;
	  if (mz>0) {
	    
	    // mzTOL = fragPPMTOL * mz / 1e6;
	    
	    foundPK = peakList_->findPeakPtr(mz, mzTOL, foundMZs);

	    theor_str_ += "M("+numberToString((*neutlosses_[type])[n],1)+"),"
	      +numberToString(iso)+","+
	      numberToString(mz,3) +";";

	    foundInt = 0;
	    if (foundPK) {
	      foundMZ = foundPK->mz;
	      if (!avoidMZ || avoidMZ->find(foundMZ)==avoidMZ->end()) {
		directEvidence += foundPK->intensity;
		foundInt = foundPK->intensity;
		
		forEvidence += foundInt;
		foundMZs->insert(make_pair(foundMZ, wt*foundInt));
		annot_tic_ += wt*foundInt;
		annot_str_ += "M("+numberToString((*neutlosses_[type])[n],1)+"),"+numberToString(iso)+","+
		  numberToString(foundMZ,3) +","+
		  numberToString(foundMZ-mz,3) +","+
		  numberToString(foundInt,1) +";";
		forCount++;

		if (iso == 0 && pep_coupled_) {
		  quantEvid[n]->quant_ = foundInt;
		  //(*neutquantvect_[type])[n]->push_back(foundInt);
		}
	
	      }
	    }
	  }
	  mz = sav_mz;
	}
	//	if (pep_coupled_) {
	//DDS TODO DIAMODE CORRECT ISOS HERE
	  
	  if (diaMode_) {
	    
	    //	    string frag = string(mpep->stripped);
	    string frag = *pep_str_;
	    
	    //frag = frag.substr(i);
	    
	    string molecform = pepToMolecForm(frag);
	    
	    mercury_->GoMercury((char*)molecform.c_str(), ch);
	    correctIsotopeErrors(quantEvid, type, mercury_);
	  }
	  
	  updateQuants(type, quantEvid);

	  //}
      }
    }
  }

  map<double, double>::iterator itr;

  gsl_vector* r;
  
  
  unsigned int max_ch = (unsigned int)mpep->charge < MAX_CH ?  (unsigned int)mpep->charge : MAX_CH;

  string::size_type pos = 0;

     
  while (pos++ < pep_unmod_str_.length()) {
    string token = pep_->getModToken(*pep_unmod_str_.substr(pos,1).c_str(),
				      pep_->mods.find(pos-nterm_)->second);
    
    foundPK = NULL;
    mz = -1;
    map<string, vector<double>>::iterator im_itr;
   
    im_itr = mpep->AAMonoisotopicImmoniumTable->find(token);
    
    int SZ = 0;

    if (im_itr != mpep->AAMonoisotopicImmoniumTable->end()) {
      SZ = im_itr->second.size();
    }

    for ( int im=0; im < SZ; im++) {
      
      mz = (*mpep->AAMonoisotopicImmoniumTable)[token][im];


      if (mz>0) {
    	foundPK = peakList_->findPeakPtr<double>(mz, mzTOL, foundMZs);

	theor_str_ += "I,-,"+
	   numberToString(mz,3) +";";

    	if (foundPK) {
    	  foundMZ = foundPK->mz;
	  if (!avoidMZ || avoidMZ->find(foundMZ)==avoidMZ->end()) {
	    directEvidence = foundPK->intensity;
	  
	    foundInt = foundPK->intensity;
	    if ((itr = foundMZs->find(foundMZ)) == foundMZs->end()) {
	      //if (1) {
	      double wt = 1;// - sqrt(fabs(foundMZ-mz)/mzTOL);
	      forEvidence += wt*foundInt;
	      //wtSum += wt;
    	    foundMZs->insert(make_pair(foundMZ, wt*foundInt));
	    annot_tic_ += wt*foundInt;
	    annot_str_ += "I,-,"+
	      numberToString(foundMZ,3) +","+
	      numberToString(foundMZ-mz,3) +","+
	      numberToString(foundInt,1) +";";
    	    forCount++;
	    }
	    else {
	      double wt = 1;// - sqrt(fabs(foundMZ-mz)/mzTOL);
	      if (itr->second <  wt*foundInt) {
		forEvidence += wt*foundInt - itr->second;
		(*foundMZs)[foundMZ] = wt*foundInt;
	      }
	    }
	  }
    	}
      }
    }
  }
  double foundEvid = 0;
  double prevEvid = 0;
  for (int i=1; i<mpep->NAA() ; i++) {
    for (unsigned int ch = 1; ch <= max_ch ; ch++) {
      for (int ni=0; ni<nions_.size(); ni++) {
	ion1 = nions_[ni];
	foundPK = NULL;
	if (ion1 == 'a' && i == 1) {
	  continue;
	}
	
	foundPK = NULL;
	
	mz = mpep->monoisotopicMZFragment(ion1, i, ch);
	if (mz>0) {
	  foundPK = peakList_->findPeakPtr<double>(mz, mzTOL, foundMZs);
	  theor_str_ += ion1+numberToString(i)+","+
	    numberToString(ch)+","+
	    numberToString(mz,3) +";"; 
	  foundInt = 0;
	  if (foundPK) {
	    foundMZ = foundPK->mz;
	    if (!avoidMZ || avoidMZ->find(foundMZ)==avoidMZ->end()) {
	      directEvidence += foundPK->intensity;
	      
	      foundInt = foundPK->intensity;
	      
	      foundEvid = 0;
	      prevEvid = 0;
	      if ((itr = foundMZs->find(foundMZ)) == foundMZs->end()) {
		//if (1) {
		double wt = 1;// - sqrt(fabs(foundMZ-mz)/mzTOL);
		forEvidence += wt*foundInt;
		foundEvid = wt*foundInt;
		//wtSum += wt;
		foundMZs->insert(make_pair(foundMZ, wt*foundInt));
		annot_tic_ += wt*foundInt;
		annot_str_ += ion1+numberToString(i)+","+numberToString(ch)+","+
		  numberToString(foundMZ,3) +","+
		  numberToString(foundMZ-mz,3) +","+
		  numberToString(foundInt,1) +";"; 
		forCount++;
	      }
	      else {
		double wt = 1;// - sqrt(fabs(foundMZ-mz)/mzTOL);
		if (itr->second <  wt*foundInt) {	      
		  forEvidence += wt*foundInt - itr->second;
		  prevEvid =  itr->second;
		  foundEvid = wt*foundInt;
		  (*foundMZs)[foundMZ] = wt*foundInt;
		}
	      }
	      
	    }
	  }   
	 
	  
	  
	  //BEGIN NEUTRAL LOSSES
	  double sav_mz = mz;
	  //int ch = charge_;
	  double wt = 1;
	  double foundEvid, prevEvid;
	  for (int type = 0; type < shift_.size() && nmods_[type] > 0 ; type++) {
	    bool have_mod = false;
	    for (int m = 0; m < mods_[type]->size(); m++) {
	      if ((*mods_[type])[m] <= i) {
		have_mod = true;
		break;
	      }
	    }
	    
	    if (have_mod && type < neutlosses_.size()) {
	      for (int n = 0; n < neutlosses_[type]->size(); n++) {
		if (n+1>quantEvid.size()) {
		  quantEvid.push_back(new NeutQuant((*neutlosses_[type])[n], 0, 0, 0));
		}
		else {
		  quantEvid[n]->quant_ = 0;
		}
		
		foundPK = NULL;
		
		mz = sav_mz+(*neutlosses_[type])[n]/ch;
		if (mz>0) {
		  
		  foundPK = peakList_->findPeakPtr(mz, mzTOL, foundMZs);

		  theor_str_ +=  ion1+numberToString(i)+
		    "("+numberToString((*neutlosses_[type])[n],1)+"),"+
		    numberToString(ch)+","+
		    numberToString(mz,3) +";";

		  foundInt = 0;
		  if (foundPK) {
		    foundMZ = foundPK->mz;
		    if (!avoidMZ || avoidMZ->find(foundMZ)==avoidMZ->end()) {
		      directEvidence += foundPK->intensity;
		      foundInt = foundPK->intensity;
		      
		      forEvidence += foundInt;
		      foundMZs->insert(make_pair(foundMZ, wt*foundInt));
		      annot_tic_ += wt*foundInt;
		      annot_str_ +=  ion1+numberToString(i)+"("+numberToString((*neutlosses_[type])[n],1)+"),"+numberToString(ch)+","+
			numberToString(foundMZ,3) +","+
			numberToString(foundMZ-mz,3) +","+
			numberToString(foundInt,1) +";";
		      forCount++;
		      //	    (*neutquants_[type])[n] += foundInt;
		      if (i > 1 && !pep_coupled_) {
			quantEvid[n]->quant_ = foundInt;
			//(*neutquantvect_[type])[n]->push_back(foundInt);
		      }
		      
		    }
		  }
		}
		mz = sav_mz;
	      }
	      if (diaMode_) {

		//DDS TODO DIAMODE CORRECT ISOS HERE
		//TODO get correct molec.form including mods
		//string frag = string(mpep->stripped);

		string frag = *pep_str_;
		string allAA = ALL_AA;
		int Nterm_ion = 0;

		for (int p=0; p < pep_str_->length(); p++) {
		  if (allAA.find((*pep_str_)[p]) != string::npos) {
		    Nterm_ion++;
		  }
		  if (Nterm_ion > i) {
		    frag = frag.substr(0,p);
		    break;
		  }

		}
		
		//frag = frag.substr(0,i);
		
		string molecform = pepToMolecForm(frag);

		mercury_->GoMercury((char*)molecform.c_str(), ch);
		correctIsotopeErrors(quantEvid, type, mercury_);
	      }
	      
	      if (!pep_coupled_)
		updateQuants(type, quantEvid);
	    }
	  }
	}
      }      
      



      for (int ci=0; ci<cions_.size(); ci++) {
	ion2 = cions_[ci];	
	foundPK = NULL;
	mz = mpep->monoisotopicMZFragment(ion2, mpep->NAA() - i, ch);
	if (mz>0) {
	  foundPK = peakList_->findPeakPtr<double>(mz, mzTOL, foundMZs);

	  theor_str_ += ion2+numberToString(mpep->NAA() -i)+","+
	    numberToString(ch)+","+
	    numberToString(mz,3) +";"; 

	  foundInt = 0;
	  if (foundPK) {
	    foundMZ = foundPK->mz;
	    if (!avoidMZ || avoidMZ->find(foundMZ)==avoidMZ->end()) {
	      directEvidence += foundPK->intensity;
	      
	      foundInt = foundPK->intensity;
	      foundEvid = 0;;
	      prevEvid = 0;
	      if ((itr = foundMZs->find(foundMZ)) == foundMZs->end()) {
		//if (1) {
		double wt = 1;// - sqrt(fabs(foundMZ-mz)/mzTOL);
		forEvidence += wt*foundInt;
		foundEvid = wt*foundInt;
		//wtSum += wt;
		foundMZs->insert(make_pair(foundMZ, wt*foundInt));
		annot_tic_ += wt*foundInt;
		annot_str_ += ion2+numberToString(mpep->NAA() -i)+","+numberToString(ch)+","+
		  numberToString(foundMZ,3) +","+
		  numberToString(foundMZ-mz,3) +","+
		  numberToString(foundInt,1) +";"; 
		forCount++;
	      }
	      else {
		double wt = 1;// - sqrt(fabs(foundMZ-mz)/mzTOL);
		if (itr->second <  wt*foundInt) {
		  forEvidence += wt*foundInt - itr->second;
		  prevEvid =  itr->second;
		  foundEvid = wt*foundInt;
		  (*foundMZs)[foundMZ] = wt*foundInt;
		  
		}
	      }
	    }
	  }
	  
	  //BEGIN NEUTRAL LOSSES
	  double sav_mz = mz;
	  //	int ch = charge_;
	  double wt = 1;
	  double foundEvid, prevEvid;
	  for (int type = 0; type < shift_.size() && nmods_[type] > 0; type++) {
	    bool have_mod = false;
	    for (int m = 0; m < mods_[type]->size(); m++) {
	      if ((*mods_[type])[m] > i) {
		have_mod = true;
		break;
	      }
	    }
	    
	    if (have_mod && type < neutlosses_.size()) {
	      for (int n = 0; n < neutlosses_[type]->size(); n++) {
		if (n+1>quantEvid.size()) {
		  quantEvid.push_back(new NeutQuant((*neutlosses_[type])[n], 0, 0, 0));
		}
		else {
		  quantEvid[n]->quant_ = 0;
		}
		
		foundPK = NULL;
		
		mz = sav_mz+(*neutlosses_[type])[n]/ch;
		if (mz>0) {
		  
		  foundPK = peakList_->findPeakPtr(mz, mzTOL, foundMZs);
		  theor_str_ +=  ion2+numberToString(mpep->NAA() -i)+
		    "("+numberToString((*neutlosses_[type])[n],1)+"),"+
		    numberToString(ch)+","+
		    numberToString(mz,3) +";";
		  foundInt = 0;
		  if (foundPK) {
		    foundMZ = foundPK->mz;
		    if (!avoidMZ || avoidMZ->find(foundMZ)==avoidMZ->end()) {
		      directEvidence += foundPK->intensity;
		      foundInt = foundPK->intensity;
		      
		      forEvidence += foundInt;
		      foundMZs->insert(make_pair(foundMZ, wt*foundInt));
		      annot_tic_ += wt*foundInt;
		      annot_str_ +=  ion2+numberToString(mpep->NAA() -i)+"("+numberToString((*neutlosses_[type])[n],1)+"),"+numberToString(ch)+","+
			numberToString(foundMZ,3) +","+
			numberToString(foundMZ-mz,3) +","+
			numberToString(foundInt,1) +";";
		      forCount++;
		      //(*neutquants_[type])[n] += foundInt;
		      if ( mpep->NAA()-i > 1 && !pep_coupled_){
			quantEvid[n]->quant_ = foundInt;
			//(*neutquantvect_[type])[n]->push_back(foundInt);
		      }
		      
		    }
		  }
		}
		mz = sav_mz;
	      }


	      if (diaMode_) {

		//DDS TODO DIAMODE CORRECT ISOS HERE
		//TODO get correct molec.form including mods
		//string frag = string(mpep->stripped);

		string frag = *pep_str_;
		string allAA = ALL_AA;
		int Cterm_ion = 0;

		for (int p=pep_str_->length()-1; p >= 0; p--) {
		  if (allAA.find((*pep_str_)[p]) != string::npos) {
		    Cterm_ion++;
		  }
		  if (Cterm_ion == mpep->NAA()-i) {
		    frag = frag.substr(p);
		    break;
		  }

		}
		

		
		//frag = frag.substr(i);
		
		string molecform = pepToMolecForm(frag);

		mercury_->GoMercury((char*)molecform.c_str(), ch);
		correctIsotopeErrors(quantEvid, type, mercury_);
	      }
	      
	      if (!pep_coupled_) 
		updateQuants(type, quantEvid);
	      
	    }
	  }
	}
      }
     
      minFor = 0;
      //      r = gsl_vector_calloc(nDECOY);
      for (idx = 0; idx < nDECOY; idx++) {
	//NULL distro using randomized peaklist
	 tmpFor = 0;

	 for (int ni=0; ni<nions_.size(); ni++) {
	  ion1 = nions_[ni];
	 
	  
	  mz = mpep->monoisotopicMZFragment(ion1, i, ch);
	  if (mz > 0) {
	    foundInt = (*decoyPeakLists_)[idx]->findPeak(mz, mzTOL);
	    tmpFor += foundInt;
	  }
	}
	
	for (int ci=0; ci<cions_.size(); ci++) {
	  ion2 = cions_[ci];
	  tmpFor = 0;
	  
	  mz = mpep->monoisotopicMZFragment(ion2, mpep->NAA() - i, ch);
	  if (mz > 0) {
	    foundInt = (*decoyPeakLists_)[idx]->findPeak(mz, mzTOL);
	    //decoyForEvidence += foundInt;
	    tmpFor += foundInt;
	  }
	}
	
	gsl_vector_set(r, idx, tmpFor);
	minFor += tmpFor;
	
      }
     
    }
  }
  //  delete foundMZs;
  
  forEvidence /= wtSum;
  
  if (decoyForEvidence < minInt_) {
    decoyForEvidence = minInt_;
  }

  if (forEvidence < minInt_) {
    forEvidence = minInt_;
  }
  

  
  double totEvidence = forEvidence;

  if (top_peaks_) {
    annot_tic_ = 0;
    
    vector<double>* foundIntens = new vector<double>();
    
    for ( itr =  foundMZs->begin(); itr != foundMZs->end(); itr++ ) {
      foundIntens->push_back(itr->second);
    }
    
    std::sort(foundIntens->begin(),  foundIntens->end(), compareIntens); 
    
    for (int i=0; i < ( top_peaks_ > foundIntens->size() ? foundIntens->size() : top_peaks_ ); i++) {
      annot_tic_ += (*foundIntens)[i];

    }
    foundIntens->clear();
    delete foundIntens;
  }
  

  if (foundMZs) {
    foundMZs->clear();
    delete foundMZs;
  }

  for (int type = 0; type < shift_.size(); type++) {
    calcQuantStats(type);
    for (int n=0; n < neutquantvect_[type]->size(); n++) {
      (*neutquantvect_[type])[n]->clear();
    }
  }
  
  quantEvid.clear();
  
  return;
}


std::string Quantic::numberToString(double num, int prec) {

  std::ostringstream streamObj;
  streamObj << std::fixed;
  streamObj << std::setprecision(prec);
  streamObj << num;

  // Get string from output string stream
  std::string strObj = streamObj.str();

  return strObj;
}


std::string Quantic::numberToString(long num) {

  std::ostringstream streamObj;

  streamObj << num;

  // Get string from output string stream
  std::string strObj = streamObj.str();

  return strObj;
}

void Quantic::correctIsotopeErrors(vector<NeutQuant*>& input, int type, CMercury8* mercury) {
  double weight = 1;
  vector<double> output;
  // vector<double> outputmeans;
  std::sort(input.begin(), input.end(), compareNeutQuants);

  if (0) { //TESTING 
    for (int n=0; n < neutlosses_[type]->size(); n++) {
      if (n==0) input[n]->quant_ = 118728440;
      if (n==1) input[n]->quant_ = 124059288;
      if (n==2) input[n]->quant_ = 58313872;
    }
  }
  
  for (int n=0; n < neutlosses_[type]->size(); n++) {
    output.push_back(0.);
    output[n] = input[n]->quant_;
    
    // outputmeans.push_back(0.);
    //outputmeans[n] = input[n]->mean_;;

    for (int nn=0; nn < n; nn++) {
      
      weight = 0;
      
      //losses are sorted smallest m/z (biggest loss / lowest peak) largest m/z

      weight += mercury->FixedData[n-nn].data/100;      
      
      //for (int nnn = 1; nnn <= n; nnn++) {

	//if (nnn == n)
	//weight += pow(mercury->FixedData[n-nnn+1].data/100, nnn);
	//else 
	//  weight -= pow(mercury->FixedData[n-nnn+1].data/100, nnn);
      // }
      //      output[n] -= weight*input[nn]->quant_;

      output[n] -= weight*output[nn];
      //      outputmeans[n] -= weight*input[nn]->mean_;
      
    }
  }

  // if (0) {
  //   double sum = 0;
  //   for (int n=0; n<output.size(); n++) {
  //     sum += output[n];
  //   }
    
  //   for (int n=0; n<output.size(); n++) {
  //     outputmeans[n] =    output[n] / sum;
  //   }
  // }
  
  for (int n=0; n < neutlosses_[type]->size(); n++) {
    input[n]->quant_ = output[n];
    //input[n]->mean_ = outputmeans[n];
    
  }
}




void Quantic::correctIsotopeErrors(int type, CMercury8* mercury) {
  double weight = 1;
  vector<double> output;
  vector<double> outputmeans;
  std::sort(neutlossquants_[type]->begin(), neutlossquants_[type]->end(), compareNeutQuants);

  if (0) { //TESTING 
    for (int n=0; n < neutlosses_[type]->size(); n++) {
      if (n==0) (*neutlossquants_[type])[n]->quant_ = 118728440;
      if (n==1) (*neutlossquants_[type])[n]->quant_ = 124059288;
      if (n==2) (*neutlossquants_[type])[n]->quant_ = 58313872;
    }
  }
  
  for (int n=0; n < neutlosses_[type]->size(); n++) {
    output.push_back(0.);
    output[n] = (*neutlossquants_[type])[n]->quant_;
    
    outputmeans.push_back(0.);
    outputmeans[n] = (*neutlossquants_[type])[n]->mean_;;

    for (int nn=0; nn < n; nn++) {
      
      weight = 0;
      
      //losses are sorted smallest m/z (biggest loss / lowest peak) largest m/z

      weight += mercury->FixedData[n-nn].data/100;      
      
      //for (int nnn = 1; nnn <= n; nnn++) {

	//if (nnn == n)
	//weight += pow(mercury->FixedData[n-nnn+1].data/100, nnn);
	//else 
	//  weight -= pow(mercury->FixedData[n-nnn+1].data/100, nnn);
      // }
      //      output[n] -= weight*(*neutlossquants_[type])[nn]->quant_;

      output[n] -= weight*output[nn];
      //      outputmeans[n] -= weight*(*neutlossquants_[type])[nn]->mean_;
      
    }
  }
  double sum = 0;
  for (int n=0; n<output.size(); n++) {
    sum += output[n];
  }

  for (int n=0; n<output.size(); n++) {
    outputmeans[n] =    output[n] / sum;
  }
  
  for (int n=0; n < neutlosses_[type]->size(); n++) {
    (*neutlossquants_[type])[n]->quant_ = output[n];
    (*neutlossquants_[type])[n]->mean_ = outputmeans[n];
    
  }
}

void Quantic::setMolForms(string& ptmMolForms) {
  string buffer = ptmMolForms;
  
  string key = "";
  string value = "";

  int sep = buffer.find("]", 0);

  
  while (sep != string::npos) {
    key = buffer.substr(0, sep+1);

    buffer = buffer.substr(sep+1);
    int comma = buffer.find(",", 0);    
    if (comma != string::npos) {
      value = buffer.substr(0, comma);
      buffer = buffer.substr(comma+1);
    }
    else {
      value = buffer.substr(0);
    }
    ptm_mol_hash_.insert(make_pair(*key.c_str(), value));

    sep = buffer.find("]", 0);

  }
  
  
}


string Quantic::pepToMolecForm(string& pep) {

  string allAA = ALL_AA;
  allAA += "nc";
  string allDIGS = "0123456789";
  float Cn = 0;
  float Hn = 2;
  float On = 1;
  float Nn = 0;
  float Pn = 0;
  float Sn = 0;

  string ptmTok = "";
  char aa;
  char elem;
  for (int i=0; i< pep.length(); i++) {
    if (allAA.find(pep[i]) != string::npos) {

      aa = pep[i];
      ptmTok = "";
      if (i+1 < pep.length()) {
	//Not at end
	if (pep[i+1] == '[') {
	  ptmTok = pep[i++];
	  while (i < pep.length() && pep[i] != ']') {
	    ptmTok += pep[i++];
	  }
	  ptmTok += "]";
	}
	
      }
      
      if (ptm_mol_hash_.find(*ptmTok.c_str()) != ptm_mol_hash_.end()) {

	string molForm = ptm_mol_hash_[*ptmTok.c_str()];

	string countStr = "";
	int count = 1;

	for (int p=0; p< molForm.length(); p++) {
	  countStr = "";
	  count = 1;
	  elem = molForm[p];
	  
	  while (p+1 < molForm.length() && allDIGS.find(molForm[p+1]) != string::npos) {
	    countStr += molForm[++p];
	  }

	  if (!countStr.empty()) {
	    count = atoi(countStr.c_str());
	  }

	  switch (elem) {
	    
	  case 'C':
	    Cn += count;
	    break;
	  case 'H':
	    Hn += count;
	    break;	    
	  case 'O':
	    On += count;
	    break;
	  case 'N':
	    Nn += count;
	    break;
	  case 'P':
	    Pn += count;
	    break;
	  case 'S':
	    Sn += count;
	    break;
	    
	  default:
	    cerr << "WARNING: unknown element " << elem << endl;
	    break;
	  }
	}
	

	if (ptmTok[0] == 'n' && Hn >= 1) {
	  Hn -= 1;
	}
	else if (ptmTok[0] == 'c' && Hn >= 1 && On >= 1) {
	  Hn -= 1;
	  On -= 1;
	}

	
	continue;
      }
    
      

      switch (aa) {

      case 'G':

	Cn += 2;
	Hn += 3;
	On += 1;
	Nn += 1;
	Pn += 0;
	Sn += 0;
	
	break;
	
      case 'A':
	
	Cn += 3;
	Hn += 5;
	On += 1;
	Nn += 1;
	Pn += 0;
	Sn += 0;
	
	break;
	
      case 'S':
	
	Cn += 3;
	Hn += 5;
	On += 2;
	Nn += 1;
	Pn += 0;
	Sn += 0;
	
	break;
	
      case 'P':
	
	Cn += 5;
	Hn += 7;
	On += 1;
	Nn += 1;
	Pn += 0;
	Sn += 0;
	
	break;
	
      case 'V':
	
	Cn += 5;
	Hn += 9;
	On += 1;
	Nn += 1;
	Pn += 0;
	Sn += 0;
	
	break;
	
      case 'T':
	
	Cn += 4;
	Hn += 7;
	On += 2;
	Nn += 1;
	Pn += 0;
	Sn += 0;
	
	break;
	
      case 'C':
	
	Cn += 3;
	Hn += 5;
	On += 1;
	Nn += 1;
	Pn += 0;
	Sn += 1;
	
	break;
	
      case 'L':
	
	Cn += 6;
	Hn += 11;
	On += 1;
	Nn += 1;
	Pn += 0;
	Sn += 0;
	
	break;
	
      case 'I':
	
	Cn += 6;
	Hn += 11;
	On += 1;
	Nn += 1;
	Pn += 0;
	Sn += 0;
	
	break;
	
      case 'N':
	
	Cn += 4;
	Hn += 6;
	On += 2;
	Nn += 2;
	Pn += 0;
	Sn += 0;
	
	break;
	
      case 'D':
	
	Cn += 4;
	Hn += 5;
	On += 3;
	Nn += 1;
	Pn += 0;
	Sn += 0;
	
	break;
	
      case 'Q':
	
	Cn += 5;
	Hn += 8;
	On += 2;
	Nn += 2;
	Pn += 0;
	Sn += 0;
	
	break;
	
      case 'K':
	
	Cn += 6;
	Hn += 12;
	On += 1;
	Nn += 2;
	Pn += 0;
	Sn += 0;
	
	break;
	
      case 'E':
	
	Cn += 5;
	Hn += 7;
	On += 3;
	Nn += 1;
	Pn += 0;
	Sn += 0;
	
	break;
	
      case 'M':
	
	Cn += 5;
	Hn += 9;
	On += 1;
	Nn += 1;
	Pn += 0;
	Sn += 1;
	
	break;
	
      case 'H':
	
	Cn += 6;
	Hn += 7;
	On += 1;
	Nn += 3;
	Pn += 0;
	Sn += 0;
	
	break;
	
      case 'F':
	
	Cn += 9;
	Hn += 9;
	On += 1;
	Nn += 1;
	Pn += 0;
	Sn += 0;
	
	break;
	
      case 'R':
	
	Cn += 6;
	Hn += 12;
	On += 1;
	Nn += 4;
	Pn += 0;
	Sn += 0;
	
	break;
	
      case 'Y':
	
	Cn += 9;
	Hn += 9;
	On += 2;
	Nn += 1;
	Pn += 0;
	Sn += 0;
	
	break;
	
      case 'W':
	
	Cn += 11;
	Hn += 10;
	On += 1;
	Nn += 2;
	Pn += 0;
	Sn += 0;
	
	break;
	
      case 'O':
	
	Cn += 5;
	Hn += 12;
	On += 2;
	Nn += 2;
	Pn += 0;
	Sn += 0;
	
	break;
	
      case 'X':
	//same as I/L
	Cn += 6;
	Hn += 11;
	On += 1;
	Nn += 1;
	Pn += 0;
	Sn += 0;
	
	break;
	
      case 'B':
	//avg N and D
	Cn += 4;
	Hn += 5.5;
	On += 2.5;
	Nn += 1.5;
	Pn += 0;
	Sn += 0;
	
	break;
	
      case 'Z':
	//avg O and E
	Cn += 5;
	Hn += 9.5;
	On += 2.5;
	Nn += 1.5;
	Pn += 0;
	Sn += 0;
	
	break;
	
      default:

	
	break;
      }
      
    }
  }

 
  std::ostringstream out;
  
  if (Cn > 0) out << "C" << (int)(Cn+0.5);
  if (Hn > 0) out << "H" << (int)(Hn+0.5);
  if (Nn > 0) out << "N" << (int)(Nn+0.5);
  if (On > 0) out << "O" << (int)(On+0.5);
  if (Pn > 0) out << "P" << (int)(Pn+0.5);
  if (Sn > 0) out << "S" << (int)(Sn+0.5);
  
  return out.str();
   
}

