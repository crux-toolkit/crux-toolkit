/***************************************************************************
*
*  FILENAME      :   pi.c
*
*  DESCRIPTION   :   Given a protein sequence, these functions return its pI.
* 
*  ARGUMENTS     :   char * seq; string containing the sequence, 
*                    all upper case letters, no leading or trailing blanks.
*        
*  AUTHORS       :   ROA from Amos' BASIC procedure
*
*  VERSION       :   1.6
*  DATE          :   1/25/95
*  
*  Copyright 1993 by Melanie/UIN/HCUG. All rights reserved.
* Performance tweaks by Insilicos LLC
*
***************************************************************************/

#include "pICalculator.h"
#include <assert.h>

pICalculator::pICalculator() { 
  common_ctor_init(-1);
}
pICalculator::pICalculator(int run_idx) {
  common_ctor_init(run_idx);
}

void pICalculator::common_ctor_init(int run_idx) {
  //DDS: pI model
  //DDS: TODO Move this stuff into the pI model class???
  run_idx_ = run_idx;
  run_pI_sum_ = 0;
  run_pI_count_ = 0;
  run_pI_used_count_ = 0;
  peps_ = new Array<char*>;
  run_pI_mean_=0;
  run_pI_stddev_=0;
}

pICalculator::~pICalculator() {
  for (int i=0; i<peps_->size(); i++) {
    delete [] (char*)(*peps_)[i];
  }
  delete peps_;
}

double pICalculator::Peptide_pI(const char* seq, ModificationInfo* modinfo) {
  return Peptide_pI(seq, modinfo, -1);
}

double pICalculator::Peptide_pI(const char* seq, ModificationInfo* modinfo, double pI) {
  double rtn = pI;
  if (rtn <= 0) {
    rtn = COMPUTE_PI(seq, (int)strlen(seq), 0, modinfo);
  }  
  pIs_.insertAtEnd(rtn);
  char* new_pep = new char[strlen(seq)+1];
  strcpy(new_pep, seq);
  //new_pep[strlen(seq)]='\0';
  peps_->insertAtEnd(new_pep);
  //DDS: pI model
  run_pI_sum_ += rtn;
  run_pI_count_++;
  return rtn;
}

static inline double sqrval(double d) {
	return d*d;
}

void pICalculator::recalc_pIstats(Array<double>* probs) {
  double numer = 0;
  double denom = 0;
  double tot_pIs = 0;
  double prob;
  run_pI_used_count_ = 0;
  // cout << "DDS: probsSize=" << probs->size() << " runpINum=" << run_pI_count_<< endl;
  //assert(probs->size() == run_pI_count_);
  int i;
  for (i=0; i<run_pI_count_; i++) {
    if ((*probs)[i] >= 0) {
      run_pI_used_count_ ++;
      prob = sqrval((*probs)[i]);
      tot_pIs += pIs_[i];
      numer += prob*pIs_[i];
      denom += prob;
    }
  }

  if (denom == 0) {
    numer = tot_pIs;
    denom = run_pI_count_;
    
    run_pI_mean_ = numer / denom;
    run_pI_stddev_ = 0;

    for (i=0; i<run_pI_count_; i++) {
      run_pI_stddev_ += sqrval((run_pI_mean_ - pIs_[i]));
    }
  }
  else {
    run_pI_mean_ = numer / denom;
    run_pI_stddev_ = 0;
    
    for (i=0; i<run_pI_count_; i++) {
      if ((*probs)[i] >= 0) {
      prob = sqrval((*probs)[i]);
      run_pI_stddev_ += prob * sqrval((run_pI_mean_ - pIs_[i]));
      }
    }
  }
  run_pI_stddev_ /= denom;

  run_pI_stddev_ = sqrt(run_pI_stddev_);  

}

void pICalculator::write_pIstats(ostream& out) {
  out <<  run_idx_ << "\t" << run_pI_mean_ << "\t" << run_pI_stddev_ << "\t" << run_pI_used_count_<< endl;
}

int pICalculator::getRunIdx() {
  return run_idx_;
}

double pICalculator::getRunPIMean() {
  return run_pI_mean_;
}

double pICalculator::getRunPIStdDev() {
  return run_pI_stddev_;
}

Boolean pICalculator::recalc_pIstats(Array<double>* probs, double min_prob, Array<int>* ntts, int min_ntt) {
  double numer = 0;
  double denom1 = 0; // BSP note using multiple varbls for denom to avoid an apparent gcc 3.4.4 
                     // optimization bug causing regression test failures 
                     // (yeah, I know, I didn't believe it either but BoundsChecker,
                     // HeapAgent, and Valgrind don't find anything, and when coded
                     // this way Cygwin, MinGW and VC8 builds now give same results)
  double denom2;
  double tot_pIs = 0;
  double all_pIs = 0;
  double prob;  
  double tmp_pI_mean = 0;
  double tmp_pI_stddev = 0;
  run_pI_used_count_ = 0;
  // cout << "DDS: probsSize=" << probs->size() << " runpINum=" << run_pI_count_<< endl;
  assert(probs->size() == run_pI_count_);
  int i;
  for (i=0; i<run_pI_count_; i++) {
    if ((*probs)[i] >= 0) {
      prob = (*probs)[i];
      all_pIs += pIs_[i];
      if (prob >= min_prob && (*ntts)[i] >= min_ntt) {
	tot_pIs += pIs_[i];
	numer += prob*pIs_[i];
	denom1 += prob;
	run_pI_used_count_ ++;
      }
    }
  }

  if (denom1 == 0) {
    numer = all_pIs;
   
    tmp_pI_mean = numer / (denom2 = run_pI_count_);
    tmp_pI_stddev = 0;
   
    for (i=0; i<run_pI_count_; i++) {
      tmp_pI_stddev += sqrval(tmp_pI_mean - pIs_[i]);
    }
  }
  else {
    tmp_pI_mean = numer / denom1;
    tmp_pI_stddev = 0;
    denom2 = 0;
    for (i=0; i<run_pI_count_; i++) {
      if ((*probs)[i] >= 0) {
	prob = (*probs)[i];
	if (prob >= min_prob && (*ntts)[i] >= min_ntt) {
	  tmp_pI_stddev += prob * sqrval(tmp_pI_mean - pIs_[i]);
	  denom2 += prob;
	}
      }
    }
  }
  tmp_pI_stddev /= denom2;

  tmp_pI_stddev = sqrt(tmp_pI_stddev); 
  
  run_pI_mean_ = tmp_pI_mean;
  run_pI_stddev_ = tmp_pI_stddev;

 
  //DDS: Remove outliers

  numer = 0;
  double denom3 = 0;
  tot_pIs = 0;
  all_pIs = 0;
  run_pI_used_count_ = 0;
  // cout << "DDS: probsSize=" << probs->size() << " runpINum=" << run_pI_count_<< endl;
  assert(probs->size() == run_pI_count_);
  for (i=0; i<run_pI_count_; i++) {
    if ((*probs)[i] >= 0) {
      prob = (*probs)[i];
      all_pIs += pIs_[i];
      if (prob >= min_prob && (*ntts)[i] >= min_ntt &&
	  pIs_[i] < tmp_pI_mean + 3*tmp_pI_stddev &&
	  pIs_[i] > tmp_pI_mean - 3*tmp_pI_stddev ) {
	tot_pIs += pIs_[i];
	numer += prob*pIs_[i];
	denom3 += prob;
	//numer += pIs_[i];
	//denom += 1;
	run_pI_used_count_ ++;
      }
    }
  }

  double denom4;
  if (denom3 == 0) {
    numer = all_pIs;
    
    run_pI_mean_ = numer / (denom4 = run_pI_count_);
    run_pI_stddev_ = 0;
   
    for (i=0; i<run_pI_count_; i++) {
      run_pI_stddev_ += sqrval((run_pI_mean_ - pIs_[i]));
    }
  }
  else {
    run_pI_mean_ = numer / denom3;
    run_pI_stddev_ = 0;
    denom4 = 0;
    for (i=0; i<run_pI_count_; i++) {
      if ((*probs)[i] >= 0) {
	prob = (*probs)[i];
	if (prob >= min_prob && (*ntts)[i] >= min_ntt &&
	  pIs_[i] < tmp_pI_mean + 3*tmp_pI_stddev &&
	  pIs_[i] > tmp_pI_mean - 3*tmp_pI_stddev ) {
	  run_pI_stddev_ += prob * sqrval((run_pI_mean_ - pIs_[i]));
	  denom4 += prob;
	}
      }
    }
  }
  run_pI_stddev_ /= denom4;

  run_pI_stddev_ = sqrt(run_pI_stddev_);  
  


  return True;
}

//DDS: pI model
void pICalculator::calc_pIstats() {
  double numer = 0;
  double denom = 0;
  run_pI_used_count_ = run_pI_count_;
  int i;
  for (i=0; i<run_pI_count_; i++) {
    numer += pIs_[i];
    denom ++;
  }
  
  run_pI_mean_ = numer/denom;
  run_pI_stddev_ = 0;
 
  for (i=0; i<run_pI_count_; i++) {
    run_pI_stddev_ += sqrval((run_pI_mean_ - pIs_[i]));
  }

  run_pI_stddev_ /= denom;
  
  run_pI_stddev_ = sqrt(run_pI_stddev_);  
}

double pICalculator::getpIScore(double pI, const char* pep) {
  //TODO: DDS hack for now should divide by average stdev over all runs
  //if (run_pI_stddev_ <= 0.05) run_pI_stddev_ = 0.5;
 
  double rtn = (pI - run_pI_mean_) / run_pI_stddev_;
  if (rtn > 5) {
    rtn = 5;
  }

  if (rtn < -5) {
    rtn = -5;
  }
  //  cout << "peptide seq: " << pep << " pI: " << pI << " Score: " << rtn <<endl;
  return rtn;
  
}

double pICalculator::ModifiedPeptide_pI(const char* seq, int offset, int len) {
  // first strip it

  char* stripped = new char[len+1];
  if(stripped == NULL) {
    cout << "error with stripped" << endl;
    exit(1);
  }
  int index = 0;
  for(int k = offset; k < offset+len && k < (int) strlen(seq); k++)
    if(seq[k] >= 'A' && seq[k] <= 'Z') {
      //cout << (k+1) << ": " << seq[k] << endl;
      stripped[index++] = seq[k];
    }
  stripped[index] = 0;
  //cout << seq << ": " << "stripped: " << stripped << endl;
  double result = Peptide_pI(stripped, NULL, -1);
  delete[] stripped;
  return result;

  //  return 0.0;
}

static double exp10negpk[PK_X][PK_Y]; // cache EXP10(-pk[][]) calcs
static bool bInit=false;


// must be stripped first
// ADD IN HANDLING OF MODIFICATIONS IN FUTURE
double pICalculator::COMPUTE_PI(const char *seq,unsigned long seq_length, int charge_increment, ModificationInfo* modinfo)
{
   double          comp[26];    /* Amino acid composition of the protein */
   register int    nterm_res,   /* N-terminal residue */
                   cterm_res;   /* C-terminal residue */
   register unsigned long i;
  
   memset((char *)comp, 0, 26 * sizeof(double));
   
   if (!bInit) { // cache some expensive calcs  Insilicos 2007
   	  for (int x=PK_X;x--;) {
        for (int y=PK_Y;y--;) {
           exp10negpk[x][y] = pow(10.0,-pk[x][y]); // cache these EXP10 calcs
        }
     }
     bInit = true;
   }

   for (i = 0; i < seq_length; i++)        /* compute the amino acid composition */
   {
     // put in adjustments for modifications 
     if(modinfo != NULL && seq[i] == 'C' && modinfo->isModifiedResidue(i) && 
	ICATMixtureDistr::isICAT(modinfo->getModifiedResidueMass(i), MOD_ERROR))
       ; // icat
     else
      comp[seq[i] - 'A']+=1.0;
   }
  
   nterm_res = seq[0] - 'A';               /* Look up N-terminal residue */
   cterm_res = seq[seq_length-1] - 'A';    /* Look up C-terminal residue */
  
   double ph_min = PH_MIN;
   double ph_max = PH_MAX;
   double ph_mid;
  
   for (i = 0; i<MAXLOOP && (ph_max - ph_min)>EPSI; i++)
   {
      ph_mid = ph_min + (ph_max - ph_min) / 2.0;

      double expnegph_mid = pow(10.0,-ph_mid);
    
      double cter = exp10negpk[cterm_res][0] / (exp10negpk[cterm_res][0] + expnegph_mid);
      double nter = expnegph_mid / (exp10negpk[nterm_res][1] + expnegph_mid);
    
      double carg = comp[R] * expnegph_mid / (exp10negpk[R][2] + expnegph_mid);
      double chis = comp[H] * expnegph_mid / (exp10negpk[H][2] + expnegph_mid);
      double clys = comp[K] * expnegph_mid / (exp10negpk[K][2] + expnegph_mid);
    
      double casp = comp[D] * exp10negpk[D][2] / (exp10negpk[D][2] + expnegph_mid);
      double cglu = comp[E] * exp10negpk[E][2] / (exp10negpk[E][2] + expnegph_mid);
      double ccys = comp[C] * exp10negpk[C][2] / (exp10negpk[C][2] + expnegph_mid);
      

      double ctyr = comp[Y] * exp10negpk[Y][2] / (exp10negpk[Y][2] + expnegph_mid);
    
      double charge = carg + clys + chis + nter + charge_increment 
         - (casp + cglu + ctyr + ccys + cter);
    
      if (charge > 0.0)
      {
         ph_min = ph_mid;
      }
      else
      {
         ph_max = ph_mid;
      }
   }
  
   return(ph_mid);
}

