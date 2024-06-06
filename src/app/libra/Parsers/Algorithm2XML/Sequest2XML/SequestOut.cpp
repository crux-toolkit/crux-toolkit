#include "SequestOut.h"

SequestOut::SequestOut() {
  // Initialize Data Structure
  szFileName[0] = '\0';
  szBaseFileName[0] = '\0';
  szMod[0] = '\0';
  szDatabase[0] = '\0';
  dAMass = 0.0;
  dMass = 0.0;
  dMass1 = 0.0;
  dMass2 = 0.0;
  dMass3 = 0.0;
  dMass4 = 0.0;
  dMass5 = 0.0;
  dMass6 = 0.0;
  dMassCT = 0.0;
  dMassNT = 0.0;
  iMassType = 0;
  bNucDb = 0;
}

SequestOut::SequestOut(const SequestOut& out) {
  // Initialize Data Structure
  strcpy(szFileName, out.szFileName);
  strcpy(szBaseFileName, out.szBaseFileName);
  strcpy(szMod, out.szMod);
  strcpy(szDatabase, out.szDatabase);
  dAMass = out.dAMass;
  dMass = out.dMass;
  dMass1 = out.dMass1;
  dMass2 = out.dMass2;
  dMass3 = out.dMass3;
  dMass4 = out.dMass4;
  dMass5 = out.dMass5;
  dMass6 = out.dMass6;
  dMassCT = out.dMassCT;
  dMassNT = out.dMassNT;
  iMassType = out.iMassType;;
  bNucDb = out.bNucDb;

}

void SequestOut::writeOutFile(FILE* outFile, SequestParams* params) {
  fprintf(outFile, "\n%s\n", szBaseFileName);
  fprintf(outFile, "TPP Sequest Out File Writer %s (c) 2006\n", szTPPVersionInfo);
  fprintf(outFile, "Seattle Proteome Center, Institute for Systems Biology, D.Shteynberg\n");
  fprintf(outFile, "Licensed under LGPL.\n");
  //TODO: Some of this bookeeping can be improved
  fprintf(outFile, "Running on %s \n", getenv("HOST"));
  fprintf(outFile, "(M+H)+ mass = %f ~ %f (+%c), fragment tol = %f, ", 
	  dAMass, params->pep_mass_tol_, szBaseFileName[strlen(szBaseFileName)-1] , params->frag_ion_tol_);
  if (params->precursor_monoisotopic_) {
    fprintf(outFile, "MONO/");
  }
  else {
    fprintf(outFile, "AVG/");
  }

  if (params->fragment_monoisotopic_) {
    fprintf(outFile, "MONO\n");
  }
  else {
    fprintf(outFile, "AVG\n");
  }

  fprintf(outFile, "total inten = ?, lowest Sp = 0.0, # matched peptides = ?\n");
  fprintf(outFile, "# amino acids = ?, # proteins = ?, %s\n", szDatabase);
  fprintf(outFile, "ion series nABY ABCDVWXYZ: %s\n", params->ion_series_);
  fprintf(outFile, "display top 20/0, ion %% = %f, CODE = 000000\n", params->ion_cutoff_percentage_);
  //Print modifications and Enzyme
  for (int k=0; k < params->modifications_->size(); k++) {
    if ((*(params->modifications_))[k]->variable && !(*(params->modifications_))[k]->terminal) {
      if ((*(params->modifications_))[k]->massdiff > 0.0) {
	fprintf(outFile, " (%c%s +%f) ", (*(params->modifications_))[k]->aa, 
		(*(params->modifications_))[k]->symbol, (*(params->modifications_))[k]->massdiff);    
      }
      else {
	fprintf(outFile, " (%c%s %f) ", (*(params->modifications_))[k]->aa, 
		(*(params->modifications_))[k]->symbol, (*(params->modifications_))[k]->massdiff);    
      }
    }
    else if ((*(params->modifications_))[k]->variable && (*(params->modifications_))[k]->terminal) {
      if ((*(params->modifications_))[k]->aa == 'c') {
	fprintf(outFile, " (ct%s %f) ", (*(params->modifications_))[k]->symbol, 
		(*(params->modifications_))[k]->massdiff);    
      }
      else if ((*(params->modifications_))[k]->aa == 'n') {
	fprintf(outFile, " (nt%s %f) ", (*(params->modifications_))[k]->symbol, 
		(*(params->modifications_))[k]->massdiff);    
      }
    }    
    else if (!(*(params->modifications_))[k]->variable && !(*(params->modifications_))[k]->terminal) {
      fprintf(outFile, " %c=%f ", (*(params->modifications_))[k]->aa, (*(params->modifications_))[k]->mass);    
    }    
    else if (!(*(params->modifications_))[k]->variable && (*(params->modifications_))[k]->terminal) {
      if ((*(params->modifications_))[k]->aa == 'c') {
	if ((*(params->modifications_))[k]->protein_terminus) {
	  fprintf(outFile, " +Cterm-prot=%f ", (*(params->modifications_))[k]->mass);    
	}
	else {
	  fprintf(outFile, " +Cterm-pep=%f ", (*(params->modifications_))[k]->mass);    
	}
      }
      else if ((*(params->modifications_))[k]->aa == 'n') {
	if ((*(params->modifications_))[k]->protein_terminus) {
	  fprintf(outFile, " +Nterm-prot=%f ", (*(params->modifications_))[k]->mass);    
	}
	else {
	  fprintf(outFile, " +Nterm-pep=%f ", (*(params->modifications_))[k]->mass);    
	}
      }
    }
  }
  fprintf(outFile, " Enzyme:%s\n\n", params->getEnzyme(params->enzyme_index_));
  fprintf(outFile, "  #   Rank/Sp    (M+H)+   deltCn   XCorr    Sp     Ions  Reference         Peptide\n");
  fprintf(outFile, " ---  -------  ---------  ------  ------   ----    ----  ---------         -------\n");
  for (int i=0; i<sequestHits_.size(); i++) {
    int idx = i + sequestHits_[i]->iDeltCnIdxDiff;
    if (idx < 0) {
      sequestHits_[i]->writeOutFile(outFile, 0.0);
    }
    else {
      sequestHits_[i]->writeOutFile(outFile, sequestHits_[idx]->dDeltCn);
    }
  }

}

void SequestOut::insertNextHit(SequestHit* hit) {
  sequestHits_.insertAtEnd(hit);
}

int SequestOut::getNumHits() {
  return  sequestHits_.size();
}

SequestHit* SequestOut::getHitByIndex(int idx) {
  if (idx < 0 || idx >= sequestHits_.size()) {
    return NULL;
  }
  else {
    return sequestHits_[idx];
  }
}

void SequestOut::getDeltaCn() {
  for (int i = 0; i < sequestHits_.size()-1; i++) {
    int minLen = (int)strlen(sequestHits_[i]->szPlainPep); 
    int pepLen = minLen;
    int noDeltaCnYet = 1;
    double dFirstDeltCn = sequestHits_[i+1]->dDeltCn - sequestHits_[i]->dDeltCn;

    if (dFirstDeltCn <= 0.0)
       dFirstDeltCn = 0.001;

    for (int j = i + 1; j < sequestHits_.size(); j++) {
      if ((int)strlen(sequestHits_[j]->szPlainPep) < minLen)
	minLen = (int)strlen(sequestHits_[j]->szPlainPep);
      int diffs = 0;
      for (int k = 0; k < minLen; k++) {
	/* K/Q and I/L don't count as differences */
	if (sequestHits_[i]->szPlainPep[k] != sequestHits_[j]->szPlainPep[k]) {
	    if (!((sequestHits_[i]->szPlainPep[k] == 'K' || sequestHits_[i]->szPlainPep[k] == 'Q')
		 && (sequestHits_[j]->szPlainPep[k] == 'K' || sequestHits_[j]->szPlainPep[k] == 'Q')) &&
		!((sequestHits_[i]->szPlainPep[k] == 'I' || sequestHits_[i]->szPlainPep[k] == 'L')
		  && (sequestHits_[j]->szPlainPep[k] == 'I' || sequestHits_[j]->szPlainPep[k] == 'L'))) {
		diffs++;
	      }
	  }
      }


      /*
       * Calculate deltCn only if sequences are less than
       * PERCENTAGE similar;  PERCENTAGE=0.75 for no good reason
       */
      if ((double) ((double) (pepLen - 3.0 - diffs) /
		    (double) (pepLen - 3.0)) < PERCENTAGE) {

	  sequestHits_[i]->dDeltCn = sequestHits_[j]->dDeltCn;
	  sequestHits_[j]->iDeltCnIdxDiff = i-j;
	  noDeltaCnYet = 0;
	  if (j - i > 1)
	    sequestHits_[i]->dSpecialDeltCn = dFirstDeltCn;
	  else
	    sequestHits_[i]->dSpecialDeltCn = 0;
	  break;
      }
    }
    if (noDeltaCnYet == 1)  {
      // Special dCn because there wasn't a second ranked score
      sequestHits_[i]->dSpecialDeltCn = sequestHits_[i]->dDeltCn;
      sequestHits_[i]->dDeltCn = 1.0;

    }
  }
}
