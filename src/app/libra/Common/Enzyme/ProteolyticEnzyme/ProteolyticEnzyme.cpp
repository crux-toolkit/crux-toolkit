
#include "ProteolyticEnzyme.h"
#include "gzstream.h"

  char ProteolyticEnzyme::CONSTANT_NAME_NONSPECIFIC_[] = "nonspecific";
  char ProteolyticEnzyme::CONSTANT_FIDELITY_SEMISPECIFIC_[] = "semispecific";
  char ProteolyticEnzyme::CONSTANT_FIDELITY_SPECIFIC_[] = "specific";

  char *ProteolyticEnzyme::default_fidelity_ = ProteolyticEnzyme::CONSTANT_FIDELITY_SPECIFIC_;
  char ProteolyticEnzyme::default_description_[] = "";
  int ProteolyticEnzyme::default_min_spacing_ = 1;
  Boolean ProteolyticEnzyme::default_independent_ = True;
  char ProteolyticEnzyme::default_nocut_[] = "";

Boolean ProteolyticEnzyme::changeDefaultMinSpacing(int default_min_spacing) {
  if (default_min_spacing < 0) return False;
  default_min_spacing_ = default_min_spacing;
  return True;
}

ProteolyticEnzyme::ProteolyticEnzyme(const char* name, const char* fidelity, Boolean indep, const char* description) {
  specified_ = False;
  cuts_ = new Array<char*>;
  no_cuts_ = new Array<char*>;
  n_sens_ = new Array<Boolean>;
  min_spacing_ = new Array<int>;


  name_ = new char[strlen(name)+1];
  strcpy(name_, name);
  if(fidelity == NULL) {
    fidelity_ = new char[strlen(default_fidelity_)+1];
    strcpy(fidelity_, default_fidelity_);
  }
  else {
    fidelity_ = new char[strlen(fidelity)+1];
    strcpy(fidelity_, fidelity);
  }
  independent_ = indep;
  if(description == NULL) {
    description_ = new char[strlen(default_description_)+1];
    strcpy(description_, default_description_);
  }
  else {
    description_ = new char[strlen(description)+1];
    strcpy(description_, description);
  }

}


ProteolyticEnzyme::ProteolyticEnzyme(Tag* tag) {
  specified_ = False;
  cuts_ = new Array<char*>;
  no_cuts_ = new Array<char*>;
  n_sens_ = new Array<Boolean>;
  min_spacing_ = new Array<int>;

  if(strcmp(tag->getName(), "sample_enzyme")) {
    cout << "error: did not pass <sample_enzyme> tag to ProteolyticEnzyme" << endl;
    exit(1);
  }

  const char* next = tag->getAttributeValue("name");
  name_ = new char[strlen(next)+1];
  strcpy(name_, next);

  cerr << "enzyme name is:" << name_ << endl;

  next = tag->getAttributeValue("fidelity");
  if(next == NULL) {
    fidelity_ = new char[strlen(default_fidelity_)+1];
    strcpy(fidelity_, default_fidelity_);
  }
  else {
    fidelity_ = new char[strlen(next)+1];
    strcpy(fidelity_, next);
  }

  next = tag->getAttributeValue("independent");
  if(next == NULL) {
    independent_ = default_independent_;
  }
  else {
    independent_ = ! strcmp(tag->getAttributeValue("independent"), "1");
  }

  next = tag->getAttributeValue("description");
  if(next == NULL) {
    description_ = new char[strlen(default_description_)+1];
    strcpy(description_, default_description_);
  }
  else {
    description_ = new char[strlen(next)+1];
    strcpy(description_, next);
  }

}

ProteolyticEnzyme::~ProteolyticEnzyme() {
  if(cuts_ != NULL) {
    for(int k = 0; k < cuts_->length(); k++)
      if((*cuts_)[k] != NULL)
	delete[] (*cuts_)[k];
    delete cuts_;
  }
  if(no_cuts_ != NULL) {
    for(int k = 0; k < no_cuts_->length(); k++)
      if((*no_cuts_)[k] != NULL)
	delete[] (*no_cuts_)[k];
    delete no_cuts_;
  }
  if(name_ != NULL)
    delete[] name_;
  if(description_ != NULL)
    delete[] description_;
  if(fidelity_ != NULL)
    delete[] fidelity_;
  if(min_spacing_ != NULL)
    delete min_spacing_;
  if(n_sens_ != NULL)
    delete n_sens_;
}


void ProteolyticEnzyme::fixSpecificity() {
  specified_ = True;
}

void ProteolyticEnzyme::enterSpecificity(const char* cut, const char* no_cut, Boolean n_sens, int min_spacing) {
  if(specified_)
    return;
  /*
  cout << "cut: " << cut << ", nocut: " << no_cut;
  if(n_sens)
    cout << " nsens";
  else 
    cout << " csens";
  cout << " min spacing: " << min_spacing << endl;
  */

  char* next = new char[strlen(cut)+1];
  strcpy(next, cut);
  cuts_->insertAtEnd(next);
  if(no_cut == NULL) {
      next = new char[strlen(default_nocut_)+1];
      strcpy(next, default_nocut_);
  }
  else {
    next = new char[strlen(no_cut)+1];
    strcpy(next, no_cut);
  }
  no_cuts_->insertAtEnd(next);

  n_sens_->insertAtEnd(n_sens);
  min_spacing_->insertAtEnd(min_spacing);


}

char* ProteolyticEnzyme::getStandardEnzymeName(const char* input_name) {
  char* name = NULL;
  if(input_name != NULL) {
    name = new char[strlen(input_name)+1];
    strcpy(name, input_name);
    for(int k = 0;input_name[k]; k++)
      name[k] = (char)(tolower((int)name[k]));

    if (strstr(name, "tryptic")) {
      strcpy(strstr(name, "tryptic"), "trypsin"); 
    }
  
  }
  return name;
}

void ProteolyticEnzyme::parseEnzymeName(const char* input_name, char **name, Boolean &isSemiSpecific) {
  isSemiSpecific = False;
  if(input_name != NULL) {
    char *lname = ProteolyticEnzyme::getStandardEnzymeName(input_name);
    
    const char *KEYWORD_semi = "semi";
    const char *KEYWORD_semi_hyphen = "semi-";
    isSemiSpecific = (strstr(lname, KEYWORD_semi)==lname || strstr(lname, KEYWORD_semi_hyphen)==lname);
    if (isSemiSpecific) {
      const char *usedKeyword = (strstr(lname, KEYWORD_semi_hyphen)==lname) ? KEYWORD_semi_hyphen : KEYWORD_semi;
      *name = new char[strlen(lname)-strlen(usedKeyword)+1];
      strcpy(*name, lname+strlen(usedKeyword));
      delete[] lname;
    }
    else {
      *name = lname;
    }
  }
}

void ProteolyticEnzyme::enterSpecificity(Tag* tag) {
  if(specified_)
    return;
  if(strcmp(tag->getName(), "specificity")) {
    tag->write(cout);
    cout << "ProteolyticEnzyme::enterSpecificity error" << endl;
    return;
  }
  //  cout << "entering...." << endl;
  //  tag->write(cout);

  const char* nextval = tag->getAttributeValue("cut");

  char* next = new char[strlen(nextval)+1];
  strcpy(next, nextval);
  cuts_->insertAtEnd(next);

  nextval = tag->getAttributeValue("no_cut");
  if(nextval == NULL) {
      next = new char[strlen(default_nocut_)+1];
      strcpy(next, default_nocut_);
  }
  else {
    next = new char[strlen(nextval)+1];
    strcpy(next, nextval);
  }
  no_cuts_->insertAtEnd(next);

  nextval = tag->getAttributeValue("sense");
  n_sens_->insertAtEnd(! strcmp(nextval, "N"));


  nextval = tag->getAttributeValue("min_spacing");
  if(nextval == NULL) 
    min_spacing_->insertAtEnd(default_min_spacing_);
  else
    min_spacing_->insertAtEnd(atoi(nextval));
}

void ProteolyticEnzyme::writeTraditionalPepXMLTags(FILE* fp) {
  Array<Tag*>* tags = getPepXMLTags();
  for(int k = 0; k < tags->length(); k++) {
    if((*tags)[k] != NULL) {
      (*tags)[k]->writeTraditional(fp);
      delete (*tags)[k];
    }
  }
  delete tags;
}
void ProteolyticEnzyme::writeTraditionalPepXMLTags(ogzstream* fp) {
  Array<Tag*>* tags = getPepXMLTags();
  for(int k = 0; k < tags->length(); k++) {
    if((*tags)[k] != NULL) {
      (*tags)[k]->writeTraditional(fp);
      delete (*tags)[k];
    }
  }
  delete tags;
}


void ProteolyticEnzyme::write(ostream& os) {
  os << "\t   " << name_;
  if (0==cuts_->length()) os << endl;
  for(int k = 0; k < cuts_->length(); k++) {
    if (k>0) {
      os << "\t   ";
      for (unsigned int i=0; i<strlen(name_); i++) os << " ";
    }
    os << " : cut(" << (*cuts_)[k] << ")";
    if(strcmp((*no_cuts_)[k], default_nocut_))
      os << " nocuts(" << (*no_cuts_)[k] << ")";
    os << " sense(" << (((*n_sens_)[k]) ? "N" : "C") << ")";
    os << endl;
  }
}


void ProteolyticEnzyme::writePepXMLTags(ostream& os) {
  Array<Tag*>* tags = getPepXMLTags();
  for(int k = 0; k < tags->length(); k++) {
    if((*tags)[k] != NULL) {
      (*tags)[k]->write(os);
      delete (*tags)[k];
    }
  }
  delete tags;
}


Array<Tag*>* ProteolyticEnzyme::getPepXMLTags() {
  Array<Tag*>* output = new Array<Tag*>;
  Tag* next = NULL;
  char text[10];
  if(cuts_->length() == 0) 
    next = new Tag("sample_enzyme", True, True);
  else
    next = new Tag("sample_enzyme", True, False);

  next->setAttributeValue("name", name_);

  if(strcmp(fidelity_, default_fidelity_))
    next->setAttributeValue("fidelity", fidelity_);
  if(strcmp(description_, default_description_))
    next->setAttributeValue("description", description_);
  if(independent_ != default_independent_) { 
    if(independent_) {
      next->setAttributeValue("independent", "1");
    }
    else { 
      next->setAttributeValue("independent", "0");
    }
  }
  output->insertAtEnd(next);

  for(int k = 0; k < cuts_->length(); k++) {
    next = new Tag("specificity", True, True);
    next->setAttributeValue("cut", (*cuts_)[k]);
    if(strcmp((*no_cuts_)[k], default_nocut_))
      next->setAttributeValue("no_cut", (*no_cuts_)[k]);
    if((*n_sens_)[k])
      next->setAttributeValue("sense", "N");
    else
      next->setAttributeValue("sense", "C");

    if((*min_spacing_)[k] != default_min_spacing_) {
      sprintf(text, "%d", (*min_spacing_)[k]);
      next->setAttributeValue("min_spacing", text);
    }
    output->insertAtEnd(next);
  }

  if(cuts_->length() != 0) {
    next = new Tag("sample_enzyme", False, True);
    output->insertAtEnd(next);
  }
  return output;
}


int ProteolyticEnzyme::getNumTolTerm(char prev, const char* pep, char next) {
  if(cuts_->length() == 0 || ! strcmp(fidelity_, "nonspecific") || pep == NULL || strlen(pep) < 1)
    return 2;

  int ntt = 0;
  int nterm = 0;
  int cterm = 0;
  if(! strcmp(fidelity_, ProteolyticEnzyme::CONSTANT_FIDELITY_SEMISPECIFIC_))
    ntt++;

  for(int k = 0; k < cuts_->length(); k++) {
    if(independent_) {
      nterm = 0;
      cterm = 0;
    }
    if(nterm == 0) {
       if(prev == '-')
	       nterm = 1;
       else if((*n_sens_)[k] && strchr((*cuts_)[k], pep[0]) != NULL &&
	       strchr((*no_cuts_)[k], prev) == NULL)
	       nterm = 1;
       else if(! (*n_sens_)[k] && strchr((*cuts_)[k], prev) != NULL &&
	       strchr((*no_cuts_)[k], pep[0]) == NULL)
	       nterm = 1;
    } // nterm 0
    if(cterm == 0) {
      if(next == '-')
	      cterm = 1;
      else if((*n_sens_)[k] && strchr((*cuts_)[k], next) != NULL &&
	      strchr((*no_cuts_)[k], pep[strlen(pep)-1]) == NULL)
	      cterm = 1;
      else if(! (*n_sens_)[k] && strchr((*cuts_)[k], pep[strlen(pep)-1]) != NULL &&
	      strchr((*no_cuts_)[k], next) == NULL)
	      cterm = 1;
    } // cterm 0
    if(nterm + cterm > ntt)
      ntt = nterm + cterm;


  } // next spec

  return ntt;


}

// specify modified positions (optionally) where cleavage cannot occur
int ProteolyticEnzyme::getNumMissedCleavages(const char* pep, ModificationInfo* modinfo) {
  int nmc = 0;
  if(cuts_->length() == 0 || pep == NULL || !*pep)
    return nmc;
  
  int lastcut = -1;
  int firstPepAA = -1;
  int lastPepAA = (int)strlen(pep)-1;
  int k;
  for(k = 0;pep[k]; k++) {
    if (pep[k] == '.' && firstPepAA == -1) {
      firstPepAA = k+1;
    }
    else if (pep[k] == '.') {
      lastPepAA = k-1;
    }
  }
  //Check if prev a.a. was tolerable in which case don't count the first amino acid in the peptide as 
  //missed cleavage in case min_spacing is greater than zero
  if (firstPepAA >= 2) {
    for(k = 0; k < cuts_->length(); k++) {
      if (!(*n_sens_)[k]) {
	//Cuts on c-term e.g. Trypsin K.xxxx
	if (strchr((*cuts_)[k], pep[firstPepAA-2]) != NULL) {
	  lastcut = firstPepAA-1;
	}
      }
      else {
	//Cuts on n-term e.g. LysN  x.Kxxx
	if (strchr((*cuts_)[k], pep[firstPepAA]) != NULL) {
	  lastcut = firstPepAA;
	}
      }
    }
  }
  if (firstPepAA == -1) firstPepAA = 0;
  for(k = firstPepAA; k <= lastPepAA; k++) {
    // check to make sure not modified position (not subject to cleavage)
    Boolean unmod = True;
    if(modinfo != NULL) {
      for(int p = 0; p < modinfo->getNumModAAs(); p++)
	if(modinfo->getModAAPos(p) == k + 1) {
	  unmod = False;
	  p = modinfo->getNumModAAs();
	}
    }
    if(unmod) {

      for(int j = 0; j < cuts_->length(); j++) {
	if(k > lastcut + (*min_spacing_)[j] && 
	   k < (int) strlen(pep) - (*min_spacing_)[j] - 1) {
	  if((*n_sens_)[j] && 
	     k > 0 && 
	     strchr((*cuts_)[j], pep[k]) != NULL && 
	     strchr((*no_cuts_)[j], pep[k-1]) == NULL) {
	    //n-term cut e.g. LysN xx<k>.K in example the <k> position is not counted as missed cleavage
	    if ((*min_spacing_)[j] == 0 || 
		pep[k+1] !=  '.' || 
		strchr((*cuts_)[j], pep[k+2]) == NULL) {
	      nmc++;
	      lastcut = k;
	      j = cuts_->length();
	    }
	  }
	  else if(! (*n_sens_)[j] && 
		  k < (int) strlen(pep) - 1 && 
		  k < lastPepAA &&
		  strchr((*cuts_)[j], pep[k]) != NULL && 
		  strchr((*no_cuts_)[j], pep[k+1]) == NULL) {
	    //c-term cut e.g. Tryps x<k>K.x in example the <k> position is not counted as missed cleavage
	    if ((*min_spacing_)[j] == 0 || 
		strchr((*cuts_)[j], pep[k+1]) == NULL || 
		pep[k+2] != '.') {
	      nmc++;
	      lastcut = k;
	      j = cuts_->length();
	    }
	  }

	} // if spacing appropriate

      } // next spec

    } // if not modified position
  } // next pep pos
  return nmc;
}

char* ProteolyticEnzyme::strip(char* pep, Boolean remove_mods) {
  int start = 0;
  int stop = (int)strlen(pep)-1;
  char* output = NULL;
  if(strlen(pep) > 4 && pep[1] == '.')
    start = 2;
  if(strlen(pep) > 4 && pep[strlen(pep)-2] == '.')
    stop = (int)strlen(pep)-3;

  if(! remove_mods) {
    output = new char[stop - start + 2];
    strncpy(output, pep+start, stop - start + 1);
    output[stop - start + 1] = 0;
    return output;
  }

  int k,length = 0;
  for(k = start; k <= stop; k++) 
    if(pep[k] >= 'A' && pep[k] <= 'Z')
      length++;
  output = new char[length+1];
  length = 0;
  for(k = start; k <= stop; k++) 
    if(pep[k] >= 'A' && pep[k] <= 'Z')
      output[length++] = pep[k];
  output[length] = 0;
  return output;
}

Boolean ProteolyticEnzyme::hasSpecificity(char aa) {
  for(int k = 0; k < cuts_->length(); k++)
    if(strchr((*cuts_)[k], aa) != NULL)
      return True;
  return False;
}

Boolean ProteolyticEnzyme::hasModAASpecificity(SearchParams* params) {
  for(int k = 0; k < params->getNumModifiedAAs(); k++)
    if(hasSpecificity(params->getModifiedAA(k)))
      return True;
      
  return False;
}

const char* ProteolyticEnzyme::getName() const {
  return name_;
}

const char* ProteolyticEnzyme::getFidelity() const {
  return fidelity_;
}

Boolean ProteolyticEnzyme::isSemiSpecific()const {
  return (!strcmp(fidelity_, ProteolyticEnzyme::CONSTANT_FIDELITY_SEMISPECIFIC_));
}

int ProteolyticEnzyme::getMinNumberTermini() const {
  if (!strcmp(name_, CONSTANT_NAME_NONSPECIFIC_)) return 0;
  return (isSemiSpecific() ? 1 : 2);
}
