#include "SearchParams.h"


SearchParams::SearchParams(const char* paramsfile) {
  paramsfile_ = new char[strlen(paramsfile)+1];
  strcpy(paramsfile_, paramsfile);
  modifications_ = new Array<Modification*>;

  sequence_constraints_ = new Array<char*>;
  enzyme_index_ = 0;
  max_num_internal_cls_ = 0;
  min_num_tol_term_ = 2;
  precursor_monoisotopic_ = False;
  fragment_monoisotopic_ = False;

  num_sequence_modifications_ = 0;
  num_terminal_modifications_ = 0;



  /*
  // standard SEQUEST modifications
  char* sequence_modification_symbols[] = {"*", "#", "@"};
  //char* sequence_modification_symbols[] = {"*", "~", "^"};
  num_sequence_modifications_ = sizeof(sequence_modification_symbols) / sizeof(char*);
  // not yet in use (SEQUEST has no variable terminal mods)
  char* terminal_modification_symbols[] = {"1", "2", "3"};
  num_terminal_modifications_ = sizeof(terminal_modification_symbols) / sizeof(char*);
  //cout << "num seq mods: " << num_sequence_modifications_ << " and term: " << num_terminal_modifications_ << endl; exit(1);

  sequence_modification_symbols_ = new char* [num_sequence_modifications_];
  for(int k = 0; k < num_sequence_modifications_; k++) {
    sequence_modification_symbols_[k] = new char[strlen(sequence_modification_symbols[k])+1];
    strcpy(sequence_modification_symbols_[k], sequence_modification_symbols[k]);
  }
  terminal_modification_symbols_ = new char* [num_terminal_modifications_];
  for(int k = 0; k < num_terminal_modifications_; k++) {
    terminal_modification_symbols_[k] = new char[strlen(terminal_modification_symbols[k])+1];
    strcpy(terminal_modification_symbols_[k], terminal_modification_symbols[k]);
  }
  */


  //num_sequence_modifications_ = sizeof(sequence_modification_symbols_) / sizeof(char*);
  sequence_modification_index_ = 0;

  //terminal_modification_symbols_ = {",", ";", ":"};
  //num_terminal_modifications_ = sizeof(terminal_modification_symbols_) / sizeof(char*);
  terminal_modification_index_ = 0;
}

SearchParams::~SearchParams() {
  if(paramsfile_ != NULL)
    delete paramsfile_;
  if(modifications_ != NULL) {
    for(int k = 0; k < modifications_->length(); k++)
      if((*modifications_)[k] != NULL) 
	delete (*modifications_)[k];
    delete modifications_;
  }
  if(sequence_constraints_ != NULL) {
    for(int k = 0; k < sequence_constraints_->length(); k++)
      if((*sequence_constraints_)[k] != NULL)
	delete (*sequence_constraints_)[k];
  }
}


void SearchParams::setMassType(Boolean monoisotopic) {
   precursor_monoisotopic_ = monoisotopic;
   fragment_monoisotopic_ = monoisotopic;
}

char* SearchParams::getStandardModifiedPeptide(const char* peptide) {
  ModificationInfo* mod_info = new ModificationInfo(peptide, modifications_);
  if(mod_info != NULL) {
    // get the stripped pep
    char* stripped = new char[strlen(peptide)+1];
    int index = 0;
    for(int k = 0; peptide[k]; k++)
      if(peptide[k] >= 'A' && peptide[k] <= 'Z')
	stripped[index++] = peptide[k];
    stripped[index] = 0;

    char* output = mod_info->getStandardModifiedPeptide(stripped, NULL, MOD_ERROR);
    delete mod_info;
    return output;
  }
  return NULL;

}

Array<Tag*>* SearchParams::getModificationInfoTags(const char* peptide) {
  ModificationInfo* mod_info = new ModificationInfo(peptide, modifications_);
  if(mod_info != NULL) {
    Array<Tag*>* output = mod_info->getModificationInfoInfoTags(NULL);
    delete mod_info;
    return output;
  }
  return NULL;
}


void SearchParams::setModificationSymbols(int num_sequence_mods, const char** sequence_modification_symbols, int num_terminal_mods, const char** terminal_modification_symbols) {
  num_sequence_modifications_ = num_sequence_mods; //sizeof(sequence_modification_symbols) / sizeof(char*);
  //  cout << "num seq: " << num_sequence_modifications_ << endl;
  num_terminal_modifications_ = num_terminal_mods; //sizeof(terminal_modification_symbols) / sizeof(char*);
  //  cout << "num term: " << num_terminal_modifications_ << endl;
  sequence_modification_symbols_ = new char*[num_sequence_modifications_];
  int k;
  for(k = 0; k < num_sequence_modifications_; k++) {
    sequence_modification_symbols_[k] = new char[strlen(sequence_modification_symbols[k])+1];
    strcpy(sequence_modification_symbols_[k], sequence_modification_symbols[k]);
  }
  terminal_modification_symbols_ = new char*[num_terminal_modifications_];
  for(k = 0; k < num_terminal_modifications_; k++) {
    terminal_modification_symbols_[k] = new char[strlen(terminal_modification_symbols[k])+1];
    strcpy(terminal_modification_symbols_[k], terminal_modification_symbols[k]);
  }

}


void SearchParams::replace(char* peptide, const char* orig, const char* replacements, Boolean termini) {
  if(strlen(orig) != strlen(replacements)) {
    cout << "error: " << orig << " has different number of symbols than " << replacements << endl;
    exit(1);
  }
  for(int k = 0; peptide[k]; k++) {
    if(! termini || k == 0 || !peptide[k+1]) {
      for(int j = 0; orig[j]; j++) {
    	if(peptide[k] == orig[j]) {
	      if (!(peptide[k] = replacements[j])) {
	    	return;
     	  }
    	}
      }
    }
  }
}

void SearchParams::updateModifications() {
  int j,k;
  for(k = 0; k < num_sequence_modifications_; k++)
    modifyAminoacidModifications(sequence_modification_symbols_[k]);
  for(k = 0; k < num_terminal_modifications_; k++)
    modifyTerminalModifications(terminal_modification_symbols_[k]);
  // now check to make sure each is unique
  for(k = 0; k < num_sequence_modifications_; k++) {
    for(j = 0; j < num_sequence_modifications_; j++)
      if(j != k && ! strcmp(sequence_modification_symbols_[k], sequence_modification_symbols_[j])) {
	cout << "error: symbol " << sequence_modification_symbols_[k] << " used more than once" << endl;
	exit(1);
      }
    for(j = 0; j < num_terminal_modifications_; j++)
      if(j != k && ! strcmp(sequence_modification_symbols_[k], terminal_modification_symbols_[j])) {
	cout << "error: symbol " << sequence_modification_symbols_[k] << " used more than once" << endl;
	exit(1);
      }
  } // next sequence mod
  for(k = 0; k < num_terminal_modifications_; k++) {
    for(j = 0; j < num_terminal_modifications_; j++)
      if(j != k && ! strcmp(terminal_modification_symbols_[k], terminal_modification_symbols_[j])) {
	cout << "error: symbol " << terminal_modification_symbols_[k] << " used more than once" << endl;
	exit(1);
      }
    for(j = 0; j < num_sequence_modifications_; j++)
      if(j != k && ! strcmp(terminal_modification_symbols_[k], sequence_modification_symbols_[j])) {
	cout << "error: symbol " << terminal_modification_symbols_[k] << " used more than once" << endl;
	exit(1);
      }
  } // next sequence mod

}

void SearchParams::setModificationSymbol(Modification* mod, Boolean advance) {
  if(mod->terminal) {
    if(terminal_modification_index_ >= num_terminal_modifications_) {
      cout << "error: exceeded maximum " << num_terminal_modifications_ << " terminal modifications" << endl;
      exit(1);
    }
    strcpy(mod->symbol, terminal_modification_symbols_[terminal_modification_index_]);
    if(advance)
      terminal_modification_index_++;
    return;
  }
  if(sequence_modification_index_ >= num_sequence_modifications_) {
    cout << "error: exceeded maximum " << num_sequence_modifications_ << " sequence modifications" << endl;
    exit(1);
  }
  strcpy(mod->symbol, sequence_modification_symbols_[sequence_modification_index_]);
  if(advance)
    sequence_modification_index_++;
}


Tag* SearchParams::getModificationTag(Modification* mod) {
  char text[100];
  
  Tag* output = NULL;
  if(mod->terminal) {
    output = new Tag("terminal_modification", True, True);
    sprintf(text, "%c", mod->aa);
    output->setAttributeValue("terminus", text);
  }
  else {
    output = new Tag("aminoacid_modification", True, True);
    sprintf(text, "%c", mod->aa);
    output->setAttributeValue("aminoacid", text);
  }
  
  //Tag* output = new Tag("search_modification", True, True);
  //if(mod->aa == 'n')
  //  output->setAttributeValue("aminoacid", "<");
  //else if(mod->aa == 'c')
  //  output->setAttributeValue("aminoacid", ">");
  //else {


  //sprintf(text, "%c", mod->aa);
  //output->setAttributeValue("aminoacid", text);
    //}
  sprintf(text, "%0.4f", mod->massdiff);
  output->setAttributeValue("massdiff", text);
  sprintf(text, "%0.4f", mod->mass);
  output->setAttributeValue("mass", text);


  if(mod->variable) {
    output->setAttributeValue("variable", "Y");
    //if(mod->aa == 'n')
    // output->setAttributeValue("symbol", "+");
    //else if(mod->aa == 'c')
    //  output->setAttributeValue("symbol", "-");
    //else 
    //sprintf(text, "%s", mod->symbol);
    //output->setAttributeValue("symbol", text);

    //#ifndef WRITE_MOD_INFO
    output->setAttributeValue("symbol", mod->symbol);
    //#endif
  }
  else
    output->setAttributeValue("variable", "N");

  if(mod->terminal) {
    if(mod->protein_terminus)
      output->setAttributeValue("protein_terminus", "Y");
    else
      output->setAttributeValue("protein_terminus", "N");
  }

  //  output->write(cout);
  return output;
}

// have a function which modifies mass diffs to be relative to naked guy (for sequest only....?)


Array<Tag*>* SearchParams::getModificationTags() {
  setModificationMasses();
  Array<Tag*>* output = new Array<Tag*>;
  // do aminoacid before terminal
  int k;
  for(k = 0; k < modifications_->length(); k++)
    if(! (*modifications_)[k]->terminal)
      output->insertAtEnd(getModificationTag((*modifications_)[k]));


  for(k = 0; k < modifications_->length(); k++)
    if((*modifications_)[k]->terminal)
    output->insertAtEnd(getModificationTag((*modifications_)[k]));
  return output;
}

char SearchParams::getModifiedAA(int index) {
  if(index < 0 || index >= modifications_->length()) {
    cout << "error in getModifiedAA with index " << index << endl;
    exit(1);
  }
  return (*modifications_)[index]->aa;
}

int SearchParams::getNumModifiedAAs() {
  return modifications_->length();
}


Tag* SearchParams::getEnzymeConstraintTag() {
  char text[100];
  if(enzyme_index_) {
    Tag* output = new Tag("enzymatic_search_constraint", True, True);
    output->setAttributeValue("enzyme", getParamsEnzyme());
    sprintf(text, "%d", max_num_internal_cls_);
    output->setAttributeValue("max_num_internal_cleavages", text);
    //    if(min_num_tol_term_ == 8)
    //      output->setAttributeValue("min_number_termini", "1-n"); // sequest values
    //   else if(min_num_tol_term_ == 9)
    //     output->setAttributeValue("min_number_termini", "1-c"); // sequest values
    //    else {
      sprintf(text, "%d", min_num_tol_term_);
      output->setAttributeValue("min_number_termini", text);  // code 8 for 1-n, 9 for 1-c
      //   }
    return output;
  }
  return NULL;
}

Tag* SearchParams::getSequenceConstraintTag(const char* seq) {
  Tag* output = new Tag("sequence_search_constraint", True, True);
  output->setAttributeValue("sequence", seq);
  return output;
}


Array<Tag*>* SearchParams::getSequenceConstraintTags() {
  Array<Tag*>* output = new Array<Tag*>;
  Tag* next = NULL;
  for(int k = 0; k < sequence_constraints_->length(); k++) 
    output->insertAtEnd(getSequenceConstraintTag((*sequence_constraints_)[k]));

  return output;
}



Array<Tag*>* SearchParams::getSearchParamTags(const char* basename, const char* engine, const char* database) {
  return getSearchParamTagsPlus(basename, engine, precursor_monoisotopic_, fragment_monoisotopic_, database);
}


Array<Tag*>* SearchParams::getSearchParamTagsPlus(const char* basename, const char* engine, Boolean parent_mono, Boolean fragment_mono, const char* database) {
  Array<Tag*>* output = new Array<Tag*>;
  Array<Tag*>* next_tags = NULL;
  Tag* next = new Tag("search_summary", True, False);
  next->setAttributeValue("base_name", basename);
  next->setAttributeValue("search_engine", engine);
  if(parent_mono)
    next->setAttributeValue("precursor_mass_type", "monoisotopic");
  else
    next->setAttributeValue("precursor_mass_type", "average");
  if(fragment_mono)
    next->setAttributeValue("fragment_mass_type", "monoisotopic");
  else
    next->setAttributeValue("fragment_mass_type", "average");
  next->setAttributeValue("out_data_type", "out");
  next->setAttributeValue("out_data", ".tgz");
  next->setAttributeValue("search_id", "1");
  output->insertAtEnd(next);

  
  next = new Tag("search_database", True, True);
  next->setAttributeValue("local_path", database);
  next->setAttributeValue("type", "AA");
  output->insertAtEnd(next);

  next = getEnzymeConstraintTag();
  if(next != NULL)
    output->insertAtEnd(next);
  next_tags = getSequenceConstraintTags();
  int k;
  for(k = 0; k < next_tags->length(); k++)
    if((*next_tags)[k] != NULL) {
      output->insertAtEnd((*next_tags)[k]);
    }
  delete next_tags;
  next_tags = getModificationTags();
  for(k = 0; k < next_tags->length(); k++)
    if((*next_tags)[k] != NULL) {
      output->insertAtEnd((*next_tags)[k]);
    }
  delete next_tags;
  next_tags = getParameterTags();
  for(k = 0; k < next_tags->length(); k++)
    if((*next_tags)[k] != NULL) {
      output->insertAtEnd((*next_tags)[k]);
    }
  delete next_tags;
  next = new Tag("search_summary", False, True);
  output->insertAtEnd(next);
  return output;

}


const char* SearchParams::getParamsEnzyme() {
  return getEnzyme(enzyme_index_);
}


double SearchParams::getMonoisotopicAAMass(char aa) {
  switch(aa) {
  case '<' : return 1.0078250321;
  case 'n' : return 1.0078250321;
  case '>' : return 17.0027396542;
  case 'c' : return 17.0027396542;
  case 'G' : return 57.0214637236;
  case 'A' : return 71.0371137878;
  case 'V' : return 99.0684139162;
  case 'L' : return 113.0840639804;
  case 'I' : return 113.0840639804;
  case 'S' : return 87.0320284099;
  case 'C' : return 103.0091844778;
  case 'T' : return 101.0476784741;
  case 'M' : return 131.0404846062;
  case 'P' : return 97.0527638520;
  case 'F' : return 147.0684139162;
  case 'Y' : return 163.0633285383;
  case 'W' : return 186.0793129535;
  case 'H' : return 137.0589118624;
  case 'K' : return 128.0949630177;
  case 'R' : return 156.1011110281;
  case 'D' : return 115.0269430320;
  case 'E' : return 129.0425930962;
  case 'N' : return 114.0429274472;
  case 'Q' : return 128.0585775114;
  case 'O' : return 114.07931;
  case 'Z' : return 128.55059;
  case 'B' : return 114.53494;
  default: return 0.0;
  } // switch
}

double SearchParams::getAverageAAMass(char aa) {
  switch(aa) {
  case '<' : return 1.0078250321;
  case 'n' : return 1.0078250321;
  case '>' : return 17.0027396542;
  case 'c' : return 17.0027396542;
  case 'G' : return 57.0519;
  case 'A' : return 71.0788;
  case 'V' : return 99.1326;
  case 'L' : return 113.1594;
  case 'I' : return 113.1594;
  case 'S' : return 87.0782;
  case 'C' : return 103.1388;
  case 'T' : return 101.1051;
  case 'M' : return 131.1926;
  case 'P' : return 97.1167;
  case 'F' : return 147.1766;
  case 'Y' : return 163.1760;
  case 'W' : return 186.2132;
  case 'H' : return 137.1411;
  case 'K' : return 128.1741;
  case 'R' : return 156.1875;
  case 'D' : return 115.0886;
  case 'E' : return 129.1155;
  case 'N' : return 114.1038;
  case 'Q' : return 128.1307;
  case 'O' : return 114.1472;
  case 'Z' : return 128.6231;
  case 'B' : return 114.5962;
  default: return 0.0;
  } // switch
}


