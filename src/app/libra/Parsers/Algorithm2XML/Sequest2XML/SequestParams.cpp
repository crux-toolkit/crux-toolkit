#include "SequestParams.h"

SequestParams::SequestParams(const char* paramsfile) : SearchParams(paramsfile) {
  cout << "params file here: " << paramsfile << endl;
  pep_mass_tol_ = 0.0;
  ion_series_[0] = 0; 

  // these are the default settings, and will no be recorded in xml unless different
  frag_ion_tol_ = 0.0;
  max_num_differential_AA_per_mod_ = 4;
  nucleotide_reading_frame_ = 0;
  remove_precursor_peak_ = 0;
  ion_cutoff_percentage_ = 0;
  match_peak_count_ = 0;
  match_peak_allowed_error_ = 1;
  match_peak_tolerance_ = 1.0;
  strcpy(protein_mass_filter_, "0 0");
  sequence_header_filter_[0] = 0;
  num_output_lines_ = 10;

  use_default_params_ = False; //True;
  print_dups_ = False;

  // standard SEQUEST modifications
  const char* sequence_modification_symbols[] = {"*", "#", "@", "^", "~", "$"};
  // not yet in use (SEQUEST has no variable terminal mods)
  const char* terminal_modification_symbols[] = {"]", "["};  // nterminal (after 1st aa), cterminal (after last aa)

  setModificationSymbols(sizeof(sequence_modification_symbols) / sizeof(char*), 
			 sequence_modification_symbols, sizeof(terminal_modification_symbols) / sizeof(char*), 
			 terminal_modification_symbols);

  // make substitutions of modification symbols here
  updateModifications();

  init();
}

void SequestParams::init() {

  ifstream fin(paramsfile_);
  if(! fin) {
    cerr << "SequestParams: error opening " << paramsfile_ << endl;
    exit(1);
  }
  const int line_width = 10000;
  char nextline[line_width];
  
  char mod_tag[] = "diff_search_options =";
  char term_mod_tag[] = "term_diff_search_options =";
  char term_mod_tag2[] = "variable_C_terminus =";
  char term_mod_tag3[] = "variable_N_terminus =";
  char first[100];
  char second[100];
  char third[100];
  char fourth[100];
  char fifth[100];
  char sixth[100];
  double massd[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double term_massd[] = {0.0, 0.0}; //c-term first, n-term second

  char static_mod_tag[] = "add_";

  char db_tag0[] = "database_name =";
  char db_tag1[] = "first_database_name =";

  char sequence_constr_tag[] = "partial_sequence =";
  //char sequence_constraint[25];
  //sequence_constraint[0] = 0;
  char enzyme_constr_tag[] = "enzyme_number =";
  //int enz_index = 0;
  char max_num_cl_tag[] = "max_num_internal_cleavage_sites =";

  char pep_mass_tol_tag[] = "peptide_mass_tolerance =";
  char frag_ion_tol_tag[] = "fragment_ion_tolerance =";
  char ion_series_tag[] = "ion_series =";

  char parent_type_tag[] = "mass_type_parent =";
  char fragment_type_tag[] = "mass_type_fragment =";

  int parent_type = 0;
  int fragment_type = 0;
  char* seq[5];
  for(int k = 0; k < 5; k++) {
    seq[k] = new char[100];
    seq[k][0] = '\0';
  }

  while(fin.getline(nextline, line_width)) {
    while(strlen(nextline) > 0 && (nextline[strlen(nextline)-1] < 32 || nextline[strlen(nextline)-1] > 126)) 
	  nextline[strlen(nextline)-1] = '\0';
    
    if(matchStart(nextline, mod_tag)) {
      first[0] = '\0';
      second[0] = '\0';
      third[0] = '\0';
      fourth[0] = '\0';
      fifth[0] = '\0';
      sixth[0] = '\0';
      sscanf(nextline + strlen(mod_tag), "%lf %s %lf %s %lf %s %lf %s %lf %s %lf %s", massd, first, massd+1, second, massd+2, third, massd+3, fourth, massd+4, fifth, massd+5, sixth);
      // now see which are more than 0
      if(massd[0] > 0.0 || massd[0] < 0.0) {

	for(int k = 0; first[k]; k++) {
	  Modification* next = new Modification();
	  next->terminal = False;
	  next->aa = first[k];
	  next->massdiff = massd[0];
	  next->variable = True;
	  next->protein_terminus = False;
	  setModificationSymbol(next, k == (int) strlen(first) - 1);
	  //next->symbol = '*';
	  modifications_->insertAtEnd(next);
	}
      }
      else if (first[0] != '\0') {
	sequence_modification_index_++;
      }

      if(massd[1] > 0.0 || massd[1] < 0.0) {

	for(int k = 0; second[k]; k++) {
	  Modification* next = new Modification();
	  next->terminal = False;
	  next->aa = second[k];
	  next->massdiff = massd[1];
	  next->variable = True;
	  next->protein_terminus = False;
	  setModificationSymbol(next, k == (int) strlen(second) - 1);
	  modifications_->insertAtEnd(next);
	}
      }
      else if (second[0] != '\0') {
	sequence_modification_index_++;
      }
     
      if(massd[2] > 0.0 || massd[2] < 0.0) {

	for(int k = 0; third[k]; k++) {
	  Modification* next = new Modification();
	  next->terminal = False;
	  next->aa = third[k];
	  next->massdiff = massd[2];
	  next->variable = True;
	  next->protein_terminus = False;
	  setModificationSymbol(next, k == (int) strlen(third) - 1);
	  modifications_->insertAtEnd(next);
	}
      }
      else if (third[0] != '\0') {
	sequence_modification_index_++;
      }

      if(massd[3] > 0.0 || massd[3] < 0.0) {

	for(int k = 0; fourth[k]; k++) {
	  Modification* next = new Modification();
	  next->terminal = False;
	  next->aa = fourth[k];
	  next->massdiff = massd[3];
	  next->variable = True;
	  next->protein_terminus = False;
	  setModificationSymbol(next, k == (int) strlen(fourth) - 1);
	  modifications_->insertAtEnd(next);
	}
      }
      else if (fourth[0] != '\0') {
	sequence_modification_index_++;
      }

      if(massd[4] > 0.0 || massd[4] < 0.0) {

	for(int k = 0; fifth[k]; k++) {
	  Modification* next = new Modification();
	  next->terminal = False;
	  next->aa = fifth[k];
	  next->massdiff = massd[4];
	  next->variable = True;
	  next->protein_terminus = False;
	  setModificationSymbol(next, k == (int) strlen(fifth) - 1);
	  modifications_->insertAtEnd(next);
	}
      }
      else if (fifth[0] != '\0') {
	sequence_modification_index_++;
      }      

      if(massd[5] > 0.0 || massd[5] < 0.0) {

	for(int k = 0; sixth[k]; k++) {
	  Modification* next = new Modification();
	  next->terminal = False;
	  next->aa = sixth[k];
	  next->massdiff = massd[5];
	  next->variable = True;
	  next->protein_terminus = False;
	  setModificationSymbol(next, k == (int) strlen(sixth) - 1);
	  modifications_->insertAtEnd(next);
	}
      }
      else if (sixth[0] != '\0') {
	sequence_modification_index_++;
      }
    } 
    else if (matchStart(nextline, term_mod_tag)) {
            sscanf(nextline + strlen(term_mod_tag), "%lf %lf", term_massd, term_massd+1);
      // now see which are more than 0
      if(term_massd[0] > 0.0 || term_massd[0] < 0.0) {
	Modification* next = new Modification();
	next->terminal = True;
	next->aa = 'c';
	next->massdiff = term_massd[0];
	next->variable = True;
	next->protein_terminus = False;
	strcpy(next->symbol, "[");
	modifications_->insertAtEnd(next);
      }

      if(term_massd[1] > 0.0 || term_massd[1] < 0.0) {
	Modification* next = new Modification();
	next->terminal = True;
	next->aa = 'n';
	next->massdiff = term_massd[1];
	next->variable = True;
	next->protein_terminus = False;
	strcpy(next->symbol, "]");
	modifications_->insertAtEnd(next);
      }

    } // have variable modification...
    else if (matchStart(nextline, term_mod_tag2)) {
       sscanf(nextline + strlen(term_mod_tag2), "%lf", term_massd);
       if(term_massd[0] > 0.0 || term_massd[0] < 0.0) {
          Modification* next = new Modification();
          next->terminal = True;
          next->aa = 'c';
          next->massdiff = term_massd[0];
          next->variable = True;
          next->protein_terminus = False;
          strcpy(next->symbol, "[");
          modifications_->insertAtEnd(next);
       }
    }
    else if (matchStart(nextline, term_mod_tag3)) {
       sscanf(nextline + strlen(term_mod_tag3), "%lf", term_massd+1);
       if(term_massd[1] > 0.0 || term_massd[1] < 0.0) {
          Modification* next = new Modification();
          next->terminal = True;
          next->aa = 'n';
          next->massdiff = term_massd[1];
          next->variable = True;
          next->protein_terminus = False;
          strcpy(next->symbol, "]");
          modifications_->insertAtEnd(next);
       }
    }

    else if(matchStart(nextline, static_mod_tag)) { // static matches
      char* result = strstr(nextline, "=");
      if(result != NULL) {
	sscanf(result+1, "%lf", massd);
	if(massd[0] > 0.0 || massd[0] < 0.0) {
	  char* orig = strstr(nextline, static_mod_tag) + strlen(static_mod_tag);
	  if(strstr(orig, "C_terminus") != NULL) {
	    Modification* next = new Modification();
	    next->terminal = True;
	    next->aa = 'c';
	    next->massdiff = massd[0];
	    next->variable = False;
	    next->protein_terminus = False;
	    modifications_->insertAtEnd(next);
	  }
	  else if(strstr(orig, "N_terminus") != NULL) {
	    Modification* next = new Modification();
	    next->terminal = True;
	    next->aa = 'n';
	    next->massdiff = massd[0];
	    next->variable = False;
	    next->protein_terminus = False;
	    modifications_->insertAtEnd(next);
	  }
	  else if(strstr(orig, "Nterm_peptide") != NULL) {
	    Modification* next = new Modification();
	    next->terminal = True;
	    next->aa = 'n';
	    next->massdiff = massd[0];
	    next->variable = False;
	    next->protein_terminus = False;
	    modifications_->insertAtEnd(next);
	  }
	  else if(strstr(orig, "Cterm_peptide") != NULL) {
	    Modification* next = new Modification();
	    next->terminal = True;
	    next->aa = 'c';
	    next->massdiff = massd[0];
	    next->variable = False;
	    next->protein_terminus = False;
	    modifications_->insertAtEnd(next);
	  }
	  else if(strstr(orig, "Nterm_protein") != NULL) {
	    Modification* next = new Modification();
	    next->terminal = True;
	    next->aa = 'n';
	    next->massdiff = massd[0];
	    next->variable = False;
	    next->protein_terminus = True;
	    modifications_->insertAtEnd(next);
	  }
	  else if(strstr(orig, "Cterm_protein") != NULL) {
	    Modification* next = new Modification();
	    next->terminal = True;
	    next->aa = 'c';
	    next->massdiff = massd[0];
	    next->variable = False;
	    next->protein_terminus = True;
	    modifications_->insertAtEnd(next);
	  }
	  else {
	    Modification* next = new Modification();
	    next->terminal = False;
	    next->aa = orig[0];
	    next->massdiff = massd[0];
	    next->variable = False;
	    next->protein_terminus = False;
	    modifications_->insertAtEnd(next);
	  } // if have match	    

	}

      }
    } // static match
    else if(matchStart(nextline, enzyme_constr_tag)) { // variable matches
      sscanf(nextline + strlen(enzyme_constr_tag), "%d", &enzyme_index_);
    }
    else if(matchStart(nextline, sequence_constr_tag)) { // variable matches
      sscanf(nextline + strlen(sequence_constr_tag), "%s %s %s %s %s", seq[0], seq[1], seq[2], seq[3], seq[4]);
      bool failed = false;
      char * buf = new char[2];
      //Check that every character in this input is an amino acid
      for(int z = 0; z < 5; z++) {  
	buf[1] = '\0';
	for (int zz = 0; seq[z][zz]; zz++) {
	  buf[0] = toupper(seq[z][zz]);
	  if (strstr(AMINO_ACIDS, buf) == NULL) {
	    failed = true;
	    break;
	  }
	}
	if (strlen(seq[z]) > 0 && !failed) {
	  char* nextz = new char[strlen(seq[z])+1];
	  strcpy(nextz, seq[z]);
	  sequence_constraints_->insertAtEnd(nextz);
	}
      }
      delete buf;
    }
    else if(matchStart(nextline, max_num_cl_tag)) {
      sscanf(nextline + strlen(max_num_cl_tag), "%d", &max_num_internal_cls_);
    }
    else if(matchStart(nextline, pep_mass_tol_tag)) {
      sscanf(nextline + strlen(pep_mass_tol_tag), "%lf", &pep_mass_tol_);
    }
    else if(matchStart(nextline, frag_ion_tol_tag)) {
      sscanf(nextline + strlen(frag_ion_tol_tag), "%lf", &frag_ion_tol_);
    }
    else if(matchStart(nextline, ion_series_tag)) {
      strcpy(ion_series_, nextline + strlen(ion_series_tag));
    }
    else if(matchStart(nextline, parent_type_tag)) {
      sscanf(nextline + strlen(parent_type_tag), "%d", &parent_type);
      precursor_monoisotopic_ = parent_type == 1; 
    }
    else if(matchStart(nextline, fragment_type_tag)) {
      sscanf(nextline + strlen(fragment_type_tag), "%d", &fragment_type);
      fragment_monoisotopic_ = fragment_type == 1; 
    }
    // all the other params here....
    else if(matchStart(nextline, "max_num_differential_AA_per_mod =")) {
      sscanf(nextline + strlen("max_num_differential_AA_per_mod ="), "%d", &max_num_differential_AA_per_mod_);
    }
    else if(matchStart(nextline, "print_duplicate_references =")) {
      int tmp;
      sscanf(nextline + strlen("print_duplicate_references ="), "%d", &tmp);
      print_dups_ = tmp == 0 ? False : True; 
    }
    else if(matchStart(nextline, "nucleotide_reading_frame =")) {
      sscanf(nextline + strlen("nucleotide_reading_frame ="), "%d", &nucleotide_reading_frame_);
    }
    else if(matchStart(nextline, "remove_precursor_peak =")) {
      sscanf(nextline + strlen("remove_precursor_peak ="), "%d", &remove_precursor_peak_);
    }
    else if(matchStart(nextline, "ion_cutoff_percentage =")) {
      sscanf(nextline + strlen("ion_cutoff_percentage ="), "%lf", &ion_cutoff_percentage_);
    }
    else if(matchStart(nextline, "match_peak_count =")) {
      sscanf(nextline + strlen("match_peak_count ="), "%d", &match_peak_count_);
    }
    else if(matchStart(nextline, "match_peak_allowed_error =")) {
      sscanf(nextline + strlen("match_peak_allowed_error ="), "%d", &match_peak_allowed_error_);
    }
    else if(matchStart(nextline, "match_peak_tolerance =")) {
      sscanf(nextline + strlen("match_peak_tolerance ="), "%lf", &match_peak_tolerance_);
    }
    else if(matchStart(nextline, "protein_mass_filter =")) {
      sscanf(nextline + strlen("protein_mass_filter ="), "%s %s", seq[0], seq[1]);
      strcpy(protein_mass_filter_, seq[0]);
      strcat(protein_mass_filter_, " ");
      strcat(protein_mass_filter_, seq[1]);
    }
    else if(matchStart(nextline, "sequence_header_filter =")) {
      sscanf(nextline + strlen("sequence_header_filter ="), "%s", sequence_header_filter_);
      //      cout << "header: " << sequence_header_filter_ << " with len: " << strlen(sequence_header_filter_) << endl;
    }
    else if(matchStart(nextline, "num_output_lines =")) {
      sscanf(nextline + strlen("num_output_lines ="), "%d", &num_output_lines_);
    }
    else if(matchStart(nextline, "NumEnzymeTermini =")) {
      sscanf(nextline + strlen("NumEnzymeTermini ="), "%d", &min_num_tol_term_);
    }
    else if(matchStart(nextline, "num_enzyme_termini =")) {
      sscanf(nextline + strlen("num_enzyme_termini ="), "%d", &min_num_tol_term_);
    }
    else if(matchStart(nextline, db_tag0)) {
      sscanf(nextline + strlen(db_tag0), "%s", database_);
    }
    else if(matchStart(nextline, db_tag1)) {
      sscanf(nextline + strlen(db_tag1), "%s", database_);
    }

  } // next line

  fin.close();

}

// DEFINE HERE WHAT AMINOACID MODIFICATION SYMBOL SUBSTITUTIONS TO MAKE
void SequestParams::modifyAminoacidModifications(char* peptide) {
  //  replace(peptide, "#", "~", 0);
  //replace(peptide, "*#@", "123", 0);
}
// DEFINE HERE WHAT TERMINAL MODIFICATION SYMBOL SUBSTITUTIONS TO MAKE
void SequestParams::modifyTerminalModifications(char* peptide) {
  replace(peptide, "", "", 1);
}

Array<Tag*>* SequestParams::getSearchParamTags(const char* basename, const char* engine) {
  return SearchParams::getSearchParamTags( basename, engine, this->database_);
}

Array<Tag*>* SequestParams::getSearchParamTags(const char* basename, const char* engine, const char* database) {
  return SearchParams::getSearchParamTags( basename, engine, database);
}

const char* SequestParams::getEnzyme(int index) {
  switch(index) {
  case 0: return "Nonspecific";
  case 1: return "Trypsin";
  case 2: return "Chymotrypsin";
  case 3: return "Clostripain";
  case 4: return "CNBr";
  case 5: return "Iodosobenzoate";
  case 6: return "Proline_endopept";
  case 7: return "Staph_protease";
  case 8: return "Trypsin_K";
  case 9: return "Trypsin_R";
  case 10: return "AspN";
  case 11: return "chymotrypic/modified";
  case 12: return "Elastase";
  case 13: return "Elastase/Trysin/Chymotrypsin";
  default: return "Nonspecific";
  }
}
/*
0.  No_Enzyme              0      -           -
1.  Trypsin                1      KR          P
2.  Chymotrypsin           1      FWY         P
3.  Clostripain            1      R           -
4.  Cyanogen_Bromide       1      M           -
5.  IodosoBenzoate         1      W           -
6.  Proline_Endopept       1      P           -
7.  Staph_Protease         1      E           -
8.  Trypsin_K              1      K           P
9.  Trypsin_R              1      R           P
10. AspN                   0      D           -
11. Cymotryp/Modified      1      FWYL        P
12. Elastase               1      ALIV        P
13. Elastase/Tryp/Chymo    1      ALIVKRWFY   P
*/


Array<Tag*>* SequestParams::getParameterTags() {
  Array<Tag*>* output = new Array<Tag*>;
  char text[100];
  Tag* next = new Tag("parameter", True, True);
  next->setAttributeValue("name", "peptide_mass_tol");
  sprintf(text, "%0.3f", pep_mass_tol_);
  next->setAttributeValue("value", text);
  output->insertAtEnd(next);
  if(! use_default_params_ || frag_ion_tol_ > 0.0) {
    next = new Tag("parameter", True, True);
    next->setAttributeValue("name", "fragment_ion_tol");
    sprintf(text, "%0.3f", frag_ion_tol_);
    next->setAttributeValue("value", text);
    output->insertAtEnd(next);
  }
  next = new Tag("parameter", True, True);
  next->setAttributeValue("name", "ion_series");
  next->setAttributeValue("value", ion_series_);
  output->insertAtEnd(next);

  if(! use_default_params_ || max_num_differential_AA_per_mod_ != 4) {
    next = new Tag("parameter", True, True);
    next->setAttributeValue("name", "max_num_differential_AA_per_mod");
    sprintf(text, "%d", max_num_differential_AA_per_mod_);
    next->setAttributeValue("value", text);
    output->insertAtEnd(next);
  }
  if(! use_default_params_ || nucleotide_reading_frame_ > 0) {
    next = new Tag("parameter", True, True);
    next->setAttributeValue("name", "nucleotide_reading_frame");
    sprintf(text, "%d", nucleotide_reading_frame_);
    next->setAttributeValue("value", text);
    output->insertAtEnd(next);
  }
  if(! use_default_params_ || num_output_lines_ != 10) {
    next = new Tag("parameter", True, True);
    next->setAttributeValue("name", "num_output_lines");
    sprintf(text, "%d", num_output_lines_);
    next->setAttributeValue("value", text);
    output->insertAtEnd(next);
  }
  if(! use_default_params_ || remove_precursor_peak_ > 0) {
    next = new Tag("parameter", True, True);
    next->setAttributeValue("name", "remove_precursor_peak");
    sprintf(text, "%d", remove_precursor_peak_);
    next->setAttributeValue("value", text);
    output->insertAtEnd(next);
  }
  if(! use_default_params_ || ion_cutoff_percentage_ > 0.0) {
    next = new Tag("parameter", True, True);
    next->setAttributeValue("name", "ion_cutoff_percentage");
    sprintf(text, "%0.1f", ion_cutoff_percentage_);
    next->setAttributeValue("value", text);
    output->insertAtEnd(next);
  }
  if(! use_default_params_ || match_peak_count_ > 0) {
    next = new Tag("parameter", True, True);
    next->setAttributeValue("name", "match_peak_count");
    sprintf(text, "%d", match_peak_count_);
    next->setAttributeValue("value", text);
    output->insertAtEnd(next);
  }
  if(! use_default_params_ || match_peak_allowed_error_ != 1) {
    next = new Tag("parameter", True, True);
    next->setAttributeValue("name", "match_peak_allowed_error");
    sprintf(text, "%d", match_peak_allowed_error_);
    next->setAttributeValue("value", text);
    output->insertAtEnd(next);
  }
  if(! use_default_params_ || match_peak_tolerance_ != 1.0) {
    next = new Tag("parameter", True, True);
    next->setAttributeValue("name", "match_peak_tolerance");
    sprintf(text, "%0.1f", match_peak_tolerance_);
    next->setAttributeValue("value", text);
    output->insertAtEnd(next);
  }
  if(! use_default_params_ || strcmp(protein_mass_filter_, "0 0")) {
    next = new Tag("parameter", True, True);
    next->setAttributeValue("name", "protein_mass_filter");
    next->setAttributeValue("value", protein_mass_filter_);
    output->insertAtEnd(next);
  }
  if(! use_default_params_ || strlen(sequence_header_filter_) > 0) {
    next = new Tag("parameter", True, True);
    next->setAttributeValue("name", "sequence_header_filter");
    next->setAttributeValue("value", sequence_header_filter_);
    output->insertAtEnd(next);
  }

  return output;
}

Boolean SequestParams::matchStart(const char* line, const char* tag) {
  const char* result = strstr(line, tag);
  return result != NULL && strlen(result) == strlen(line);
}

void SequestParams::writeParams(FILE* paramFile) {
  fprintf(paramFile, "[SEQUEST CombineOut]\n");
  fprintf(paramFile, "database_name = %s\n", database_);
  fprintf(paramFile, "peptide_mass_tolerance = %f\n", pep_mass_tol_);
  fprintf(paramFile, "email_address = \n");
  fprintf(paramFile, "create_output_files = 1\n");
  fprintf(paramFile, "create_output_files = 1\n");
  fprintf(paramFile, "ion_series = %s\n", ion_series_);
  fprintf(paramFile, "fragment_ion_tolerance = %f\n", frag_ion_tol_);
  fprintf(paramFile, "num_output_lines = %d\n", num_output_lines_);
  fprintf(paramFile, "num_description_lines = 5\n");
  fprintf(paramFile, "show_fragment_ions = 0\n");
  fprintf(paramFile, "print_duplicate_references = 0\n");
  fprintf(paramFile, "enzyme_number = 0\n");
  fprintf(paramFile, "max_num_differential_AA_per_mod = %d\n", max_num_differential_AA_per_mod_);
  fprintf(paramFile, "diff_search_options =");
  
  double nterm_diff = 0;
  double cterm_diff = 0;
  int k;
 
  for(k = 0; k < modifications_->length(); k++) {
    if((*modifications_)[k]->variable && !(*modifications_)[k]->terminal) {
      fprintf(paramFile, " %f %c", (*modifications_)[k]->massdiff, (*modifications_)[k]->aa);
    }
    else if ((*modifications_)[k]->variable) {
      if ((*modifications_)[k]->aa == 'n') {
	nterm_diff = (*modifications_)[k]->massdiff;	
      }
      else if ((*modifications_)[k]->aa == 'c') {
	cterm_diff = (*modifications_)[k]->massdiff;	
      }
    }
  }
  while (k+1 < 3) { //min number of mods required for sequest
    fprintf(paramFile, " 0.0 X");
    k++;
  }
  fprintf(paramFile, "\n");
  fprintf(paramFile, "term_diff_search_options = %f %f\n", cterm_diff, nterm_diff);
  fprintf(paramFile, "mass_type_parent = %d\n",  precursor_monoisotopic_);
  fprintf(paramFile, "mass_type_fragment = %d\n",  fragment_monoisotopic_);
  fprintf(paramFile, "remove_precursor_peak = %d\n", remove_precursor_peak_);
  fprintf(paramFile, "ion_cutoff_percentage = %f\n", ion_cutoff_percentage_);
  fprintf(paramFile, "max_num_internal_cleavage_sites = %d\n", max_num_internal_cls_);
  fprintf(paramFile, "protein_mass_filter = %s\n",  protein_mass_filter_);
  fprintf(paramFile, "match_peak_count = %d\n", match_peak_count_);
  fprintf(paramFile, "match_peak_allowed_error = %d\n", match_peak_allowed_error_);
  fprintf(paramFile, "match_peak_tolerance = %f\n", match_peak_tolerance_);
  fprintf(paramFile, "create_output_files = 1\n");
  fprintf(paramFile, "partial_sequence = 1\n");
  fprintf(paramFile, "sequence_header_filter = %s\n",  sequence_header_filter_);
  
  for(k = 0; k < modifications_->length(); k++) {
    if(!(*modifications_)[k]->variable && !(*modifications_)[k]->terminal) {
      switch ((*modifications_)[k]->aa)
	 {
	 case 'A':
	   fprintf(paramFile, "add_A_Alanine = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'S':
	   fprintf(paramFile, "add_S_Serine = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'P':
	   fprintf(paramFile, "add_P_Proline = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'V':
	   fprintf(paramFile, "add_V_Valine = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'T':
	   fprintf(paramFile, "add_T_Threonine = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'C':
	   fprintf(paramFile, "add_C_Cysteine = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'L':
	   fprintf(paramFile, "add_L_Leucine = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'I':
	   fprintf(paramFile, "add_I_Isoleucine = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'X':
	   fprintf(paramFile, "add_X_LorI = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'N':
	   fprintf(paramFile, "add_N_Asparagine = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'O':
	   fprintf(paramFile, "add_O_Ornithine = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'B':
	   fprintf(paramFile, "add_B_avg_NandD = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'D':
	   fprintf(paramFile, "add_D_Aspartic_Acid = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'Q':
	   fprintf(paramFile, "add_Q_Glutamine = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'K':
	   fprintf(paramFile, "add_K_Lysine = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'Z':
	   fprintf(paramFile, "add_Z_avg_QandE = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'E':
	   fprintf(paramFile, "add_E_Glutamic_Acid = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'M':
	   fprintf(paramFile, "add_M_Methionine = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'H':
	   fprintf(paramFile, "add_H_Histidine = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'F':
	   fprintf(paramFile, "add_F_Phenylalanine = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'R':
	   fprintf(paramFile, "add_R_Arginine = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'Y':
	   fprintf(paramFile, "add_Y_Tyrosine = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 case 'W':
	   fprintf(paramFile, "add_W_Tryptophan = %f\n", (*modifications_)[k]->massdiff);
	   break;
	 default:
	   break;
      }
    }
    else if (!(*modifications_)[k]->variable) {
      if ((*modifications_)[k]->aa == 'n') {
	if ((*modifications_)[k]->protein_terminus) {
	  fprintf(paramFile, "add_Nterm_protein = %f\n", (*modifications_)[k]->massdiff);
	}
	else {
	  fprintf(paramFile, "add_Nterm_peptide = %f\n", (*modifications_)[k]->massdiff); 
	}
      }
      else if ((*modifications_)[k]->aa == 'c') {
	if ((*modifications_)[k]->protein_terminus) {
	  fprintf(paramFile, "add_Cterm_protein = %f\n", (*modifications_)[k]->massdiff);
	}
	else {
	  fprintf(paramFile, "add_Cterm_peptide = %f\n", (*modifications_)[k]->massdiff); 
	}
      }
    }
  }

  fprintf(paramFile, "[SEQUEST_ENZYME_INFO]\n");
  fprintf(paramFile, "0.	No_Enzyme				0	-		-\n");
  fprintf(paramFile, "1.	Trypsin				1	KR		-\n");
  fprintf(paramFile, "2.	Trypsin(KRLNH)				1	KRLNH		-\n");
  fprintf(paramFile, "3.	Chymotrypsin				1	FWYL		-\n");
  fprintf(paramFile, "4.	Chymotrypsin(FWY)				1	FWY		P\n");
  fprintf(paramFile, "5.	Clostripain				1	R		-\n");
  fprintf(paramFile, "6.	Cyanogen_Bromide				1	M		-\n");
  fprintf(paramFile, "7.	IodosoBenzoate				1	W		-\n");
  fprintf(paramFile, "8.	Proline_Endopept				1	P		-\n");
  fprintf(paramFile, "9.	Staph_Protease				1	E		-\n");
  fprintf(paramFile, "10.	Trypsin_K				1	K		P\n");
  fprintf(paramFile, "11.	Trypsin_R				1	R		P\n");
  fprintf(paramFile, "12.	GluC				1	ED		-\n");
  fprintf(paramFile, "13.	LysC				1	K		-\n");
  fprintf(paramFile, "14.	AspN				0	D		-\n");
  fprintf(paramFile, "15.	Elastase				1	ALIV		P\n");
  fprintf(paramFile, "16.	Elastase/Tryp/Chymo				1	ALIVKRWFY		P\n");
  
}

void SequestParams::setModificationMasses() {
  // modify masses to be with respect to unmodified
  // first modify static by adding unmodified
  // then modify variable by adding (modified) static
  Boolean mono = precursor_monoisotopic_ && fragment_monoisotopic_;
  int k;
  for(k = 0; k < modifications_->length(); k++){ 
    if(! (*modifications_)[k]->variable){
      if(mono){
	(*modifications_)[k]->mass = (*modifications_)[k]->massdiff + getMonoisotopicAAMass((*modifications_)[k]->aa);
      }else{
	(*modifications_)[k]->mass = (*modifications_)[k]->massdiff + getAverageAAMass((*modifications_)[k]->aa);
      }
    }
  }

  //now the variables (with respect to statics)
  for(k = 0; k < modifications_->length(); k++) 
    if((*modifications_)[k]->variable) {
      double base = 0.0;
      if(mono)
	base += getMonoisotopicAAMass((*modifications_)[k]->aa);
      else
	base += getAverageAAMass((*modifications_)[k]->aa);
      // now see if must use static baseline instead
      for(int j = 0; j < modifications_->length(); j++)
	if(! (*modifications_)[j]->variable && (*modifications_)[j]->aa == (*modifications_)[k]->aa) {
	  base = (*modifications_)[j]->mass;
	  j = modifications_->length();
	}

      (*modifications_)[k]->mass = (*modifications_)[k]->massdiff + base;


    } // if variable
  /*
  for(int k = 0; k < modifications_->length(); k++) { 
    cout << (*modifications_)[k]->aa;
    if((*modifications_)[k]->variable)
      cout << (*modifications_)[k]->symbol;
    cout << " : " << (*modifications_)[k]->mass << endl;
  }
  */
}

