/*
Program       : ModificationInfo
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: ModificationInfo.cpp 7582 2017-05-08 19:50:49Z real_procopio $

Copyright (C) 2003 Andrew Keller

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

Andrew Keller
Institute for Systems Biology
401 Terry Avenue North
Seattle, WA  98109  USA
akeller@systemsbiology.org
*/

#include "ModificationInfo.h"


ModificationInfo::ModificationInfo(double nterm, double cterm) {
  mod_nterm_mass_ = nterm;
  mod_cterm_mass_ = cterm;
  mod_aa_positions_ = new Array<int>;
  mod_aa_masses_ = new Array<double>;
  mod_peptide_ = NULL;
}

ModificationInfo::ModificationInfo(const char* embedded_pep, Array<Modification*>* mods)
{
  mod_nterm_mass_ = 0.0;
  mod_cterm_mass_ = 0.0;
  mod_aa_positions_ = new Array<int>;
  mod_aa_masses_ = new Array<double>;
  mod_peptide_ = NULL;

  mod_peptide_ = new char[strlen(embedded_pep)+1];
  strcpy(mod_peptide_, embedded_pep);

  double nextmass;

  // now start going through peptide
  int counter = 0; // for stripped pep pos
  int k = 0;

  // first deal with n terminus

  // legal: M*]HLSTR[
  int j;
  for(j = 0; j < mods->length(); j++) {
    // make sure symbol is present somewhere (located after first amino acid and its modification
    if((*mods)[j]->aa == 'n' && strstr(embedded_pep, (*mods)[j]->symbol) != NULL && (*mods)[j]->symbol[0] ) {
      mod_nterm_mass_ = (*mods)[j]->mass;
      j = mods->length();
      //k++; // advance past this one
    } // if
  } // next mod

  if(mod_nterm_mass_ == 0.0) { // check for static
    for(int j = 0; j < mods->length(); j++) {
      if((*mods)[j]->aa == 'n' && strlen((*mods)[j]->symbol) == 0) { // static
        mod_nterm_mass_ = (*mods)[j]->mass;
        j = mods->length();
      } // if
    }
  } // static


  int elen = (int)strlen(embedded_pep);
  for(; k < elen; k++) {

    if(embedded_pep[k] >= 'A' && embedded_pep[k] <= 'Z')
      counter++;

    nextmass = 0.0;
    for(int j = 0; j < mods->length(); j++) {
      if(! (*mods)[j]->terminal && (*mods)[j]->aa == embedded_pep[k]
          && k < elen - 1 && (*mods)[j]->symbol[0]
          && embedded_pep[k+1] == (*mods)[j]->symbol[0])
      {
          nextmass = (*mods)[j]->mass;
          j = mods->length();
          k++; //advance past symbol
      } // if
    } // next mod

    if(nextmass == 0.0) { // check for static
      for(int j = 0; j < mods->length(); j++) {
        if((*mods)[j]->aa == embedded_pep[k] && strlen((*mods)[j]->symbol) == 0) { // static
            nextmass = (*mods)[j]->mass;
            j = mods->length();
        } // if
      }
    } // static

    // check for normal aa at this position
    if(nextmass > 0.0) {
      enterAAMod(counter, nextmass);
    }
  } // next pos

  // finally deal with c terminus
  for(j = 0; j < mods->length(); j++) {
    if((*mods)[j]->aa == 'c' && embedded_pep[strlen(embedded_pep)-1] == (*mods)[j]->symbol[0]) {
      mod_cterm_mass_ = (*mods)[j]->mass;
      j = mods->length();
      k++; // advance past this one
    } // if
  } // next mod

  if(mod_cterm_mass_ == 0.0) { // check for static
    for(int j = 0; j < mods->length(); j++) {
      if((*mods)[j]->aa == 'c' && strlen((*mods)[j]->symbol) == 0) { // static
        mod_cterm_mass_ = (*mods)[j]->mass;
        j = mods->length();
      } // if
    }
  } // static

}


ModificationInfo::ModificationInfo(Array<Tag*>* tags) {
  mod_nterm_mass_ = 0.0;
  mod_cterm_mass_ = 0.0;
  mod_aa_positions_ = new Array<int>;
  mod_aa_masses_ = new Array<double>;
  mod_peptide_ = NULL;

  for(int k = 0; k < tags->length(); k++) {
    Tag* tag = (*tags)[k];
    if(tag->isStart() && ! strcmp(tag->getName(), "modification_info")) {
      if(tag->getAttributeValue("mod_nterm_mass") != NULL)
	mod_nterm_mass_ = atof(tag->getAttributeValue("mod_nterm_mass"));
      if(tag->getAttributeValue("mod_cterm_mass") != NULL)
	mod_cterm_mass_ = atof(tag->getAttributeValue("mod_cterm_mass"));
      if(tag->getAttributeValue("modified_peptide") != NULL) {
	mod_peptide_ = new char[strlen(tag->getAttributeValue("modified_peptide"))+1];
	strcpy(mod_peptide_, tag->getAttributeValue("modified_peptide"));
      }
    }
    else if(tag->isStart() && ! strcmp(tag->getName(), "mod_aminoacid_mass"))
      enterAAMod(atoi(tag->getAttributeValue("position")), atof(tag->getAttributeValue("mass")));

  } //  next tag
  /*
  if(mod_cterm_mass_ > 0.0)
    cout << "MOD C TERM" << endl;
  if(mod_nterm_mass_ > 0.0)
    cout << "MOD N TERM" << endl;
  */
}

ModificationInfo::ModificationInfo(const ModificationInfo* modinfo) {
  mod_nterm_mass_ = 0.0;
  mod_cterm_mass_ = 0.0;
  mod_aa_positions_ = new Array<int>;
  mod_aa_masses_ = new Array<double>;
  mod_peptide_ = NULL;
  if(modinfo != NULL) {
    mod_nterm_mass_ = modinfo->getNtermModMass();
    mod_cterm_mass_ = modinfo->getCtermModMass();
    for(int k = 0; k < modinfo->getNumModAAs(); k++)
      enterAAMod(modinfo->getModAAPos(k), modinfo->getModAAMass(k));
    if(modinfo->getModifiedPeptide() != NULL) {
      mod_peptide_ = new char[strlen(modinfo->getModifiedPeptide())+1];
      strcpy(mod_peptide_, modinfo->getModifiedPeptide());
    }
  } // not null
}

ModificationInfo::~ModificationInfo() {
  if(mod_aa_positions_ != NULL)
    delete mod_aa_positions_;
  if(mod_aa_masses_ != NULL)
    delete mod_aa_masses_;
  if(mod_peptide_ != NULL)
    delete[] mod_peptide_;
}

/*
double ModificationInfo::getModifiedPeptideMass(char* peptide, Boolean monoisotopic) {
  double mass = 0.0;
  if(mod_nterm_mass_)
    mass += mod_nterm_mass_;
  else
    mass += ResidueMass::getMass('n', monoisotopic);
  if(mod_cterm_mass_)
    mass += mod_cterm_mass_;
  else
    mass += ResidueMass::getMass('c', monoisotopic);
  int index = 0;
  for(int k = 0; k < strlen(peptide); k++)
    if(index < mod_aa_positions_->length() && (*mod_aa_positions_)[index] == k + 1)
      mass += (*mod_aa_masses_)[index++];
    else
      mass += ResidueMass::getMass(peptide[k], monoisotopic);
  return mass;
}
*/
Boolean ModificationInfo::equivalentModification(const ModificationInfo* modinfo, double error, const char* peptide, const char* quant_labels) const {

  // make an array
  //DDS: original  if(modinfo == NULL || peptide == NULL)
  //TODO: Is this the right thing?"
  if(modinfo == NULL || peptide == NULL || mod_peptide_ == NULL)
    return False;
  if(quant_labels == NULL || strchr(quant_labels, 'n') == NULL)

    if(modinfo->getNtermModMass() - mod_nterm_mass_ > error ||
       mod_nterm_mass_ - modinfo->getNtermModMass() > error)
      return False;

  if(quant_labels == NULL || strchr(quant_labels, 'c') == NULL)
    if(modinfo->getCtermModMass() - mod_cterm_mass_ > error ||
       mod_cterm_mass_ - modinfo->getCtermModMass() > error)
      return False;
  
  if(modinfo->getNumModAAs() != getNumModAAs())
    return False;
  
  for(int k = 0; k < modinfo->getNumModAAs(); k++) {
    if(quant_labels == NULL || strchr(quant_labels, peptide[modinfo->getModAAPos(k)-1]) == NULL)

      if(modinfo->getModAAPos(k) != getModAAPos(k) ||
	 fabs(modinfo->getModAAMass(k) - getModAAMass(k)) > error)
	return False;

  }
  return True;
}

void ModificationInfo::enterAAMod(int pos, double mass) {
  mod_aa_positions_->insertAtEnd(pos);
  mod_aa_masses_->insertAtEnd(mass);
}

double ModificationInfo::getNtermModMass() const {
  return mod_nterm_mass_;
}

double ModificationInfo::getCtermModMass() const {
  return mod_cterm_mass_;
}

Boolean ModificationInfo::isModified() const {
  return mod_nterm_mass_ > 0.0 || mod_cterm_mass_ > 0.0 || getNumModAAs() > 0;
}

Boolean ModificationInfo::isModifiedResidue(int pos) const {
  for(int k = 0; k < mod_aa_positions_->length(); k++)
    if((*mod_aa_positions_)[k] == pos + 1)
      return True;
  return False;
}

char* ModificationInfo::getModifiedPeptide() const {
  return mod_peptide_;
}

double ModificationInfo::getModifiedResidueMass(int pos) const {
  for(int k = 0; k < mod_aa_positions_->length(); k++)
    if((*mod_aa_positions_)[k] == pos + 1)
      return (*mod_aa_masses_)[k];
  return 0.0;
}

int ModificationInfo::getNumModAAs() const {
  return mod_aa_positions_->length();
}

int ModificationInfo::getModAAPos(int index) const {
  if(index < 0 || index >= mod_aa_positions_->length()) {
    cout << "WARNING: modification info at index " << index << " cannot be found for peptide " <<  mod_peptide_ << endl;
    return -1;
  }
  return (*mod_aa_positions_)[index];
}

double ModificationInfo::getModAAMass(int index) const {
  if(index < 0 || index >= mod_aa_masses_->length()) {
    cout << "WARNING: modification info at index " << index << " cannot be found for peptide " <<  mod_peptide_ << endl;
    return -1;
  }
  return (*mod_aa_masses_)[index];
}

Array<Tag*>* ModificationInfo::getModificationInfoInfoTags(char* std_peptide_name) const {
  char text[20];
  char modification_tag[] = "modification_info";
  char mod_aminoacid_tag[] = "mod_aminoacid_mass";
  Array<Tag*>* output = new Array<Tag*>;

  if(! isModified())
    return output;

  Tag* first = new Tag(modification_tag, True, False);
  if(std_peptide_name != NULL)
    first->setAttributeValue("modified_peptide", std_peptide_name);
  if(mod_nterm_mass_ > 0.0) {
    sprintf(text, "%f", mod_nterm_mass_);
    first->setAttributeValue("mod_nterm_mass", text);
  }
  if(mod_cterm_mass_ > 0.0) {
    sprintf(text, "%f", mod_cterm_mass_);
    first->setAttributeValue("mod_cterm_mass", text);
  }
  output->insertAtEnd(first);
  for(int k = 0; k < getNumModAAs(); k++) {
    Tag* next = new Tag(mod_aminoacid_tag, True, True);
    sprintf(text, "%d", (*mod_aa_positions_)[k]);
    next->setAttributeValue("position", text);
    sprintf(text, "%f", (*mod_aa_masses_)[k]);
    next->setAttributeValue("mass", text);
    output->insertAtEnd(next);
  }
  output->insertAtEnd(new Tag(modification_tag, False, True));
  return output;
}

char* ModificationInfo::getStandardModifiedPeptide(const char* peptide, Array<StaticModificationCount>* static_consts, double error, const char* starttag, const char* endtag) const {
  return getStandardModifiedPeptide(peptide, static_consts, error, starttag, endtag, false);
}

char* ModificationInfo::getStandardModifiedPeptide(const char* peptide, Array<StaticModificationCount>* static_consts, double error, const char* starttag, const char* endtag, bool out_static) const {
  int next_mod = 0;
  int len = (int)strlen(peptide) + 1;
  // approximate length is 8 per aa mod, 9 per terminal,  increased to 20 to allow for longer mods....
  if(mod_nterm_mass_ > 0.0)
    len += (int)(20 + strlen(starttag) + strlen(endtag));
  if(mod_cterm_mass_ > 0.0)
    len += (int)(20 + strlen(starttag) + strlen(endtag));
  len += getNumModAAs() * (int)(20 + strlen(starttag) + strlen(endtag));

  char* output = new char[len];

  if(! isModified()) {
    strcpy(output, peptide);
    return output;
  }

  output[0] = 0;

  char text[20];
  int curr_len = 0;

  //  cout << "Peptide: " << peptide << endl;
  if(mod_nterm_mass_ > 0.0) {
    // check whether static const to be ignored
    //    cout << "nterm mass: " << mod_nterm_mass_ << " ";
    Boolean found = False;
    if(!out_static && static_consts != NULL) {
      for(int k = 0; k < static_consts->length(); k++) {
	//	cout << (*static_consts)[k].mod << ":" << (*static_consts)[k].mass << " ";
	if((*static_consts)[k].mod == 'n' && (*static_consts)[k].mass - mod_nterm_mass_ <= error &&
	   mod_nterm_mass_ - (*static_consts)[k].mass <= error) {
	  found = True;
	  k = static_consts->length();
	}
	//	if(found)
	//	  cout << "found! ";
      }
      //      cout << endl;
    }
    if(! found) {
      sprintf(text, "%0.0f", mod_nterm_mass_);
      strcat(output, "n");
      strcat(output, starttag);
      strcat(output, text);
      strcat(output, endtag);
      //      strcat(output, "]");
    }
  }
  for(int k = 0; peptide[k]; k++) {
    curr_len = (int)strlen(output);
    output[curr_len++] = peptide[k];
    output[curr_len] = 0;
    
    for (next_mod = 0; next_mod < getNumModAAs(); next_mod++) {
      // now check for mod
      if(next_mod < getNumModAAs() && getModAAPos(next_mod) == k + 1) {
	// check whether static const to be ignored
	Boolean found = False;
	double nextmass = getModAAMass(next_mod);
	if(!out_static && static_consts != NULL) {
	  for(int j = 0; j < static_consts->length(); j++) {
	    //cout << (*static_consts)[j].mod << "-" << (*static_consts)[j].mass << " vs " << peptide[j] << "-" << nextmass << endl;
	    if((*static_consts)[j].mod == peptide[k] && (*static_consts)[j].mass - nextmass <= error &&
	       nextmass - (*static_consts)[j].mass <= error) {
	      found = True;
	      j = static_consts->length();
	    }
	  }
	}
	if(! found) {
	  sprintf(text, "%0.0f", nextmass);
	  //	strcat(output, "[");
	  strcat(output, starttag);
	  strcat(output, text);
	  strcat(output, endtag);
	  //	strcat(output, "]");
	}
      
	//	next_mod++;
	
      } // if
    }


  } // next aa
  if(mod_cterm_mass_ > 0.0) {

    // check whether static const to be ignored
    Boolean found = False;
    if(!out_static && static_consts != NULL) {
      for(int k = 0; k < static_consts->length(); k++)
	if((*static_consts)[k].mod == 'c' && (*static_consts)[k].mass - mod_cterm_mass_ <= error &&
	   mod_cterm_mass_ - (*static_consts)[k].mass <= error) {
	  found = True;
	  k = static_consts->length();
	}
    }
    if(! found) {
      sprintf(text, "%0.0f", mod_cterm_mass_);
      strcat(output, "c");
	strcat(output, starttag);
      strcat(output, text);
	strcat(output, endtag);
	//      strcat(output, "]");
    }
  }

  return output;
}

// prints out query string for cgi
Boolean ModificationInfo::setQueryString(char* query, int maxlength) {
  int len = 0;
  char prefix[] = "Mod";
  char mid[] = "=";
  char suffix[] = "&amp;";
  query[0] = 0;
  char text[500];
  if(getNtermModMass() > 0.0) {
    sprintf(text, "%sN%s%0.4f%s", prefix, mid, getNtermModMass(), suffix);
    if(len + (int) strlen(text) <= maxlength) {
      strcat(query, text);
      len += (int)strlen(text);
    }
    else
      return False;
  }
  if(getCtermModMass() > 0.0) {
    sprintf(text, "%sC%s%0.4f%s", prefix, mid, getCtermModMass(), suffix);
    if(len + (int) strlen(text) <= maxlength) {
      strcat(query, text);
      len += (int)strlen(text);
    }
    else
      return False;
  }
  for(int k = 0; k < getNumModAAs(); k++) {
    sprintf(text, "%s%d%s%0.4f%s", prefix, getModAAPos(k), mid, getModAAMass(k), suffix);
    if(len + (int) strlen(text) <= maxlength) {
      strcat(query, text);
      len += (int)strlen(text);
    }
    else
      return False;
  } // next mod pos
  return True;
}

char* ModificationInfo::getStandardModifiedPeptide(const char* peptide, Array<StaticModificationCount>* static_consts, double error) const {
  return getStandardModifiedPeptide(peptide, static_consts, error, "[", "]", false);
}


char* ModificationInfo::getStandardModifiedPeptide(const char* peptide, Array<StaticModificationCount>* static_consts, double error, bool out_static) const {
  return getStandardModifiedPeptide(peptide, static_consts, error, "[", "]", out_static);
}
