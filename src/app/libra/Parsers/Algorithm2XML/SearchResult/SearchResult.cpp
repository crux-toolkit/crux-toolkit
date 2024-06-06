#include "SearchResult.h"

/*

Program       : SequestResult for discr_calc of PeptideProphet 
Author        : Andrew Keller <akeller@systemsbiology.org>                                                       
Date          : 11.27.02 


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


/*
format -1 uncertain
format 0 no mw column
format 1 has mw column
*/

#include <stdlib.h>

SearchResult::SearchResult() { 
 init();
 varmod_count_ = 0;
}

void SearchResult::setRunIdx(int run_idx) {
  run_idx_ = run_idx;
}

void SearchResult::setRunName(string* run_name) {
  run_name_ = run_name;
}

SearchResult::SearchResult(Array<Tag*>* tags) {
  char result_tag[] = "spectrum_query";
  //  char result_tag[] = "search_result";
  char hit_tag[] = "search_hit";
  //char score_tag[] = "search_score";
  char pepproph_tag[] = "peptideprophet_result";
 char iproph_tag[] = "interprophet_result";
  double neutral_prec_mass = 0;
  
  init();
  proteins_ = new Array<char*>();
  exp_label_ = "";
  /*
  pI_ = -1.0;
  adjusted_prob_ = False;
  incomplete_prob_ = False;
  */
  const int num_nec_fields = 2;
  Boolean found[num_nec_fields];
  int k;
  for(k = 0; k < num_nec_fields; k++)
    found[k] = False;
 
  processed_ = True;

  //  for(int k = 0; k < tags->size(); k++) 
  //    (*tags)[k]->write(cout);
  varmod_count_ = 0;
  adduct_ = 0.;
  
#ifdef USE_STD_MODS
  mod_info_ = NULL;
  Array<Tag*>* modifications = NULL;
  Boolean mod_on = False;
#endif

  bool gotprob = false;
  
  Tag* tag;
  for(k = 0; k < tags->size(); k++) {
    tag = (*tags)[k];

    //    tag->write(cout);

    if(! strcmp(tag->getName(), result_tag) && tag->isStart()) {
      gotprob = false;
      found[0] = True;
      if (tag->getAttributeValue("experiment_label")) {
	exp_label_ = tag->getAttributeValue("experiment_label");
      }
      spectrum_ = strCopy(tag->getAttributeValue("spectrum"));
      //cout << "spectrum: " << spectrum_ << endl;
      charge_ = atoi(tag->getAttributeValue("assumed_charge"));
      //cout << "ch: " << charge_ << endl;
      neutral_prec_mass = atof(tag->getAttributeValue("precursor_neutral_mass"));
      //cout << "mass: " << neutral_prec_mass << endl;
      scan_ = (atoi(tag->getAttributeValue("start_scan")) + atoi(tag->getAttributeValue("end_scan"))) / 2;

      RT_ = tag->getAttributeValue("retention_time_sec") != NULL ? atof(tag->getAttributeValue("retention_time_sec")) : 0.;
      CV_ = tag->getAttributeValue("compensation_voltage") != NULL ? -1*atof(tag->getAttributeValue("compensation_voltage")) : 0.;

    }
    else if(! strcmp(tag->getName(), hit_tag) && tag->isEnd()) {
      found[1] = True;
    }
    else if(!found[1] && ! strcmp(tag->getName(), "alternative_protein")) {
      proteins_->insertAtEnd(strCopy(tag->getAttributeValue("protein")));
    }
    else if(! strcmp(tag->getName(), hit_tag) && tag->isStart()) {
      if(!found[1] && !strcmp(tag->getAttributeValue("hit_rank"), "1")) {
	// only want the top hit
	// DDS: Sometimes there can be multiple hits with hit_rank==1 (e.g. Myrimatch) , TODO: how to deal? Disallow in schema or filter out or use first appearing? 
	// FOR NOW: Using First Appearing hit_rank==1
	proteins_->insertAtEnd(strCopy(tag->getAttributeValue("protein")));
	protein_ =  (*proteins_)[0];
	//cout << "prot: " << protein_ << endl;
	// peptide, add on prev and following aa's
	
	// change if MOD_INFO
#ifdef USE_STD_MODS
	modified_peptide_ = strCopy(tag->getAttributeValue("peptide"));
	peptide_ = new char[strlen(tag->getAttributeValue("peptide")) + 1]; //
	strip_peptide_ = string(tag->getAttributeValue("peptide"));
	strcpy(peptide_, tag->getAttributeValue("peptide"));
	if(tag->getAttributeValue("peptide_prev_aa") != NULL && strlen(tag->getAttributeValue("peptide_prev_aa")) > 0) 
	  prev_aa_ = (tag->getAttributeValue("peptide_prev_aa"))[0];
	else
	  prev_aa_ = '-';
	if(tag->getAttributeValue("peptide_next_aa") != NULL && strlen(tag->getAttributeValue("peptide_next_aa")) > 0) 
	  next_aa_ = (tag->getAttributeValue("peptide_next_aa"))[0];
	else
	  next_aa_ = '-';
	//      cout << tag->getAttributeValue("peptide_prev_aa") << "." << tag->getAttributeValue("peptide") << "." << tag->getAttributeValue("peptide_next_aa") << endl;
	//      cout << prev_aa_ << "." << peptide_ << "." << next_aa_ << endl;
#else
	peptide_ = new char[strlen(tag->getAttributeValue("peptide")) + 5]; // room for prev and foll
	if(tag->getAttributeValue("peptide_prev_aa") != NULL && strlen(tag->getAttributeValue("peptide_prev_aa")) > 0) {
	  strcpy(peptide_, tag->getAttributeValue("peptide_prev_aa"));
	  strcat(peptide_, ".");
	}
	else
	  peptide_[0] = 0;
	
	strcat(peptide_, tag->getAttributeValue("peptide"));
	if(tag->getAttributeValue("peptide_next_aa") != NULL && strlen(tag->getAttributeValue("peptide_next_aa")) > 0) {
	  strcat(peptide_, ".");
	  strcat(peptide_, tag->getAttributeValue("peptide_next_aa"));
	}
#endif
	//cout << "pep: " << peptide_ << endl;
	
	if(tag->getAttributeValue("calc_pI") != NULL) {
	  pI_ = atof(tag->getAttributeValue("calc_pI"));
	}
	else {
	  pI_ = -1;
	}
	if (tag->getAttributeValue("num_tot_proteins") == NULL ) {
	  cerr << "ERROR parsing search result for: " << spectrum_ << endl;
	  exit(1);
	}
      
	
	degen_ = strcmp(tag->getAttributeValue("num_tot_proteins"), "1");
	massdiff_ = atof(tag->getAttributeValue("massdiff"));
	neutral_mass_ = neutral_prec_mass + massdiff_;
	xlink_type_ = na;
	if (tag->getAttributeValue("xlink_type") ) { 
	  if (!strcmp(tag->getAttributeValue("xlink_type"), "xl")) {
	    xlink_type_ = xl;
	  }
	  else   if (!strcmp(tag->getAttributeValue("xlink_type"), "loop")) {
	    xlink_type_ = loop;
	  }
	}
	if (strcmp(tag->getName(), "linked_peptide") && !tag->getAttributeValue("xlink_type") ) { 
	  num_matched_ions_ = tag->getAttributeValue("num_matched_ions") == NULL ? -1 : atoi(tag->getAttributeValue("num_matched_ions"));
	  tot_num_ions_ = tag->getAttributeValue("tot_num_ions") == NULL ? -1 : atoi(tag->getAttributeValue("tot_num_ions"));
	  num_tol_term_ = tag->getAttributeValue("num_tol_term") == NULL ? -1 : atoi(tag->getAttributeValue("num_tol_term"));
	  num_missed_cl_ = tag->getAttributeValue("num_missed_cleavages") == NULL ? -1 : atoi(tag->getAttributeValue("num_missed_cleavages"));
	}
	else {
	  num_matched_ions_ = 0;
	  tot_num_ions_ = 0;
	  num_tol_term_ = tag->getAttributeValue("num_tol_term") == NULL ? -1 : atoi(tag->getAttributeValue("num_tol_term"));
	  num_missed_cl_ = tag->getAttributeValue("num_missed_cleavages") == NULL ? -1 : atoi(tag->getAttributeValue("num_missed_cleavages"));

	}
	unweight_spec_entropy_ = tag->getAttributeValue("unweight_spec_entropy") == NULL ? -1 : atof(tag->getAttributeValue("unweight_spec_entropy"));
	delta_RT_loess_ = tag->getAttributeValue("delta_RT_loess") == NULL ? -1 : atof(tag->getAttributeValue("delta_RT_loess"));
      }	

      
    /*
    else if(! strcmp(tag->getName(), score_tag) && tag->isStart()) {
      if(! strcmp(tag->getAttributeValue("name"), "sequest_xcorr"))
	xcorr_ = atof(tag->getAttributeValue("value"));
      else if(! strcmp(tag->getAttributeValue("name"), "sequest_deltacn"))
	delta_ = atof(tag->getAttributeValue("value"));
      else if(! strcmp(tag->getAttributeValue("name"), "sequest_deltacnstar"))
	deltastar_ = atoi(tag->getAttributeValue("value"));
      else if(! strcmp(tag->getAttributeValue("name"), "sequest_sprank"))
	rank_ = atoi(tag->getAttributeValue("value"));
      else if(! strcmp(tag->getAttributeValue("name"), "sequest_spscore"))
	sp_score_ = atof(tag->getAttributeValue("value"));

    }
    */
    } // end if top hit
    else if(!gotprob && ! strcmp(tag->getName(), pepproph_tag) && tag->isStart()) {
      probability_ = atof(tag->getAttributeValue("probability"));
      //      cout << "probability for " << spectrum_ << ": " << probability_ << endl;
      if(tag->getAttributeValue("analysis") != NULL) {
	//	cout << "analysis: " << tag->getAttributeValue("analysis") << endl;
	adjusted_prob_ = ! strcmp(tag->getAttributeValue("analysis"), "adjusted");
	incomplete_prob_ = ! strcmp(tag->getAttributeValue("analysis"), "incomplete");
      }
    }
    
    else if(! strcmp(tag->getName(), iproph_tag) && tag->isStart()) {
      probability_ = atof(tag->getAttributeValue("probability"));
      //      cout << "probability for " << spectrum_ << ": " << probability_ << endl;
      if(tag->getAttributeValue("analysis") != NULL) {
	//	cout << "analysis: " << tag->getAttributeValue("analysis") << endl;
	adjusted_prob_ = ! strcmp(tag->getAttributeValue("analysis"), "adjusted");
	incomplete_prob_ = ! strcmp(tag->getAttributeValue("analysis"), "incomplete");
      }
      gotprob = true;
    }
#ifdef USE_STD_MODS
    else if(!found[1] && tag->isStart() && ! strcmp(tag->getName(), "mod_aminoacid_mass")) {
      	    if (tag->getAttributeValue("source") != NULL &&
		! strcmp(tag->getAttributeValue("source"), "adduct")) 
	      adduct_ = atof(tag->getAttributeValue("variable"));

	    if ( (tag->getAttributeValue("variable") != NULL &&
		  tag->getAttributeValue("param") != NULL &&
		  ! strcmp(tag->getAttributeValue("variable"), "param") ) ||
		 (tag->getAttributeValue("variable") != NULL &&
		  tag->getAttributeValue("param") == NULL) ) {
	      char aa;
	      int pos = atoi(tag->getAttributeValue("position"))-1;
	      if (pos < 0) {
		aa = 'n';
	      }
	      else if (pos >= strip_peptide_.length()) {
		aa = 'c';
	      }
	      else {
		aa = strip_peptide_[pos];
	      }
	      if ((*var_mods_hash_).find(aa)==(*var_mods_hash_).end()) {
		(*var_mods_hash_).insert(make_pair(aa, new vector<double>()));
	      }
	      (*var_mods_hash_)[aa]->push_back(atof(tag->getAttributeValue("variable")));
	    }
	    else if  (tag->getAttributeValue("mass") != NULL) {
	      double mass = atof(tag->getAttributeValue("mass"));
	      char aa;
	      int pos = atoi(tag->getAttributeValue("position"))-1;
	      if (pos < 0) {
		aa = 'n';
	      }
	      else if (pos >= strip_peptide_.length()) {
		aa = 'c';
	      }
	      else {
		aa = strip_peptide_[pos];
	      }

	      
	      if ((*var_mods_hash_).find(aa)==(*var_mods_hash_).end()) {
		(*var_mods_hash_).insert(make_pair(aa, new vector<double>()));
	      }
	      
	      mass = mass - ResidueMass::getMass(aa,true);
	      (*var_mods_hash_)[aa]->push_back(mass);
	    
	    }
	    
    }    
    else if(!found[1] && tag->isStart() && ! strcmp(tag->getName(), "modification_info")) {

      modified_peptide_ = strCopy(tag->getAttributeValue("modified_peptide"));



      //IMPORTANT: Relying on fact that variable mods will be in bracket notation
      //            and static mods will NOT be in bracket notation
      //
      // <modification_info modified_peptide="KLGCSVFLLPEDIVEVN[115]QKM">
      //  <mod_aminoacid_mass position="4" mass="160.030649" static="57.021464"/>
      //  <mod_aminoacid_mass position="17" mass="115.026943" variable="0.984016" source="param"/>
      // </modification_info>

      for (int i=0; modified_peptide_[i]!='\0'; i++) {
	if (modified_peptide_[i]=='[')
	  varmod_count_++;
      }
      
      
      //cout << "here with that!" << endl;
      modifications = new Array<Tag*>;
      //      if(modifications == NULL)
      //	cout << "NULL" << endl;
      //      else
      //	cout << "go the memory!" << endl;
      modifications->insertAtEnd(tag);
      mod_on = !tag->isEnd();
      //cout << "done with that!" << endl;
      if (!mod_on) { // tag already closed, process it now
	mod_info_ = new ModificationInfo(modifications);
	delete modifications;
      }
      if  (tag->getAttributeValue("mod_nterm_mass") != NULL) {
	double mass = atof(tag->getAttributeValue("mod_nterm_mass"));
	char aa = 'n';
	
	if ((*var_mods_hash_).find(aa)==(*var_mods_hash_).end()) {
	  (*var_mods_hash_).insert(make_pair(aa, new vector<double>()));
	}
	mass = mass - ResidueMass::getMass(aa,true);
	(*var_mods_hash_)[aa]->push_back(mass);
      }
       if  (tag->getAttributeValue("mod_cterm_mass") != NULL) {
	double mass = atof(tag->getAttributeValue("mod_cterm_mass"));
	char aa = 'c';
	
	if ((*var_mods_hash_).find(aa)==(*var_mods_hash_).end()) {
	  (*var_mods_hash_).insert(make_pair(aa, new vector<double>()));
	}
	mass = mass - ResidueMass::getMass(aa,true);
	(*var_mods_hash_)[aa]->push_back(mass);
      }
    }
    else if(!found[1] && mod_on && tag->isEnd() && ! strcmp(tag->getName(), "modification_info")) {
      modifications->insertAtEnd(tag);
      mod_on = False;
      mod_info_ = new ModificationInfo(modifications);
      delete modifications;
    }
    else if(!found[1] && mod_on) {
      modifications->insertAtEnd(tag);
    }
    else if (!strcmp(tag->getName(), "linked_peptide") && tag->getAttributeValue("peptide")) {
      strip_peptide_ = string(tag->getAttributeValue("peptide"));
    }
     
    
#endif



  } // next tag

  for(k = 0; k < num_nec_fields; k++)
    if(! found[k])
      processed_ = False;

  if (fabs(massdiff_) < 1e-12 && fabs(adduct_) > 1e-12)
    massdiff_ = adduct_;
  
  }



SearchResult::~SearchResult() {
  if(spectrum_ != NULL) {
    delete [] spectrum_;
  }
  if (proteins_ != NULL) {
    for (int i = 0; i<proteins_->length(); i++)
      delete [] (*proteins_)[i];
    delete proteins_;
  }

  if(peptide_ != NULL) {
    delete [] peptide_;
  }
  for (TPP_HASHMAP_T<char, vector<double>*>::iterator itr = (*var_mods_hash_).begin(); itr != (*var_mods_hash_).end(); itr++) {
    itr->second->clear();
    delete itr->second;
  }
  (*var_mods_hash_).clear();
  delete var_mods_hash_;
#ifdef USE_STD_MODS
  if(mod_info_ != NULL)
    delete mod_info_;

  if (modified_peptide_ != NULL) {
    delete [] modified_peptide_;
  }

#endif
}


void SearchResult::init() {
  maldi_ = False;
  run_idx_ = -1;
  processed_ = False;
  spectrum_ = NULL;
  protein_ = NULL;
  proteins_ = NULL;
  peptide_ = NULL;
  modified_peptide_ = NULL;
  probability_ = -1.0;
  pI_ = -1.0;
  run_pI_diff_ = 15;
  RT_ = 0.0;
  run_RT_diff_ = 15;
  adjusted_prob_ = False;
  incomplete_prob_ = False;

  var_mods_hash_ = new TPP_HASHMAP_T<char, vector<double>*>();
  #ifdef USE_STD_MODS
  mod_info_ = NULL;
  #endif
}

Boolean SearchResult::isProcessed() { return processed_; }



char* SearchResult::strCopy(const char* orig) {
  if (orig == NULL) 
    return NULL;
  char* output = new char[strlen(orig)+1];
  strcpy(output, orig);
  output[strlen(orig)] = 0;
  return output;
}

#ifdef USE_STD_MODS
Boolean SearchResult::setModificationEncoding(char* buffer, int buffer_len) {
  char text[100];
  if(buffer == NULL || ! buffer_len)
    return False;
  buffer[0] = 0;
  if(mod_info_ == NULL)
    return True; // nothing to report
  if(mod_info_->getNtermModMass() > 0.0) {
    sprintf(text, "n[%0.2f]", mod_info_->getNtermModMass());
    if((int) (strlen(buffer) + strlen(text)) < buffer_len)
      strcat(buffer, text);
    else return False;
  }
  if(mod_info_->getCtermModMass() > 0.0) {
    sprintf(text, "c[%0.2f]", mod_info_->getCtermModMass());
    if((int) (strlen(buffer) + strlen(text)) < buffer_len)
      strcat(buffer, text);
    else return False;
  }
  for(int k = 0; k < mod_info_->getNumModAAs(); k++) {
    sprintf(text, "%d[%0.2f]", mod_info_->getModAAPos(k), mod_info_->getModAAMass(k));
    if((int) (strlen(buffer) + strlen(text)) < buffer_len)
      strcat(buffer, text);
    else return False;
  } // next mod pos
  return True;
}
#endif

Boolean SearchResult::isMaldi(char* spec) {
  return False;
  // need to find a way in future....
}

ostream& SearchResult::print(ostream& os) {
  os << charge_ << " " << protein_ << " " << peptide_ << " " << spectrum_ << endl;
  return os;
}

char* SearchResult::stripHTML(const char* orig) {

  int iBufCt=0;
  int iNewBufCt=0;
  char* output = new char[strlen(orig)+1];
  size_t iLen = strlen(orig);
  int iSlashCount=0;
  while (iBufCt < (int)iLen)
    {
      if (orig[iBufCt]=='<')
	{
	  while (orig[iBufCt++] != '>');
	}
      else
	{
	  output[iNewBufCt]=orig[iBufCt];
          if (iNewBufCt-1 >= 0 && output[iNewBufCt]=='/' && output[iNewBufCt-1]!='.' && iSlashCount<1)
          {
            output[iNewBufCt]=' ';
            iSlashCount++;
          }

	  iBufCt++;
	  iNewBufCt++;
	}
    }
  output[iNewBufCt]='\0';
  return output;
}

char* SearchResult::extractDatabaseWithTags(const char* html, const char* start_tag, const char* end_tag) {
  char* output = NULL;
  const char* result = strstr(html, start_tag);
  if(result == NULL)
    return output;
  const char* second = strstr(result + strlen(start_tag), end_tag);
  if(second == NULL || strlen(result) == strlen(second))
    return output;
  size_t length = strlen(result + strlen(start_tag)) - strlen(second);
  output = new char[length+1];
  strncpy(output, result + strlen(start_tag), length);
  output[length] = 0;
  unCygwinify(output); // normalize path seps etc - no effect in cygwin builds
  // do we need to fix up the path at all?
  struct stat statbuf;
  if (stat(output,&statbuf)) { // didn't find it
	 std::string buf = resolve_root(output); // is it in the wwwroot area?
     if (!stat(buf.c_str(),&statbuf)) {
        // that worked
        delete[] output;
        output = new char[strlen(buf.c_str())+1];
        strcpy(output,buf.c_str());
     }
  }
  if (!stat(output,&statbuf)) {
     // use the full path
     char *dbfullpath = makeFullPath(output); 
     output = new char[strlen(dbfullpath)+1];
     strcpy(output,dbfullpath);
     free(dbfullpath);
  }
  return output;
}
