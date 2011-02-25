#include "xhhc.h"
#include "XLinkBondMap.h"


#include "parameter.h"
#include "peptide.h"
#include "DelimitedFile.h"

#include "ProteinPeptideIterator.h"

using namespace std;

FLOAT_T LinkedPeptide::linker_mass;


map<string, vector<PEPTIDE_T*> > sequence_peptide_map; //hack to keep track of peptides.

void get_linear_peptides(set<string>& peptides,
			 DATABASE_PROTEIN_ITERATOR_T* protein_iterator,
			 PEPTIDE_CONSTRAINT_T* peptide_constraint) {

  ProteinPeptideIterator* peptide_iterator = NULL;
  Protein* protein;
  PEPTIDE_T* peptide;

  string sequence = "";
  string last_sequence = "zz";
  //bool missed_cleavage = false;
  // keep track of whether the next peptide contains the previous one or not
  //size_t index;
  while (database_protein_iterator_has_next(protein_iterator)) {
    protein = database_protein_iterator_next(protein_iterator);
    peptide_iterator = new ProteinPeptideIterator(protein, peptide_constraint);
    // missed_cleavages must be TRUE in protein.c for this to work
    peptide_iterator->prepareMc(true);
    while (peptide_iterator->hasNext()) {
      //peptide = database_peptide_iterator_next(peptide_iterator);
      peptide = peptide_iterator->next();
      sequence = get_peptide_sequence(peptide); 
      carp(CARP_INFO,"Adding linear peptide:%s",get_peptide_sequence(peptide));
      peptides.insert(sequence);

      map<string, vector<PEPTIDE_T*> >::iterator find_iter;

      find_iter = sequence_peptide_map.find(sequence);

      carp(CARP_DEBUG,"Adding to map:%s,",sequence.c_str());

      if (find_iter == sequence_peptide_map.end()) {
        vector<PEPTIDE_T*> peptide_vector;
        peptide_vector.push_back(peptide);
        sequence_peptide_map.insert(make_pair(sequence, peptide_vector));
      } else {
        find_iter -> second.push_back(peptide);
      }
      

    }
  } 
}

vector<PEPTIDE_T*>& get_peptides_from_sequence(string& sequence) {
  return sequence_peptide_map[sequence];
}

void free_peptides() {
  map<string, vector<PEPTIDE_T*> >::iterator map_iter;

  for (map_iter = sequence_peptide_map.begin();
       map_iter != sequence_peptide_map.end();
       ++map_iter) {
    vector<PEPTIDE_T*>& peptides = map_iter -> second;
    for (unsigned int i=0;i<peptides.size();i++) {
      free_peptide(peptides[i]);
    }
  }
  sequence_peptide_map.clear();
}


// a hack, works for EDC linker only
void get_linkable_peptides(set<string>& peptides, 
	DATABASE_PROTEIN_ITERATOR_T* protein_iterator,
	PEPTIDE_CONSTRAINT_T* peptide_constraint) 
{
  ProteinPeptideIterator* peptide_iterator = NULL;
  Protein* protein;
  PEPTIDE_T* peptide;
  string sequence = "";
  string last_sequence = "zz";
  bool missed_cleavage = false;
  // keep track of whether the next peptide contains the previous one or not
  size_t index;
  while (database_protein_iterator_has_next(protein_iterator)) {
    protein = database_protein_iterator_next(protein_iterator);
    peptide_iterator = new ProteinPeptideIterator(protein, peptide_constraint);
    // missed_cleavages must be TRUE in protein.c for this to work
    peptide_iterator->prepareMc(true); 
    while (peptide_iterator->hasNext()) {
      //peptide = database_peptide_iterator_next(peptide_iterator);
      peptide = peptide_iterator->next();
      sequence = get_peptide_sequence(peptide); 

      map<string, vector<PEPTIDE_T*> >::iterator find_iter;

      find_iter = sequence_peptide_map.find(sequence);

      carp(CARP_DEBUG,"Adding to map:%s",sequence.c_str());

      if (find_iter == sequence_peptide_map.end()) {
        vector<PEPTIDE_T*> peptide_vector;
        peptide_vector.push_back(peptide);
        sequence_peptide_map.insert(make_pair(sequence, peptide_vector));
      } else {
        find_iter -> second.push_back(peptide);
      }
      



      index = sequence.find(last_sequence);
      // if doesn't contain last peptide
      if (sequence[0] == 'R' && sequence[1] != 'P') { continue;}
      if (index == string::npos || missed_cleavage) {
        missed_cleavage = !missed_cleavage;
        if (!missed_cleavage && last_sequence[last_sequence.size()-1] != 'K') {
	  carp(CARP_DETAILED_DEBUG, "skipping1 %s", get_peptide_sequence(peptide));
	  continue;
	}
	carp(CARP_DETAILED_DEBUG, "peptide %s", get_peptide_sequence(peptide));
        peptides.insert(string(get_peptide_sequence(peptide)));
      } else {
	carp(CARP_DETAILED_DEBUG, "skipping2 %s", get_peptide_sequence(peptide));
	missed_cleavage = false;
      }
      last_sequence = string(get_peptide_sequence(peptide));
    }
  } 
}

void print_precursor_count(vector<LinkedPeptide>& all_ions) {
  int linear_peptides = 0;
  int self_loop_products = 0;
  int dead_end_products = 0;
  int xlink_products = 0;

  for (unsigned int idx=0;idx<all_ions.size();idx++) {
    LinkedPeptide& current = all_ions[idx];
    if (current.is_single())
      linear_peptides++;
    else if (current.is_dead_end())
      dead_end_products++;
    else if (current.is_self_loop())
      self_loop_products++;
    else if (current.is_linked())
      xlink_products++;
    else
      carp(CARP_ERROR,"Unknown product");
  }
  
  carp(CARP_INFO,"Linear Peptides:%d",linear_peptides);
  carp(CARP_INFO,"Inter/Intra Xlinks:%d",xlink_products);
  carp(CARP_INFO,"Dead End Products:%d",dead_end_products);
  carp(CARP_INFO,"Self Loop Products:%d",self_loop_products);

}

// creates an index of all linked peptides from a fasta file
void find_all_precursor_ions(vector<LinkedPeptide>& all_ions, 
			     const char* links, 
			     const char* missed_link_cleavage,
		             const char* database_file,
			     int charge)
{

  carp(CARP_DEBUG,"missed link cleavage:%s", missed_link_cleavage);
  carp(CARP_DEBUG,"find_all_precursor_ions: start()");
  DATABASE_T* db = new_database(database_file, FALSE);
  carp(CARP_DEBUG,"peptide constraint");
  PEPTIDE_CONSTRAINT_T* peptide_constraint = new_peptide_constraint_from_parameters();
  // add 
  set_peptide_constraint_num_mis_cleavage(peptide_constraint, 1);
  //set_verbosity_level(CARP_INFO);
  //Protein* protein = NULL;
  carp(CARP_DEBUG,"protein iterator");
  DATABASE_PROTEIN_ITERATOR_T* protein_iterator = new_database_protein_iterator(db);
  //PROTEIN_PEPTIDE_ITERATOR_T* peptide_iterator = NULL;
  string bonds_string = string(links);
  set<string> peptides;
  carp(CARP_DEBUG,"get_linkable_peptides");
  get_linkable_peptides(peptides, protein_iterator, peptide_constraint);
  carp(CARP_DEBUG,"add_linked_peptides");
  add_linked_peptides(all_ions, peptides, bonds_string, charge);
  
  /*
  if (get_boolean_parameter("xlink-include-linears")) {
    free_database_protein_iterator(protein_iterator);
    protein_iterator = new_database_protein_iterator(db);
    peptides.clear();    
    get_linear_peptides(peptides, protein_iterator, peptide_constraint);
    
    set<string>::iterator iter;
    for (iter = peptides.begin();
	 iter != peptides.end();
	 ++iter) {
      LinkedPeptide lp = LinkedPeptide(charge);
      Peptide p = Peptide(*iter);
      lp.add_peptide(p);
      all_ions.push_back(lp);
    }
  }
  */

  free_database_protein_iterator(protein_iterator);


  carp(CARP_DEBUG,"find_all_precursor_ions: done()");

  //sort by increasing mass.

  //LinkedPeptide::sortByMass(all_ions);

  IF_CARP(CARP_DEBUG,print_precursor_count(all_ions));

 
}



// modified version of crux's estimate_weibull_parameters_from_xcorrs
BOOLEAN_T hhc_estimate_weibull_parameters_from_xcorrs(
  FLOAT_T* scores,
  int num_scores,
  FLOAT_T* eta,
  FLOAT_T* beta,
  FLOAT_T* shift,
  FLOAT_T* correlation,
  Spectrum* spectrum,
  int charge
  ){

  // check that we have the minimum number of matches
  if( num_scores < MIN_WEIBULL_MATCHES ){
    carp(CARP_DETAILED_INFO, "Too few psms (%i) to estimate "
         "p-value parameters for spectrum %i, charge %i",
         num_scores, spectrum->getFirstScan(), charge);
    // set eta, beta, and shift to something???
    return FALSE;
  }

  // reverse sort the first num_samples of them
  sort(scores, scores + num_scores, compareDescending());

  // use only a fraction of the samples, the high-scoring tail
  // this parameter is hidden from the user
  double fraction_to_fit = get_double_parameter("fraction-top-scores-to-fit");
  assert( fraction_to_fit >= 0 && fraction_to_fit <= 1 );
  int num_tail_samples = (int)(num_scores * fraction_to_fit);
  carp(CARP_DEBUG, "Estimating Weibull params with %d psms (%.2f of %i)",
       num_tail_samples, fraction_to_fit, num_scores);

  // do the estimation
  fit_three_parameter_weibull(scores, num_tail_samples, num_scores,
      MIN_XCORR_SHIFT, MAX_XCORR_SHIFT, XCORR_SHIFT, 0,  /*CORR_THRESHOLD*/
      eta, beta, shift, correlation);

  carp(CARP_DETAILED_DEBUG,
      "Corr: %.6f Eta: %.6f Beta: %.6f Shift: %.6f",
      *correlation, *eta, *beta, *shift);

  return TRUE;
}

// helper function for get_linkable_peptides
// creates single and cross-linked peptides, considering every pair
// of peptides and every possible link site

void add_linked_peptides(vector<LinkedPeptide>& all_ions, set<string>& peptides, string links, int charge) {
  XLinkBondMap bonds(links); 
  vector<LinkedPeptide> ions;

  // iterate over both sequences, adding linked peptides with correct links
  for (set<string>::iterator pepA = peptides.begin(); pepA != peptides.end(); ++pepA) {
    char* sequenceA = (char*) pepA->c_str();
    // add unlinked precursor
    LinkedPeptide lp = LinkedPeptide(charge);
    Peptide p = Peptide(sequenceA);
    lp.add_peptide(p);

    //TODO separate linears from xlinking stuff.
    if (get_boolean_parameter("xlink-include-linears")) {
      ions.push_back(lp);
    }

    string seqA = *pepA;
    vector<PEPTIDE_T*>& crux_peptides = get_peptides_from_sequence(seqA);
    
    if (get_boolean_parameter("xlink-include-deadends")) {

      for (unsigned int seq_idx=0;seq_idx < pepA->length();seq_idx++) {
        bool canLink = false;

        for (unsigned int pep_idx=0;pep_idx < crux_peptides.size();pep_idx++) {
          canLink |= bonds.canLink(crux_peptides[pep_idx], seq_idx);
        }
        if (canLink) {
          if (seq_idx == pepA->length()-1 && pepA->at(seq_idx) == 'K') {
            continue;
          } else {
            ions.push_back(LinkedPeptide(sequenceA, NULL, seq_idx, -1, charge));
          }
        }  
      }
    } /* xlink-include-dead-ends */

    if (get_boolean_parameter("xlink-include-selfloops")) {
      
      for (unsigned int seq_idx1 = 0;seq_idx1 < pepA->length()-1; seq_idx1++) {
        for (unsigned int seq_idx2=seq_idx1+1;seq_idx2 < pepA->length()-1; seq_idx2++) {
          bool canLink = false;
          for (unsigned int pep_idx1=0;pep_idx1 < crux_peptides.size();pep_idx1++) {
            canLink |= bonds.canLink(crux_peptides[pep_idx1],seq_idx1,seq_idx2);
          }
          if (canLink) {
            if (seq_idx1 == pepA->length()-1 && pepA->at(seq_idx1) == 'K') continue;
            else if (seq_idx2 == pepA->length()-1 && pepA->at(seq_idx2) == 'K') continue;
            else {
              ions.push_back(LinkedPeptide(sequenceA, NULL, seq_idx1, seq_idx2, charge));
            }
          }
        }
      }
    } /* xlink-include-selfloops */

    {
      set<string>::iterator pepB = pepA;
      pepB++;
      for (;pepB != peptides.end(); ++pepB) {
        string seqB = *pepB;
        vector<PEPTIDE_T*>& crux_peptides2 = get_peptides_from_sequence(seqB);
        
        for (unsigned int seq_idx1 = 0; seq_idx1 < pepA->length();seq_idx1++) {
          for (unsigned int seq_idx2 = 0; seq_idx2 < pepB->length();seq_idx2++) {
            bool canLink = false;
            for (unsigned int pep_idx1=0;pep_idx1 < crux_peptides.size();pep_idx1++) {
              for (unsigned int pep_idx2=0;pep_idx2 < crux_peptides2.size();pep_idx2++) {
                canLink |= bonds.canLink(crux_peptides[pep_idx1],crux_peptides2[pep_idx2],seq_idx1,seq_idx2);
              }
            }
            if (canLink) {
              if (seq_idx1 == pepA->length()-1 && pepA->at(seq_idx1) == 'K') continue;
              else if (seq_idx2 == pepB->length()-1 && pepB->at(seq_idx2) == 'K') continue;
              else {
                char* sequenceB = (char*) pepB->c_str();
                ions.push_back(LinkedPeptide(sequenceA, sequenceB, seq_idx1, seq_idx2, charge));
              }
            }
          }
        }
      }
    } /* xlink-inter/intra links */
  } // get next pepA 
   all_ions.insert(all_ions.end(), ions.begin(), ions.end());
}

// append one shuffled decoy to decoys vector
void add_decoys(vector<LinkedPeptide>& decoys, LinkedPeptide& lp) {
  vector<Peptide> peptides = lp.peptides();
  LinkedPeptide decoy = LinkedPeptide(lp.charge()); 
  Peptide pepA_shuffled = shuffle(peptides[0]);
  decoy.add_peptide(pepA_shuffled);
  if (lp.size() == 2) {
    Peptide pepB_shuffled = shuffle(peptides[1]);
    decoy.add_peptide(pepB_shuffled);
  }
  decoy.set_decoy();
  decoy.calculate_mass(AVERAGE);
  decoys.push_back(decoy);
}

// return a shuffled peptide, preserving any links
Peptide shuffle(Peptide peptide) {
  string shuffled = string(peptide.sequence());
  Peptide shuffled_peptide = Peptide(peptide.sequence());
  for (size_t i = 0; i < shuffled.length(); ++i) {
    if (peptide.has_link_at(i)) shuffled_peptide.add_link(i);
  }

  int start_idx = 1;
  int end_idx = peptide.length() - 2;
  int switch_idx = 0;
  char temp_char = 0;
  while(start_idx <= end_idx){
    switch_idx = get_random_number_interval(start_idx, end_idx);
    temp_char = shuffled[start_idx];
    shuffled[start_idx] = shuffled[switch_idx];
    shuffled[switch_idx] = temp_char;
    if (shuffled_peptide.has_link_at(switch_idx)) {
      //if not a self loop
      if (!shuffled_peptide.has_link_at(start_idx)) {
        shuffled_peptide.remove_link(switch_idx);
        shuffled_peptide.add_link(start_idx);
      }
    } else if (shuffled_peptide.has_link_at(start_idx)) {
      shuffled_peptide.remove_link(start_idx);
      shuffled_peptide.add_link(switch_idx);
    }
    ++start_idx;
  }
  shuffled_peptide.set_sequence(shuffled);
  return shuffled_peptide;
}

//
//
// Peptide and LinkedPeptide class method definitions below
//
//////////////////////////////////////////////////////////


void Peptide::set_sequence(string sequence) {
  sequence_ = sequence;
  length_ = sequence.length();  
}

void Peptide::remove_link(int index) {
  links[index] = false;
  num_links--;
}

bool LinkedPeptide::is_single() {
  return (peptides_.size() == 1 && peptides_[0].link_site() == -1);
} 

bool LinkedPeptide::is_dead_end() {
  return (peptides_.size() == 1 && peptides_[0].get_num_links() == 1);
}

bool LinkedPeptide::is_self_loop() {
  return (peptides_.size() == 1 && peptides_[0].get_num_links() == 2);
}

LinkedPeptide::LinkedPeptide(char* peptide_A, char* peptide_B, int posA, int posB, int charge) {
  mass_calculated[MONO] = false;
  mass_calculated[AVERAGE] = false;
  charge_ = charge;
  decoy_ = false;
  type_ = (ION_TYPE_T)NULL;
  Peptide pepA = Peptide(peptide_A);
  // if a self link or dead end
  if (peptide_B == NULL) {
     pepA.add_link(posA);
    if (posB >= 0)
      pepA.add_link(posB);
    peptides_.push_back(pepA);
  } else {
    carp(CARP_DETAILED_DEBUG, "adding links at %d and %d", posA, posB);
    Peptide pepB = Peptide(peptide_B);
    pepA.add_link(posA);
    pepB.add_link(posB);
    peptides_.push_back(pepA);
    peptides_.push_back(pepB);
  }
}




int Peptide::link_site() {
  for (int i = 0; i < length_; ++i) {
    if (has_link_at(i)) 
      return i;
  }
  return -1;
}

FLOAT_T Peptide::mass(MASS_TYPE_T mass_type) {
  if (mass_calculated[mass_type]) 
    return mass_[mass_type];
  else {
    mass_[mass_type] = calc_sequence_mass((char*)sequence_.c_str(), mass_type);
    mass_calculated[mass_type] = true;
    return mass_[mass_type];
  }
}

// calculates mass of linked peptide,
// remove H2O from mass if it's a b-ion

FLOAT_T LinkedPeptide::mass(MASS_TYPE_T mass_type) {
  if (!mass_calculated[mass_type])
    calculate_mass(mass_type);
  return mass_[mass_type];
   
}

void LinkedPeptide::calculate_mass(MASS_TYPE_T mass_type) {

  if (mass_type == MONO) {
    //cout <<"MONO ";
  } else {
    //cout <<"AVERAGE ";
  }

  mass_[mass_type] = calc_sequence_mass((char*)peptides_[0].sequence().c_str(), mass_type);   
  //cout << "Mass of "<<peptides_[0].sequence().c_str()<<":"
  //    <<calc_sequence_mass((char*)peptides_[0].sequence().c_str(), mass_type)<<endl;

  if (peptides_[0].get_num_links() > 0) {
    mass_[mass_type] += linker_mass;
  }

  if (size() == 2) {
    mass_[mass_type] += calc_sequence_mass((char*)peptides_[1].sequence().c_str(), mass_type);
     //cout << "Mass of "<<peptides_[1].sequence().c_str() <<":"
       //<<calc_sequence_mass((char*)peptides_[1].sequence().c_str(), mass_type)<<endl;
  }

  //cout <<"Total Mass of" << (*this) << ":" <<mass_[mass_type]<<endl;

  //mass += MASS_H_MONO*charge_;
  if (type_ == B_ION) {
    if (mass_type == MONO) {
      mass_[mass_type] -= MASS_H2O_MONO;
    }
    else {
      mass_[mass_type] -= MASS_H2O_AVERAGE;
    }
  } 
  mass_calculated[mass_type] = true;
}

// calculate mass to charge 
FLOAT_T LinkedPeptide::get_mz(MASS_TYPE_T mass_type) {
  //if (mz < 5)
  //mz = mass_ / charge_;
  if (mass_type == MONO) {
    mz[MONO] = (mass(MONO) + MASS_PROTON * charge_) / charge_;
  } else {
    mz[AVERAGE] = (mass(AVERAGE) + MASS_H_AVERAGE * charge_) / charge_;
  }
  return mz[mass_type];
}

void Peptide::add_link(int index) {
  links[index] = true;
  num_links++;
}

// split between every amino acid on every peptide in the
// linked peptide.
void LinkedPeptide::split(vector<pair<LinkedPeptide, LinkedPeptide> >& ion_pairs) {
  bool is_loop = false;
  Peptide peptideA = peptides_[0];  
  Peptide peptideB = peptides_[0];
  // dead end
  if (is_dead_end()) {
   peptideB.set_sequence("");
  }
  // if a loop
  if (peptideA.get_num_links() == 2) {
    is_loop = true;
  } else if (is_linked()) {
    peptideB = peptides_[1];
  }
  for (int i = 1; i < peptideA.length(); ++i) {
    peptideA.split_at(i, ion_pairs, charge_, peptideB, is_loop);
  }
 
  if (is_linked()) {
    for (int i = 1; i < peptideB.length(); ++i) {
      peptideB.split_at(i, ion_pairs, charge_, peptideA, is_loop);
    }
  } 
}

void LinkedPeptide::splitA(vector<pair<LinkedPeptide, LinkedPeptide> >& ion_pairs) {
  Peptide peptideA = peptides_[0];
  Peptide peptideB = peptides_[1];

  for (int i = 1; i < peptideA.length(); ++i) {
    peptideA.split_at(i, ion_pairs, charge_, peptideB, false);
  }
}

void LinkedPeptide::splitB(vector<pair<LinkedPeptide, LinkedPeptide> >& ion_pairs) {
  Peptide peptideA = peptides_[0];
  Peptide peptideB = peptides_[1];

  for (int i = 1; i < peptideB.length(); ++i) {
    peptideB.split_at(i, ion_pairs, charge_, peptideA, false);
  }
}


// temporary
std::ostream &operator<< (std::ostream& os, LinkedPeptide& lp) {
  vector<Peptide> peptides = lp.peptides();
  ostringstream link_positions;
  link_positions << "(";
  for (int i = 0; i < peptides[0].length(); ++i) {
	if (peptides[0].has_link_at(i)) link_positions << (i+1) << "," ;
  }
  if (peptides.size() == 2) {
    for (int i = 0; i < peptides[1].length(); ++i) {
	if (peptides[1].has_link_at(i)) link_positions << (i+1) << ")";
    }
    return os << peptides[0].sequence() << ", " << peptides[1].sequence() << " " << link_positions.str();// << " +" << lp.charge();
  }
  string final = link_positions.str();
  if (final.length() > 1) final.erase(final.end()-1);
  return os << peptides[0].sequence() << " " << final << ")";// +" << lp.charge();
}

// this needs to change
void Peptide::split_at(int index, vector<pair<LinkedPeptide, LinkedPeptide> >& pairs, int charge, Peptide& other, bool is_loop) {
  bool self_flag = false;
  // for every charge state
  for (int c = 0; c <= charge; ++c) {
    Peptide pepA = Peptide(sequence_.substr(0, index));
    Peptide pepB = Peptide(sequence_.substr(index, length_ - index));
    LinkedPeptide linkedA = LinkedPeptide(c);
    LinkedPeptide linkedB = LinkedPeptide(charge - c);
    self_flag = true;
    // for every position on pepA
    for (int i = 0; i < index; i++) {
      if (has_link_at(i)) {
        pepA.add_link(i); 
        if (!other.empty() && !is_loop) linkedA.add_peptide(other);
        // if a loop, skip cleavages in between link sites (same mass as precursor)
	if (is_loop) self_flag = !self_flag;
      }
    }
    if (!self_flag) continue;
    // for every position on pepB
    for (int i = index; i < length_; i++) {
      if (has_link_at(i)) {
        pepB.add_link(i - index);
	if (!other.empty() && !is_loop) linkedB.add_peptide(other);
	//else (self_flag = !self_flag);
      }
    } 
    linkedA.add_peptide(pepA);
    linkedB.add_peptide(pepB);
    pairs.push_back(pair<LinkedPeptide, LinkedPeptide> (linkedA, linkedB));
  }
}


bool compareMassAverage(const LinkedPeptide& lp1, const LinkedPeptide& lp2) {
  return lp1.mass_[AVERAGE] < lp2.mass_[AVERAGE];
}

bool compareMassMono(const LinkedPeptide& lp1, const LinkedPeptide& lp2) {
  return lp1.mass_[MONO] < lp2.mass_[MONO];
}


void LinkedPeptide::sortByMass(std::vector<LinkedPeptide>& linked_peptides, MASS_TYPE_T mass_type) {
  //TODO : should we put code here to make sure that we have
  //calculated the all of the masses for a particular mass type?

  if (mass_type == MONO) {
    sort(linked_peptides.begin(), linked_peptides.end(), compareMassMono);
  }
  else {
    sort(linked_peptides.begin(), linked_peptides.end(), compareMassAverage);
  }
}
