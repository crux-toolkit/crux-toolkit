#include "xhhc.h"
#include "XLinkBondMap.h"
#include "LinkedPeptide.h"
#include "XHHC_Peptide.h"

#include "parameter.h"
#include "Peptide.h"
#include "DelimitedFile.h"
#include "DatabaseProteinIterator.h"
#include "ProteinPeptideIterator.h"

using namespace std;




map<string, vector<Peptide*> > sequence_peptide_map; //hack to keep track of peptides.

void get_linear_peptides(set<string>& peptides,
			 DatabaseProteinIterator* protein_iterator,
			 PeptideConstraint* peptide_constraint) {

  ProteinPeptideIterator* peptide_iterator = NULL;
  Protein* protein;
  Peptide* peptide;

  string sequence = "";
  string last_sequence = "zz";
  //bool missed_cleavage = false;
  // keep track of whether the next peptide contains the previous one or not
  //size_t index;
  while (protein_iterator->hasNext()) {
    protein = protein_iterator->next();
    peptide_iterator = new ProteinPeptideIterator(protein, peptide_constraint);
    // missed_cleavages must be all in protein.c for this to work
    peptide_iterator->prepareMc(500);
    while (peptide_iterator->hasNext()) {
      //peptide = database_peptide_iterator_next(peptide_iterator);
      peptide = peptide_iterator->next();
      sequence = peptide->getSequence(); 
      carp(CARP_INFO,"Adding linear peptide:%s",peptide->getSequence());
      peptides.insert(sequence);

      map<string, vector<Peptide*> >::iterator find_iter;

      find_iter = sequence_peptide_map.find(sequence);

      carp(CARP_DEBUG,"Adding to map:%s,",sequence.c_str());

      if (find_iter == sequence_peptide_map.end()) {
        vector<Peptide*> peptide_vector;
        peptide_vector.push_back(peptide);
        sequence_peptide_map.insert(make_pair(sequence, peptide_vector));
      } else {
        find_iter -> second.push_back(peptide);
      }
      

    }
  } 
}

vector<Peptide*>& get_peptides_from_sequence(string& sequence) {
  return sequence_peptide_map[sequence];
}

void free_peptides() {
  map<string, vector<Peptide*> >::iterator map_iter;

  for (map_iter = sequence_peptide_map.begin();
       map_iter != sequence_peptide_map.end();
       ++map_iter) {
    vector<Peptide*>& peptides = map_iter -> second;
    for (unsigned int i=0;i<peptides.size();i++) {
      delete peptides[i];
    }
  }
  sequence_peptide_map.clear();
}


// a hack, works for EDC linker only
void get_linkable_peptides(set<string>& peptides, 
	DatabaseProteinIterator* protein_iterator,
	PeptideConstraint* peptide_constraint) 
{
  ProteinPeptideIterator* peptide_iterator = NULL;
  Protein* protein;
  Peptide* peptide;
  string sequence = "";
  string last_sequence = "zz";
  bool missed_cleavage = false;
  // keep track of whether the next peptide contains the previous one or not
  size_t index;
  while (protein_iterator->hasNext()) {
    protein = protein_iterator->next();
    peptide_iterator = new ProteinPeptideIterator(protein, peptide_constraint);
    // missed_cleavages must be all in protein.c for this to work
    //TODO either fix this code or make an option to allow all missed cleavages...
    peptide_iterator->prepareMc(500); 
    while (peptide_iterator->hasNext()) {
      //peptide = database_peptide_iterator_next(peptide_iterator);
      peptide = peptide_iterator->next();
      sequence = peptide->getSequence(); 

      map<string, vector<Peptide*> >::iterator find_iter;

      find_iter = sequence_peptide_map.find(sequence);

      carp(CARP_DEBUG,"Adding to map:%s",sequence.c_str());

      if (find_iter == sequence_peptide_map.end()) {
        vector<Peptide*> peptide_vector;
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
	  carp(CARP_DETAILED_DEBUG, "skipping1 %s", peptide->getSequence());
	  continue;
	}
	carp(CARP_DETAILED_DEBUG, "peptide %s", peptide->getSequence());
        peptides.insert(string(peptide->getSequence()));
      } else {
	carp(CARP_DETAILED_DEBUG, "skipping2 %s", peptide->getSequence());
	missed_cleavage = false;
      }
      last_sequence = string(peptide->getSequence());
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
    if (current.isLinear())
      linear_peptides++;
    else if (current.isDeadEnd())
      dead_end_products++;
    else if (current.isSelfLoop())
      self_loop_products++;
    else if (current.isCrossLinked())
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
  Database* db = new Database(database_file, false);
  carp(CARP_DEBUG,"peptide constraint");
  PeptideConstraint* peptide_constraint = 
    PeptideConstraint::newFromParameters();
  // add 
  peptide_constraint->setNumMisCleavage(get_int_parameter("missed-cleavages") + 1);
  //set_verbosity_level(CARP_INFO);
  //Protein* protein = NULL;
  carp(CARP_DEBUG,"protein iterator");
  DatabaseProteinIterator* protein_iterator = new DatabaseProteinIterator(db);
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

  delete protein_iterator;


  carp(CARP_DEBUG,"find_all_precursor_ions: done()");

  //TODO - sort by increasing mass.
  IF_CARP(CARP_DEBUG,print_precursor_count(all_ions));
}



// modified version of crux's estimate_weibull_parameters_from_xcorrs
bool hhc_estimate_weibull_parameters_from_xcorrs(
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
    return false;
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

  return true;
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
    XHHC_Peptide p = XHHC_Peptide(sequenceA);
    lp.addPeptide(p);

    //TODO separate linears from xlinking stuff.
    if (get_boolean_parameter("xlink-include-linears")) {
      ions.push_back(lp);
    }

    string seqA = *pepA;
    vector<Peptide*>& crux_peptides = get_peptides_from_sequence(seqA);
    
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

      for (set<string>::iterator pepB = pepA ;pepB != peptides.end(); ++pepB) {
        string seqB = *pepB;
        vector<Peptide*>& crux_peptides2 = get_peptides_from_sequence(seqB);
        
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
  vector<XHHC_Peptide> peptides = lp.getPeptides();
  LinkedPeptide decoy = LinkedPeptide(lp.getCharge()); 
  XHHC_Peptide pepA_shuffled = peptides[0].shuffle();
  decoy.addPeptide(pepA_shuffled);
  if (lp.size() == 2) {
    XHHC_Peptide pepB_shuffled = peptides[1].shuffle();
    decoy.addPeptide(pepB_shuffled);
  }
  decoy.setDecoy();
  decoy.calculateMass(AVERAGE);
  decoys.push_back(decoy);
}




