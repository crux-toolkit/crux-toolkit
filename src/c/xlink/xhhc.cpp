#include "xhhc.h"
#include "XLinkBondMap.h"
#include "LinkedPeptide.h"
#include "XHHC_Peptide.h"
#include "XLinkablePeptide.h"


#include "parameter.h"
#include "Peptide.h"
#include "DelimitedFile.h"
#include "DatabaseProteinIterator.h"
#include "ProteinPeptideIterator.h"

using namespace std;
using namespace Crux;

/**
 * Global Variable
 */
map<string, vector<Peptide*> > sequence_peptide_map; //hack to keep track of peptides.


void get_linear_peptides(
  set<string>& peptides, 
  DatabaseProteinIterator* protein_iterator,
  PeptideConstraint* peptide_constraint
  ) {

  int max_missed_cleavages = get_int_parameter("missed-cleavages");

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
      char* seq = peptide->getSequence();
      sequence = seq; 
      free(seq);

      if (peptide->getMissedCleavageSites() <= max_missed_cleavages) {
        carp(CARP_DEBUG, "Adding linear peptide:%s",peptide->getSequence());
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
      } else {
        delete peptide;
      }

    }
    delete peptide_iterator;
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


void get_linkable_peptides(
        set<XLinkablePeptide>& peptides, 
        XLinkBondMap& bondmap,
  DatabaseProteinIterator* protein_iterator,
  PeptideConstraint* peptide_constraint) 
{
  ProteinPeptideIterator* peptide_iterator = NULL;
  Protein* protein;
  Peptide* peptide;
  string sequence = "";
  string last_sequence = "zz";
  // keep track of whether the next peptide contains the previous one or not
  while (protein_iterator->hasNext()) {
    protein = protein_iterator->next();
    peptide_iterator = new ProteinPeptideIterator(protein, peptide_constraint);
    // missed_cleavages must be all in protein.c for this to work
    //TODO either fix this code or make an option to allow all missed cleavages...
    peptide_iterator->prepareMc(500); 
    while (peptide_iterator->hasNext()) {
      //peptide = database_peptide_iterator_next(peptide_iterator);
      peptide = peptide_iterator->next();
      char* seq = peptide->getSequence();
      sequence = seq;
      free(seq);
 
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

      XLinkablePeptide xlinkable(peptide, bondmap);
  
      if (xlinkable.isLinkable()) {
        peptides.insert(xlinkable);
      }

    }
    delete peptide_iterator;
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
void find_all_precursor_ions(
  vector<LinkedPeptide>& all_ions
  ) {

  char* database_file = get_string_parameter("protein-database");
  
  carp(CARP_DEBUG,"find_all_precursor_ions: start()");
  Database* db = new Database(database_file, false);
  carp(CARP_DEBUG,"peptide constraint");
  PeptideConstraint* peptide_constraint = 
    PeptideConstraint::newFromParameters();
  // add two to account for a self loop that prevents two cleavages
  peptide_constraint->setNumMisCleavage(get_int_parameter("missed-cleavages") + 2);

  carp(CARP_DEBUG,"protein iterator");
  DatabaseProteinIterator* protein_iterator = new DatabaseProteinIterator(db);

  set<XLinkablePeptide> peptides;
  carp(CARP_DEBUG,"get_linkable_peptides");
  XLinkBondMap bondmap;
  get_linkable_peptides(peptides, bondmap, protein_iterator, peptide_constraint);
  carp(CARP_DEBUG, "add_linked_peptides");
  add_linked_peptides(all_ions, peptides, bondmap, 1);
  
  if (get_boolean_parameter("xlink-include-linears")) {
    delete protein_iterator;
    protein_iterator = new DatabaseProteinIterator(db);
    set<string> linear_peptides;
    get_linear_peptides(linear_peptides, protein_iterator, peptide_constraint);
    
    set<string>::iterator iter;
    for (iter = linear_peptides.begin();
      iter != linear_peptides.end();
      ++iter) {

      LinkedPeptide lp = LinkedPeptide(1);
      XHHC_Peptide p(*iter);
      lp.addPeptide(p);
      all_ions.push_back(lp);
    }
  }
  

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
  sort(scores, scores + num_scores, greater<FLOAT_T>());

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

void add_linked_peptides(
  vector<LinkedPeptide>& all_ions, 
  set<XLinkablePeptide>& peptides, 
  XLinkBondMap& bondmap, 
  int charge) {

  vector<LinkedPeptide> ions;

  // iterate over both sequences, adding linked peptides with correct links
  for (set<XLinkablePeptide>::iterator iterA = peptides.begin(); iterA != peptides.end(); ++iterA) {

    XLinkablePeptide pepA = *iterA;
    string seqA = pepA.getModifiedSequenceString();
    
    // add unlinked precursor
    LinkedPeptide lp = LinkedPeptide(charge);
    XHHC_Peptide p = XHHC_Peptide(seqA);
    lp.addPeptide(p);

    //vector<Peptide*>& crux_peptides = get_peptides_from_sequence(seqA);

    int max_missed_cleavages = get_int_parameter("missed-cleavages");
    int seqA_missed_cleavages = pepA.getMissedCleavageSites();
    carp(CARP_DEBUG,"seqA:%s %i %i",seqA.c_str(), pepA.numLinkSites(), seqA_missed_cleavages);
    if (get_boolean_parameter("xlink-include-deadends") && 
      seqA_missed_cleavages <= (max_missed_cleavages+1)) {
      
      for (size_t link_idx = 0;link_idx < pepA.numLinkSites();link_idx++) {
        int seq_idx = pepA.getLinkSite(link_idx);
        int total_cleavages = seqA_missed_cleavages;
        if (pepA.linkSitePreventsCleavage(link_idx)) {
          total_cleavages--;
        }
        if (total_cleavages <= max_missed_cleavages) {
          ions.push_back(LinkedPeptide((char*)seqA.c_str(), NULL, seq_idx, -1, charge));
        }
      }
    } /* xlink-include-dead-ends */

    if (get_boolean_parameter("xlink-include-selfloops") && 
        (seqA_missed_cleavages <= (max_missed_cleavages+2))) {
      
      for (size_t link_idx1 = 0;link_idx1 < pepA.numLinkSites()-1;link_idx1++) {
        for (size_t link_idx2 = link_idx1+1;link_idx2 < pepA.numLinkSites();link_idx2++) {
          int seq_idx1 = pepA.getLinkSite(link_idx1);
          int seq_idx2 = pepA.getLinkSite(link_idx2);

          if (bondmap.canLink(pepA.getPeptide(), seq_idx1, seq_idx2)) {
            int total_missed_cleavages = seqA_missed_cleavages;
            if (pepA.linkSitePreventsCleavage(link_idx1)) {
              total_missed_cleavages--;
            }
            if (pepA.linkSitePreventsCleavage(link_idx2)) {
              total_missed_cleavages--;
            }
            if (total_missed_cleavages <= max_missed_cleavages) {
              ions.push_back(LinkedPeptide((char*)seqA.c_str(), NULL, seq_idx1, seq_idx2, charge));
            }
          }
        }
      }
    } /* xlink-include-selfloops */

    /* xlink-inter/intra links */
    if (seqA_missed_cleavages <= max_missed_cleavages+1) {
    
      for (set<XLinkablePeptide>::iterator iterB = iterA; iterB != peptides.end();++iterB) {
        XLinkablePeptide pepB = *iterB;
        string seqB = pepB.getModifiedSequenceString();
        int seqB_missed_cleavages = pepB.getMissedCleavageSites();
        if (seqB_missed_cleavages <= (max_missed_cleavages+1)) {
          for (size_t link_idx1 = 0;link_idx1 < pepA.numLinkSites();link_idx1++) {
            for (size_t link_idx2 = 0; link_idx2 < pepB.numLinkSites();link_idx2++) {
              int seq_idx1 = pepA.getLinkSite(link_idx1);
              int seq_idx2 = pepB.getLinkSite(link_idx2);
              if (bondmap.canLink(pepA.getPeptide(), pepB.getPeptide(), seq_idx1, seq_idx2)) {

                int totalA = seqA_missed_cleavages;
                if (pepA.linkSitePreventsCleavage(link_idx1)) {
                  totalA--;
                }

                int totalB = seqB_missed_cleavages;
                if (pepB.linkSitePreventsCleavage(link_idx2)) {
                  totalB--;
                }

                if (totalA <= max_missed_cleavages && 
                    totalB <= max_missed_cleavages) {
                  ions.push_back(LinkedPeptide((char*)seqA.c_str(), (char*)seqB.c_str(), seq_idx1, seq_idx2, charge));
                }
              }
            }
          }
        }
      }
    }
  } // get next pepA 
  all_ions.insert(all_ions.end(), ions.begin(), ions.end());
}

// append one shuffled decoy to decoys vector
void add_decoy(vector<LinkedPeptide>& decoys, LinkedPeptide& lp) {
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

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */


