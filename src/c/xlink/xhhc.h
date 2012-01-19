#ifndef XHHC_H
#define XHHC_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <fstream>

/*Crux Includes*/
#include "utils.h"
#include "Spectrum.h"
#include "carp.h"
#include "PeptideConstraint.h"
#include "Peptide.h"
#include "PeptideSrc.h"
#include "Protein.h"
#include "Database.h"
#include "DatabaseProteinIterator.h"
#include "parse_arguments.h"
#include "parameter.h"
#include "objects.h"
#include "crux-utils.h"
#include "Ion.h"

// get rid of these
//#define PARAM_ESTIMATION_SAMPLE_COUNT 500
#define MIN_WEIBULL_MATCHES 40 
#define MIN_XCORR_SHIFT -5.0
#define MAX_XCORR_SHIFT  5.0
#define XCORR_SHIFT 0.05

// mine
#define BONFERRONI_CUT_OFF_P 0.0001
#define BONFERRONI_CUT_OFF_NP 0.01
#define MIN_WEIBULL_SAMPLES 750 
#define MIN_PRECURSORS 3


//HACK for getting protein ids.
std::vector<Peptide*>& get_peptides_from_sequence(std::string& sequence);
void free_peptides();

class XHHC_Peptide;

class XHHC_Peptide {
 public:
  XHHC_Peptide() {mass_calculated[MONO] = false; mass_calculated[AVERAGE] = false;}
  XHHC_Peptide(std::string sequence) : num_links(0), sequence_(sequence), length_(sequence.length())  {
    mass_calculated[MONO] = false;
    mass_calculated[AVERAGE] = false;
    int length = sequence.length();
    for (int i = 0; i < length; ++i)
      links.push_back(false);
  }
    // cleave the peptide at an index, adding b and y ions
    // skips if a cleave site on self loop between the link
    void split_at(int index, std::vector<std::pair<LinkedPeptide, LinkedPeptide> >& pairs, int charge, XHHC_Peptide& other, bool is_loop);
    bool has_link_at(int index)		{ return (links[index]); }
    int length()			{ return length_; }
    bool has_link() 			{ return num_links > 0; }
    std::string sequence() 		{ return sequence_;}	
    bool empty() 			{ return length_ == 0; }
    int get_num_links() 		{ return num_links; }
    void add_link(int index);
    void remove_link(int index);
    int link_site(); 
    void set_sequence(std::string sequence);
    FLOAT_T mass(MASS_TYPE_T mass_type);
private:
    bool mass_calculated[NUMBER_MASS_TYPES];
    int num_links;
    std::vector<bool> links;
    std::string sequence_;
    int length_;
    FLOAT_T mass_[NUMBER_MASS_TYPES];
};


// methods used by xhhc_search.cpp and xhhc_make_histogram.cpp
// defined in xhhc.cpp
//////////////////////////////////////////

XHHC_Peptide shuffle(XHHC_Peptide peptide);

bool hhc_estimate_weibull_parameters_from_xcorrs(
  FLOAT_T* scores,
  int num_scores,
  FLOAT_T* eta,
  FLOAT_T* beta,
  FLOAT_T* shift,
  FLOAT_T* correlation,
  Spectrum* spectrum,
  int charge
  );

void get_linkable_peptides(std::set<std::string>& peptides, 
	DatabaseProteinIterator* protein_iterator,
	PeptideConstraint* peptide_constraint); 

void add_linked_peptides(std::vector<LinkedPeptide>& all_ions, std::set<std::string>& peptides, std::string links, int charge);

//void add_decoys(vector<LinkedPeptide>& decoys, LinkedPeptide& lp, char* links, int charge, FLOAT_T linker_mass);

void add_decoys(std::vector<LinkedPeptide>& decoys, LinkedPeptide& lp);


void find_all_precursor_ions(std::vector<LinkedPeptide>& all_ions, 
	const char* links, 
	const char* missed_link_cleavage, 
	const char* database_file,
	int charge);


#endif
