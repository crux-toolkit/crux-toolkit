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
#include "peptide_constraint.h"
#include "peptide.h"
#include "peptide_src.h"
#include "Protein.h"
#include "database.h"
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
std::vector<PEPTIDE_T*>& get_peptides_from_sequence(std::string& sequence);
void free_peptides();


/*
void prepare_protein_peptide_iterator(
    PROTEIN_PEPTIDE_ITERATOR_T* iterator
  );
*/

//using namespace std;

//typedef std::map<char, std::set<char> > BondMap;

class BondMap: public std::map<char, std::set<char> > {

public:
  BondMap();
  virtual ~BondMap();
  BondMap(std::string links_string);

};


class Peptide;
class LinkedPeptide {
 public:
  // constructor for a linked peptide. If A or B null, then
  // a self link will be created. If an index is -1, a link to nothing
  // will be created.
  LinkedPeptide() {mass_calculated[MONO] = false; mass_calculated[AVERAGE]=false;}
    LinkedPeptide(char* peptide_A, char* peptide_B, int posA, int posB, int charge);
    LinkedPeptide(int charge) : charge_(charge), decoy_(false) {mass_calculated[MONO]=false; mass_calculated[AVERAGE]=false;}
    std::vector<Peptide>& peptides() 	{ return peptides_;}
    int charge() 			{ return charge_;}
    void set_charge(int charge)		{ charge_ = charge; }
    void set_type(ION_TYPE_T type)   	{ type_ = type; }
    ION_TYPE_T type()			{ return type_;}
    int size()				{ return peptides_.size(); }

    void set_decoy()			{ decoy_ = true; }
    bool is_decoy()			{ return decoy_;}
    void add_peptide(Peptide& peptide)  { peptides_.push_back(peptide); } 
    bool is_linked() 			{ return size() == 2; }
    bool is_dead_end();
    bool is_self_loop();
    bool is_single();
    FLOAT_T get_mz(MASS_TYPE_T mass_type);
    void calculate_mass(MASS_TYPE_T mass_type);
    FLOAT_T mass(MASS_TYPE_T mass_type);
  
    void splitA(std::vector<std::pair<LinkedPeptide, LinkedPeptide> >& ion_pairs);
    void splitB(std::vector<std::pair<LinkedPeptide, LinkedPeptide> >& ion_pairs);

    void split(std::vector<std::pair<LinkedPeptide, LinkedPeptide> >& ion_pairs);
    friend std::ostream &operator<< (std::ostream& os, LinkedPeptide& lp); 

    // for sorting LinkedPeptides by mass
    friend bool operator< (const LinkedPeptide &lp1, const LinkedPeptide &lp2) {
      return lp1.mass_ < lp2.mass_;
      //return lp1.mz < lp2.mz;
    }

    friend bool compareMassMono(const LinkedPeptide& lp1, const LinkedPeptide& lp2);
    friend bool compareMassAverage(const LinkedPeptide& lp1, const LinkedPeptide& lp2);

    static void sortByMass(std::vector<LinkedPeptide>& linked_peptides, MASS_TYPE_T mass_type=MONO);
    


    static FLOAT_T linker_mass;
   private:
    bool mass_calculated[NUMBER_MASS_TYPES]; //MONO or AVERAGE.
    int charge_;
    bool decoy_;
    ION_TYPE_T type_; //B_ION or Y_ION
      FLOAT_T mass_[NUMBER_MASS_TYPES];
      FLOAT_T mz[NUMBER_MASS_TYPES];
      std::vector<Peptide> peptides_;
};

//FLOAT_T LinkedPeptide::linker_mass;

class Peptide {
 public:
  Peptide() {mass_calculated[MONO] = false; mass_calculated[AVERAGE] = false;}
  Peptide(std::string sequence) : num_links(0), sequence_(sequence), length_(sequence.length())  {
    mass_calculated[MONO] = false;
    mass_calculated[AVERAGE] = false;
    int length = sequence.length();
    for (int i = 0; i < length; ++i)
      links.push_back(false);
  }
    // cleave the peptide at an index, adding b and y ions
    // skips if a cleave site on self loop between the link
    void split_at(int index, std::vector<std::pair<LinkedPeptide, LinkedPeptide> >& pairs, int charge, Peptide& other, bool is_loop);
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
    //std::map<int, Peptide> links_;
    std::vector<bool> links;
    std::string sequence_;
    int length_;
    FLOAT_T mass_[NUMBER_MASS_TYPES];
};


// methods used by xhhc_search.cpp and xhhc_make_histogram.cpp
// defined in xhhc.cpp
//////////////////////////////////////////

Peptide shuffle(Peptide peptide);

BOOLEAN_T hhc_estimate_weibull_parameters_from_xcorrs(
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
	DATABASE_PROTEIN_ITERATOR_T* protein_iterator,
	PEPTIDE_CONSTRAINT_T* peptide_constraint); 

void add_linked_peptides(std::vector<LinkedPeptide>& all_ions, std::set<std::string>& peptides, std::string links, int charge);

//void add_decoys(vector<LinkedPeptide>& decoys, LinkedPeptide& lp, char* links, int charge, FLOAT_T linker_mass);

void add_decoys(std::vector<LinkedPeptide>& decoys, LinkedPeptide& lp);


void find_all_precursor_ions(std::vector<LinkedPeptide>& all_ions, 
	const char* links, 
	const char* missed_link_cleavage, 
	const char* database_file,
	int charge);


#endif
