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
