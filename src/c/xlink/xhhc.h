/**
 * \file xhhc.h
 * $Revision: 1.00$
 * \brief utility routines for xlinking programs
 */
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

#include "XLinkBondMap.h"
#include "XLinkablePeptide.h"

// TODO: get rid of these
#define MIN_WEIBULL_MATCHES 40 
#define MIN_XCORR_SHIFT -5.0
#define MAX_XCORR_SHIFT  5.0
#define XCORR_SHIFT 0.05

// mine
#define BONFERRONI_CUT_OFF_P 0.0001
#define BONFERRONI_CUT_OFF_NP 0.01
#define MIN_WEIBULL_SAMPLES 750 
#define MIN_PRECURSORS 3

/**
 * \returns the peptide objects from a sequence string
 * TODO - eliminate this from xlink (this is a HACK)
 */
std::vector<Crux::Peptide*>& get_peptides_from_sequence(
  std::string& sequence ///< the peptide sequence string
  );

/**
 * frees allocated peptides
 */
void free_peptides();


bool hhc_estimate_weibull_parameters_from_xcorrs(
  FLOAT_T* scores,
  int num_scores,
  FLOAT_T* eta,
  FLOAT_T* beta,
  FLOAT_T* shift,
  FLOAT_T* correlation,
  Crux::Spectrum* spectrum,
  int charge
  );

void get_linkable_peptides(
  std::set<XLinkablePeptide>& peptides, 
  XLinkBondMap& bondmap,
  DatabaseProteinIterator* protein_iterator,
  PeptideConstraint* peptide_constraint
  ); 

void add_linked_peptides(
  std::vector<LinkedPeptide>& all_ions, 
  std::set<XLinkablePeptide>& peptides, 
  XLinkBondMap& bondmap, 
  int charge
  );

/**
 * appends one shuffled decoy to a vector
 */
void add_decoy(
  std::vector<LinkedPeptide>& decoys, ///< decoy list to add to
  LinkedPeptide& lp ///< linked peptide to generate decoy from
  );


/**
 * Generates all of the product precursors for the
 * search.
 */
void find_all_precursor_ions(
  std::vector<LinkedPeptide>& all_ions ///< all of the precursors found -out 
  );

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
