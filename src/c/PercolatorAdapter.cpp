/**
 * \file PercolatorAdapter.cpp
 * $Revision$
 * \brief Converts Percolator result objects to Crux result objects.
 */

#include "PercolatorAdapter.h"

using namespace Crux;

PercolatorAdapter::PercolatorAdapter() {
}

PercolatorAdapter::~PercolatorAdapter() {
  // TODO Auto-generated destructor stub
}

MatchCollection* PercolatorAdapter::psmScoresToMatchCollection(Scores* scores) {

  // Create new MatchCollection object that will be the converted Percolator Scores
  MatchCollection* converted = new MatchCollection();

  // Iterate over each ScoreHolder in Scores object
  for (
    vector<ScoreHolder>::iterator score_itr = scores->begin();
    score_itr != scores->end();
    score_itr++
    ) {

    PSMDescription* psm = score_itr->pPSM;
    int charge_state = parseChargeState(psm->id);

      //const char*         id, ///< The protein sequence id.
      //const char*   sequence, ///< The protein sequence.
      //unsigned int length, ///< The length of the protein sequence.
      //const char* annotation,  ///< Optional protein annotation.  -in
      //unsigned long int offset,///< The file location in the source file in the database -in
      //unsigned int protein_idx,///< The index of the protein in its database. -in
      //Database* database ///< the database of its origin
    PostProcessProtein* parent_protein = new PostProcessProtein();
    parent_protein->setId((*psm->proteinIds.begin()).c_str());
    int start_idx = parent_protein->findStart(psm->getPeptideSequence(), "", "");
      //unsigned char length,     ///< The length of the peptide -in
      //FLOAT_T peptide_mass,       ///< The neutral mass of the peptide -in
      //Protein* parent_protein, ///< The parent_protein of this peptide -in
      //int start_idx ///< Start index of peptide in the protein sequence -in
    FLOAT_T peptide_mass = (psm->expMass * charge_state) - ((charge_state - 1) * MASS_PROTON);
    Crux::Peptide* peptide = new Crux::Peptide(psm->getPeptideSequence().size(), peptide_mass, parent_protein, start_idx);

      //int               first_scan,         ///< number of the first scan -in
      //int               last_scan,          ///< number of the last scan -in
      //FLOAT_T           precursor_mz,       ///< m/z of the precursor
      //const std::vector<int>& possible_z,   ///< possible charge states
      //const char*       filename
    vector<int> zStates;
    zStates.push_back(charge_state);
    Crux::Spectrum* spectrum = new Crux::Spectrum(psm->scan, psm->scan, psm->expMass, zStates, "");

      //FLOAT_T neutral_mass,
      //int charge
    SpectrumZState zState(peptide_mass, charge_state);

      //Peptide* peptide, ///< the peptide for this match
      //Spectrum* spectrum, ///< the spectrum for this match
      //SpectrumZState& zstate, ///< the charge/mass of the spectrum
      //bool is_decoy);///< is the peptide a decoy or not
    Match* match = new Match(peptide, spectrum, zState, score_itr->isDecoy());
    match->setScore(PERCOLATOR_SCORE, score_itr->score);
    match->setScore(PERCOLATOR_QVALUE, psm->q);
    match->setScore(PERCOLATOR_PEP, psm->pep);

    converted->addMatch(match);
  }

  converted->forceScoredBy(PERCOLATOR_SCORE);
  converted->forceScoredBy(PERCOLATOR_QVALUE);
  converted->forceScoredBy(PERCOLATOR_PEP);
  converted->populateMatchRank(PERCOLATOR_SCORE);
  return converted;

}

void PercolatorAdapter::addPsmScores(ProteinMatchCollection* collection, Scores* scores) {

  MatchCollection* matches = psmScoresToMatchCollection(scores);
  collection->addMatches(matches);
  delete matches;

}

void PercolatorAdapter::addProteinScores(ProteinMatchCollection* collection, Scores* scores) {

  // Iterate over each ScoreHolder in Scores object
  for (
    vector<ScoreHolder>::iterator score_itr = scores->begin();
    score_itr != scores->end();
    score_itr++
    ) {

    PSMDescription* psm = score_itr->pPSM;

    // Set scores
    ProteinMatch* protein_match = collection->getProteinMatch(*psm->proteinIds.begin());
    protein_match->setScore(PERCOLATOR_SCORE, score_itr->score);
    protein_match->setScore(PERCOLATOR_QVALUE, psm->q);
    protein_match->setScore(PERCOLATOR_PEP, psm->pep);
  }

}

void PercolatorAdapter::addPeptideScores(ProteinMatchCollection* collection, Scores* scores) {

  // Iterate over each ScoreHolder in Scores object
  for (
    vector<ScoreHolder>::iterator score_itr = scores->begin();
    score_itr != scores->end();
    score_itr++
    ) {

    PSMDescription* psm = score_itr->pPSM;
    int charge_state = parseChargeState(psm->id);

    // Set scores
    PeptideMatch* peptide_match = collection->getPeptideMatch(psm->getPeptideSequence());
    peptide_match->setScore(PERCOLATOR_SCORE, score_itr->score);
    peptide_match->setScore(PERCOLATOR_QVALUE, psm->q);
    peptide_match->setScore(PERCOLATOR_PEP, psm->pep);
  }

}

/**
 * Given a Percolator psm_id in the form ".*_([0-9]+)_[^_]*",
 * find the charge state (matching group)
 */
int PercolatorAdapter::parseChargeState(string psm_id)
{
  size_t charge_begin, charge_end;

  charge_end = psm_id.rfind("_");
  if (charge_end < 0)
  {
    return -1;
  }

  charge_begin = psm_id.rfind("_", charge_end - 1) + 1;
  if (charge_begin < 0)
  {
    return -1;
  }

  string charge_str = psm_id.substr(charge_begin, charge_end - charge_begin);
  return atoi(charge_str.c_str());
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
