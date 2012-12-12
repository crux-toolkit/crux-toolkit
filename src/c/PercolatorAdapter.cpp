/**
 * \file PercolatorAdapter.cpp
 * $Revision$
 * \brief Converts Percolator result objects to Crux result objects.
 */

#include "PercolatorAdapter.h"

PercolatorAdapter::PercolatorAdapter() {
  scores_ = NULL;
}

PercolatorAdapter::PercolatorAdapter(Scores* scores) {
  scores_ = scores;
}

PercolatorAdapter::~PercolatorAdapter() {
  // TODO Auto-generated destructor stub
}

void PercolatorAdapter::setScores(Scores* scores) {
  scores_ = scores;
}

MatchCollection* PercolatorAdapter::convertFromPsms() {

  // Create new MatchCollection object that will be the converted Percolator Scores
  MatchCollection* converted = new MatchCollection();

  // Iterate over each ScoreHolder in Scores object
  for (
    vector<ScoreHolder>::iterator score_itr = scores_->begin();
    score_itr != scores_->end();
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

      //unsigned char length,     ///< The length of the peptide -in
      //FLOAT_T peptide_mass,       ///< The neutral mass of the peptide -in
      //Protein* parent_protein, ///< The parent_protein of this peptide -in
      //int start_idx ///< Start index of peptide in the protein sequence -in
    FLOAT_T peptide_mass = (psm->expMass * charge_state) - ((charge_state - 1) * MASS_PROTON);
    Peptide peptide(psm->getPeptideSequence().size(), peptide_mass, parent_protein, 0);

      //int               first_scan,         ///< number of the first scan -in
      //int               last_scan,          ///< number of the last scan -in
      //FLOAT_T           precursor_mz,       ///< m/z of the precursor
      //const std::vector<int>& possible_z,   ///< possible charge states
      //const char*       filename
    vector<int> zStates;
    zStates.push_back(charge_state);
    Spectrum spectrum(psm->scan, psm->scan, psm->expMass, zStates, "");

      //FLOAT_T neutral_mass,
      //int charge
    SpectrumZState zState(peptide_mass, charge_state);

      //Peptide* peptide, ///< the peptide for this match
      //Spectrum* spectrum, ///< the spectrum for this match
      //SpectrumZState& zstate, ///< the charge/mass of the spectrum
      //bool is_decoy);///< is the peptide a decoy or not
    Match* match = new Match(&peptide, &spectrum, zState, score_itr->isDecoy());
    match->setScore(PERCOLATOR_SCORE, score_itr->score);
    match->setScore(PERCOLATOR_QVALUE, psm->q);
    match->setScore(PERCOLATOR_PEP, psm->pep);

    converted->addMatch(match);
  }

  return converted;

}

MatchCollection* PercolatorAdapter::convertFromPsms(Scores* scores) {

  PercolatorAdapter* adapter = new PercolatorAdapter(scores);
  MatchCollection* collection = adapter->convertFromPsms();

  return collection;

}

MatchCollection* PercolatorAdapter::convertFromProteins() {
  // TODO
  return 0;
}

MatchCollection* PercolatorAdapter::convertFromProteins(Scores* scores) {
  PercolatorAdapter* adapter = new PercolatorAdapter(scores);
  MatchCollection* collection = adapter->convertFromProteins();

  return collection;
}

MatchCollection* PercolatorAdapter::convertFromPeptides() {
  // TODO
  return 0;
}

MatchCollection* PercolatorAdapter::convertFromPeptides(Scores* scores) {
  PercolatorAdapter* adapter = new PercolatorAdapter(scores);
  MatchCollection* collection = adapter->convertFromPeptides();

  return collection;
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
