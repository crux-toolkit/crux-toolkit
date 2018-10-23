/**
 * \file PercolatorAdapter.cpp
 * $Revision$
 * \brief Converts Percolator result objects to Crux result objects.
 */

#include "PercolatorAdapter.h"
#include "DataSet.h"
#include "util/AminoAcidUtil.h"
#include "util/StringUtils.h"
#include "FeatureNames.h"

#include <map>

using namespace std;

/**
 * Constructor for PercolatorAdapter. 
 */
PercolatorAdapter::PercolatorAdapter() : Caller() {
  collection_ = new ProteinMatchCollection();
  decoy_collection_ = new ProteinMatchCollection();
}

/**
 * Destructor for PercolatorAdapter
 */
PercolatorAdapter::~PercolatorAdapter() {
  // delete match collections created by this adapter
  int collectionsDeleted = 0;
  for (vector<MatchCollection*>::iterator iter = match_collections_made_.begin();
       iter != match_collections_made_.end();
       ++iter) {
    delete *iter;
    ++collectionsDeleted;
  }
  // delete proteins created by this adapter
  for (vector<PostProcessProtein*>::iterator iter = proteins_made_.begin();
       iter != proteins_made_.end();
       ++iter) {
    delete *iter;
  }

  deleteCollections();
  carp(CARP_DEBUG, "PercolatorAdapter::~PercolatorAdapter - done. %d "
       "MatchCollections deleted.", collectionsDeleted);
}

void PercolatorAdapter::deleteCollections() {
  if (collection_) {
    delete collection_;
    collection_ = NULL;
  }
  if (decoy_collection_) {
    delete decoy_collection_;
    decoy_collection_ = NULL;
  }
}

/**
 * Adds PSM scores from Percolator objects into a ProteinMatchCollection
 */
void PercolatorAdapter::processPsmScores(Scores& allScores) {
  if (collection_ == NULL || decoy_collection_ == NULL) {
    return;
  }

  // Create new MatchCollection object that will be the converted Percolator Scores
  MatchCollection* targets = new MatchCollection();
  MatchCollection* decoys = new MatchCollection();
  match_collections_made_.push_back(targets);
  match_collections_made_.push_back(decoys);

  // Find out which feature is lnNumSP and get indices of charge state features
  bool lnNumDSP = false;
  int lnNumSPIndex = findFeatureIndex("lnnumsp");
  if (lnNumSPIndex == -1) {
    if ((lnNumSPIndex = findFeatureIndex("lnnumdsp")) != -1) {
      lnNumDSP = true;
    }
  }

  map<int, int> chargeStates = mapChargeFeatures();

  Normalizer* normalizer = Normalizer::getNormalizer();
  double* normSub = normalizer->getSub();
  double* normDiv = normalizer->getDiv();

  // Iterate over each ScoreHolder in Scores object
  for (vector<ScoreHolder>::iterator score_itr = allScores.begin();
       score_itr != allScores.end();
       score_itr++) {
    bool is_decoy = score_itr->isDecoy();

    PSMDescription* psm = score_itr->pPSM;

    int psm_file_idx = -1, psm_charge;
    parsePSMId(psm->id_, psm_file_idx, psm_charge);

    // Try to look up charge state in map
    int charge_state = -1;
    for (map<int, int>::const_iterator i = chargeStates.begin();
         i != chargeStates.end();
         ++i) {
      if (psm->features[i->first] > 0) {
        charge_state = i->second;
        break;
      }
    }

    if (charge_state == -1) {
      carp_once(CARP_WARNING, "Could not determine charge state of PSM");
    }

    Crux::Peptide* peptide = extractPeptide(psm, charge_state, is_decoy);
    if (peptide == NULL) {
      deleteCollections();
      return;
    }

    SpectrumZState zState;
    zState.setSinglyChargedMass(psm->expMass, charge_state);
    // calcMass/expMass = singly charged mass
    Crux::Spectrum* spectrum = new Crux::Spectrum(
      psm->scan, psm->scan, zState.getMZ(), vector<int>(1, charge_state),
      Crux::Match::getFilePath(psm_file_idx));

    Crux::Match* match = new Crux::Match(peptide, spectrum, zState, is_decoy);
    match->setScore(PERCOLATOR_SCORE, score_itr->score);
    match->setScore(PERCOLATOR_QVALUE, score_itr->q);
    match->setScore(PERCOLATOR_PEP, score_itr->pep);

    match->setFileIndex(psm_file_idx);

    // Get matches/spectrum
    if (lnNumSPIndex < 0) {
      match->setLnExperimentSize(-1);
    } else {
      double lnNumSP = unnormalize(psm, lnNumSPIndex, normDiv, normSub);
      match->setLnExperimentSize(lnNumSP);
    }

    if (!is_decoy) {
      targets->addMatch(match);
    } else {
      decoys->addMatch(match);
    }
    match->setPostProcess(true); // so spectra get deleted when match does
    Crux::Match::freeMatch(match); // so match gets deleted when collection does
  }

  if (lnNumSPIndex >= 0 && lnNumDSP) {
    targets->setHasDistinctMatches(true);
    decoys->setHasDistinctMatches(true);
  }

  targets->setScoredType(PERCOLATOR_SCORE, true);
  targets->setScoredType(PERCOLATOR_QVALUE, true);
  targets->setScoredType(PERCOLATOR_PEP, true);
  targets->populateMatchRank(PERCOLATOR_SCORE);

  decoys->setScoredType(PERCOLATOR_SCORE, true);
  decoys->setScoredType(PERCOLATOR_QVALUE, true);
  decoys->setScoredType(PERCOLATOR_PEP, true);
  decoys->populateMatchRank(PERCOLATOR_SCORE);

  // sort by q-value
  targets->sort(PERCOLATOR_QVALUE);
  decoys->sort(PERCOLATOR_QVALUE);

  collection_->addMatches(targets);
  decoy_collection_->addMatches(decoys);
}

/**
 * Adds peptide scores from Percolator objects into a ProteinMatchCollection
 */
void PercolatorAdapter::processPeptideScores(Scores& allScores) {
  if (collection_ == NULL || decoy_collection_ == NULL) {
    return;
  }

  carp(CARP_DEBUG, "Setting peptide scores");

  // Iterate over each ScoreHolder in Scores object
  for (vector<ScoreHolder>::iterator score_itr = allScores.begin();
       score_itr != allScores.end();
       score_itr++) {

    PSMDescription* psm = score_itr->pPSM;
    string sequence;
    vector<Crux::Modification> mods;
    Crux::Modification::FromSeq(psm->getFullPeptideSequence(), &sequence, &mods);
    string peptideId = Crux::Peptide::getId(sequence, mods);

    // Set scores
    PeptideMatch* peptide_match = !score_itr->isDecoy()
      ? collection_->getPeptideMatch(peptideId)
      : decoy_collection_->getPeptideMatch(peptideId);
    if (peptide_match == NULL) {
      deleteCollections();
      return;
    }
    peptide_match->setScore(PERCOLATOR_SCORE, score_itr->score);
    peptide_match->setScore(PERCOLATOR_QVALUE, score_itr->q);
    peptide_match->setScore(PERCOLATOR_PEP, score_itr->pep);
  }
}
  
/**
 * Adds protein scores from Percolator objects into a ProteinMatchCollection
 */
void PercolatorAdapter::processProteinScores(ProteinProbEstimator* protEstimator) {
  if (collection_ == NULL || decoy_collection_ == NULL) {
    return;
  }

  vector<ProteinMatch*> matches;
  vector<ProteinMatch*> decoy_matches;
  const vector<ProteinScoreHolder>& protein_scores = protEstimator->getProteins();
  
  for (vector<ProteinScoreHolder>::const_iterator score_iter = protein_scores.begin();
       score_iter != protein_scores.end();
       score_iter++) {
    // Set scores
    bool decoy = score_iter->isDecoy();
    ProteinMatch* protein_match = !decoy
      ? collection_->getProteinMatch(score_iter->getName())
      : decoy_collection_->getProteinMatch(score_iter->getName());
    if (!protein_match) {
      PostProcessProtein* protein = makeProtein(score_iter->getName());
      protein_match = !decoy
        ? collection_->getProteinMatch(protein, true)
        : decoy_collection_->getProteinMatch(protein, true);
    }
    if (!decoy) {
      matches.push_back(protein_match);
    } else {
      decoy_matches.push_back(protein_match);
    }
    protein_match->setScore(PERCOLATOR_SCORE, -log(score_iter->getP()));
    protein_match->setScore(PERCOLATOR_QVALUE, score_iter->getQ());
    protein_match->setScore(PERCOLATOR_PEP, score_iter->getPEP());
  }

  // set percolator score ranks
  std::sort(matches.begin(), matches.end(),
            PercolatorAdapter::comparePercolatorScores);
  std::sort(decoy_matches.begin(), decoy_matches.end(),
            PercolatorAdapter::comparePercolatorScores);
  int cur_rank = 0;
  for (vector<ProteinMatch*>::iterator iter = matches.begin();
       iter != matches.end();
       ++iter) {
    (*iter)->setRank(PERCOLATOR_SCORE, ++cur_rank);
  }
  cur_rank = 0;
  for (vector<ProteinMatch*>::iterator iter = decoy_matches.begin();
       iter != decoy_matches.end();
       ++iter) {
    (*iter)->setRank(PERCOLATOR_SCORE, ++cur_rank);
  }
}

/*
 *\returns the ProteinMatchCollection, to be called after Caller::run() is finished
 */
ProteinMatchCollection* PercolatorAdapter::getProteinMatchCollection() {
  return collection_;
}

/*
 *\returns the decoy ProteinMatchCollection, to be called after Caller::run() is finished
 */
ProteinMatchCollection* PercolatorAdapter::getDecoyProteinMatchCollection() {
  return decoy_collection_;
}

bool PercolatorAdapter::parsePSMId(
  const string& psm_id, ///< psm id to parse information from
  int& file_idx, ///< file index of psm
  int& charge ///< charge of psm
) {
  // <target|decoy>_<fileindex>_<scan>_<charge>_<rank> OR
  // <filestem>_<scan>_<charge>_<rank>
  vector<string> tokens = StringUtils::Split(psm_id, '_');
  if (tokens.size() < 4) {
    return false;
  }
  try {
    if (tokens.size() == 5 && (tokens[0] == "target" || tokens[0] == "decoy")) {
      // Parse as <target|decoy>_<fileindex>_<scan>_<charge>_<rank>
      file_idx = StringUtils::FromString<int>(tokens[1]);
      // tokens[2] = scan
      charge = StringUtils::FromString<int>(tokens[3]);
      // tokens[4] = rank
    } else {
      // Try to parse as <filestem>_<scan>_<charge>_<rank>
      tokens.pop_back(); // pop rank
      charge = StringUtils::FromString<int>(tokens.back());
      tokens.pop_back(); // pop charge
      string scanStr = tokens.back();
      tokens.pop_back(); // pop scan
      // just <filestem> remaining
      if (tokens.size() >= 4 &&
          tokens[tokens.size() - 3] == "SII" &&
          tokens[tokens.size() - 2] == scanStr) {
        // handle msgf+ case: <filestem>_SII_<scan>_<?>
        tokens.erase(tokens.end() - 3, tokens.end());
      }
      stringstream ss;
      for (vector<string>::const_iterator i = tokens.begin(); i != tokens.end(); i++) {
        if (i != tokens.begin()) {
          ss << '_' << *i;
        } else {
          ss << *i;
        }
      }
      file_idx = Crux::Match::addUniqueFilePath(ss.str(), true);
    }
    return true;
  } catch (...) {
    return false;
  }
}

/**
 * Compares two AbstractMatches by Percolator score
 */
bool PercolatorAdapter::comparePercolatorScores(
  AbstractMatch* lhs, ///< first match with Percolator scores to compare
  AbstractMatch* rhs ///< second match with Percolator scores to compare
) {
  if (!lhs->hasScore(PERCOLATOR_SCORE) || !rhs->hasScore(PERCOLATOR_SCORE)) {
    carp(CARP_FATAL, "Could not compare matches by Percolator score.");
  }
  return lhs->getScore(PERCOLATOR_SCORE) < rhs->getScore(PERCOLATOR_SCORE);
}

/**
* \returns a Crux peptide from the PSM
*/
Crux::Peptide* PercolatorAdapter::extractPeptide(
  PSMDescription* psm, ///< psm
  int charge_state, ///< charge state
  bool is_decoy ///< is psm a decoy?
) {
  string sequence = psm->getFullPeptideSequence();
  string n_term = "", c_term = "";
  // Get flanking AA if they exist
  if (sequence.length() >= 5 &&
      sequence[1] == '.' && sequence[sequence.length() - 2] == '.') {
    n_term += sequence[0];
    c_term += sequence[sequence.length() - 1];
  }

  string unmodifiedSeq;
  vector<Crux::Modification> mods;
  Crux::Modification::FromSeq(sequence, &unmodifiedSeq, &mods);

  // add proteins
  Crux::Peptide* peptide = NULL;
  for (vector<string>::const_iterator i = psm->proteinIds.begin();
       i != psm->proteinIds.end();
       ++i) {
    PostProcessProtein* protein = makeProtein(i->c_str());
    int start_idx = protein->findStart(unmodifiedSeq, n_term, c_term);
    if (peptide == NULL) {
      peptide = new Crux::Peptide(unmodifiedSeq.length(), protein, start_idx);
      peptide->setMods(mods);
      // TODO peptide->setModifiedAASequence(mod_seq, is_decoy);
    } else {
      peptide->addPeptideSrc(new PeptideSrc(NON_SPECIFIC_DIGEST, protein, start_idx));
    }
  }

  return peptide;
}

PostProcessProtein* PercolatorAdapter::makeProtein(const string& name) {
  PostProcessProtein* protein = new PostProcessProtein();
  proteins_made_.push_back(protein);
  protein->setId(name);
  return protein;
}

// Finds the index of the given feature name (case insensitive).
int PercolatorAdapter::findFeatureIndex(string feature) {
  feature = StringUtils::ToLower(feature);
  vector<string> features = StringUtils::Split(DataSet::getFeatureNames().getFeatureNames(), '\t');
  for (int i = 0; i < features.size(); ++i) {
    if (StringUtils::ToLower(features[i]) == feature) {
      return i;
    }
  }
  return -1;
}

// Generate a map where the keys are feature indices and their values are the
// charge state of that feature.
map<int, int> PercolatorAdapter::mapChargeFeatures() {
  map<int, int> charges;
  vector<string> features = StringUtils::Split(DataSet::getFeatureNames().getFeatureNames(), '\t');
  for (int i = 0; i < features.size(); ++i) {
    if (StringUtils::StartsWith(StringUtils::ToLower(features[i]), "charge")) {
      size_t charge = StringUtils::FromString<size_t>(features[i].substr(6));
      charges[i] = charge;
    }
  }
  return charges;
}

// Unnormalize a PSM feature to obtain its original value.
// normDiv and normSub are optional; if they are NULL, they will be
// automatically retrieved.
double PercolatorAdapter::unnormalize(
  const PSMDescription* psm,
  int featureIndex,
  double* normDiv,
  double* normSub
) {
  if (normDiv == NULL || normSub == NULL) {
    Normalizer* normalizer = Normalizer::getNormalizer();
    normDiv = normalizer->getDiv();
    normSub = normalizer->getSub();
  }
  return psm->features[featureIndex] * normDiv[featureIndex] + normSub[featureIndex];
}

void PercolatorAdapter::printScores(Scores* scores, int label, ostream& os) {
  std::vector<ScoreHolder>::iterator scoreIt = scores->begin();

  bool lnNumDSP = false;
  int lnNumSPIdx = PercolatorAdapter::findFeatureIndex("lnnumsp");
  if (lnNumSPIdx == -1) {
    if ((lnNumSPIdx = PercolatorAdapter::findFeatureIndex("lnnumdsp")) != -1) {
      lnNumDSP = true;
    }
  }

  os //<< "file\t"
     << "file_idx\t"
     << "scan\t"
     << "charge\t"
     << "spectrum precursor m/z\t"
     << "spectrum neutral mass\t"
     << "peptide mass\t"
     << "percolator score\t"
     << "percolator q-value\t"
     << "percolator PEP\t"
     << (lnNumDSP ? "distinct matches/spectrum" : "total matches/spectrum") << '\t'
     << "sequence\t"
     //<< "modifications\t"
     << "protein id\t"
     << "flanking aa\n";

  Normalizer* normalizer = Normalizer::getNormalizer();
  double* nSub = normalizer->getSub();
  double* nDiv = normalizer->getDiv();

  for ( ; scoreIt != scores->end(); ++scoreIt) {
    if (scoreIt->label != label) {
      continue;
    }
    int fileIdx, charge;
    if (!parsePSMId(scoreIt->pPSM->getId(), fileIdx, charge)) {
      fileIdx = -1;
      charge = -1;
    }

    std::string flankingStr = "XX";
    std::string seq = scoreIt->pPSM->peptide;
    if (seq.length() >= 5 && seq[1] == '.' && seq[seq.length() - 2] == '.') {
      flankingStr[0] = seq[0];
      flankingStr[1] = seq[seq.length() - 1];
      seq = seq.substr(2, seq.length() - 4);
    }

    double neutralMass = scoreIt->pPSM->expMass - MASS_PROTON;
    double peptideMass = MASS_H2O_MONO; // Reported as 0 if a problem occurs
    for (size_t i = 0; i < seq.length(); i++) {
      if (seq[i] == '[') {
        double modMass = 0.0;
        size_t j = seq.find(']', ++i);
        if (j == string::npos || !StringUtils::TryFromString(seq.substr(i, j - i), &modMass)) {
          peptideMass = 0;
          break;
        }
        peptideMass += modMass;
        i = j;
      } else {
        try {
          peptideMass += AminoAcidUtil::GetMass(seq[i], true);
        } catch (...) {
          peptideMass = 0;
          break;
        }
      }
    }

    int precision = Params::GetInt("precision");
    int massPrecision = Params::GetInt("mass-precision");
    os //<< "" << '\t' // file
       << fileIdx << '\t'
       << scoreIt->pPSM->scan << '\t'
       << charge << '\t'
       << StringUtils::ToString((charge > 0 ? neutralMass/charge + MASS_PROTON : 0), massPrecision) << '\t'
       << StringUtils::ToString(neutralMass, massPrecision) << '\t'
       << StringUtils::ToString(peptideMass, massPrecision) << '\t'
       << StringUtils::ToString(scoreIt->score, precision) << '\t'
       << StringUtils::ToString(scoreIt->q, precision, false) << '\t'
       << StringUtils::ToString(scoreIt->pep, precision, false) << '\t'
       << (lnNumSPIdx >= 0 ? exp(PercolatorAdapter::unnormalize(scoreIt->pPSM, lnNumSPIdx, nDiv, nSub)) : 0) << '\t'
       << seq << '\t'
       //<< "" << '\t' // mods
       << StringUtils::Join(scoreIt->pPSM->proteinIds, ',') << '\t'
       << flankingStr << std::endl;
  }
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

