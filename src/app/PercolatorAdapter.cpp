/**
 * \file PercolatorAdapter.cpp
 * $Revision$
 * \brief Converts Percolator result objects to Crux result objects.
 */

#include "PercolatorAdapter.h"
#include "build/src/percolator/src/DataSet.h"
#include "FeatureNames.h"

#include <map>

using namespace std;

/**
 * Constructor for PercolatorAdapter. 
 */
PercolatorAdapter::PercolatorAdapter() : Caller() {
  carp(CARP_DEBUG, "PercolatorAdapter::PercolatorAdapter");
  collection_ = new ProteinMatchCollection();
  decoy_collection_ = new ProteinMatchCollection();
}

/**
 * Destructor for PercolatorAdapter
 */
PercolatorAdapter::~PercolatorAdapter() {
  carp(CARP_DEBUG, "PercolatorAdapter::~PercolatorAdapter");

  // delete match collections created by this adapter
  int collectionsDeleted = 0;
  for (vector<MatchCollection*>::iterator iter = match_collections_made_.begin();
       iter != match_collections_made_.end();
       ++iter) {
    delete *iter;
    ++collectionsDeleted;
  }
  // delete proteins created by this adapter
  int proteinsDeleted = 0;
  for (vector<PostProcessProtein*>::iterator iter = proteins_made_.begin();
       iter != proteins_made_.end();
       ++iter) {
    delete *iter;
    ++proteinsDeleted;
  }

  //delete collection_;  //TODO find out why this crashes
  carp(CARP_DEBUG, "PercolatorAdapter::~PercolatorAdapter - done. %d "
       "MatchCollections deleted.", collectionsDeleted);
}

/**
 * Converts a set of Percolator scores into a Crux MatchCollection
 */
void PercolatorAdapter::psmScoresToMatchCollection(
  MatchCollection** match_collection,  ///< out parameter for targets
  MatchCollection** decoy_match_collection ///< out parameter for decoys
) {

  // Create new MatchCollection object that will be the converted Percolator Scores
  *match_collection = new MatchCollection();
  match_collections_made_.push_back(*match_collection);
  *decoy_match_collection = new MatchCollection();
  match_collections_made_.push_back(*decoy_match_collection);

  // Find out which feature is lnNumSP and get indices of charge state features
  Normalizer* normalizer = Normalizer::getNormalizer();
  double* normSub = normalizer->getSub();
  double* normDiv = normalizer->getDiv();
  FeatureNames& features = DataSet::getFeatureNames();
  string featureNames = features.getFeatureNames();
  vector<string> featureTokens;
  tokenize(featureNames, featureTokens);
  int lnNumSPIndex = -1, massIndex = -1;
  map<int, int> chargeStates; // index of feature -> charge
  for (int i = 0; i < featureTokens.size(); ++i) {
    string featureName = featureTokens[i];
    transform(featureName.begin(), featureName.end(),
              featureName.begin(), ::tolower);
    if (featureName == "lnnumsp") {
      lnNumSPIndex = i;
    } else if (featureName.find("charge") == 0) {
      size_t chargeNum = atoi(featureName.substr(6).c_str());
      chargeStates[i] = chargeNum;
    } else if (featureName == "mass") {
      massIndex = i;
    }
  }

  // Iterate over each ScoreHolder in Scores object
  for (vector<ScoreHolder>::iterator score_itr = fullset.begin();
       score_itr != fullset.end();
       score_itr++) {

    bool is_decoy = score_itr->isDecoy();

    PSMDescription* psm = score_itr->pPSM;
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
    int psm_file_idx;
    int psm_scan;
    int psm_charge;
    int psm_rank;
    parsePSMId(psm->id, 
      psm_file_idx,
      psm_scan,
      psm_charge,
      psm_rank);
    if (charge_state == -1) {
      // Failed, try to parse charge state from id
      charge_state = psm_charge;
      if (charge_state == -1) {
        carp_once(CARP_WARNING, "Could not determine charge state of PSM");
      }
    }
    Crux::Peptide* peptide = extractPeptide(psm, charge_state, is_decoy);

    FLOAT_T obsMass = (massIndex < 0) ?
      0 : psm->getFeatures()[massIndex] * normDiv[massIndex] + normSub[massIndex];
    SpectrumZState zState;
    zState.setSinglyChargedMass(obsMass, charge_state);
    // calcMass/expMass = singly charged mass
    Crux::Spectrum* spectrum = new Crux::Spectrum(
      psm->scan, psm->scan, zState.getMZ(), vector<int>(1, charge_state), ""
    );

    Crux::Match* match = new Crux::Match(peptide, spectrum, zState, is_decoy);
    match->setScore(PERCOLATOR_SCORE, score_itr->score);
    match->setScore(PERCOLATOR_QVALUE, psm->q);
    match->setScore(PERCOLATOR_PEP, psm->pep);

    match->setFileIndex(psm_file_idx);

    // Get matches/spectrum
    if (lnNumSPIndex < 0) {
      match->setLnExperimentSize(-1);
    } else {
      double lnNumSP = psm->getFeatures()[lnNumSPIndex]
        * normDiv[lnNumSPIndex] + normSub[lnNumSPIndex];
      match->setLnExperimentSize(lnNumSP);
    }

    if (!is_decoy) {
      (*match_collection)->addMatch(match);
    } else {
      (*decoy_match_collection)->addMatch(match);
    }
    match->setPostProcess(true); // so spectra get deleted when match does
    Crux::Match::freeMatch(match); // so match gets deleted when collection does
  }

  (*match_collection)->forceScoredBy(PERCOLATOR_SCORE);
  (*match_collection)->forceScoredBy(PERCOLATOR_QVALUE);
  (*match_collection)->forceScoredBy(PERCOLATOR_PEP);
  (*match_collection)->populateMatchRank(PERCOLATOR_SCORE);

  (*decoy_match_collection)->forceScoredBy(PERCOLATOR_SCORE);
  (*decoy_match_collection)->forceScoredBy(PERCOLATOR_QVALUE);
  (*decoy_match_collection)->forceScoredBy(PERCOLATOR_PEP);
  (*decoy_match_collection)->populateMatchRank(PERCOLATOR_SCORE);

  // sort by q-value
  (*match_collection)->sort(PERCOLATOR_QVALUE);
  (*decoy_match_collection)->sort(PERCOLATOR_QVALUE);

}

/**
 * Adds PSM scores from Percolator objects into a ProteinMatchCollection
 */
void PercolatorAdapter::addPsmScores() {
  MatchCollection* targets;
  MatchCollection* decoys;
  psmScoresToMatchCollection(&targets, &decoys);
  collection_->addMatches(targets);
  decoy_collection_->addMatches(decoys);
}

/**
 * Adds protein scores from Percolator objects into a ProteinMatchCollection
 */
void PercolatorAdapter::addProteinScores() {

  vector<ProteinMatch*> matches;
  vector<ProteinMatch*> decoy_matches;
  map<const string,Protein*> protein_scores = protEstimator->getProteins();
  
  for (map<const string,Protein*>::iterator score_iter = protein_scores.begin();
       score_iter != protein_scores.end();
       score_iter++) {
    if (score_iter->second == NULL) {
      continue;
    }
    // Set scores
    ProteinMatch* protein_match;
    if (!score_iter->second->getIsDecoy()) {
      protein_match = collection_->getProteinMatch(score_iter->second->getName());
      matches.push_back(protein_match);
    } else {
      protein_match = decoy_collection_->getProteinMatch(score_iter->second->getName());
      decoy_matches.push_back(protein_match);
    }
      protein_match->setScore(PERCOLATOR_SCORE, -log(score_iter->second->getP()));
      protein_match->setScore(PERCOLATOR_QVALUE, score_iter->second->getQ());
      protein_match->setScore(PERCOLATOR_PEP, score_iter->second->getPEP());
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
    ProteinMatch* match = *iter;
    match->setRank(PERCOLATOR_SCORE, ++cur_rank);
  }
  cur_rank = 0;
  for (vector<ProteinMatch*>::iterator iter = decoy_matches.begin();
       iter != decoy_matches.end();
       ++iter) {
    ProteinMatch* match = *iter;
    match->setRank(PERCOLATOR_SCORE, ++cur_rank);
  }

}

/**
 * Adds peptide scores from Percolator objects into a ProteinMatchCollection
 */
void PercolatorAdapter::addPeptideScores() {

  carp(CARP_DEBUG, "Setting peptide scores");

  // Iterate over each ScoreHolder in Scores object
  for (vector<ScoreHolder>::iterator score_itr = fullset.begin();
       score_itr != fullset.end();
       score_itr++) {

    PSMDescription* psm = score_itr->pPSM;
    string sequence;
    FLOAT_T peptide_mass;
    MODIFIED_AA_T* mod_seq = getModifiedAASequence(psm, sequence, peptide_mass);

    // Set scores
    PeptideMatch* peptide_match;
    if (!score_itr->isDecoy()) {
      peptide_match = collection_->getPeptideMatch(mod_seq);
    } else {
      peptide_match = decoy_collection_->getPeptideMatch(mod_seq);
    }
    if (peptide_match == NULL) {
      carp(CARP_FATAL, "Cannot find peptide %s %i",
                       psm->getFullPeptideSequence().c_str(), score_itr->isDecoy());
    }
    peptide_match->setScore(PERCOLATOR_SCORE, score_itr->score);
    peptide_match->setScore(PERCOLATOR_QVALUE, psm->q);
    peptide_match->setScore(PERCOLATOR_PEP, psm->pep);

    free(mod_seq);

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

void PercolatorAdapter::parsePSMId(
  const string& psm_id, ///< psm id to parse information from
  int& file_idx, ///< file index of psm
  int& scan, ///< scan number of psm
  int& charge, ///< charge of psm
  int& rank ///< rank of psm
) {
  // <target|decoy>_<fileindex>_<scan>_<charge>_<rank> OR
  // <filestem>_<scan>_<charge>_<rank>
  vector<string> tokens;
  tokenize(psm_id, tokens, '_');
  if (tokens.size() < 4) {
    carp(CARP_FATAL, "PSMID should be (((target|decoy)_fileidx)|filestem)_"
                     "scan_charge_rank, but was %s", psm_id.c_str());
  }
  if (tokens.size() == 5 && (tokens[0] == "target" || tokens[0] == "decoy")) {
    // Parse as <target|decoy>_<fileindex>_<scan>_<charge>_<rank>
    from_string<int>(file_idx, tokens[1]);
    from_string<int>(scan, tokens[2]);
    from_string<int>(charge, tokens[3]);
    from_string<int>(rank, tokens[4]);
  } else {
    // Try to parse as <filestem>_<scan>_<charge>_<rank>
    from_string<int>(rank, tokens.back());
    tokens.pop_back();
    from_string<int>(charge, tokens.back());
    tokens.pop_back();
    from_string<int>(scan, tokens.back());
    tokens.pop_back();
    stringstream ss;
    for (vector<string>::const_iterator i = tokens.begin(); i != tokens.end(); i++) {
      if (i != tokens.begin()) {
        ss << '_' << *i;
      } else {
        ss << *i;
      }
    }
    file_idx = Crux::Match::findFileIndex(ss.str(), true);
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

  string seq;
  FLOAT_T peptide_mass;
  
  MODIFIED_AA_T* mod_seq = getModifiedAASequence(psm, seq, peptide_mass);

  carp(CARP_DEBUG, "seq:%s",seq.c_str());
  char* mod_seq_str = modified_aa_string_to_string_with_masses(mod_seq, seq.length(), MOD_MASS_ONLY);
  carp(CARP_DEBUG, "mod seq:%s", mod_seq_str);
  free(mod_seq_str);
  
  string full_peptide(psm->getFullPeptide());
  carp(CARP_DEBUG, "full peptide:%s", full_peptide.c_str());
  carp(CARP_DEBUG, "=======================");
  string n_term = "";
  string c_term = "";
  if (!full_peptide.empty()) {
    n_term += full_peptide[0];
    c_term += full_peptide[full_peptide.length() - 1];
  }

  PostProcessProtein* parent_protein = new PostProcessProtein();
  proteins_made_.push_back(parent_protein);
  parent_protein->setId((*psm->proteinIds.begin()).c_str());
  int start_idx = parent_protein->findStart(seq, n_term, c_term);

  Crux::Peptide* peptide = new Crux::Peptide(seq.length(), peptide_mass, parent_protein, start_idx);

  // add other proteins
  bool skip_one = true;
  for (set<string>::iterator iter = psm->proteinIds.begin();
       iter != psm->proteinIds.end();
       ++iter) {
    if (skip_one) {
      skip_one = false;
      continue;
    }
    PostProcessProtein* secondary_protein = new PostProcessProtein();
    proteins_made_.push_back(secondary_protein);
    secondary_protein->setId(iter->c_str());
    int secondary_idx = secondary_protein->findStart(seq, n_term, c_term);
    peptide->addPeptideSrc(
      new PeptideSrc(NON_SPECIFIC_DIGEST, secondary_protein, secondary_idx)
    );
  }

  peptide->setModifiedAASequence(mod_seq, is_decoy);

  free(mod_seq);
  return peptide;
}

/**
 * \returns the modified and unmodified peptide sequence
 * for the psm
 */
MODIFIED_AA_T* PercolatorAdapter::getModifiedAASequence(
  PSMDescription* psm, ///< psm -in
  string& seq, ///< sequence -out
  FLOAT_T& peptide_mass ///< calculated mass of peptide with modifications -out
  ) {

  std::stringstream ss_seq;
  string perc_seq = psm->getFullPeptideSequence();
  size_t perc_seq_len = perc_seq.length();
  if (perc_seq_len >= 5 &&
      perc_seq[1] == '.' && perc_seq[perc_seq_len - 2] == '.') {
    // Trim off flanking AA if they exist
    perc_seq = perc_seq.substr(2, perc_seq_len - 4);
    perc_seq_len -= 4;
  }
  peptide_mass = 0.0;
  
  if (perc_seq.find("UNIMOD") != string::npos) {
    carp(CARP_FATAL, 
      "UNIMOD modifications currently not supported:%s", 
      perc_seq.c_str());
  }

  MODIFIED_AA_T* mod_seq = NULL;
  carp(CARP_DEBUG, "PercolatorAdapter::getModifiedAASequence(): seq:%s", perc_seq.c_str());
  int mod_len = convert_to_mod_aa_seq(perc_seq.c_str(), &mod_seq, MOD_MASS_ONLY);
  peptide_mass = get_mod_aa_seq_mass(mod_seq, get_mass_type_parameter("isotopic-mass"));
  char* cseq = modified_aa_to_unmodified_string(mod_seq, mod_len);
  seq = cseq;
  free(cseq);
  
  carp(CARP_DEBUG, "Peptide mass:%lf", peptide_mass);
  
  return mod_seq;
}

/** 
 * Executes the flow of the percolator process:
 * 1. reads in the input file
 * 2. trains the SVM
 * 3. calculate PSM probabilities
 * 4. (optional) calculate peptide probabilities
 * 5. (optional) calculate protein probabilities
 */
int PercolatorAdapter::run() {  

  time(&startTime);
  startClock = clock();
  if (VERB > 0) {
    cerr << extendedGreeter();
  }
  // populate tmp input file with cin information if option is enabled
  if(readStdIn){
    ofstream tmpInputFile;
    tmpInputFile.open(xmlInterface.getXmlInputFN().c_str());
    while(cin) {
      char buffer[1000];
      cin.getline(buffer, 1000);
      tmpInputFile << buffer << endl;
    }
    tmpInputFile.close();
  }
  
  // Reading input files (pin or temporary file)
  if(!readFiles()) {
    throw MyException("ERROR: Failed to read in file, check if the correct file-format was used.");
  }
  // Copy feature data to Scores object
  fillFeatureSets();
  
  // delete temporary file if reading from stdin
  if(readStdIn) {
    remove(xmlInterface.getXmlInputFN().c_str());
  }
  if(VERB > 2){
    std::cerr << "FeatureNames::getNumFeatures(): "<< FeatureNames::getNumFeatures() << endl;
  }
  int firstNumberOfPositives = crossValidation.preIterationSetup(fullset, pCheck, pNorm);
  if (VERB > 0) {
    cerr << "Estimating " << firstNumberOfPositives << " over q="
        << test_fdr << " in initial direction" << endl;
  }
  
  time_t procStart;
  clock_t procStartClock = clock();
  time(&procStart);
  double diff = difftime(procStart, startTime);
  if (VERB > 1) cerr << "Reading in data and feature calculation took "
      << ((double)(procStartClock - startClock)) / (double)CLOCKS_PER_SEC
      << " cpu seconds or " << diff << " seconds wall time" << endl;
  
  // Do the SVM training
  crossValidation.train(pNorm);
  crossValidation.postIterationProcessing(fullset, pCheck);
  // calculate psms level probabilities
  
  //PSM probabilities TDA or TDC
  calculatePSMProb(false, &fullset, procStart, procStartClock, diff, target_decoy_competition);
  addPsmScores();
  if (xmlInterface.getXmlOutputFN().size() > 0){
    xmlInterface.writeXML_PSMs(fullset);
  }
  
  // calculate unique peptides level probabilities WOTE
  if(reportUniquePeptides){
    calculatePSMProb(true, &fullset, procStart, procStartClock, diff, target_decoy_competition);
    addPeptideScores();
    if (xmlInterface.getXmlOutputFN().size() > 0){
      xmlInterface.writeXML_Peptides(fullset);
    }
  }
  // calculate protein level probabilities with FIDO
  if(ProteinProbEstimator::getCalcProteinLevelProb()){
    calculateProteinProbabilitiesFido();
    addProteinScores();
  }
  // write output to file
  xmlInterface.writeXML(fullset, protEstimator, call);  
  //return 1; don't know why percolator returns 1 here, but we'll return 0
  return 0;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

