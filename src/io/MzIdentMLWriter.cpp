/**
 * \file MzIdentMLWriter.cpp
 * \brief Writes search results in the MzIdentML (mzid) format.
 */
#include <iostream>
#include <fstream>
#include <string>
#include "MzIdentMLWriter.h"
#include "util/crux-utils.h"
#include "util/MathUtil.h"
#include "util/Params.h"
#include "util/StringUtils.h"
#include "model/MatchCollection.h"
#include "model/Modification.h"
#include "pwiz/data/identdata/Serializer_mzid.hpp"
#include "model/ProteinMatch.h"
#include "model/ProteinMatchCollection.h"
#include "model/PeptideMatch.h"
#include "model/SpectrumMatch.h"

using namespace std;
using namespace Crux;
using namespace pwiz;
using namespace identdata;

#define calculateMassToCharge(peptide_mass, charge) (FLOAT_T) ((peptide_mass + (charge*MASS_PROTON))/charge)

MzIdentMLWriter::MzIdentMLWriter() : PSMWriter() {
  //data_ = NULL;
  fout_ = NULL;
  sir_idx_ = 0;
  sil_idx_ = 0;
  sii_idx_ = 0;
  peptide_idx_ = 0;
  peptide_evidence_idx_ = 0;
  dbs_idx_ = 0;
  pag_idx_ = 0;
  pdh_idx_ = 0;
}

MzIdentMLWriter::~MzIdentMLWriter() {
  closeFile();
}

void MzIdentMLWriter::openFile(
  const char* filename,
  bool overwrite) {
  fout_ = create_stream_in_path(filename, NULL, overwrite);  
  mzid_ = IdentDataPtr(new IdentData());
}

void MzIdentMLWriter::openFile(
  const string& filename,
  bool overwrite) {
  openFile(filename.c_str(), overwrite);
}

void MzIdentMLWriter::openFile(
  CruxApplication* application,
  string filename,
  MATCH_FILE_TYPE type) {
  openFile(filename.c_str(), Params::GetBool("overwrite"));
}

/**
 * Close the file, if open.
 */
void MzIdentMLWriter::closeFile() {
  if (mzid_ != NULL && fout_ != NULL) {
    // Add modification information to AnalysisProtocolCollection -> ModificationParams
    AnalysisProtocolCollection& apc = mzid_->analysisProtocolCollection;
    SpectrumIdentificationProtocolPtr sipp(new SpectrumIdentificationProtocol());
    apc.spectrumIdentificationProtocol.push_back(sipp);

    vector<const ModificationDefinition*> mods = ModificationDefinition::AllMods();
    for (vector<const ModificationDefinition*>::const_iterator i = mods.begin();
         i != mods.end();
         i++) {
      SearchModificationPtr smp(new SearchModification());
      smp->fixedMod = (*i)->Static();
      smp->massDelta = MathUtil::Round((*i)->DeltaMass(), Params::GetInt("mod-precision"));
      smp->residues = vector<char>((*i)->AminoAcids().begin(), (*i)->AminoAcids().end());
      switch ((*i)->Position()) {
        case PEPTIDE_N:
          smp->specificityRules = CVParam(MS_modification_specificity_peptide_N_term);
          break;
        case PEPTIDE_C:
          smp->specificityRules = CVParam(MS_modification_specificity_peptide_C_term);
          break;
        case PROTEIN_N:
          smp->specificityRules = CVParam(MS_modification_specificity_protein_N_term);
          break;
        case PROTEIN_C:
          smp->specificityRules = CVParam(MS_modification_specificity_protein_C_term);
          break;
      }
      sipp->modificationParams.push_back(smp);
    }

    Serializer_mzIdentML serializer;
    serializer.write(*fout_, *mzid_);
    fout_->close();
  }
  delete fout_;
  fout_ = NULL;
}

/**
 * \returns the MzIdentML Peptide object, creating it if
 * it doesn't exist already
 */
PeptidePtr MzIdentMLWriter::getPeptide(
  Crux::Peptide* peptide ///< Peptide -in
  ) {
  char* seqTmp = peptide->getSequence();
  string sequence(seqTmp);
  free(seqTmp);

  int modPrecision = Params::GetInt("mod-precision");

  vector<Crux::Modification> mods = peptide->getVarMods();
  for (vector<PeptidePtr>::const_iterator i = mzid_->sequenceCollection.peptides.begin();
       i != mzid_->sequenceCollection.peptides.end();
       ++i) {
    if ((*i)->peptideSequence == sequence && (*i)->modification.size() == mods.size()) {
      if (mods.empty()) {
        return *i;
      }
      bool match = true;
      for (vector<ModificationPtr>::const_iterator j = (*i)->modification.begin();
           j != (*i)->modification.end();
           ++j) {
        bool matchCur = false;
        for (vector<Crux::Modification>::const_iterator k = mods.begin();
             k != mods.end();
             k++) {
          if ((*j)->location == k->Index() &&
              MathUtil::AlmostEqual((*j)->monoisotopicMassDelta, k->DeltaMass(), modPrecision)) {
            matchCur = true;
            break;
          }
        }
        if (!matchCur) {
          match = false;
          break;
        }
      }
      //all modifications match, return peptide.
      if (match) {
        return *i;
      }
    }
  }

  //Okay, we didn't find a match, so create a new peptide object.
  PeptidePtr peptide_p(new pwiz::identdata::Peptide(
    "PEP_" + StringUtils::ToString(peptide_idx_++)));
  //peptide_p->id = sequence_str;
  peptide_p->peptideSequence = sequence;

  //add the modifications.
  for (vector<Crux::Modification>::const_iterator i = mods.begin(); i != mods.end(); i++) {
    ModificationPtr mod_p(new pwiz::identdata::Modification());
    switch (i->Position()) {
      case PEPTIDE_N:
      case PROTEIN_N:
        mod_p->location = 0;
        break;
      case PEPTIDE_C:
      case PROTEIN_C:
        mod_p->location = i->Index() + 2;
        break;
      default:
        mod_p->location = i->Index() + 1;
        break;
    }
    mod_p->monoisotopicMassDelta = MathUtil::Round(i->DeltaMass(), modPrecision);
    mod_p->residues.push_back(sequence[i->Index()]);
    mod_p->set(MS_unknown_modification);
    peptide_p->modification.push_back(mod_p);
  }

  mzid_->sequenceCollection.peptides.push_back(peptide_p);
  return peptide_p;
}

DBSequencePtr MzIdentMLWriter::getDBSequence(std::string& protein_id) {
  vector<DBSequencePtr>::iterator dbs_iter;
  for (dbs_iter = mzid_->sequenceCollection.dbSequences.begin();
       dbs_iter != mzid_->sequenceCollection.dbSequences.end();
       ++dbs_iter) {
    DBSequencePtr dbs_ptr = *dbs_iter;
    if (protein_id == dbs_ptr->accession) {
      return dbs_ptr;
    }
  }

  //I don't know what to do here, there is a problem if we don't have the
  //protein database, we can't really assign the sequence.
  //we can fake it as below, but then we have multiple DBSequences for
  //the same protein.  Here we might not have the full sequence.
  DBSequencePtr dbs_ptr(new DBSequence("DBS_"+boost::lexical_cast<string>(dbs_idx_++)));
  dbs_ptr->accession = protein_id;
  
  return (dbs_ptr);
}

/**
 * \returns DBSequence for the protein source.  If it doesn't exist, 
 * then first create the object in the mzid object
 */
DBSequencePtr MzIdentMLWriter::getDBSequence(
  Crux::Peptide* peptide,  ///< peptide -in
  PeptideSrc* src ///< Source of the peptide -in
  ) {
  string protein_id = src->getParentProtein()->getIdPointer();
  bool is_post_process = src->getParentProtein()->isPostProcess();
  vector<DBSequencePtr>::iterator dbs_iter;
 
  string sequence_str;
  if (is_post_process) {
    char* seq = peptide->getSequence();
    sequence_str = seq;
    free(seq);
  }

  for (dbs_iter = mzid_->sequenceCollection.dbSequences.begin();
       dbs_iter != mzid_->sequenceCollection.dbSequences.end();
       ++dbs_iter) {
    DBSequencePtr dbs_ptr = *dbs_iter;
    if (protein_id == dbs_ptr->accession) {
      if (is_post_process) {
        //sequence str should be the peptide
        if (sequence_str == dbs_ptr->seq) {
          return dbs_ptr;
        }
      } else {
        //we have the full sequence.
        return dbs_ptr;
      }
    }
  }

  DBSequencePtr dbs_ptr(new DBSequence("DBS_"+boost::lexical_cast<string>(dbs_idx_++)));
  dbs_ptr->accession = protein_id;
  if (is_post_process) {
    dbs_ptr->length = sequence_str.length();
    dbs_ptr->seq = sequence_str;
  } else {
    dbs_ptr->length = src->getParentProtein()->getLength();
    dbs_ptr->seq = src->getParentProtein()->getSequencePointer();
    //TODO add description
  }
  
  mzid_->sequenceCollection.dbSequences.push_back(dbs_ptr);

  return dbs_ptr;
}

PeptideEvidencePtr MzIdentMLWriter::getPeptideEvidence(
  Crux::Peptide* peptide,
  bool is_decoy,
  string& protein_id) {

  char* seq = peptide->getSequence();
  string sequence_str = seq;
  free(seq);

  //Is there already a peptide evidence ptr?
  for (vector<PeptideEvidencePtr>::iterator i = mzid_->sequenceCollection.peptideEvidence.begin();
       i != mzid_->sequenceCollection.peptideEvidence.end();
       ++i) {
    PeptideEvidencePtr pe_ptr = *i;
    if (pe_ptr->peptidePtr->peptideSequence == sequence_str &&
        pe_ptr->dbSequencePtr->accession == protein_id) {
      return pe_ptr;
    }
  }

  carp(CARP_FATAL, "Couldn't find %s in %s", sequence_str.c_str(), protein_id.c_str());
  return PeptideEvidencePtr(); // Avoid compiler warning.
}

/**
 * \returns PeptideEvidence for the peptide and src.
 * creates it if it doesn't exist
 */
PeptideEvidencePtr MzIdentMLWriter::getPeptideEvidence(
  Crux::Peptide* peptide, ///< peptide -in
  bool is_decoy,  ///< is this peptide a decoy? -in
  PeptideSrc* src ///< where to peptide comes from -in
  ) {

  char* seq = peptide->getSequence();
  string sequence_str = seq;
  free(seq);

  string protein_id = src->getParentProtein()->getId();

  //Is there already a peptide evidence ptr?
  vector<PeptideEvidencePtr>::iterator pe_iter;

  for (vector<PeptideEvidencePtr>::iterator i = mzid_->sequenceCollection.peptideEvidence.begin();
       i != mzid_->sequenceCollection.peptideEvidence.end();
       ++i) {
    PeptideEvidencePtr pe_ptr = *i;
    if (pe_ptr->peptidePtr->peptideSequence == sequence_str &&
        pe_ptr->dbSequencePtr->accession == protein_id) {
      return pe_ptr;
    }
  }

  DBSequencePtr dbs_ptr = getDBSequence(peptide, src);
  PeptideEvidencePtr pe_ptr(new PeptideEvidence(
    "PE_" + StringUtils::ToString(peptide_evidence_idx_++)));

  if (src->getParentProtein()->isPostProcess()) {
    pe_ptr->start = 0;
    pe_ptr->end = peptide->getLength();
  } else {
    pe_ptr->start = src->getStartIdx();
    pe_ptr->end = src->getStartIdx()+peptide->getLength();
  }
  pe_ptr->dbSequencePtr = dbs_ptr;
  pe_ptr->peptidePtr = getPeptide(peptide);
  pe_ptr->isDecoy = is_decoy; 

  mzid_->sequenceCollection.peptideEvidence.push_back(pe_ptr);
  return pe_ptr;  
}

/**
 * \returns the SpectrumIdentificationList, creating it if 
 * it doesn't exist yet.
 */
SpectrumIdentificationListPtr MzIdentMLWriter::getSpectrumIdentificationList() {
  SpectrumIdentificationListPtr silp;

  if (mzid_->dataCollection.analysisData.spectrumIdentificationList.size() > 0) {
    silp = mzid_->dataCollection.analysisData.spectrumIdentificationList.back();
  } else {
    silp = SpectrumIdentificationListPtr(new SpectrumIdentificationList(
      "SIL_" + StringUtils::ToString(sil_idx_++)));
    mzid_->dataCollection.analysisData.spectrumIdentificationList.push_back(silp);   
  }
  return silp;
}

ProteinDetectionListPtr MzIdentMLWriter::getProteinIdentificationList() {
  ProteinDetectionListPtr pdlp;
  if (mzid_ -> dataCollection.analysisData.proteinDetectionListPtr != NULL) {
    pdlp = mzid_->dataCollection.analysisData.proteinDetectionListPtr;
  } else {
    pdlp = ProteinDetectionListPtr(
      new ProteinDetectionList("PDL_1"));
    mzid_->dataCollection.analysisData.proteinDetectionListPtr = pdlp;
  }
  return pdlp;
}

ProteinAmbiguityGroupPtr MzIdentMLWriter::getProteinAmbiguityGroup(
  std::string& protein_id
) {
  ProteinDetectionListPtr pdlp = getProteinIdentificationList();
  vector<ProteinAmbiguityGroupPtr>& pag_vec = pdlp->proteinAmbiguityGroup;

  for (size_t pa_idx = 0; pa_idx < pag_vec.size() ;pa_idx++) {
    std::vector<ProteinDetectionHypothesisPtr>& pdh_vec = pag_vec[pa_idx]->proteinDetectionHypothesis;  
    for (size_t pdh_idx = 0 ; pdh_idx < pdh_vec.size() ; pdh_idx++) {
      if (pdh_vec[pdh_idx]->dbSequencePtr->accession == protein_id) {

        return pag_vec[pa_idx]; 
      }
    }
  }

  //not found, create it
  ProteinAmbiguityGroupPtr pagp = ProteinAmbiguityGroupPtr(
    new ProteinAmbiguityGroup("PAG_"+boost::lexical_cast<string>(pag_idx_++)));

  ProteinDetectionHypothesisPtr pdhp = getProteinDetectionHypothesis(pagp, protein_id);
  pag_vec.push_back(pagp);
  return pagp;
}
  
ProteinDetectionHypothesisPtr MzIdentMLWriter::getProteinDetectionHypothesis(
  ProteinAmbiguityGroupPtr pagp,
  std::string& protein_id
  ) {
  std::vector<ProteinDetectionHypothesisPtr>& pdh_vec = pagp->proteinDetectionHypothesis;
  
  for (size_t pdh_idx = 0; pdh_idx < pdh_vec.size() ; pdh_idx++) {
    if (pdh_vec[pdh_idx]->dbSequencePtr->accession == protein_id) {
      return pdh_vec[pdh_idx];
    }
  }
 
  ProteinDetectionHypothesisPtr pdhp = ProteinDetectionHypothesisPtr(
    new ProteinDetectionHypothesis("PDH_"+boost::lexical_cast<string>(pdh_idx_++)));

  pdhp->dbSequencePtr = getDBSequence(protein_id);
  pdh_vec.push_back(pdhp);
  return pdhp;
}

ProteinDetectionHypothesisPtr MzIdentMLWriter::getProteinDetectionHypothesis(
  string& protein_id 
  ) {
  ProteinAmbiguityGroupPtr pagp = getProteinAmbiguityGroup(protein_id);
  return getProteinDetectionHypothesis(pagp, protein_id);
}

PeptideHypothesis& MzIdentMLWriter::getPeptideHypothesis(
  ProteinMatch* protein_match,
  PeptideMatch* peptide_match) {

  string protein_id = protein_match->getId();
  PeptideSrc* src = peptide_match->getSrc(protein_match);
  ProteinDetectionHypothesisPtr pdhp = getProteinDetectionHypothesis(protein_id);
  PeptideEvidencePtr pe_ptr = getPeptideEvidence(peptide_match->getPeptide(), true, src);

  for (size_t ph_idx = 0;ph_idx < pdhp->peptideHypothesis.size();ph_idx++) {
    if (pdhp->peptideHypothesis[ph_idx].peptideEvidencePtr == pe_ptr) {
      return pdhp->peptideHypothesis[ph_idx];
    }
  }

  PeptideHypothesis pep_hyp;
  pdhp->peptideHypothesis.push_back(pep_hyp);
  pdhp->peptideHypothesis.back().peptideEvidencePtr = pe_ptr;

  return pdhp->peptideHypothesis.back();  
}

/**
 * \returns the SpectrumIdentificationResult for the spectrum.
 * creating it first if it doesn't exist
 */
SpectrumIdentificationResultPtr MzIdentMLWriter::getSpectrumIdentificationResult(
  Crux::Spectrum* spectrum ///< Crux spectrum object -in
  ) {

  string spectrum_idStr = 
    StringUtils::ToString(spectrum->getFirstScan()) + 
    "-" +
    StringUtils::ToString(spectrum->getLastScan());

  SpectrumIdentificationListPtr silp = getSpectrumIdentificationList();

  vector<SpectrumIdentificationResultPtr>::iterator sirp_iter;

  for (sirp_iter = silp->spectrumIdentificationResult.begin();
       sirp_iter != silp->spectrumIdentificationResult.end();
       ++sirp_iter) {
  
    SpectrumIdentificationResultPtr sirp = *sirp_iter;
    if (sirp->spectrumID == spectrum_idStr) {
      return sirp;
    }
  }

  SpectrumIdentificationResultPtr sirp(new SpectrumIdentificationResult());
  
  sirp->id = "SIR_"+boost::lexical_cast<string>(sir_idx_++);
  sirp->spectrumID = spectrum_idStr;
  silp->spectrumIdentificationResult.push_back(sirp);
  return sirp;
}

SpectrumIdentificationItemPtr MzIdentMLWriter::getSpectrumIdentificationItem(
  SpectrumMatch* spectrum_match
  ) {
  Crux::Spectrum* spectrum = spectrum_match->getSpectrum();
  Crux::Peptide* crux_peptide = spectrum_match->getPeptideMatch()->getPeptide();
  SpectrumZState& zstate = spectrum_match->getZState();
  int charge_state = zstate.getCharge();
  double exp_mz = zstate.getMZ();
   
  PeptidePtr peptide_ptr = getPeptide(crux_peptide);
  SpectrumIdentificationResultPtr sirp = getSpectrumIdentificationResult(spectrum_match->getSpectrum());
  vector<SpectrumIdentificationItemPtr>::iterator siip_iter;

  for (siip_iter = sirp->spectrumIdentificationItem.begin(); siip_iter != sirp->spectrumIdentificationItem.end(); ++siip_iter) {
    SpectrumIdentificationItemPtr siip = *siip_iter;
    if (siip->peptidePtr == peptide_ptr && siip->chargeState == charge_state && fabs(siip->experimentalMassToCharge - exp_mz) < 0.0001) {
      return siip;
    }
  }

  //wasn't found, create it.
  SpectrumIdentificationItemPtr siip(new SpectrumIdentificationItem(
    "SII_"+boost::lexical_cast<string>(sii_idx_++)));
  siip->chargeState = zstate.getCharge();
  siip->experimentalMassToCharge = zstate.getMZ();

  siip->calculatedMassToCharge = calculateMassToCharge(crux_peptide->calcModifiedMass(), (FLOAT_T) zstate.getCharge());

  addSpectrumScores(spectrum_match, siip);
  siip->passThreshold = true;
  siip->peptidePtr = peptide_ptr;
  addPeptideEvidences(crux_peptide, false, siip);
  sirp->spectrumIdentificationItem.push_back(siip);
  return siip;
}

void MzIdentMLWriter::addSpectrumScores(
  SpectrumMatch* spectrum_match,
  SpectrumIdentificationItemPtr siip) {

  for (ScoreMapIterator iter = spectrum_match->scoresBegin();
       iter != spectrum_match->scoresEnd();
       ++iter) {
    SCORER_TYPE_T score_type = iter->first;
    FLOAT_T score = iter->second;
    CVID cvparam_type = getScoreCVID(score_type);
    if (cvparam_type != CVID_Unknown) {
      CVParam cvparam(cvparam_type, score);
      siip->cvParams.push_back(cvparam);
    } else {
      carp_once(CARP_WARNING, "Unknown CVPARAM %d", (int)score_type);
    }

  }
}

/**
 * Adds all the peptide evidences to the SpectrumIdentificationItem
 * using the peptide's protein sources
 */
void MzIdentMLWriter::addPeptideEvidences(
  Crux::Peptide* peptide, ///< peptide to add evidence for
  bool is_decoy, ///< is peptide a decoy?
  SpectrumIdentificationItemPtr siip ///< item to add evidences to
  ) {

  for (PeptideSrcIterator src_iter = peptide->getPeptideSrcBegin();
    src_iter != peptide->getPeptideSrcEnd();
    ++src_iter) {

    PeptideEvidencePtr peptide_evidence = getPeptideEvidence(peptide, is_decoy, *src_iter);
    siip->peptideEvidencePtr.push_back(peptide_evidence);
  }
}

/**
 * \returns the mapping of SCORER_TYPE_T to the cvParam id
 */
CVID MzIdentMLWriter::getScoreCVID(
  SCORER_TYPE_T type ///< type to convert
  ) {
  switch(type) {
    case XCORR:
      return MS_SEQUEST_xcorr;
    case SP:
      return MS_SEQUEST_PeptideSp;
    case TIDE_SEARCH_EXACT_PVAL:
      return MS_p_value;
    case PERCOLATOR_SCORE:
      return MS_percolator_score;
    case PERCOLATOR_QVALUE:
      return MS_percolator_Q_value;
    case PERCOLATOR_PEP:
      return MS_percolator_PEP;
    default:
      return CVID_Unknown;  
  }
}

/**
 * \returns the mapping of SCORER_TYPE_T to cvParam id for ranks
 */
CVID MzIdentMLWriter::getRankCVID(
  SCORER_TYPE_T type ///< type to convert
  ) {
  switch(type) {
    case SP:
      return MS_SEQUEST_PeptideRankSp;
    default:
      return CVID_Unknown;
  }
}

/**
 * Adds the Match scores to the SpectrumIdentificationItem
 */
void MzIdentMLWriter::addScores(
  MatchCollection* match_collection, ///< Parent collection of match
  Match* match, ///< Match to add
  SpectrumIdentificationItemPtr item ///< item to add the scores to
  ) {

  for(int score_idx = (int)SP; 
      score_idx < (int)NUMBER_SCORER_TYPES; 
      score_idx++) {
    SCORER_TYPE_T score_type = (SCORER_TYPE_T)score_idx;
    if (match_collection->getScoredType(score_type)) {
      if (score_type == XCORR && match_collection->exact_pval_search_) {
        CVParam exactPval(MS_p_value,
                          match->getScore(TIDE_SEARCH_EXACT_PVAL));
        item->cvParams.push_back(exactPval);
        CVParam refactXCORR(MS_SEQUEST_xcorr,
                            match->getScore(TIDE_SEARCH_REFACTORED_XCORR));
        item->cvParams.push_back(refactXCORR);
      } else {
        CVID cvparam_type = getScoreCVID(score_type);
        if (cvparam_type != CVID_Unknown) {
          CVParam cvparam(cvparam_type, match->getScore(score_type));
          item->cvParams.push_back(cvparam);   
        }
      }
    }
  }

  if (match_collection->getScoredType(XCORR)) {
    CVParam delta_cn(MS_SEQUEST_deltacn, match->getScore(DELTA_CN));
    item->cvParams.push_back(delta_cn);
  }

  if (match_collection->getScoredType(BY_IONS_MATCHED)) {
    CVParam matched_ions(MS_SEQUEST_matched_ions, match->getScore(BY_IONS_MATCHED));
    item->cvParams.push_back(matched_ions);
  }
  if (match_collection->getScoredType(BY_IONS_TOTAL)) {
    CVParam total_ions(MS_SEQUEST_total_ions, match->getScore(BY_IONS_TOTAL)); 
    item->cvParams.push_back(total_ions);
  }
}

/**
 * Adds the match ranks to the SpectrumIdentificationItem
 */
void MzIdentMLWriter::addRanks(
  MatchCollection* match_collection, ///< Parent collection of the match
  Match* match, ///< Match to add
  SpectrumIdentificationItemPtr item ///< item to add the ranks to
  ) {

  if (match_collection->getScoredType(SP)) {
    CVParam sp_rank(getRankCVID(SP), match->getRank(SP));
    item->cvParams.push_back(sp_rank);
  }
}
  
void MzIdentMLWriter::write(
  MatchCollection* collection,
  string database) {
  addMatches(collection);
}

/**
 * Adds the matches in the match collection to
 * the mzid objects
 */
void MzIdentMLWriter::addMatches(
  MatchCollection* collection  ///< matches to add
  ) {
  MatchIterator match_iter(collection);
  while (match_iter.hasNext()) {
    addMatch(collection, match_iter.next());
  }
}

/**
 * Adds the match to the mzIdentML object
 */
void MzIdentMLWriter::addMatch(
  MatchCollection* collection, ///< parent collection
  Match* match ///< match to add
  ) {
  Crux::Spectrum* spectrum = match->getSpectrum();
  Crux::Peptide* peptide = match->getPeptide();
  SpectrumZState zstate = match->getZState();
   
  SpectrumIdentificationItemPtr siip(new SpectrumIdentificationItem(
    "SII_"+boost::lexical_cast<string>(sii_idx_++)));

  siip->chargeState = zstate.getCharge();
  siip->experimentalMassToCharge = zstate.getMZ();

  siip->calculatedMassToCharge = calculateMassToCharge(peptide->calcModifiedMass(), (FLOAT_T) zstate.getCharge());

  if (collection->getScoredType(PERCOLATOR_SCORE)) {
    siip->rank = match->getRank(PERCOLATOR_SCORE);
  } else if (collection->getScoredType(XCORR)) { 
    siip->rank = match->getRank(XCORR);
  }
  addScores(collection, match, siip);
  addRanks(collection, match, siip);
  siip->passThreshold = true;
  siip->peptidePtr = getPeptide(peptide);
  addPeptideEvidences(peptide, match->isDecoy(), siip);

  SpectrumIdentificationResultPtr sirp =
    getSpectrumIdentificationResult(spectrum);
  sirp->spectrumIdentificationItem.push_back(siip);
}

void MzIdentMLWriter::addProteinScores(
  ProteinDetectionHypothesisPtr pdhp,
  ProteinMatch* protein_match) {

  //TODO should we use these instead?
    /// q-value for peptides: Peptide identification confidence metric q-value.
  //    MS_q_value_for_peptides = 1001868,

    /// q-value for proteins: Protein identification confidence metric q-value.
  //   MS_q_value_for_proteins = 1001869,

  for (ScoreMapIterator iter = protein_match->scoresBegin();
    iter != protein_match->scoresEnd();
    ++iter) {

    SCORER_TYPE_T score_type = iter->first;
    if (score_type == DELTA_CN || score_type == DELTA_LCN ||
        score_type == BY_IONS_MATCHED || score_type == BY_IONS_TOTAL) {
      continue;
    }
    FLOAT_T score = iter->second;
    CVID cvparam_type = getScoreCVID(score_type);
    if (cvparam_type != CVID_Unknown) {
      CVParam cvparam(cvparam_type, score);
      pdhp->cvParams.push_back(cvparam);
    } else {
      carp(CARP_WARNING, "Unknown parameter type for score type:%d", score_type);
      //TODO create a user param.
    }
  }
}

void MzIdentMLWriter::addPeptideScores(
  PeptideMatch* peptide_match
  ) {
  //So PeptideHypothesis doesn't let you add cvparams and they
  //are not unique,  the Peptide tag requires uniqueness, so
  //I'll add the q-values and scores to those elements (SJM).
  PeptidePtr peptidep = getPeptide(peptide_match->getPeptide());

  //we already added the scores for this peptide.
  if (peptidep->cvParams.size() != 0) {
    return;
  }

  for (ScoreMapIterator iter = peptide_match->scoresBegin();
    iter != peptide_match->scoresEnd();
    ++iter) {

    SCORER_TYPE_T score_type = iter->first;
    if (score_type == DELTA_CN || score_type == DELTA_LCN ||
        score_type == BY_IONS_MATCHED || score_type == BY_IONS_TOTAL) {
      continue;
    }
    FLOAT_T score = iter->second;
    CVID cvparam_type = getScoreCVID(score_type);
    if (cvparam_type != CVID_Unknown) {
      CVParam cvparam(cvparam_type, score);
      peptidep->cvParams.push_back(cvparam);
    } else {
      carp(CARP_WARNING, "Unknown parameter type for score type:%d", score_type);
      //TODO create a user param.
    }
  }
}

void MzIdentMLWriter::addSpectrumMatches(
  ProteinMatch* protein_match,
  PeptideMatch* peptide_match
  ) {
  PeptideHypothesis& peptide_hypothesis = getPeptideHypothesis(protein_match, peptide_match);
  for (SpectrumMatchIterator spectrum_iter = peptide_match->spectrumMatchBegin();
       spectrum_iter != peptide_match->spectrumMatchEnd();
       ++spectrum_iter) {
    SpectrumMatch* spectrum_match = *spectrum_iter;
    SpectrumIdentificationItemPtr sip = getSpectrumIdentificationItem(spectrum_match);
    peptide_hypothesis.spectrumIdentificationItemPtr.push_back(sip);
  }
}

void MzIdentMLWriter::addPeptideMatches(
  ProteinMatch* protein_match
  ) {
  for (PeptideMatchIterator peptide_iter = protein_match->peptideMatchBegin();
       peptide_iter != protein_match->peptideMatchEnd();
       ++peptide_iter) {
    PeptideMatch* peptide_match = *peptide_iter;
    addPeptideScores(peptide_match);
    addSpectrumMatches(protein_match, peptide_match);
  }
}

/**
 * Adds the protein matches to the mzid object
 */
void MzIdentMLWriter::addProteinMatches(
  ProteinMatchCollection* protein_match_collection
) {
  if (protein_match_collection == NULL) {
    carp(CARP_FATAL, "ProteinMatchCollection was null");
  }
  //now add the protein matches.
  for (ProteinMatchIterator match_iter = protein_match_collection->proteinMatchBegin();
       match_iter != protein_match_collection->proteinMatchEnd();
       match_iter++) {
    addProteinMatch(*match_iter);
  }
}

/**
 * Adds a protein match to the mzIdentML object
 */
void MzIdentMLWriter::addProteinMatch(
  ProteinMatch* protein_match
  ) {
  string protein_id = protein_match->getId();

  //For now, each protein should get its own protein ambiguity group
  ProteinDetectionHypothesisPtr pdhp = getProteinDetectionHypothesis(protein_id);
  addProteinScores(pdhp, protein_match);
  addPeptideMatches(protein_match);
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */


