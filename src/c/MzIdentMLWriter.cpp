/**
 * \file MzIdentMLWriter.cpp
 * \brief Writes search results in the MzIdentML (mzid) format.
 */
#include <iostream>
#include <fstream>
#include <string>
#include "MzIdentMLWriter.h"
#include "crux-utils.h"
#include "MatchCollection.h"
#include "pwiz/data/identdata/Serializer_mzid.hpp"


using namespace std;
using namespace pwiz;
using namespace identdata;

MzIdentMLWriter::MzIdentMLWriter() {
  //data_ = NULL;
  fout_ = NULL;
  sir_idx_ = 0;
  sil_idx_ = 0;
  sii_idx_ = 0;
  peptide_idx_ = 0;
  peptide_evidence_idx_ = 0;
  dbs_idx_ = 0;
}

MzIdentMLWriter::~MzIdentMLWriter()
{
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

/**
 * Close the file, if open.
 */
void MzIdentMLWriter::closeFile(){
  if( mzid_ != NULL && fout_ != NULL){
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
 
  char* sequence = peptide->getSequence();
  string sequence_str = sequence;
  free(sequence);

  vector<PeptidePtr>::iterator peptide_iter;

  for (peptide_iter = mzid_->sequenceCollection.peptides.begin();
       peptide_iter != mzid_->sequenceCollection.peptides.end();
       ++peptide_iter) {
    if ((*peptide_iter)->peptideSequence == sequence_str) {
      return (*peptide_iter);
    }
  }

  PeptidePtr peptide_p(new pwiz::identdata::Peptide("PEP_"+boost::lexical_cast<string>(peptide_idx_++)));
  //peptide_p->id = sequence_str;
  peptide_p->peptideSequence = sequence_str;
  mzid_->sequenceCollection.peptides.push_back(peptide_p);
  //TODO handle peptide modifications.

  return peptide_p;

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
    dbs_ptr->length=sequence_str.length();
    dbs_ptr->seq=sequence_str;
  } else {
    dbs_ptr->length=src->getParentProtein()->getLength();
    dbs_ptr->seq=src->getParentProtein()->getSequencePointer();
    //TODO add description
  }
  
  mzid_->sequenceCollection.dbSequences.push_back(dbs_ptr);

  return dbs_ptr;
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

  for (pe_iter = mzid_->sequenceCollection.peptideEvidence.begin();
       pe_iter != mzid_->sequenceCollection.peptideEvidence.end();
	 ++pe_iter) {
    PeptideEvidencePtr pe_ptr = *pe_iter;
    if ((pe_ptr->peptidePtr->peptideSequence == sequence_str) && (pe_ptr->dbSequencePtr->accession==protein_id)) {
      return pe_ptr;
    }

  }

  DBSequencePtr dbs_ptr = getDBSequence(peptide, src);
  PeptideEvidencePtr pe_ptr(new PeptideEvidence(
    "PE_"+boost::lexical_cast<string>(peptide_evidence_idx_++)));

  if (src->getParentProtein()->isPostProcess()) {
    pe_ptr->start=0;
    pe_ptr->end=peptide->getLength();
  } else {
    pe_ptr->start=src->getStartIdx();
    pe_ptr->end=src->getStartIdx()+peptide->getLength();
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
    silp = SpectrumIdentificationListPtr(
        new SpectrumIdentificationList("SIL_"+boost::lexical_cast<string>(sil_idx_++)));
    mzid_->dataCollection.analysisData.spectrumIdentificationList.push_back(silp);   
  }
  return silp;

}

/*
ProteinDetectionListPtr MzIdentMLWriter::getProteinIdentificationList() {

  ProteinDetectionListPtr pdlp;
  if (mzid_ -> dataCollection.analysisData.proteinDetectionListPtr != NULL) {
    pdlp = mzid_->dataCollection.analysisData.proteinDetectionListPtr;
  } else {
    pdlp = ProteinDetectionListPtr(
      new ProteinDetectionListPtr("PDL_1"));
    mzid_->dataCollection.analysisData.proteinDetectionListPtr = pdlp;
  }
  return pdlp;
}

ProteinAmbiguityGroupPtr MzIdentMLWriter::getProteinAmbiguityGroup() {

  ProteinDetectionListPtr pdlp = getProteinIdentificationList();

  ProteinAmbiguityGroupPtr pagp;  

  if (pdlp->proteinAmbiguityGroup.size() > 0) {
    pagp = pdlp->proteinAmbiguityGroup[0];
  } else {
    pagp = ProteinAmbiguityGroupPtr(
      new ProteinAmbiguityGroup("PAG_1"));
    pdlp->proteinAmbiguityGroup.push_back(pagp);
  }
  return pagp;
}

ProteinDetectionHypothesisPtr MzIdentMLWriter::getProteinDetectionHypothesis() {

  ProteinAmbiguityGroupPtr pagp = getProteinAmbiguityGroup();
  ProteinDetectionHypothesis pdhp;

  return pdhp;
}


PeptideHypothesis& MzIdentMLWriter::getPeptideHypothesis() {

  ProteinDetectionHypothesisPtr pdhp = getProteinDetectionHypothesis();
  return pdhp->peptideHypothesis[0];


}
*/

/**
 * \returns the SpectrumIdentificationResult for the spectrum.
 * creating it first if it doesn't exist
 */
SpectrumIdentificationResultPtr MzIdentMLWriter::getSpectrumIdentificationResult(
  Crux::Spectrum* spectrum ///< Crux spectrum object -in
  ) {

  string spectrum_idStr = DelimitedFileWriter::to_string(spectrum->getFirstScan());

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

  SpectrumIdentificationResultPtr sirp(
   new SpectrumIdentificationResult());
  
  sirp->id = "SIR_"+boost::lexical_cast<string>(sir_idx_++);
  sirp->spectrumID = spectrum_idStr;
  silp->spectrumIdentificationResult.push_back(sirp);
  return sirp;
  

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
      return MS_Sequest_xcorr;
    case SP:
      return MS_Sequest_PeptideSp;
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
      return MS_Sequest_PeptideRankSp;
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
      CVID cvparam_type = getScoreCVID(score_type);
      if (cvparam_type != CVID_Unknown) {
        CVParam cvparam(cvparam_type, match->getScore(score_type));
        item->cvParams.push_back(cvparam);   
      }
    }
  }

  CVParam delta_cn(MS_Sequest_deltacn, match->getDeltaCn());
  item->cvParams.push_back(delta_cn);

  if (match_collection->getScoredType(SP)) {
    CVParam matched_ions(MS_Sequest_matched_ions, match->getBYIonMatched());
    item->cvParams.push_back(matched_ions);
    CVParam total_ions(MS_Sequest_total_ions, match->getBYIonPossible()); 
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
  siip->calculatedMassToCharge = FLOAT_T((peptide->getPeptideMass()+MASS_PROTON)/(double)zstate.getCharge());
  siip->rank = match->getRank(XCORR);
  addScores(collection, match, siip);
  addRanks(collection, match, siip);
  siip->passThreshold = true;
  siip->peptidePtr = getPeptide(peptide);
  addPeptideEvidences(peptide, match->isDecoy(), siip);

  SpectrumIdentificationResultPtr sirp =
    getSpectrumIdentificationResult(spectrum);
  sirp->spectrumIdentificationItem.push_back(siip);
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */


