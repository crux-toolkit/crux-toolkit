/*************************************************************************//**
 * \file MzIdentMLReader.cpp
 * \brief Object for parsing pepxml files
 ****************************************************************************/

#include "MzIdentMLReader.h"
#include "mass.h"
#include "expat.h"


#include "Protein.h"
#include "Peptide.h"


#include <cstdio>
#include <cstring>

#include <iostream>

#include "DelimitedFile.h"
#include "parameter.h"
#include "MatchCollectionParser.h"
#include "Peptide.h"


using namespace std;
using namespace Crux;
using namespace pwiz;
using namespace identdata;


/**
 * Initializes the object
 */
void MzIdentMLReader::init() {
  match_collection_ = NULL;
  use_pass_threshold_ = get_boolean_parameter("mzid-use-pass-threshold");
}

/**
 * \returns an initialized object
 */
MzIdentMLReader::MzIdentMLReader() {
  init();
}

/**
 * \returns an object initialized with the file_path
 */
MzIdentMLReader::MzIdentMLReader(
  const string& file_path ///< the path of the pep.xml file
  ) {
  
  init();
  file_path_ = file_path;
}

/**
 * \returns an object initialized with the xml path, and the target,decoy databases
 */
MzIdentMLReader::MzIdentMLReader(
  const string& file_path, ///< the path of the pep.xml
  Database* database, ///< the protein database
  Database* decoy_database ///< the decoy protein database (can be null)
  ) {

  file_path_ = file_path;
  database_ = database;
  decoy_database_ = decoy_database;

}

/**
 * sets the target database for the parser
 */
void MzIdentMLReader::setDatabase(
  Database* database ///< the target protein database
  ) {

  database_ = database;

}

/**
 * sets the decoy protein database for the parser
 */
void MzIdentMLReader::setDecoyDatabase(
  Database* decoy_database ///<  the decoy protein database
  ) {

  decoy_database_ = decoy_database;

}

/**
 * default destructor
 */
MzIdentMLReader::~MzIdentMLReader() {

}

/**
 * \returns the MatchCollection resulting from the parsed xml file
 */

void MzIdentMLReader::addScores(
  const SpectrumIdentificationItem& item, ///<proteowizard psm
  Match* match ///<our psm
) {
  vector<CVParam>::const_iterator iter = item.cvParams.begin();
  
  FLOAT_T fvalue;
  int ivalue;

  for (; iter != item.cvParams.end(); ++iter) {
    switch (iter->cvid) {
      case MS_SEQUEST_xcorr:
        from_string(fvalue, iter->value);
        match->setScore(XCORR, fvalue);
        match_collection_->setScoredType(XCORR, true);
        break;
      case MS_SEQUEST_PeptideSp:
        from_string(fvalue, iter->value);
        match->setScore(SP, fvalue);
        match_collection_->setScoredType(SP, true);
        break;
      case MS_SEQUEST_PeptideRankSp:
        from_string(ivalue, iter->value);
        match->setRank(SP, ivalue);
        break;
      case MS_SEQUEST_deltacn:
        from_string(fvalue, iter->value);
        match->setDeltaCn(fvalue);
        break;
      case MS_SEQUEST_matched_ions:
        from_string(ivalue, iter->value);
        match->setBYIonMatched(ivalue);
        break;
      case MS_SEQUEST_total_ions:
        from_string(ivalue, iter->value);
        match->setBYIonPossible(ivalue);
        break;
      default:
        carp(CARP_DEBUG, "Unknown score type, will be set in custom scores");
    }
    //go ahead and set all custom scores to the cvParam names.
    string name = cvTermInfo((*iter).cvid).name;
    from_string(fvalue, iter->value);
    match->setCustomScore(name, fvalue);
  }


  vector<UserParam>::const_iterator iter2 = item.userParams.begin();

  for (; iter2 != item.userParams.end(); ++iter2) {
    string name = iter2->name;

    bool success = from_string(fvalue, iter2->value);
    if (success) {
      match->setCustomScore(name, fvalue);
    }
  }
}

/**
 * \returns the MatchCollection resulting from the parsed xml file
 */
MatchCollection* MzIdentMLReader::parse(
    const char* path, ///< path of the xml file
    Database* database, ///< target protein database
    Database* decoy_database ///< decoy protein database (can be null)
  ) {

  MzIdentMLReader* reader = new MzIdentMLReader(path);
  reader->setDatabase(database);
  reader->setDecoyDatabase(decoy_database);

  MatchCollection* collection = reader->parse();

  delete reader;

  return collection;

}
 
/**
 * \returns the MatchCollection resulting from the parsed xml file
 */
MatchCollection* MzIdentMLReader::parse() {

  match_collection_ = new MatchCollection();
  match_collection_ -> preparePostProcess();
  carp(CARP_DEBUG, "MzIdentMLReader::opening file:%s", file_path_.c_str());
  pwiz_reader_ = new IdentDataFile(file_path_);

  parseDatabaseSequences();

  parsePSMs();

  return match_collection_;
}

/**
 * parses the database sequences from the mzid file
 */
void MzIdentMLReader::parseDatabaseSequences() {
  vector<DBSequencePtr>::const_iterator db_iter;
  vector<DBSequencePtr>::const_iterator db_end;

  db_iter = pwiz_reader_->sequenceCollection.dbSequences.begin();
  db_end = pwiz_reader_->sequenceCollection.dbSequences.end();

  for (; db_iter != db_end; ++db_iter) {

    DBSequence& db_sequence = (**db_iter);
    int length = db_sequence.length;
    string accession = db_sequence.accession;
    string seq = db_sequence.seq;

    //cerr <<"protein accession:"<<accession<<endl;
    //cerr <<"protein length:"<<length<<endl;
    //cerr <<"protein sequence:"<<seq<<endl;
    //cerr <<"sequence length:"<<seq.length()<<endl;
    Protein* protein = NULL;
    bool is_decoy;

    if (seq.length() == 0) {
      protein = MatchCollectionParser::getProtein(
        database_, decoy_database_, accession, is_decoy);
    } else {
      protein = MatchCollectionParser::getProtein(
        database_, decoy_database_, accession, seq, is_decoy);
    }
  }
  


}

/**
 * parses the psms from the mzid file
 */
void MzIdentMLReader::parsePSMs() {

  vector<SpectrumIdentificationListPtr>::const_iterator sil_iter;
  vector<SpectrumIdentificationListPtr>::const_iterator sil_end;
  vector<SpectrumIdentificationResultPtr>::const_iterator sir_iter;
  vector<SpectrumIdentificationItemPtr>::const_iterator sii_iter;

  sil_iter = pwiz_reader_->dataCollection.analysisData.spectrumIdentificationList.begin();
  sil_end = pwiz_reader_->dataCollection.analysisData.spectrumIdentificationList.end();

  int count = 0;

  for (; sil_iter != sil_end; ++sil_iter) {
    for (sir_iter = (**sil_iter).spectrumIdentificationResult.begin();
      sir_iter != (**sil_iter).spectrumIdentificationResult.end();
      ++sir_iter) {

      SpectrumIdentificationResult& result = **sir_iter;
      string idStr = result.spectrumID;
      carp(CARP_DEBUG, "idStr:%s", idStr.c_str());
      //TODO crux requires scan numbers to be integer, where mzid can have
      //them be strings. I don't know what the appropiate thing to do here.
      //For now, lets parse as a first_scan-last_scan string, since our
      //writer outputs it in that format.  SJM. 
      int first_scan;
      int last_scan;
      if (!get_first_last_scan_from_string(idStr, first_scan, last_scan)) {
        carp(CARP_ERROR, "Cannot find first,last scan from spectrumID:%s",
             idStr.c_str());    
      }

      //We are not using this yet, but we might.
      string filename = "";
      if (result.spectraDataPtr != NULL) {
        string filename = result.spectraDataPtr->location;
        carp(CARP_DEBUG, "filename:%s", filename.c_str());
      }
      for (sii_iter = result.spectrumIdentificationItem.begin();
        sii_iter != result.spectrumIdentificationItem.end();
        ++sii_iter) {
        SpectrumIdentificationItem& item = **sii_iter;
        if (!use_pass_threshold_ || item.passThreshold) {
          int charge = item.chargeState;
          FLOAT_T obs_mz = item.experimentalMassToCharge;
          
          SpectrumZState zstate;
          zstate.setMZ(obs_mz, charge);
          vector<int> charge_vec;
          charge_vec.push_back(charge);


          Spectrum* spectrum = 
            new Spectrum(first_scan,last_scan,obs_mz, charge_vec, "");

          FLOAT_T calc_mz = item.calculatedMassToCharge;
          FLOAT_T calc_mass = (calc_mz - MASS_PROTON ) * (FLOAT_T)charge;
          int rank = item.rank;


          PeptidePtr peptide_ptr = item.peptidePtr;
          string sequence = peptide_ptr->peptideSequence;
        
          vector<PeptideEvidencePtr>& peptide_evidences = item.peptideEvidencePtr;
          PeptideEvidencePtr peptide_evidence_ptr = peptide_evidences.front();
          string protein_id = peptide_evidence_ptr->dbSequencePtr->accession;
          int start_idx = peptide_evidence_ptr->start;

          bool is_decoy = false;
          bool is_decoy_test;

          carp(CARP_DEBUG,"getting protein %s",protein_id.c_str());

          Protein* protein = MatchCollectionParser::getProtein(
            database_, decoy_database_, protein_id, is_decoy_test);
       
          if (peptide_evidence_ptr->isDecoy != is_decoy_test) {
            carp_once(CARP_WARNING, "Protein 0) %s : mzid says isdecoy: %d, "
              "but database says isdecoy: %d, did you set decoy-prefix?", 
              protein_id.c_str(), 
              peptide_evidence_ptr->isDecoy, is_decoy_test);
          }
          //If there is one PeptideEvidence object that is a decoy
          //then mark the match as a decoy.
          is_decoy = is_decoy || peptide_evidence_ptr->isDecoy || is_decoy_test;

          start_idx = protein->findStart(sequence, "", "");
          if (start_idx == -1) {
            carp(CARP_FATAL, "can't find sequence %s in first protein %s",sequence.c_str(), protein->getIdPointer());
          }
          int length = sequence.length();
        
          carp(CARP_DEBUG, "creating peptide %s %f %i",sequence.c_str(), calc_mass, start_idx);

          Crux::Peptide* peptide = 
            new Crux::Peptide(length, calc_mass, protein, start_idx);
        
          for (int pe_idx = 1; pe_idx < peptide_evidences.size();pe_idx++) {
            PeptideEvidencePtr peptide_evidence_ptr = peptide_evidences[pe_idx];
            int start = peptide_evidence_ptr->start;
            int end = peptide_evidence_ptr->end;
            protein_id = peptide_evidence_ptr->dbSequencePtr->accession; 
            carp(CARP_DEBUG, "id: %s start:%i end: %i decoy: %i", protein_id.c_str(),
             start, end, is_decoy);

            protein = MatchCollectionParser::getProtein(
              database_, decoy_database_, protein_id, is_decoy_test);
            
            if (peptide_evidence_ptr->isDecoy != is_decoy_test) {
              carp_once(CARP_WARNING, "Protein %i) %s: mzid says isdecoy: %i, "
                "but database says isdecoy: %i, did you set \"decoy-prefix\"?", 
                pe_idx, protein_id.c_str(), peptide_evidence_ptr->isDecoy, is_decoy_test);
            }
            is_decoy = is_decoy || peptide_evidence_ptr->isDecoy || is_decoy_test;
            start_idx = protein->findStart(sequence, "", "");
            if (start_idx != -1) {
              PeptideSrc* src = new PeptideSrc((DIGEST_T)0, protein, start_idx);
              peptide->addPeptideSrc(src);
            }
          }
          Match* match = new Match(peptide, spectrum, zstate, is_decoy);  
          match_collection_->addMatchToPostMatchCollection(match);

          match->setRank(XCORR, rank); // Is it safe to assume this?
          carp(CARP_DEBUG, "charge: %i obs mz: %f calc mass: %f sequence: %s", charge, obs_mz, calc_mass, sequence.c_str());
          addScores(item, match);
        }
      }
    }


    count++;
  }
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
