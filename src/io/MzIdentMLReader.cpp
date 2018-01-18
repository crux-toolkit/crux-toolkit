/*************************************************************************
 * \file MzIdentMLReader.cpp
 * \brief Object for parsing pepxml files
 *************************************************************************/

#include "MzIdentMLReader.h"
#include "util/mass.h"
#include "util/Params.h"
#include "util/StringUtils.h"
#include "expat.h"
#include "model/Protein.h"
#include "model/Peptide.h"

#include <cstdio>
#include <cstring>
#include <iostream>

#include "DelimitedFile.h"
#include "parameter.h"
#include "MatchCollectionParser.h"

using namespace std;
using namespace Crux;
using namespace pwiz;
using namespace identdata;

/**
 * Initializes the object
 */
void MzIdentMLReader::init() {
  match_collection_ = NULL;
  use_pass_threshold_ = Params::GetBool("mzid-use-pass-threshold");
}

/**
 * \returns an initialized object
 */
MzIdentMLReader::MzIdentMLReader() : PSMReader() {
  init();
}

/**
 * \returns an object initialized with the file_path
 */
MzIdentMLReader::MzIdentMLReader(
  const string& file_path ///< the path of the pep.xml file
  ) : PSMReader(file_path) {
  init();
}

/**
 * \returns an object initialized with the xml path, and the target,decoy databases
 */
MzIdentMLReader::MzIdentMLReader(
  const string& file_path, ///< the path of the pep.xml
  Database* database, ///< the protein database
  Database* decoy_database ///< the decoy protein database (can be null)
  ) : PSMReader(file_path, database, decoy_database) {
  init();
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
  FLOAT_T fvalue;
  int ivalue;

  for (vector<CVParam>::const_iterator i = item.cvParams.begin();
       i != item.cvParams.end();
       ++i) {
    switch (i->cvid) {
      case MS_SEQUEST_xcorr:
        fvalue = StringUtils::FromString<FLOAT_T>(i->value);
        match->setScore(XCORR, fvalue);
        match_collection_->setScoredType(XCORR, true);
        break;
      case MS_SEQUEST_PeptideSp:
        fvalue = StringUtils::FromString<FLOAT_T>(i->value);
        match->setScore(SP, fvalue);
        match_collection_->setScoredType(SP, true);
        break;
      case MS_SEQUEST_PeptideRankSp:
        ivalue = StringUtils::FromString<int>(i->value);
        match->setRank(SP, ivalue);
        break;
      case MS_SEQUEST_deltacn:
        fvalue = StringUtils::FromString<FLOAT_T>(i->value);
        match->setScore(DELTA_CN, fvalue);
        break;
      case MS_SEQUEST_matched_ions:
        ivalue = StringUtils::FromString<int>(i->value);
        match->setScore(BY_IONS_MATCHED, ivalue);
        break;
      case MS_SEQUEST_total_ions:
        ivalue = StringUtils::FromString<int>(i->value);
        match->setScore(BY_IONS_TOTAL, ivalue);
        break;
      case MS_p_value: // exact p-value
        fvalue = StringUtils::FromString<FLOAT_T>(i->value);
        match->setScore(TIDE_SEARCH_EXACT_PVAL, fvalue);
        match_collection_->setScoredType(TIDE_SEARCH_EXACT_PVAL, true);
        break;
      default:
        carp(CARP_DEBUG, "Unknown score type, will be set in custom scores");
    }
    //go ahead and set all custom scores to the cvParam names.
    string name = cvTermInfo(i->cvid).name;
    fvalue = StringUtils::FromString<FLOAT_T>(i->value);
    match->setCustomScore(name, fvalue);
  }

  for (vector<UserParam>::const_iterator i = item.userParams.begin();
       i != item.userParams.end();
       ++i) {
    if (StringUtils::TryFromString(i->value, &fvalue)) {
      match->setCustomScore(i->name, fvalue);
    }
  }
}

/**
 * \returns the MatchCollection resulting from the parsed xml file
 */
MatchCollection* MzIdentMLReader::parse(
  const string& path, ///< path of the xml file
  Database* database, ///< target protein database
  Database* decoy_database ///< decoy protein database (can be null)
  ) {
  MzIdentMLReader reader(path);
  reader.setDatabase(database);
  reader.setDecoyDatabase(decoy_database);
  MatchCollection* collection = reader.parse();
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
  parseMods();
  parseDatabaseSequences();
  parsePSMs();
  return match_collection_;
}

void MzIdentMLReader::parseMods() {
  AnalysisProtocolCollection& apc = pwiz_reader_->analysisProtocolCollection;
  for (vector<SpectrumIdentificationProtocolPtr>::const_iterator i = apc.spectrumIdentificationProtocol.begin();
       i != apc.spectrumIdentificationProtocol.end();
       i++) {
    for (vector<SearchModificationPtr>::const_iterator j = (*i)->modificationParams.begin();
         j != (*i)->modificationParams.end();
         j++) {
      const SearchModification& sm = **j;
      ModPosition position;
      switch (sm.specificityRules.cvid) {
        case MS_modification_specificity_peptide_N_term: position = PEPTIDE_N; break;
        case MS_modification_specificity_peptide_C_term: position = PEPTIDE_C; break;
        case MS_modification_specificity_protein_N_term: position = PROTEIN_N; break;
        case MS_modification_specificity_protein_C_term: position = PROTEIN_C; break;
        default:                                         position = ANY; break;
      }
      ModificationDefinition::New(
        StringUtils::Join(sm.residues), sm.massDelta, position, sm.fixedMod);
    }
  }
}

/**
 * parses the database sequences from the mzid file
 */
void MzIdentMLReader::parseDatabaseSequences() {
  for (vector<DBSequencePtr>::const_iterator i = pwiz_reader_->sequenceCollection.dbSequences.begin();
       i != pwiz_reader_->sequenceCollection.dbSequences.end();
       i++) {
    DBSequence& db_sequence = (**i);
    int length = db_sequence.length;
    string accession = db_sequence.accession;
    string seq = db_sequence.seq;

    Protein* protein = NULL;
    bool is_decoy;

    if (seq.empty()) {
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
  int maxRank = Params::GetInt("top-match-in");
  for (vector<SpectrumIdentificationListPtr>::const_iterator i =
         pwiz_reader_->dataCollection.analysisData.spectrumIdentificationList.begin();
       i != pwiz_reader_->dataCollection.analysisData.spectrumIdentificationList.end();
       i++) {
    for (vector<SpectrumIdentificationResultPtr>::const_iterator j =
           (*i)->spectrumIdentificationResult.begin();
         j != (*i)->spectrumIdentificationResult.end();
         ++j) {
      SpectrumIdentificationResult& result = **j;
      string idStr = result.spectrumID;
      carp(CARP_DEBUG, "idStr:%s", idStr.c_str());
      //TODO crux requires scan numbers to be integer, where mzid can have
      //them be strings. I don't know what the appropiate thing to do here.
      //For now, lets parse as a first_scan-last_scan string, since our
      //writer outputs it in that format.  SJM.
      int first_scan, last_scan;
      if (!get_first_last_scan_from_string(idStr, first_scan, last_scan)) {
        carp(CARP_ERROR, "Cannot find first,last scan from spectrumID:%s",
             idStr.c_str());
      }

      //We are not using this yet, but we might.
      string filename;
      if (result.spectraDataPtr != NULL) {
        string filename = result.spectraDataPtr->location;
        carp(CARP_DEBUG, "filename:%s", filename.c_str());
      }
      for (vector<SpectrumIdentificationItemPtr>::const_iterator k =
             result.spectrumIdentificationItem.begin();
           k != result.spectrumIdentificationItem.end();
           k++) {
        SpectrumIdentificationItem& item = **k;
        if (use_pass_threshold_ && !item.passThreshold) {
          continue;
        } else if (maxRank != 0 && item.rank > maxRank) {
          continue;
        }
        int charge = item.chargeState;
        FLOAT_T obs_mz = item.experimentalMassToCharge;

        SpectrumZState zstate;
        zstate.setMZ(obs_mz, charge);
        vector<int> charge_vec;
        charge_vec.push_back(charge);

        Spectrum* spectrum = new Spectrum(first_scan, last_scan, obs_mz, charge_vec, "");
        FLOAT_T calc_mass = (item.calculatedMassToCharge - MASS_PROTON) * (FLOAT_T)charge;
        string sequence = item.peptidePtr->peptideSequence;

        vector<PeptideEvidencePtr>& peptide_evidences = item.peptideEvidencePtr;
        PeptideEvidencePtr peptide_evidence_ptr = peptide_evidences.front();

        string protein_id = peptide_evidence_ptr->dbSequencePtr->accession;
        int start_idx = peptide_evidence_ptr->start;

        carp(CARP_DEBUG, "getting protein %s", protein_id.c_str());
        bool is_decoy_test;
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
        bool is_decoy = is_decoy || peptide_evidence_ptr->isDecoy || is_decoy_test;

        start_idx = protein->findStart(sequence, "", "");

        if (start_idx == -1) {
          carp(CARP_FATAL, "can't find sequence %s in first protein %s", sequence.c_str(), protein->getIdPointer().c_str());
        }

        int length = sequence.length();
        carp(CARP_DEBUG, "creating peptide %s %f %i", sequence.c_str(), calc_mass, start_idx);

        Crux::Peptide* peptide = new Crux::Peptide(length, protein, start_idx);

        // Add mods to the peptide
        PeptidePtr pep_ptr = item.peptidePtr;
        for (vector<ModificationPtr>::const_iterator l = pep_ptr->modification.begin();
             l != pep_ptr->modification.end();
             l++) {
          const ModificationDefinition* mod = ModificationDefinition::Find(
            (*l)->monoisotopicMassDelta, false);
          if (mod != NULL) {
            unsigned char index = (*l)->location;
            switch (mod->Position()) {
              case PEPTIDE_N:
              case PROTEIN_N:
                index = 0;
                break;
              case PEPTIDE_C:
              case PROTEIN_C:
                index = length - 1;
                break;
              default:
                index = (*l)->location - 1;
            }
            peptide->addMod(mod, index);
          } else {
            carp(CARP_ERROR, "Unknown modification: %f", (*l)->monoisotopicMassDelta);
          }
        }

        for (int pe_idx = 1; pe_idx < peptide_evidences.size(); pe_idx++) {
          PeptideEvidencePtr peptide_evidence_ptr = peptide_evidences[pe_idx];
          int start = peptide_evidence_ptr->start;
          int end = peptide_evidence_ptr->end;
          protein_id = peptide_evidence_ptr->dbSequencePtr->accession;
          carp(CARP_DEBUG, "id: %s start:%i end: %i decoy: %i",
               protein_id.c_str(), start, end, is_decoy);

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
        // trying to set lnExperimentSize to -1 -> representing unavailable information on distinct matches/spectrum
        match->setLnExperimentSize(-1);
        match_collection_->addMatchToPostMatchCollection(match);
        match->setRank(XCORR, item.rank); // Is it safe to assume this?
        carp(CARP_DEBUG, "charge: %i obs mz: %f calc mass: %f sequence: %s",
             charge, obs_mz, calc_mass, sequence.c_str());
        addScores(item, match);
      }
    }
  }
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
