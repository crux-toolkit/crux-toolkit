/**
 * \file MzIdentMLWriter.h
 * \brief Writes search results in the MzIdentML (mzid) format.
 */
#ifndef MZIDENTMLWRITER_H
#define MZIDENTMLWRITER_H

#include <string>
#include <vector>
#include "objects.h"
#include "pwiz/data/identdata/IdentData.hpp"

class MzIdentMLWriter{

 protected:
  std::ofstream* fout_; ///< stream to write the mzid to
  pwiz::identdata::IdentDataPtr mzid_; ///< structure of mzidentml 
  size_t sir_idx_; ///< counter for SpectrumIdentification results
  size_t sil_idx_; ///<counter for SpectrumIdentificationList
  size_t sii_idx_; ///<counter for SpectrumIdentificationItem
  size_t peptide_idx_; ///<counter for pwiz Peptide
  size_t peptide_evidence_idx_; ///<counter for PeptideEvidence
  size_t dbs_idx_; ///< counter for DBSequence

  /**
   * \returns DBSequence for the protein source.  If it doesn't exist, 
   * then first create the object in the mzid object
   */
  pwiz::identdata::DBSequencePtr getDBSequence(
    Crux::Peptide* peptide, ///< peptide -in 
    PeptideSrc* src ///< Source of the peptide -in
  );

  /**
   * \returns PeptideEvidence for the peptide and src.
   * creates it if it doesn't exist
   */
  pwiz::identdata::PeptideEvidencePtr getPeptideEvidence(
    Crux::Peptide* peptide, ///< peptide -in
    bool is_decoy, ///< is this peptide a decoy? -in
    PeptideSrc* src ///< where the peptide is coming from -in
  );

  /**
   * \returns the MzIdentML Peptide object, creating it if
   * it doesn't exist already
   */
  pwiz::identdata::PeptidePtr getPeptide(
    Crux::Peptide* peptide ///< peptide -in
  );

  /**
   * \returns the SpectrumIdentificationList, creating it if 
   * it doesn't exist yet.
   */
  pwiz::identdata::SpectrumIdentificationListPtr getSpectrumIdentificationList();

  /**
   * \returns the SpectrumIdentificationResult for the spectrum.
   * creating it first if it doesn't exist
   */
  pwiz::identdata::SpectrumIdentificationResultPtr getSpectrumIdentificationResult(
    Crux::Spectrum* spectrum ///< Crux spectrum object -in
  );

  /**
   * Adds all the peptide evidences to the SpectrumIdentificationItem
   * using the peptide's protein sources
   */
  void addPeptideEvidences(
    Crux::Peptide* peptide, ///< peptide to add evidence for
    bool is_decoy, ///< is peptide a decoy?
    pwiz::identdata::SpectrumIdentificationItemPtr siip ///<item to add evidences to.
  );

  /**
   * \returns the mapping of SCORER_TYPE_T to the cvParam id
   */
  pwiz::cv::CVID getScoreCVID(
    SCORER_TYPE_T type ///< type to convert
  );

  /**
   * \returns the mapping of SCORER_TYPE_T to the cvPARAM id for ranks
   */
  pwiz::cv::CVID getRankCVID(
    SCORER_TYPE_T type ///< type to convert
  );

  /**
   * Adds the Match scores to the SpectrumIdentificationItem
   */
  void addScores(
    MatchCollection* match_collection, ///< Parent collection of match
    Crux::Match* match, ///< Match to add
    pwiz::identdata::SpectrumIdentificationItemPtr item ///< item to add the scores to
  );

  /**
   * Adds the match ranks to the SpectrumIdentificationItem
   */
  void addRanks(
    MatchCollection* match_collection, ///< Parent collection of the match
    Crux::Match* match, ///< Match to add
    pwiz::identdata::SpectrumIdentificationItemPtr item ///< item to add the ranks to
  );

 public:

  /**
   * Basic constructor
   */
  MzIdentMLWriter();

  /**
   * Basic destructor
   */
  virtual ~MzIdentMLWriter();

  /**
   * Open a file of the given name.  Replace an existing file if
   * overwrite is true, else exit if an existing file is found.
   */
  void openFile(
    const std::string& filename,
    bool overwrite
    );

  void openFile(
    const char* filename, 
    bool overwrite
  );

  /**
   * Writes out the mzid and frees the memory
   */
  void closeFile();

  /**
   * Adds the matches in the match collection to
   * the mzid objects
   */
  void addMatches(
    MatchCollection* collection ///< matches to add
  );

  /**
   * Adds the match to the mzIdentML object
   */
  void addMatch(
    MatchCollection* collection, ///< parent collection
    Crux::Match* match ///< match to add
  );

};


#endif // MZIDENTMLWRITER_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

