/**
 * \file MzIdentMLReader.h
 * $Revision: 1.00 $ 
 * DATE: July 11th, 2012
 * AUTHOR: Sean McIlwain
 * \brief Object for reading mzidentml.xml.  This object will read a mzidentml file,
 * creating a matchcollection object.  Use proteowizard.
 * 
 **************************************************************************/

#ifndef MZIDENTMLREADER_H
#define MZIDENTMLREADER_H
#include <string>
#include <vector>

#include "pwiz/data/identdata/IdentDataFile.hpp"
#include "objects.h"

class MzIdentMLReader {

 protected:

  Database* database_; ///< target database of proteins
  Database* decoy_database_; ///< decoy database of proteins
  std::string file_path_; ///< path of the mzidentml file

  pwiz::identdata::IdentDataFile* pwiz_reader_; ///< proteowizard's reader for mzid.

  MatchCollection* match_collection_; ///<resulting match collection
  bool use_pass_threshold_; ///<indicator of whether to use the passThreshold attribute


  /**
   * Initializes the object
   */
  void init();

  /**
   * parses the psms from the mzid file
   */
  void parsePSMs();

  /**
   * parses the database sequences from the mzid file
   */
  void parseDatabaseSequences();

  /**
   * adds the scores for a particular psm from the mzid file to a match object
   */
  void addScores(
    const pwiz::identdata::SpectrumIdentificationItem& item, ///<proteowizard psm
    Crux::Match* match ///<our psm
  );

 public:  

  /**
   * \returns an initialized object
   */
  MzIdentMLReader();


  /**
   * \returns an object initialized with the file_path
   */
  MzIdentMLReader(
    const std::string& file_path_ ///< the path of the pep.xml file
  );

  /**
   * \returns an object initialized with the xml path, and the target,decoy databases
   */
  MzIdentMLReader(
    const std::string& file_path_, ///< the path of the pep.xml
    Database* database, ///< the protein database
    Database* decoy_database=NULL ///< the decoy protein database (can be null)
    );

  /**
   * sets the target protein database
   */
  void setDatabase(
    Database* database ///< the target protein database
  );

  /**
   * sets the decoy protein database
   */
  void setDecoyDatabase(
    Database* decoy_database ///<  the decoy protein database
  );

  /**
   * \returns the MatchCollection resulting from the parsed xml file
   */
  MatchCollection* parse();

  /**
   * \returns the MatchCollection resulting from the parsed xml file
   */
  static MatchCollection* parse(
    const char* path, ///< path of the xml file
    Database* database, ///< target protein database
    Database* decoy_database ///< decoy protein database (can be null)
  );

  /**
   * default destructor
   */
  virtual ~MzIdentMLReader();
};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
