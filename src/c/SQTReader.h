/**
 * \file SQTReader.h
 * $Revision: 1.00 $ 
 * DATE: July 11th, 2012
 * AUTHOR: Sean McIlwain
 * \brief Object for reading pep.xml.  This object will read a pepxml file,
 * creating a matchcollection object.  Use the expat library in MSToolkit.
 * 
 **************************************************************************/

#ifndef SQTREADER_H
#define SQTREADER_H
#include <string>
#include <vector>

#include "Spectrum.h"
#include "Match.h"
#include "MatchCollection.h"

enum SQT_LINE_T {

  SQT_LINE_NONE,
  SQT_LINE_HEADER,
  SQT_LINE_SPECTRUM,
  SQT_LINE_MATCH,
  SQT_LINE_LOCUS
};


class SQTReader {

 protected:

  Database* database_; ///< target database of proteins
  Database* decoy_database_; ///< decoy database of proteins

  std::vector<std::string> headers_;

  SpectrumZState current_zstate_; ///< keeps track of the current zstate
  std::string file_path_; ///< path of the xml file
  Crux::Spectrum* current_spectrum_; ///< Keeps track of the current spectrum object
  int current_num_matches_; ///< Keeps track of the number of matches assigned to the current spectrum.
  FLOAT_T current_ln_experiment_size_;
  FLOAT_T ln_experiment_size_; 
  Crux::Match* current_match_; ///< keeps track of the current match object
  std::string current_peptide_sequence_; ///< keeps track of the current peptide sequence
  std::string current_prev_aa_; ///< keeps track of the current previous amino acid
  std::string current_next_aa_; ///< keeps track of the current next amino acid

  MatchCollection* current_match_collection_; ///< keeps track of the current match collection object

  SQT_LINE_T last_parsed_;

  /**
   * /returns the start position of the peptide sequence within the protein
   */
  int findStart(
    Crux::Protein* protein, ///< the protein to find the sequence 
    std::string peptide_sequence, ///< the peptide sequence to find
    std::string prev_aa, ///< the amino acid before the sequence in the protein
    std::string next_aa ///< the next amino acid after the sequence in the protein
  );

  /**
   * Initializes the object
   */
  void init();

  void parseHeader(std::string& line);
  void parseSpectrum(std::string& line);
  void parseMatch(std::string& line);
  void parseLocus(std::string& line);


 public:  

  /**
   * \returns an initialized object
   */
  SQTReader();


  /**
   * \returns an object initialized with the file_path
   */
  SQTReader(
    const std::string& file_path_ ///< the path of the pep.xml file
  );

  /**
   * \returns an object initialized with the xml path, and the target,decoy databases
   */
  SQTReader(
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
    Database* decoy_database ///< sets the decoy protein database
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


  /*open/close handlers*/
  friend void open_handler(void* data, const char* element, const char** attr);
  friend void close_handler(void* data, const char* el);

  /**
   * default destructor
   */
  virtual ~SQTReader();
};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
