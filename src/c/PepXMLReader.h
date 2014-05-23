/**
 * \file PepXMLReader.h
 * $Revision: 1.00 $ 
 * DATE: July 11th, 2012
 * AUTHOR: Sean McIlwain
 * \brief Object for reading pep.xml.  This object will read a pepxml file,
 * creating a matchcollection object.  Use the expat library in MSToolkit.
 * 
 **************************************************************************/

#ifndef PEPXMLREADER_H
#define PEPXMLREADER_H
#include <string>
#include <vector>

#include "Spectrum.h"
#include "Match.h"
#include "MatchCollection.h"

class PepXMLReader {

 protected:

  Database* database_; ///< target database of proteins
  Database* decoy_database_; ///< decoy database of proteins
  SpectrumZState current_zstate_; ///< keeps track of the current zstate
  std::string file_path_; ///< path of the xml file
  Crux::Spectrum* current_spectrum_; ///< Keeps track of the current spectrum object
  Crux::Match* current_match_; ///< keeps track of the current match object
  std::string current_peptide_sequence_; ///< keeps track of the current peptide sequence
  MatchCollection* current_match_collection_; ///< keeps track of the current match collection object

  /*State variable for element tags */
  bool spectrum_query_open_; ///< are we within a spectrum_query?
  bool search_result_open_; ///< are we within a search_result element
  bool search_hit_open_;  ///< are we within a search 
  bool modification_info_open_; ///< are we within a modification_info element?
  bool mod_aminoacid_mass_open_; ///< are we within a mod_aminoacid_mass element?
  bool alternative_protein_open_; ///< are we within a alternative_protein element
  bool search_score_open_; ///< are we within a search_score element
  bool peptideprophet_result_open_; ///< are we within a peptideprophet_result element?
  
  /**
   * /returns the start position of the peptide sequence within the protein
   */
  int findStart(
    Crux::Protein* protein, ///< the protein to find the sequence 
    std::string peptide_sequence, ///< the peptide sequence to find
    std::string prev_aa, ///< the amino acid before the sequence in the protein
    std::string next_aa ///< the next amino acid after the sequence in the protein
  );

  /* Specific element tag handler functions */

  /**
   * Handles the spectrum_query open tag event 
   */
  void spectrumQueryOpen(
    const char** attr ///< attribute array for element
  );

  /**
   * Handles the spectrum_query close tag event
   */
  void spectrumQueryClose();

  /**
   * Handles the search_result open tag event
   */
  void searchResultOpen();
  
  /**
   * Handles the search_result close tag event
   */
  void searchResultClose();

  /**
   * Handles the search_hit open tag event
   */
  void searchHitOpen(
    const char** attr ///< atttribute array for element
  );

  /**
   * Handles the search_hit close tag event
   */
  void searchHitClose();

  /**
   * Handles the modification_info open tag event
   */
  void modificationInfoOpen(
    const char** attr ///< attribute array for element
  );
  
  /**
   * Handles the modification_info close tag event
   */
  void modificationInfoClose();
  
  /**
   * Handles the mod_aminoacid_mass open tag event
   */
  void modAminoAcidMassOpen(
    const char** attr ///< attribute array for element
  );
  
  /**
   * Handles the mod_aminoacid_mass close tag event
   */
  void modAminoAcidMassClose();
  
  /**
   * Handles the alternative_protein open tag event
   */
  void alternativeProteinOpen(
    const char** attr ///< attribute array for element
  );

  /**
   * Handles the alternative_protein close tag event
   */
  void alternativeProteinClose();

  /**
   * Handles the search_score open tag event
   */
  void searchScoreOpen(const char** attr);

  /**
   * Handles the search_score close tag event
   */
  void searchScoreClose();

  /**
   * Handles the peptideprophet_result open tag event
   */
  void peptideProphetResultOpen(
    const char** attr ///< attribute array for element
    );  

  /**
   * Handles the peptideprophet_result close tag event
   */
  void peptideProphetResultClose();  

  /**
   * Initializes the object
   */
  void init();


 public:  

  /**
   * \returns an initialized object
   */
  PepXMLReader();


  /**
   * \returns an object initialized with the file_path
   */
  PepXMLReader(
    const std::string& file_path_ ///< the path of the pep.xml file
  );

  /**
   * \returns an object initialized with the xml path, and the target,decoy databases
   */
  PepXMLReader(
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
  virtual ~PepXMLReader();
};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
