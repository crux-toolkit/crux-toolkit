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

#include "PSMReader.h"
#include "model/Spectrum.h"
#include "model/Match.h"
#include "model/MatchCollection.h"

class PepXMLReader : public PSMReader {

 protected:
  int maxRank_;

  SpectrumZState current_zstate_; ///< keeps track of the current zstate
  Crux::Spectrum* current_spectrum_; ///< Keeps track of the current spectrum object
  Crux::Match* current_match_; ///< keeps track of the current match object
  std::string current_peptide_sequence_; ///< keeps track of the current peptide sequence
  MatchCollection* current_match_collection_; ///< keeps track of the current match collection object

  /*State variable for element tags */
  bool aminoacid_modification_open_;
  bool spectrum_query_open_; ///< are we within a spectrum_query?
  bool search_result_open_; ///< are we within a search_result element
  bool search_hit_open_;  ///< are we within a search 
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

  void aminoacidModificationOpen(const char** attr);
  void aminoacidModificationClose();

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

  void modificationInfoOpen(const char** attr);

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
    Database* decoy_database = NULL ///< the decoy protein database (can be null)
    );

  /**
   * \returns the MatchCollection resulting from the parsed xml file
   */
  MatchCollection* parse();

  /**
   * \returns the MatchCollection resulting from the parsed xml file
   */
  static MatchCollection* parse(
    const std::string& path, ///< path of the xml file
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
