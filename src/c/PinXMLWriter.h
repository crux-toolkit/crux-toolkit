/**
 * \file PinXMLWriter.h
 * \brief Writes search results in the TPP-compliant .pin.xml format.
 */
#ifndef PINXMLWRITER_H
#define PINXMLWRITER_H

#include <set>
#include <string>
#include <vector>
#include "objects.h"
#include "Spectrum.h"
#include "mass.h"
#include "Match.h"
#include "MatchCollection.h"
#include "SpectrumZState.h"
#include "Spectrum.h"
#include "Peptide.h"
#include <limits>

//#include "MatchReader.h"

class PinXMLWriter{

 public:
  PinXMLWriter();
  ~PinXMLWriter();
  PinXMLWriter(const char* output_file);
 
  void write(
    MatchIterator* iterator,
    Crux::Spectrum* spectrum,
    int top_rank
  );
 
  void write (
    MatchCollection* target_collection, 
    std::vector<MatchCollection*>& decoy_matches,
    Crux::Spectrum* spectrum, 
    int top_rank
  );
 
  //write pin.xml file without spectrum, when reading matches form sqt or txt file
  void write (
    MatchCollection* target_collection,
    std::vector<MatchCollection*>& decoy_psms,
    int top_rank
  );
  
  //header 
  void printHeader(
  );

  //get input file pathes
  std::string getFilePath(char* file); 

  void closeFile();
  void openFile(
    const char* filename, 
    const char* ouput_directory,
    bool overwrite
  ); 

  //feature Description 
  void printFeatureDescription(); 
  
  //footer
  void printFooter();   

  bool setProcessInfo( 
    const char* target_file_path,
    const char* decoy_file_path
  );    
 std::string absolutPath(const char* filename);  
 /*
 template<typename T>
 inline bool isInfinite(T value){
  return std::numeric_limits<T>::has_infinity &&
  value == std::numeric_limits<T>::infinity();
 } 
 */ 
 protected:

  FILE* output_file_;
  std::string enzyme_; 
  std::string decoy_file_path_;
  std::string target_file_path_;
  std::string output_file_path_; 
  std::string directory_;
  
  int precision_;
  int mass_precision_;
  int scan_number_;
  bool is_sp_; 
  bool is_decoy_; 
  std::string decoy_prefix_;
  std::set<int> charges_; 
  void init(); 

  //proccess information 
  void printProcessInfo();
 
  void printFeatures(
    Crux::Match* match, 
    bool is_sp
  );

  //write PSM
  void printPSM(
    Crux::Match* match, 
    Crux::Spectrum* spectrum, 
    bool is_decoy 
  );

  //write Peptide sequence 
  void printPeptideSequence(
    Crux::Peptide* pepitde
  ); 
  
  //write features 
  void printFeatures(
    Crux::Match* match
  );

  //write occurance 
  void printOccurence(Crux::Peptide* peptide); 

  //write FooterPSM
  void printPSMsFooter();

  //write fragSpectrumMatch 
  void printHeaderFragSpectrumMatch(int scan_number);

  //footer frag
  void printFragSpecFooter();
  
  //set decoy for sqt files 
  bool isDecoy(
    Crux::Match* match
  );

  //return id for PSM 
  string getId(
    int charge,  
    bool is_decoy,
    int scan_number,
    int rank
  ); 

  //get peptide mass with modifications
  FLOAT_T calcMassOfMods(Crux::Peptide* peptide);
 
  //calculating deltaCn and deltaLCn
  void calculateDeltaCN(map<pair<int, int>, vector<Crux::Match*> >& scan_charge_to_matches);
  void calculateDeltaCN(vector<Crux::Match*>& collection);
  void calculateDeltaCN(MatchCollection* target_collection, std::vector<MatchCollection*>& decoys);

};
#endif // PINXMLWRITER_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

