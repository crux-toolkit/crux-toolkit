/**
 * \file PinXMLWriter.h
 * \brief Writes search results in the TPP-compliant .pin.xml format.
 */
#ifndef PINXMLWRITER_H
#define PINXMLWRITER_H

#include <string>
#include <vector>
#include "objects.h"
#include "Spectrum.h"
#include "Match.h"
#include "MatchCollection.h"
#include "SpectrumZState.h"
#include "Spectrum.h"
#include "Peptide.h"

//#include "MatchReader.h"

class PinXMLWriter{

 public:
  PinXMLWriter();
  ~PinXMLWriter();
  PinXMLWriter(const char* output_file);
 
  void write(
    MatchIterator* iterator,
    Spectrum* spectrum,
    int top_rank
  );
 
  void write (
    MatchCollection* target_collection, 
    std::vector<MatchCollection*>& decoy_matches,
    Spectrum* spectrum, 
    int top_rank
  );
 
  
  void write (
    MatchCollection* target_collection,
    std::vector<MatchCollection*>& decoy_psms 
  );
  
  //header 
  void printHeader(
    //std::string output_file_path, 
    //std::string target_file_path, 
    //std::string decoy_file_path
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
 protected:

  FILE* output_file_;
  std::string enzyme_; 
  std::string decoy_file_path_;
  std::string target_file_path_;
  
  int precision_;
  int mass_precision_;
  int scan_number_;
  bool is_sp_; 
   
  void init(); 

  //proccess information 
  void printProcessInfo(
    std::string output_file_path, 
    std::string target_file_path, 
    std::string decoy_file_path 
  );
  void printFeatures(Match* match, bool is_sp);

  //write PSM
  void printPSM(
    Match* match, 
    Spectrum* spectrum, 
    int is_decoy 
  );

  //write Peptide sequence 
  void printPeptideSequence(Peptide* pepitde); 
  
  //write features 
  void printFeatures(Match* match);

  //write occurance 
  void printOccurence(char flankC, char flankN, Peptide* peptide); 

  //write FooterPSM
  void printPSMsFooter();

  //write fragSpectrumMatch 
  void printHeaderFragSpectrumMatch(int scan_number);

  //footer frag
  void printFragSpecFooter();
  


  //return id for PSM 
  string getId(
    int charge,  
    bool is_decoy,
    int scan_number,
    int rank
  ); 
 //set charges for features 
 void setCharges(
   int charge_state, 
   bool charge1,  
   bool charge2, 
   bool charge3
 );

  //TODO: calculate deltLCn and deltCn
  //DeltaLCn=(xcorr(i)-xcorr(last))/xcorr(i)
  //DeltaCn=(xcorr(i)-xcorr(1))/xcorr(i)

};
#endif // PINXMLWRITER_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

