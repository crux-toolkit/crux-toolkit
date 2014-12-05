/**
 * \file PinWriter.h
 * \brief Writes search results in the .pin format.
 */
#ifndef PINWRITER_H
#define PINWRITER_H

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

class PinWriter{
 public:
  PinWriter();
  ~PinWriter();
  PinWriter(const char* output_file);
 
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
 
  //write pin file without spectrum, when reading matches form sqt or txt file
  void write (
    MatchCollection* target_collection,
    std::vector<MatchCollection*>& decoy_psms,
    int top_rank
  );
  
  void printHeader();

  void closeFile();
  void openFile(
    const char* filename, 
    const char* ouput_directory,
    bool overwrite
  ); 

 protected:
  FILE* output_file_;
  std::string decoy_file_path_;
  std::string target_file_path_;
  std::string output_file_path_; 
  std::string directory_;
  
  ENZYME_T enzyme_; 
  MASS_TYPE_T isotopic_mass_;
  int precision_;
  int mass_precision_;
  int scan_number_;
  bool is_sp_; 
  std::string decoy_prefix_;
  std::set<int> charges_; 

  void printFeatures(
    Crux::Match* match, 
    bool is_sp
  );

  void printPSM(
    Crux::Match* match, 
    Crux::Spectrum* spectrum
  );

  std::string getPeptide(Crux::Peptide* peptide);
  std::string getProteins(Crux::Peptide* peptide);
  bool isDecoy(Crux::Match* match);
  bool isInfinite(FLOAT_T x);

  std::string getId(
    Crux::Match* match,
    int scan_number
  ); 

  FLOAT_T calcMassOfMods(Crux::Peptide* peptide);
 
  void calculateDeltaCN(map<pair<int, int>, vector<Crux::Match*> >& scan_charge_to_matches);
  void calculateDeltaCN(vector<Crux::Match*>& collection);
  void calculateDeltaCN(MatchCollection* target_collection, std::vector<MatchCollection*>& decoys);
  void calculateDeltaCN(MatchCollection* matchCollection);
};
#endif // PINWRITER_H

