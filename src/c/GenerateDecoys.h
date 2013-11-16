#ifndef GENERATE_DECOYS_H
#define GENERATE_DECOYS_H

#include <fstream>

#include "CruxApplication.h"
#include "Peptide.h"

class GenerateDecoys : public CruxApplication {

protected:
  struct massCompare {
    // Sort peptides by mass (descending)
    bool operator() (const std::string& lhs, const std::string& rhs) const {
      return Crux::Peptide::calcSequenceMass(lhs.c_str(), massType_) >
             Crux::Peptide::calcSequenceMass(rhs.c_str(), massType_);
    }
  };
  static MASS_TYPE_T massType_;

public:

  /**
   * Constructor
   */
  GenerateDecoys();

  /**
   * Destructor
   */
  ~GenerateDecoys();

  /**
   * Main method
   */
  virtual int main(int argc, char** argv);

  /**
   * Check if we can generate decoy proteins with the current settings.
   */
  static bool canGenerateDecoyProteins();

  /**
   * Given a FASTA file, read in all protein IDs/sequences and cleave them.
   * Return a map of protein IDs to digested peptides from that protein
   */
  static void readFasta(
    const std::string& fastaName,  ///< FASTA file name
    std::map< std::string, std::vector<std::string> >& outProteins, ///< map to store proteins
    std::set<std::string>& outPeptides,  ///< set of unique peptides
    std::ofstream* reversedFasta, ///< optional stream to write reversed proteins
    std::set<std::string>* outReversedPeptides  ///< optional set of peptides from rev fasta
  );

  /**
   * Reads the next protein ID and corresponding sequence from the FASTA stream
   * Returns false if no more proteins in stream
   */
  static bool getNextProtein(
    std::ifstream& fasta,  ///< FASTA stream
    std::string& outId,  ///< string to store protein ID
    std::string& outSequence ///< string to store sequence
  );

  /**
   * Cleave protein sequence using specified enzyme and store results in vector
   */
  static void cleaveProtein(
    const std::string& sequence, ///< Protein sequence to cleave
    ENZYME_T enzyme,  ///< Enzyme to use for cleavage
    DIGEST_T digest,  ///< Digestion to use for cleavage
    int missedCleavages,  ///< Maximum allowed missed cleavages
    int minLength,  //< Min length of peptides to return
    int maxLength,  //< Max length of peptides to return
    std::vector<std::string>& outPeptides ///< vector to store peptides
  );

  /**
   * Cleave protein sequence using specified enzyme and store results in vector
   * Vector also contains start location of each peptide within the protein
   */
  static void cleaveProtein(
    const std::string& sequence, ///< Protein sequence to cleave
    ENZYME_T enzyme,  ///< Enzyme to use for cleavage
    DIGEST_T digest,  ///< Digestion to use for cleavage
    int missedCleavages,  ///< Maximum allowed missed cleavages
    int minLength,  //< Min length of peptides to return
    int maxLength,  //< Max length of peptides to return
    std::vector< std::pair<std::string, int> >& outPeptides ///< vector to store peptides
  );

  /**
   * Makes a decoy from the sequence.
   * Returns false on failure, and decoyOut will be the same as seq.
   */
  static bool makeDecoy(
    const std::string& seq, ///< sequence to make decoy from
    const std::set<std::string>& targetSeqs,  ///< targets to check against
    const std::set<std::string>& decoySeqs,  ///< decoys to check against
    bool shuffle, ///< shuffle (if false, reverse)
    std::string& decoyOut ///< string to store decoy
  );

  /**
   * Shuffles the peptide sequence.
   * Returns false if no different sequence was generated
   */
  static bool shufflePeptide(
    std::string& seq  ///< Peptide sequence to shuffle
  );

  /**
   * Reverses the peptide sequence.
   * Returns false if no different sequence was generated
   */
  static bool reversePeptide(
    std::string& seq ///< Peptide sequence to reverse
  );

  /**
   * Returns the command name
   */
  virtual std::string getName();

  /**
   * Returns the command description
   */
  virtual std::string getDescription();

  /**
   * Returns whether the application needs the output directory or not. (default false)
   */
  virtual bool needsOutputDirectory();

  virtual COMMAND_T getCommand();
  
};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
