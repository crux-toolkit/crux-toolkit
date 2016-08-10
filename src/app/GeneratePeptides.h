#ifndef GENERATE_PEPTIDES_H
#define GENERATE_PEPTIDES_H

#include <fstream>

#include "CruxApplication.h"
#include "model/Peptide.h"

class GeneratePeptides : public CruxApplication {

 protected:
  static MASS_TYPE_T massType_;

 public:

  class OrderedPeptide {
   public:
    explicit OrderedPeptide(const std::string& sequence):
      sequence_(sequence), sequencePtr_(NULL),
      mass_(Crux::Peptide::calcSequenceMass(sequence, massType_)) {}
    explicit OrderedPeptide(const std::string* sequence):
      sequence_(""), sequencePtr_(sequence),
      mass_(Crux::Peptide::calcSequenceMass(*sequence, massType_)) {}

    std::string Sequence() const { return sequencePtr_ ? *sequencePtr_ : sequence_; }
    unsigned int Length() const { return Sequence().length(); }
    FLOAT_T Mass() const { return mass_; }
    bool operator <(const OrderedPeptide& rhs) const { return Sequence() < rhs.Sequence(); }
   protected:
    std::string sequence_;
    const std::string* sequencePtr_;
    FLOAT_T mass_;
  };

  class CleavedPeptide : public OrderedPeptide {
   public:
    CleavedPeptide(const std::string& sequence, unsigned int position):
      OrderedPeptide(sequence), position_(position) {}
    unsigned int Position() const { return position_; }
   private:
    unsigned int position_;
  };

  /**
   * Constructor
   */
  GeneratePeptides();

  /**
   * Destructor
   */
  ~GeneratePeptides();

  /**
   * Main method
   */
  virtual int main(int argc, char** argv);

  void processFasta(
    const std::string& fastaPath,
    std::ofstream* targetList,
    const std::string& decoyFastaPath,
    std::ofstream* decoyList,
    DECOY_TYPE_T decoyType
  );

  /**
   * Check if we can generate decoy proteins with the current settings.
   */
  static bool canGenerateDecoyProteins();

  /**
   * Reads the next protein ID and corresponding sequence from the FASTA stream
   * Returns false if no more proteins in stream
   */
  static bool getNextProtein(
    std::ifstream& fasta,  ///< FASTA stream
    std::string* outId,  ///< string to store protein ID
    std::string* outSequence ///< string to store sequence
  );

  /**
   * Cleave protein sequence using specified enzyme and store results in vector
   * Vector also contains start location of each peptide within the protein
   */
  static std::vector<CleavedPeptide> cleaveProtein(
    const std::string& sequence, ///< Protein sequence to cleave
    ENZYME_T enzyme,  ///< Enzyme to use for cleavage
    DIGEST_T digest,  ///< Digestion to use for cleavage
    int missedCleavages,  ///< Maximum allowed missed cleavages
    int minLength,  //< Min length of peptides to return
    int maxLength  //< Max length of peptides to return
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
    std::string& seq,  ///< Peptide sequence to shuffle
    unsigned int maxShuffleAttempts = 6 ///< Maximum number of shuffle attempts
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
  virtual std::string getName() const;

  /**
   * Returns the command description
   */
  virtual std::string getDescription() const;

  virtual std::vector<std::string> getArgs() const;
  virtual std::vector<std::string> getOptions() const;
  virtual std::vector< std::pair<std::string, std::string> > getOutputs() const;

  /**
   * Returns whether the application needs the output directory or not. (default false)
   */
  virtual bool needsOutputDirectory() const;

  virtual COMMAND_T getCommand() const;

};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
