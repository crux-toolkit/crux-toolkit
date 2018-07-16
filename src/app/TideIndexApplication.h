#ifndef TIDEINDEXAPPLICATION_H
#define TIDEINDEXAPPLICATION_H

#include "CruxApplication.h"

#include <sys/stat.h>

#ifndef _MSC_VER
#include <unistd.h>
#endif
#include <errno.h>
#include <gflags/gflags.h>
#include "header.pb.h"
#include "tide/records.h"
#include "tide/peptide.h"
#include "tide/theoretical_peak_set.h"
#include "tide/abspath.h"
#include "TideSearchApplication.h"
#include "util/crux-utils.h"

using namespace std;

std::string getModifiedPeptideSeq(const pb::Peptide* peptide, const ProteinVec* proteins);

class TideIndexApplication : public CruxApplication {

  friend class TideSearchApplication;

 public:

  /**
   * Constructor
   */
  TideIndexApplication();

  /**
   * Destructor
   */
  ~TideIndexApplication();

  /**
   * Main method
   */
  virtual int main(int argc, char** argv);

  int main(const string& fasta, const string& index, string cmd_line = "");

  /**
   * Returns the command name
   */
  virtual string getName() const;

  /**
   * Returns the command description
   */
  virtual string getDescription() const;

  /**
   * Returns the command arguments
   */
  virtual vector<string> getArgs() const;

  /**
   * Returns the command options
   */
  virtual vector<string> getOptions() const;

  /**
   * Returns the command outputs
   */
  virtual vector< pair<string, string> > getOutputs() const;

  /**
   * Returns whether the application needs the output directory or not. (default false)
   */
  virtual bool needsOutputDirectory() const;

  virtual COMMAND_T getCommand() const;

 protected:

  class TideIndexPeptide {
   private:
    double mass_;
    int length_;
    int proteinId_;
    int proteinPos_;
    const char* residues_;  // points at protein sequence
    int decoyIdx_; // -1 if not a decoy
   public:
    TideIndexPeptide() {}
    TideIndexPeptide(double mass, int length, string* proteinSeq,
                     int proteinId, int proteinPos, int decoyIdx = -1) {
      mass_ = mass;
      length_ = length;
      proteinId_ = proteinId;
      proteinPos_ = proteinPos;
      if (decoyIdx == -1) { // It's a target
        residues_ = proteinSeq->data() + proteinPos;
      } else { // It's a decoy.
        residues_ = proteinSeq->data();
      }
      decoyIdx_ = decoyIdx;
    }
    TideIndexPeptide(const TideIndexPeptide& other) {
      mass_ = other.mass_;
      length_ = other.length_;
      proteinId_ = other.proteinId_;
      proteinPos_ = other.proteinPos_;
      residues_ = other.residues_;
      decoyIdx_ = other.decoyIdx_;
    }
    double getMass() const { return mass_; }
    int getLength() const { return length_; }
    int getProteinId() const { return proteinId_; }
    int getProteinPos() const { return proteinPos_; }
    string getSequence() const { return string(residues_, length_); }
    bool isDecoy() const { return decoyIdx_ >= 0; }
    int decoyIdx() const { return decoyIdx_; }

    friend bool operator >(
      const TideIndexPeptide& lhs, const TideIndexPeptide& rhs) {
      if (&lhs == &rhs) {
        return false;
      } else if (lhs.mass_ != rhs.mass_) {
        return lhs.mass_ > rhs.mass_;
      } else if (lhs.length_ != rhs.length_) {
        return lhs.length_ > rhs.length_;
      } else {
        int strncmpResult = strncmp(lhs.residues_, rhs.residues_, lhs.length_);
        if (strncmpResult != 0) {
          return strncmpResult > 0;
        }
      }
      if (lhs.decoyIdx_ != rhs.decoyIdx_) {
        return lhs.decoyIdx_ > rhs.decoyIdx_;
      }
      return false;
    }
    friend bool operator ==(
      const TideIndexPeptide& lhs, const TideIndexPeptide& rhs) {
      return (lhs.mass_ == rhs.mass_ && lhs.length_ == rhs.length_ &&
              strncmp(lhs.residues_, rhs.residues_, lhs.length_) == 0) &&
              lhs.decoyIdx_ == rhs.decoyIdx_;
    }
  };

  struct ProteinInfo {
    string name;
    const string* sequence;
    ProteinInfo(const string& proteinName, const string* proteinSequence)
      : name(proteinName), sequence(proteinSequence) {}
  };

  struct TargetInfo {
    ProteinInfo proteinInfo;
    int start;
    FLOAT_T mass;
    TargetInfo(const ProteinInfo& protein, int startLoc, FLOAT_T pepMass)
      : proteinInfo(protein), start(startLoc), mass(pepMass) {}
  };

  static void fastaToPb(
    const std::string& commandLine,
    const ENZYME_T enzyme,
    const DIGEST_T digestion,
    int missedCleavages,
    FLOAT_T minMass,
    FLOAT_T maxMass,
    int minLength,
    int maxLength,
    bool dups,
    MASS_TYPE_T massType,
    DECOY_TYPE_T decoyType,
    const std::string& fasta,
    const std::string& proteinPbFile,
    pb::Header& outProteinPbHeader,
    std::vector<TideIndexPeptide>& outPeptideHeap,
    std::vector<string*>& outProteinSequences,
    std::ofstream* decoyFasta
  );

  static void writePeptidesAndAuxLocs(
    std::vector<TideIndexPeptide>& peptideHeap, // will be destroyed.
    const std::string& peptidePbFile,
    const std::string& auxLocsPbFile,
    pb::Header& pbHeader
  );

  static FLOAT_T calcPepMassTide(
    const std::string& sequence,
    MASS_TYPE_T massType
  );

  static void writePbProtein(
    HeadedRecordWriter& writer,
    int id,
    const std::string& name,
    const std::string& residues,
    int targetPos = -1 // -1 if not a decoy
  );

  static void writeDecoyPbProtein(
    int id,
    const ProteinInfo& targetProteinInfo,
    std::string decoyPeptideSequence,
    int startLoc,
    HeadedRecordWriter& proteinWriter
  );

  static void getPbPeptide(
    int id,
    const TideIndexPeptide& peptide,
    pb::Peptide& outPbPeptide
  );

  static void addAuxLoc(
    int proteinId,
    int proteinPos,
    pb::AuxLocation& outAuxLoc
  );

  /**
   * Generates decoy for the target peptide, writes the decoy protein to pbProtein
   * and adds decoy to the heap.
   */
  static void generateDecoys(
    int numDecoys,
    const string& setTarget,
    std::map< const string, std::vector<const string*> >& targetToDecoy,
    set<string>* setTargets,
    set<string>* setDecoys,
    DECOY_TYPE_T decoyType,
    bool allowDups,
    unsigned int& failedDecoyCnt,
    unsigned int& decoysGenerated,
    int& curProtein,
    const ProteinInfo& proteinInfo,
    const int startLoc,
    HeadedRecordWriter& proteinWriter,
    FLOAT_T pepMass,
    vector<TideIndexPeptide>& outPeptideHeap,
    vector<string*>& outProteinSequences
  );

  virtual void processParams();
};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
