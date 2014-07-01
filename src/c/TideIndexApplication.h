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
#include "crux-utils.h"

using namespace std;

class TideIndexApplication : public CruxApplication {

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

  /**
   * Returns the command name
   */
  virtual string getName();

  /**
   * Returns the command description
   */
  virtual string getDescription();

  /**
   * Returns whether the application needs the output directory or not. (default false)
   */
  virtual bool needsOutputDirectory();

  virtual COMMAND_T getCommand();

protected:

  class TideIndexPeptide {
  private:
    double mass_;
    int length_;
    int proteinId_;
    int proteinPos_;
    const char* residues_;  // points at protein sequence
    bool decoy_;
  public:
    TideIndexPeptide() {}
    TideIndexPeptide(double mass, int length, string* proteinSeq,
                     int proteinId, int proteinPos, bool decoy) {
      mass_ = mass;
      length_ = length;
      proteinId_ = proteinId;
      proteinPos_ = proteinPos;
      residues_ = proteinSeq->data() + proteinPos;
      decoy_ = decoy;
    }
    TideIndexPeptide(const TideIndexPeptide& other) {
      mass_ = other.mass_;
      length_ = other.length_;
      proteinId_ = other.proteinId_;
      proteinPos_ = other.proteinPos_;
      residues_ = other.residues_;
      decoy_ = other.decoy_;
    }
    double getMass() const { return mass_; }
    int getLength() const { return length_; }
    int getProteinId() const { return proteinId_; }
    int getProteinPos() const { return proteinPos_; }
    string getSequence() const { return string(residues_, length_); }
    bool isDecoy() const { return decoy_; }

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
      return false;
    }
    friend bool operator ==(
      const TideIndexPeptide& lhs, const TideIndexPeptide& rhs) {
      return (lhs.mass_ == rhs.mass_ && lhs.length_ == rhs.length_ &&
              strncmp(lhs.residues_, rhs.residues_, lhs.length_) == 0);
    }
  };

  typedef pair<string, int> PeptideInfo;  // sequence, start location

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

  static void getPbProtein(
    int id,
    const std::string& name,
    const std::string& residues,
    pb::Protein& outPbProtein
  );

  static void getDecoyPbProtein(
    int id,
    const ProteinInfo& targetProteinInfo,
    std::string decoyPeptideSequence,
    int startLoc,
    pb::Protein& outPbProtein
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

  virtual void writeParamFile();

};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
