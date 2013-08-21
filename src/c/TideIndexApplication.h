#ifndef TIDEINDEXAPPLICATION_H
#define TIDEINDEXAPPLICATION_H

#include "CruxApplication.h"

#include <sys/stat.h>
#include <unistd.h>
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

  enum DECOY_TYPE {
    NONE,
    SHUFFLE,
    REVERSE
  };

  class TideIndexPeptide {
  private:
    double mass_;
    int length_;
    int proteinId_;
    int proteinPos_;
    const char* residues_;  // points at protein sequence
  public:
    TideIndexPeptide() {}
    TideIndexPeptide(double mass, int length,
                     string* proteinSeq, int proteinId, int proteinPos) {
      mass_ = mass;
      length_ = length;
      proteinId_ = proteinId;
      proteinPos_ = proteinPos;
      residues_ = proteinSeq->data() + proteinPos;
    }
    TideIndexPeptide(const TideIndexPeptide& other) {
      mass_ = other.mass_;
      length_ = other.length_;
      proteinId_ = other.proteinId_;
      proteinPos_ = other.proteinPos_;
      residues_ = other.residues_;
    }
    double getMass() const { return mass_; }
    int getLength() const { return length_; }
    int getProteinId() const { return proteinId_; }
    int getProteinPos() const { return proteinPos_; }
    string getSequence() const { return string(residues_, length_); }

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
    DECOY_TYPE decoyType,
    const std::string& fasta,
    const std::string& proteinPbFile,
    pb::Header& outProteinPbHeader,
    std::vector<TideIndexPeptide>& outPeptideHeap,
    std::vector<string*>& outProteinSequences
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

};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
