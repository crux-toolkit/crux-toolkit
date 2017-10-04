#ifndef LOCALIZEMODIFICATION_H
#define LOCALIZEMODIFICATION_H

#include "CruxApplication.h"
#include "model/Modification.h"
#include "tide/modifications.h"
#include "raw_proteins.pb.h"

class LocalizeModificationApplication : public CruxApplication {
 public:
  LocalizeModificationApplication();
  virtual ~LocalizeModificationApplication();

  virtual int main(int argc, char** argv);
  virtual std::string getName() const;
  virtual std::string getDescription() const;
  virtual std::vector<std::string> getArgs() const;
  virtual std::vector<std::string> getOptions() const;
  virtual std::vector< std::pair<std::string, std::string> > getOutputs() const;
  virtual bool needsOutputDirectory() const;
  virtual bool hidden() const;

 private:
  void reportProgress(int curTarget, int numTargets);
  void writeIndexProteins(
    const std::string& proteinsFile,
    Crux::Peptide* peptide
  ) const;
  void writeIndexPeptides(
    const std::string& peptidesFile,
    const std::string& auxLocsFile,
    Crux::Match* match
  ) const;
  static MODS_SPEC_TYPE_T modTypeToTide(ModPosition position);
  static double calcModMass(Crux::Match* match);
  VariableModTable* getModTable(
    Crux::Match* match
  ) const;

  std::set<int> progress_;
};

#endif

