#ifndef LOCALIZEMODIFICATION_H
#define LOCALIZEMODIFICATION_H

#include "CruxApplication.h"
#include "app/TideMatchSet.h"
#include "model/Peptide.h"
#include "tide/modifications.h"
#include "raw_proteins.pb.h"
#include "peptides.pb.h"
#include "tide/peptide.h"

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
  virtual COMMAND_T getCommand() const;
  virtual bool hidden() const;

 private:
  class Results {
   public:
    Results(VariableModTable* modTable) {
      TideMatchSet::initModMap(*modTable->ParsedModTable(), ANY);
      TideMatchSet::initModMap(*modTable->ParsedNtpepModTable(), PEPTIDE_N);
      TideMatchSet::initModMap(*modTable->ParsedCtpepModTable(), PEPTIDE_C);
    }
    ~Results() {
      for (std::vector< std::pair<Crux::Peptide*, FLOAT_T> >::const_iterator i = results_.begin();
           i != results_.end();
           i++) {
        delete i->first;
      }
      for (std::set<const ModificationDefinition*>::iterator i = mods_.begin(); i != mods_.end(); i++) {
        ModificationDefinition::Remove(*i);
      }
    }
    void Add(Crux::Peptide* cruxPeptide, const Peptide* tidePeptide, FLOAT_T xcorr) {
      Crux::Peptide* peptide = new Crux::Peptide(cruxPeptide);
      vector<Crux::Modification> mods = TideMatchSet::getMods(tidePeptide);
      peptide->setMods(mods);
      for (vector<Crux::Modification>::const_iterator i = mods.begin(); i != mods.end(); i++) {
        mods_.insert(i->Definition());
      }
      results_.push_back(std::make_pair(peptide, xcorr));
    }
    void Sort() { std::sort(results_.begin(), results_.end(), Sorter()); }
    size_t Size() const { return results_.size(); }
    Crux::Peptide* Peptide(size_t i) { return results_[i].first; }
    FLOAT_T XCorr(size_t i) { return results_[i].second; }
   private:
    struct Sorter {
      bool operator() (const std::pair<Crux::Peptide*, FLOAT_T>& x,
                       const std::pair<Crux::Peptide*, FLOAT_T>& y) {
        return y.second < x.second;
      }
    };
    std::vector< std::pair<Crux::Peptide*, FLOAT_T> > results_;
    std::set<const ModificationDefinition*> mods_;
  };

  void reportProgress(uint64_t curTarget, uint64_t numTargets);
  std::vector<const pb::Protein*> createPbProteins(Crux::Peptide* peptide) const;
  std::vector<pb::Peptide> createPbPeptides(
    Crux::Match* match,
    VariableModTable* modTable,
    std::vector<pb::AuxLocation>* outAuxLocs
  ) const;
  static MODS_SPEC_TYPE_T modTypeToTide(ModPosition position);
  static double calcModMass(Crux::Match* match);
  VariableModTable* getModTable(
    Crux::Match* match
  ) const;

  std::set<int> progress_;
};

#endif

