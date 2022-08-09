// Benjamin Diament
//

#include <stdio.h>
#ifndef _MSC_VER
#include <unistd.h>
#else
#include <process.h>
#include <windows.h>
#endif
#include <string>
#include <set>
#include <vector>
#include <algorithm>
#include <numeric>
#include <gflags/gflags.h>
#include "abspath.h"
#include "records.h"
#include "records_to_vector-inl.h"
#include "header.pb.h"
#include "raw_proteins.pb.h"
#include "peptides.pb.h"
#include "mass_constants.h"
#include "modifications.h"
#include "util/FileUtils.h"
#include "util/MathUtil.h"
#include "io/carp.h"
#include "app/tide/peptide.h"

using namespace std;

#define CHECK(x) GOOGLE_CHECK(x)

DEFINE_int32(buf_size, 1024,
  "Buffer size for files, in KBytes.");
DEFINE_int32(max_mods, 255,
  "Maximum number of modifications that can be applied to a single peptide.");
DEFINE_int32(min_mods, 0,
  "Minimum number of modifications that can be applied to a single peptide.");


static string GetTempName(const string& tempDir, int filenum) {
  char buf[64];
  sprintf(buf, "modified_peptides_partial_%d_%d", getpid(), filenum);
  if (!tempDir.empty()) {
    return FileUtils::Join(tempDir, buf);
  }
#ifdef _MSC_VER
  char buf2[261];
  GetTempPath(261, buf2);
  return FileUtils::Join(string(buf2), buf);
#else
  return FileUtils::Join(string("/tmp/"), buf);
#endif
}

class IModsOutputter {
 public:
  virtual void Output(pb::Peptide* peptide) = 0;
  virtual uint64_t Total() const = 0;
  virtual bool GetTempFileNames(vector<string>& filenames) = 0;
};

// Alternative class to generate modified peptides. Writes to temporary files
// before merging them. Bins modified peptides by mass, then writes them to the
// temporary file associated with that bin; the number of temporary files
// required is therefore bounded by the variance in peptide masses, rather than
// the number of possible modifications.
class ModsOutputterAlt : public IModsOutputter {
 public:
  ModsOutputterAlt(string tmpDir,
                   const vector<const pb::Protein*>& proteins,
                   VariableModTable* vmt,
                   HeadedRecordWriter* final_writer,
                   unsigned long long memory_limit)
    : tempDir_(tmpDir), proteins_(proteins), modTable_(vmt),
      maxMods_(0), writer_(final_writer), totalWritten_(0), 
      temp_file_cnt_(0), memory_limit_(memory_limit) {
    modMaxCounts_.clear();
    const vector<int>* maxCounts = vmt->MaxCounts();
    for (char c = 'A'; c <= 'Z'; c++) {
      for (int i = 0; i < vmt->NumPoss(c, MOD_SPEC); i++) {
        int deltIx = vmt->PossDeltIx(c, i, MOD_SPEC);
        int maxCt = maxCounts->at(deltIx);
        modMaxCounts_[deltIx] = maxCt;
        maxMods_ += maxCt;
      }
      for (int i = 0; i < vmt->NumPoss(c, NTPEP); i++) {
        int deltIx = vmt->PossDeltIx(c, i, NTPEP);
        int maxCt = maxCounts->at(deltIx);
        modMaxCounts_[deltIx] = maxCt;
        maxMods_ += maxCt;
      }
      for (int i = 0; i < vmt->NumPoss(c, NTPRO); i++) {
        int deltIx = vmt->PossDeltIx(c, i, NTPRO);
        int maxCt = maxCounts->at(deltIx);
        modMaxCounts_[deltIx] = maxCt;
        maxMods_ += maxCt;
      }
      for (int i = 0; i < vmt->NumPoss(c, CTPRO); i++) {
        int deltIx = vmt->PossDeltIx(c, i, CTPRO);
        int maxCt = maxCounts->at(deltIx);
        modMaxCounts_[deltIx] = maxCt;
        maxMods_ += maxCt;
      }
      for (int i = 0; i < vmt->NumPoss(c, CTPEP); i++) {
        int deltIx = vmt->PossDeltIx(c, i, CTPEP);
        int maxCt = maxCounts->at(deltIx);
        modMaxCounts_[deltIx] = maxCt;
        maxMods_ += maxCt;
      }
    }
    if (maxMods_ > FLAGS_max_mods) {
      maxMods_ = FLAGS_max_mods;
    }
  }

  ~ModsOutputterAlt() {
  }

  void Output(pb::Peptide* peptide) {
    if (maxMods_ < 0) {
      return;
    } else if (FLAGS_min_mods < 1) {
      WritePeptide(peptide); // write unmodified peptide
    }
    if (maxMods_ == 0) {
      return;
    }

    ResultMods resultMods(modTable_, modMaxCounts_, maxMods_, peptide, proteins_);
    double unmodifiedMass = peptide->mass();

    while (resultMods.Next()) {
      resultMods.ModifyPeptide();
      WritePeptide(peptide); 
      if (totalWritten_ % 10000000 == 0) {
        carp(CARP_INFO, "Wrote %lu peptides to temp files", totalWritten_);
      }
    }
  }

  // Return the total number of peptides written
  uint64_t Total() const { return totalWritten_; }

 private:
  unsigned long long memory_limit_;
  vector<pb::Peptide> pb_peptide_list_;
  unsigned long long temp_file_cnt_;
  vector<string> temp_file_names_;
  class ResultMods {
   private:
    class ModState { // mod state for a single residue
     public:
      ModState(vector<int>::const_iterator begin, vector<int>::const_iterator end)
        : begin_(begin), end_(end), current_(begin) {}
      vector<int>::const_iterator begin_, end_, current_;
    };

    bool NextCombination() {
      do {
        if (modifiedIndices_.Size() < FLAGS_min_mods) {
          modifiedIndices_ = MathUtil::Combination(modifiedIndices_.N(), FLAGS_min_mods);
        } else if (modifiedIndices_.HasNext()) {
          modifiedIndices_.Next();
        } else {
          size_t newK = modifiedIndices_.Size() + 1;
          if (newK > maxMods_ || newK > modifiedIndices_.N()) {
            return false;
          }
          modifiedIndices_ = MathUtil::Combination(modifiedIndices_.N(), newK);
        }
        // check that all residues in this combination can be modified
        bool valid = true;
        const vector<size_t>& modIndices = modifiedIndices_.Values();
        for (vector<size_t>::const_iterator i = modIndices.begin(); i != modIndices.end(); i++) {
          if (modIndexStates_[*i].empty()) {
            valid = false;
            break;
          }
        }
        if (valid) {
          return true;
        }
      } while (true); // compiler doesn't optimize tail recursive call here
    }

    const VariableModTable* modTable_;
    const map<int, int>& modMaxCounts_;
    int maxMods_;
    pb::Peptide* peptide_;
    double peptideUnmodifiedMass_;
    vector< vector<int> > modIndexStates_;
    vector<ModState> modStates_;

    MathUtil::Combination modifiedIndices_;
    bool needNewCombination_;

   public:
    ResultMods(const VariableModTable* modTable, const map<int, int>& modMaxCounts, int maxMods,
      pb::Peptide* peptide, const vector<const pb::Protein*>& proteins)
      : modTable_(modTable), modMaxCounts_(modMaxCounts), maxMods_(maxMods),
        peptide_(peptide), peptideUnmodifiedMass_(peptide->mass()),
        modifiedIndices_(MathUtil::Combination(peptide->length(), 0)), needNewCombination_(true) {
      const pb::Location& loc = peptide->first_location();
      string seq(proteins[loc.protein_id()]->residues().data() + loc.pos(), peptide->length());
      // generate possible mod index states
      modIndexStates_ = vector< vector<int> >(seq.length(), vector<int>());
      int prot_len = proteins[loc.protein_id()]->residues().length();
      for (size_t i = 0; i < seq.length(); i++) {
        char aa = seq[i];
        if (i == 0) {
          int numPossNtPep = modTable->NumPoss(aa, NTPEP);
          for (int j = 0; j < numPossNtPep; j++) {
            int possDeltIx = modTable->PossDeltIx(aa, j, NTPEP);
            modIndexStates_[i].push_back(possDeltIx);
          }
        }
        if (i == 0 && loc.pos() == 0) {
          int numPossNtPro = modTable->NumPoss(aa, NTPRO);
          for (int j = 0; j < numPossNtPro; j++) {
            int possDeltIx = modTable->PossDeltIx(aa, j, NTPRO);
            modIndexStates_[i].push_back(possDeltIx);
          }
        }
        if (i == seq.length() - 1) {
          int numPossCtPep = modTable->NumPoss(aa, CTPEP);
          for (int j = 0; j < numPossCtPep; j++) {
            int possDeltIx = modTable->PossDeltIx(aa, j, CTPEP);
            modIndexStates_[i].push_back(possDeltIx);
          }
        }
        if (i == seq.length() - 1 && i+loc.pos() == prot_len - 1) {
          int numPossCtPep = modTable->NumPoss(aa, CTPRO);
          for (int j = 0; j < numPossCtPep; j++) {
            int possDeltIx = modTable->PossDeltIx(aa, j, CTPRO);
            modIndexStates_[i].push_back(possDeltIx);
          }
        }

        // for(auto aa = aas.begin(); aa < aas.end(); aa++) {
          // //if this is protein's C_terminal
          // if((peptide_->first_location().pos() + peptide_->length()) == prot_len)
            // PlaceCTermVariableMod(pos, *aa, CTPRO, counts);
        
        int numPoss = modTable->NumPoss(aa, MOD_SPEC);
        for (int j = 0; j < numPoss; j++) {
          int possDeltIx = modTable->PossDeltIx(aa, j, MOD_SPEC);
          modIndexStates_[i].push_back(possDeltIx);
        }
      }

      // initialize mod state trackers
      for (vector< vector<int> >::const_iterator i = modIndexStates_.begin();
           i != modIndexStates_.end();
           i++) {
        modStates_.push_back(ModState(i->begin(), i->end()));
      }
    }

    bool Next() {
      do {
        const vector<size_t>& modIndices = modifiedIndices_.Values();
        // check if we need a new combination
        if (needNewCombination_) {
          if (!NextCombination()) {
            return false;
          }
          needNewCombination_ = false;
          // return first modified peptide for this combination
          for (vector<size_t>::const_iterator i = modIndices.begin(); i != modIndices.end(); i++) {
            ModState& modState = modStates_[*i];
            modState.current_ = modState.begin_;
          }
          if (ValidIndividualLimits()) {
            return true;
          }
        }
        for (vector<size_t>::const_iterator i = modIndices.begin(); i != modIndices.end(); i++) {
          ModState& modState = modStates_[*i];
          if (++(modState.current_) != modState.end_) {
            if (ValidIndividualLimits()) {
              return true;
            }
          } else if (i < modIndices.end() - 1) {
            modState.current_ = modState.begin_;
          }
        }
        needNewCombination_ = true;
      } while (true); // compiler doesn't optimize tail recursive call here
    }

    bool ValidIndividualLimits() const {
      map<int, int> curModCounts;
      const vector<size_t>& modIndices = modifiedIndices_.Values();
      for (vector<size_t>::const_iterator i = modIndices.begin(); i != modIndices.end(); i++) {
        const ModState& modState = modStates_[*i];
        int modIndex = *(modState.current_);
        map<int, int>::iterator j = curModCounts.find(modIndex);
        map<int, int>::const_iterator modLimit = modMaxCounts_.find(modIndex);
        if (j == curModCounts.end()) {
          curModCounts[modIndex] = 1;
        } else {
          if (++(j->second) > modLimit->second) {
            return false;
          }
        }
      }
      return true;
    }

    void ModifyPeptide() {
      peptide_->clear_modifications();
      double mass = peptideUnmodifiedMass_;
      const vector<size_t>& modIndices = modifiedIndices_.Values();
      for (vector<size_t>::const_iterator i = modIndices.begin(); i != modIndices.end(); i++) {
        int curIndex = *(modStates_[*i].current_);
        peptide_->add_modifications(modTable_->EncodeMod(*i, curIndex));
        mass += modTable_->PossDelta(curIndex);
      }
      peptide_->set_mass(mass);
    }
  };

  static double PbPeptideMass(VariableModTable* modTable, const pb::Peptide& pep) {
    int aaIndex, uniqueDeltaIndex;
    double mass = pep.mass();
    for (int i = 0; i < pep.modifications_size(); i++) {
      modTable->DecodeMod(pep.modifications(i), &aaIndex, &uniqueDeltaIndex);
      mass += modTable->PossDelta(uniqueDeltaIndex);
    }
    return mass;
  }

  struct PbPeptideSort {
    PbPeptideSort() {}
    inline bool operator() (const pb::Peptide& x, const pb::Peptide& y) {
      return x.mass() < y.mass();
    }
  };

  void WritePeptide(const pb::Peptide* peptide) {
    
    pb_peptide_list_.push_back(*peptide);
    ++totalWritten_;
    if (pb_peptide_list_.size() >= memory_limit_) {
      DumpPeptides();
    }
  }

  RecordWriter* GetTempWriter(string file) {
    RecordWriter* writer = new RecordWriter(file, FLAGS_buf_size << 10);
    if (!writer->OK()) {
      DeleteTempFiles();
      CHECK(writer->OK());
      carp(CARP_FATAL, "Failed to create temp files for modified peptides. Check free disk space.");
    }
    return writer;
  }

  // Write peptides to disk
  void DumpPeptides() {
    
    if (pb_peptide_list_.size() == 0)
      return;

    std::sort(pb_peptide_list_.begin(), pb_peptide_list_.end(), PbPeptideSort());

    string temp_file = GetTempName(tempDir_, temp_file_cnt_++);
    temp_file_names_.push_back(temp_file);
    RecordWriter* writer = GetTempWriter(temp_file); 
  
    for (vector<pb::Peptide>::iterator pept=pb_peptide_list_.begin(); pept != pb_peptide_list_.end(); ++pept) {
      if (!writer->Write(&(*pept))) {
        DeleteTempFiles();
        carp(CARP_FATAL, "I/O error writing modified peptide. Check free disk space.");
      }
    }
    delete writer;
    pb_peptide_list_.clear();
    vector<pb::Peptide> tmp;
    pb_peptide_list_.swap(tmp);
  }
  
  bool GetTempFileNames(vector<string>& filenames){
    DumpPeptides();
    for (vector<string>::iterator tf = temp_file_names_.begin(); tf != temp_file_names_.end(); ++tf) { 
      filenames.push_back(*tf);
    }
    return true;
  }
  
  void DeleteTempFiles() {
    for (vector<string>::iterator tf = temp_file_names_.begin(); tf != temp_file_names_.end(); ++tf) { 
      unlink((*tf).c_str());
    }
  }

  string tempDir_;
  const vector<const pb::Protein*>& proteins_;
  VariableModTable* modTable_;
  map<int, int> modMaxCounts_;
  int maxMods_;
  HeadedRecordWriter* writer_;

  map< int, pair<string, RecordWriter*> > tempFiles_; // id -> file, writer (ids must be in ascending order of mass)
  uint64_t totalWritten_;
};

unsigned long long AddMods(HeadedRecordReader* reader,
             string out_file,
             string tmpDir,
             const pb::Header& header,
             const vector<const pb::Protein*>& proteins,
             vector<string>& temp_file_name,
             unsigned long long memory_limit,
             VariableModTable* var_mod_table) {
  VariableModTable tempTable;
 
  HeadedRecordWriter writer(out_file, header, FLAGS_buf_size << 10);
  CHECK(writer.OK());

  memory_limit = memory_limit*1000000000/(sizeof(pb::Peptide)*2);
  ModsOutputterAlt outputAlt(tmpDir, proteins, var_mod_table, &writer, memory_limit);
  IModsOutputter* outputter;

  outputter = &outputAlt;

  pb::Peptide peptide;
  while (!reader->Done()) {
    CHECK(reader->Read(&peptide));
    outputter->Output(&peptide);
  }

  CHECK(reader->OK());
  unsigned long long peptide_num = outputter->Total();
  outputter->GetTempFileNames(temp_file_name);
  return peptide_num;
  
}

