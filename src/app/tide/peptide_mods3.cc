// Benjamin Diament
//

#include <stdio.h>
#ifndef _MSC_VER
#include <unistd.h>
#else
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

DEFINE_int32(buf_size, 1024, "Buffer size for files, in KBytes.");
DEFINE_int32(max_mods, 255, "Maximum number of modifications that can be applied "
                            "to a single peptide.");
DEFINE_int32(min_mods, 0, "Minimum number of modifications that can be applied "
                          "to a single peptide.");

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
    for (size_t i = 0; i < seq.length(); i++) {
      char aa = seq[i];
      if (i == 0) {
        int numPossNtPep = modTable->NumPoss(aa, NTPEP);
        for (int j = 0; j < numPossNtPep; j++) {
          int possDeltIx = modTable->PossDeltIx(aa, j, NTPEP);
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

class ModsOutputter {
 public:
  ModsOutputter(string tmpDir,
                const vector<const pb::Protein*>& proteins,
                VariableModTable* vmt,
                const map<int, int>& modMaxCounts,
                int maxMods,
                HeadedRecordWriter* final_writer)
    : tempDir_(tmpDir), proteins_(proteins), modTable_(vmt),
      modMaxCounts_(modMaxCounts), maxMods_(maxMods),
      writer_(final_writer), totalWritten_(0) {
  }

  ~ModsOutputter() { DeleteTempFiles(); }

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
      if (totalWritten_ % 10000 == 0) {
        carp(CARP_INFO, "Wrote %d peptides to temp files", totalWritten_);
      }
    }
  }

  // Return the total number of peptides written
  int64_t Total() const { return totalWritten_; }

  // Combine all temp files into the final file
  void Merge() {
    // delete all temp file writers, since destructor writes end-of-records marker
    for (map< int, pair<string, RecordWriter*> >::iterator i = tempFiles_.begin();
         i != tempFiles_.end();
         i++) {
      delete i->second.second;
      i->second.second = NULL;
    }
    while (!tempFiles_.empty()) {
      map< int, pair<string, RecordWriter*> >::iterator i = tempFiles_.begin();
      carp(CARP_DEBUG, "Reading temp file %s (id: %d)", i->second.first.c_str(), i->first);
      RecordReader reader(i->second.first, FLAGS_buf_size << 10);
      CHECK(reader.OK());
      vector<pb::Peptide> peptides;
      while (!reader.Done()) {
        peptides.push_back(pb::Peptide());
        reader.Read(&peptides[peptides.size() - 1]);
        CHECK(reader.OK());
      }
      DeleteTempFile(i);
      carp(CARP_DEBUG, "Read %d peptides from temp file, sorting...", peptides.size());
      std::sort(peptides.begin(), peptides.end(), PbPeptideSort());
      int64_t id = 0;
      for (vector<pb::Peptide>::iterator j = peptides.begin(); j != peptides.end(); j++) {
        j->set_id(id++);
        writer_->Write(&*j);
        CHECK(writer_->OK());
      }
    }
  }

 private:
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
    const int factor = 50; // 50 mass range per file
    RecordWriter* writer = GetTempWriter(int(peptide->mass()) / factor);
    if (!writer->Write(peptide)) {
      DeleteTempFiles();
      carp(CARP_FATAL, "I/O error writing modified peptide");
    }
    ++totalWritten_;
    if (get_verbosity_level() >= CARP_DETAILED_DEBUG) {
      Peptide pep(*peptide, proteins_);
      carp(CARP_DETAILED_DEBUG, "Wrote to temp: %s (mass %f)", pep.SeqWithMods().c_str(), peptide->mass());
    }
  }

  RecordWriter* GetTempWriter(int writerId) {
    map< int, pair<string, RecordWriter*> >::const_iterator i = tempFiles_.find(writerId);
    if (i != tempFiles_.end()) {
      return i->second.second;
    }

    string file;
    char buf[64];
    sprintf(buf, "modified_peptides_partial_%d", tempFiles_.size());
    if (!tempDir_.empty()) {
      file = FileUtils::Join(tempDir_, buf);
    } else {
#ifdef _MSC_VER
      char buf2[261];
      GetTempPath(261, buf2);
      file = FileUtils::Join(string(buf2), buf);
#else
      file = FileUtils::Join(string("/tmp/"), buf);
#endif
    }
    RecordWriter* writer = new RecordWriter(file, FLAGS_buf_size << 10);
    if (!writer->OK()) {
      DeleteTempFiles();
      CHECK(writer->OK());
    }
    tempFiles_[writerId] = make_pair(file, writer);
    return writer;
  }

  void DeleteTempFile(map< int, pair<string, RecordWriter*> >::iterator i) {
    carp(CARP_DEBUG, "Deleting temp file %s", i->second.first.c_str());
    unlink(i->second.first.c_str());
    if (i->second.second) {
      delete i->second.second;
    }
    tempFiles_.erase(i);
  }

  void DeleteTempFiles() {
    while (!tempFiles_.empty()) {
      DeleteTempFile(tempFiles_.begin());
    }
  }

  string tempDir_;
  const vector<const pb::Protein*>& proteins_;
  VariableModTable* modTable_;
  map<int, int> modMaxCounts_;
  int maxMods_;
  HeadedRecordWriter* writer_;

  map< int, pair<string, RecordWriter*> > tempFiles_; // id -> file, writer (ids must be in ascending order of mass)
  int64_t totalWritten_;
};

void AddMods(HeadedRecordReader* reader,
             string out_file,
             string tmpDir,
             const pb::Header& header,
             const vector<const pb::Protein*>& proteins,
             VariableModTable* var_mod_table = NULL) {
  VariableModTable tempTable;
  if (!var_mod_table) {
    tempTable.Init(header.peptides_header().mods());
    var_mod_table = &tempTable;
  }

  map<int, int> modMaxCounts;
  int maxMods = 0;
  const vector<int>* maxCounts = var_mod_table->MaxCounts();
  for (char c = 'A'; c <= 'Z'; c++) {
    for (int i = 0; i < var_mod_table->NumPoss(c, MOD_SPEC); i++) {
      int deltIx = var_mod_table->PossDeltIx(c, i, MOD_SPEC);
      int maxCt = maxCounts->at(deltIx);
      modMaxCounts[deltIx] = maxCt;
      maxMods += maxCt;
    }
    for (int i = 0; i < var_mod_table->NumPoss(c, NTPEP); i++) {
      int deltIx = var_mod_table->PossDeltIx(c, i, NTPEP);
      int maxCt = maxCounts->at(deltIx);
      modMaxCounts[deltIx] = maxCt;
      maxMods += maxCt;
    }
    for (int i = 0; i < var_mod_table->NumPoss(c, CTPEP); i++) {
      int deltIx = var_mod_table->PossDeltIx(c, i, CTPEP);
      int maxCt = maxCounts->at(deltIx);
      modMaxCounts[deltIx] = maxCt;
      maxMods += maxCt;
    }
  }
  if (maxMods > FLAGS_max_mods) {
    maxMods = FLAGS_max_mods;
  }

  CHECK(reader->OK());
  HeadedRecordWriter writer(out_file, header, FLAGS_buf_size << 10);
  CHECK(writer.OK());
  ModsOutputter outputter(tmpDir, proteins, var_mod_table, modMaxCounts, maxMods, &writer);
  pb::Peptide peptide;
  while (!reader->Done()) {
    CHECK(reader->Read(&peptide));
    outputter.Output(&peptide);
  }
  carp(CARP_INFO, "Created %d peptides.", outputter.Total());
  CHECK(reader->OK());
  outputter.Merge();
}

