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

DEFINE_int32(buf_size, 1024,
  "Buffer size for files, in KBytes.");
DEFINE_int32(max_mods, 255,
  "Maximum number of modifications that can be applied to a single peptide.");
DEFINE_int32(min_mods, 0,
  "Minimum number of modifications that can be applied to a single peptide.");
DEFINE_int32(modsoutputter_file_threshold, 1000,
  "Maximum number of temporary files that would be opened by ModsOutputter "
  "before switching to ModsOutputterAlt.");

static string GetTempName(const string& tempDir, int filenum) {
  char buf[64];
  sprintf(buf, "modified_peptides_partial_%d", filenum);
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
  virtual int64_t Total() const = 0;
};

// Original class to generate modified peptides. Writes to temporary files
// before merging them. As the number of possible modifications increases, the
// number of required temporary files can grow extremely large.
class ModsOutputter : public IModsOutputter {
 public:
  ModsOutputter(string tmpDir,
                const vector<const pb::Protein*>& proteins,
                VariableModTable* var_mod_table,
                HeadedRecordWriter* final_writer)
    : tmpDir_(tmpDir),
      modPeptideCnt_(0),
      proteins_(proteins),
      mod_table_(var_mod_table),
      max_counts_(*mod_table_->MaxCounts()),
      counts_mapper_vec_(max_counts_.size(), 0),
      final_writer_(final_writer),
      count_(0) {
    numFiles_ = 1;
    for (int i = 0; i < max_counts_.size(); ++i) {
      counts_mapper_vec_[i] = numFiles_;
      if (max_counts_[i] == 0)
        numFiles_ *= (max_counts_[i]+2);
      else
        numFiles_ *= (max_counts_[i]+1);
    }
  }

  ~ModsOutputter() {
    for (int i = 0; i < writers_.size(); ++i) {
      delete writers_[i];
    }
    if (modPeptideCnt_ > 0) {
      Merge();
    }
  }

  int NumFiles() const {
    return numFiles_;
  }

  int64_t Total() const {
    return modPeptideCnt_;
  }

  void InitCountsMapper() {
    writers_.resize(numFiles_);
    if (numFiles_ > 100) {
      carp(CARP_INFO, "Opening %d files for modifications.", numFiles_);
    }

    for (int i = 0; i < numFiles_; ++i) {
      writers_[i] = new RecordWriter(GetTempName(tmpDir_, i), FLAGS_buf_size << 10);
      if (!writers_[i]->OK()) {
        // delete temporary files
        for (int j = 0; j < i; ++j)
          unlink(GetTempName(tmpDir_, j).c_str());
        CHECK(writers_[i]->OK());
      }
    }

    const vector<double>& deltas = *mod_table_->OriginalDeltas();
    delta_by_file_.resize(numFiles_);
    for (int i = 0; i < numFiles_; ++i) {
      double total_delta = 0;
      int x = i;
      for (int j = max_counts_.size() - 1; j >= 0; --j) {
        int digit = x / counts_mapper_vec_[j];
        x %= counts_mapper_vec_[j];
        total_delta += deltas[j] * digit;
      }
      delta_by_file_[i] = total_delta;
    }
  }

  void Output(pb::Peptide* peptide) {
    peptide_ = peptide;
    const pb::Location& loc = peptide->first_location();
    residues_ = proteins_[loc.protein_id()]->residues().data() + loc.pos();
    vector<int> counts(max_counts_.size(), 0);
    OutputNtermMods(0, counts);
  }

 private:
  string tmpDir_;
  int numFiles_;
  int64_t modPeptideCnt_;

  class PepReader {
   public:
    PepReader(const string& filename)
      : reader_(filename, FLAGS_buf_size << 10) {
      CHECK(reader_.OK());
    }

    bool operator<(const PepReader& other) {
      double mass = current_.mass();
      double other_mass = other.current_.mass();
      if (mass < other_mass)
        return true;
      if (mass > other_mass)
        return false;
      return current_.id() < other.current_.id();
    }

    bool Advance() {
      if (reader_.Done())
        return false;
      reader_.Read(&current_);
      CHECK(reader_.OK());
      return true;
    }

    pb::Peptide* Current() { return &current_; }

   private:
    RecordReader reader_;
    pb::Peptide current_;
  };

  struct greater_pepreader : public binary_function<PepReader*, PepReader*, bool> {
    bool operator()(PepReader* x, PepReader* y) {
      return x != y && *y < *x;
    }
  };

  //terminal modifications count as a modification and hence
  //it is taken into account the modification limit.
  void OutputMods(int pos, vector<int>& counts) {
    if (TotalMods(counts) > FLAGS_max_mods) {
      return;
    }
    if (pos == peptide_->length()) {
      OutputCtermMods(pos-1, counts);
    } else {
      if (pos == peptide_->length()-1){
        OutputCtermMods(pos, counts);
      } else {
        char aa = residues_[pos];
        int num_poss = mod_table_->NumPoss(aa);
        for (int i = 0; i < num_poss; ++i) {
          int poss_max_ct = mod_table_->PossMaxCt(aa, i);
          if (counts[poss_max_ct] < max_counts_[poss_max_ct]) {
            ++counts[poss_max_ct];
            int delta_index = mod_table_->PossDeltIx(aa, i);
            peptide_->add_modifications(mod_table_->EncodeMod(pos, delta_index));
            OutputMods(pos+1, counts);
            peptide_->mutable_modifications()->RemoveLast();
            --counts[poss_max_ct];
          }
        }
        // Having this call to OutputMods come last is, in fact, correct, but it's
        // tricky to see why. When modified peptides have equal mass, we want
        // modified positions toward the front of the peptide to appear before those
        // that come toward the end of the peptide. Having this call at the end
        // achieves that.
        OutputMods(pos+1, counts); // without further mods
      }
    }
  }

  void OutputNtermMods(int pos, vector<int>& counts) {
    if (TotalMods(counts) > FLAGS_max_mods) {
      return;
    }
    bool any_term_modification = false;

    //add static N-terminal modifications
    char aa = residues_[0];
    int num_poss = mod_table_->NumPoss(aa, NTPEP);
    for (int i = 0; i < num_poss; ++i) {
      int poss_max_ct = mod_table_->PossMaxCt(aa, i, NTPEP);
      if (max_counts_[poss_max_ct] == 0) {
        int delta_index = mod_table_->PossDeltIx(aa, i, NTPEP);
        peptide_->add_modifications(mod_table_->EncodeMod(pos, delta_index));
        OutputMods(1, counts);
        peptide_->mutable_modifications()->RemoveLast();
        any_term_modification = true;
      }
    }
    aa = 'X';
    num_poss = mod_table_->NumPoss(aa, NTPEP);
    for (int i = 0; i < num_poss; ++i) {
      int poss_max_ct = mod_table_->PossMaxCt(aa, i, NTPEP);
      if (max_counts_[poss_max_ct] == 0) {
        int delta_index = mod_table_->PossDeltIx(aa, i, NTPEP);
        peptide_->add_modifications(mod_table_->EncodeMod(pos, delta_index));
        OutputMods(1, counts);
        peptide_->mutable_modifications()->RemoveLast();
        any_term_modification = true;
      }
    }
    if (!any_term_modification) {
      //if there were no static modificatinos add variable terminal modifications
      aa = residues_[0];
      num_poss = mod_table_->NumPoss(aa, NTPEP);
      for (int i = 0; i < num_poss; ++i) {
        int poss_max_ct = mod_table_->PossMaxCt(aa, i, NTPEP);
        if (max_counts_[poss_max_ct] == 1) {
          ++counts[poss_max_ct];
          int delta_index = mod_table_->PossDeltIx(aa, i, NTPEP);
          peptide_->add_modifications(mod_table_->EncodeMod(pos, delta_index));
          OutputMods(1, counts);
          peptide_->mutable_modifications()->RemoveLast();
          --counts[poss_max_ct];
          any_term_modification = true;
        }
      }
      aa = 'X';
      num_poss = mod_table_->NumPoss(aa, NTPEP);
      for (int i = 0; i < num_poss; ++i) {
        int poss_max_ct = mod_table_->PossMaxCt(aa, i, NTPEP);
        if (max_counts_[poss_max_ct] == 1) {
          ++counts[poss_max_ct];
          int delta_index = mod_table_->PossDeltIx(aa, i, NTPEP);
          peptide_->add_modifications(mod_table_->EncodeMod(pos, delta_index));
          OutputMods(1, counts);
          peptide_->mutable_modifications()->RemoveLast();
          --counts[poss_max_ct];
          any_term_modification = true;
        }
      }
      OutputMods(0, counts);
    }
  }

  void OutputCtermMods(int pos, vector<int>& counts) {
    int total = TotalMods(counts);
    if (total > FLAGS_max_mods) {
      return;
    } else if (total == FLAGS_max_mods) {
      if (total >= FLAGS_min_mods) {
        peptide_->set_id(count_++);
        Write(counts);
      }
      return;
    }

    bool any_term_modification = false;
    char aa = residues_[pos];
    int num_poss = mod_table_->NumPoss(aa, CTPEP);
    for (int i = 0; i < num_poss; ++i) {
      int poss_max_ct = mod_table_->PossMaxCt(aa, i, CTPEP);
      if (max_counts_[poss_max_ct] == 0) {
        int delta_index = mod_table_->PossDeltIx(aa, i, CTPEP);
        peptide_->add_modifications(mod_table_->EncodeMod(pos, delta_index));

        if (TotalMods(counts) >= FLAGS_min_mods) {
          peptide_->set_id(count_++);
          Write(counts);
        }

        peptide_->mutable_modifications()->RemoveLast();
        any_term_modification = true;
      }
    }
    aa = 'X';
    num_poss = mod_table_->NumPoss(aa, CTPEP);
    for (int i = 0; i < num_poss; ++i) {
      int poss_max_ct = mod_table_->PossMaxCt(aa, i, CTPEP);
      if (max_counts_[poss_max_ct] == 0) {
        int delta_index = mod_table_->PossDeltIx(aa, i, CTPEP);
        peptide_->add_modifications(mod_table_->EncodeMod(pos, delta_index));

        if (TotalMods(counts) >= FLAGS_min_mods) {
          peptide_->set_id(count_++);
          Write(counts);
        }

        peptide_->mutable_modifications()->RemoveLast();
        any_term_modification = true;
      }
    }
    if (!any_term_modification) {
      //if there were no static modifications add amino acid mods
      char aa = residues_[pos];
      int num_poss = mod_table_->NumPoss(aa);
      for (int i = 0; i < num_poss; ++i) {
        int poss_max_ct = mod_table_->PossMaxCt(aa, i);
        if (counts[poss_max_ct] < max_counts_[poss_max_ct]) {
          ++counts[poss_max_ct];
          int delta_index = mod_table_->PossDeltIx(aa, i);
          peptide_->add_modifications(mod_table_->EncodeMod(pos, delta_index));

          if (TotalMods(counts) >= FLAGS_min_mods) {
            peptide_->set_id(count_++);
            Write(counts);
          }

          peptide_->mutable_modifications()->RemoveLast();
          --counts[poss_max_ct];
        }
      }

      //add variable c-terminal mods
      aa = residues_[pos];
      num_poss = mod_table_->NumPoss(aa, CTPEP);
      for (int i = 0; i < num_poss; ++i) {
        int poss_max_ct = mod_table_->PossMaxCt(aa, i, CTPEP);
        if (max_counts_[poss_max_ct] == 1) {
          ++counts[poss_max_ct];
          int delta_index = mod_table_->PossDeltIx(aa, i, CTPEP);
          peptide_->add_modifications(mod_table_->EncodeMod(pos, delta_index));

          if (TotalMods(counts) >= FLAGS_min_mods) {
            peptide_->set_id(count_++);
            Write(counts);
          }

          peptide_->mutable_modifications()->RemoveLast();
          --counts[poss_max_ct];
          any_term_modification = true;
        }
      }
      aa = 'X';
      num_poss = mod_table_->NumPoss(aa, CTPEP);
      for (int i = 0; i < num_poss; ++i) {
        int poss_max_ct = mod_table_->PossMaxCt(aa, i, CTPEP);
        if (max_counts_[poss_max_ct] == 1) {
          ++counts[poss_max_ct];
          int delta_index = mod_table_->PossDeltIx(aa, i, CTPEP);
          peptide_->add_modifications(mod_table_->EncodeMod(pos, delta_index));

          if (TotalMods(counts) >= FLAGS_min_mods) {
            peptide_->set_id(count_++);
            Write(counts);
          }

          peptide_->mutable_modifications()->RemoveLast();
          --counts[poss_max_ct];
          any_term_modification = true;
        }
      }
      if (TotalMods(counts) >= FLAGS_min_mods) {
        peptide_->set_id(count_++);
        Write(counts);
      }
    }
  }

  void Merge() {
    int num_files = writers_.size();
    vector<PepReader*> readers(num_files);
    for (int i = 0; i < num_files; ++i)
      readers[i] = new PepReader(GetTempName(tmpDir_, i));

    // initialize heap
    PepReader** heap_end = &(readers[0]) + num_files;
    for (PepReader** reader = &(readers[0]); reader < heap_end; ++reader)
      if (!(*reader)->Advance())
        swap(*reader--, *--heap_end);
    make_heap(&(readers[0]), heap_end, greater_pepreader());

    // do heap merge
    int id = 0;
#ifndef NDEBUG
    double last_mass = 0.0;
#endif
    while (heap_end > &(readers[0])) {
      pop_heap(&(readers[0]), heap_end, greater_pepreader());
      pb::Peptide* current = (*(heap_end-1))->Current();
      current->set_id(id++);
#ifndef NDEBUG
      assert(current->mass() >= last_mass);
      last_mass = current->mass();
#endif
      final_writer_->Write(current);
      CHECK(final_writer_->OK());
      if ((*(heap_end-1))->Advance()) {
        push_heap(&(readers[0]), heap_end, greater_pepreader());
      } else {
        --heap_end;
      }
    }

    // delete temporary files
    for (int i = 0; i < num_files; ++i) {
      delete readers[i];
      unlink(GetTempName(tmpDir_, i).c_str());
    }
  }

  int TotalMods(const vector<int>& counts) {
    return accumulate(counts.begin(), counts.end(), 0);
  }

  int DotProd(const vector<int>& counts) {
    int dot = 0;
    for (int i = 0; i < counts.size(); ++i)
      dot += counts_mapper_vec_[i] * counts[i];
    return dot;
  }

  RecordWriter* Write(const vector<int>& counts) {
    ++modPeptideCnt_;
    int index = DotProd(counts);
    double mass = peptide_->mass();
    peptide_->set_mass(delta_by_file_[index] + mass);
    if (!writers_[index]->Write(peptide_)) {
      carp(CARP_FATAL, "I/O error writing modifications");
    }
    peptide_->set_mass(mass);
    return writers_[index];
  }

  const vector<const pb::Protein*>& proteins_;
  VariableModTable* mod_table_;
  const vector<int>& max_counts_;
  vector<int> counts_mapper_vec_;
  vector<RecordWriter*> writers_;
  vector<double> delta_by_file_;
  HeadedRecordWriter* final_writer_;
  int count_;

  pb::Peptide* peptide_;
  const char* residues_;
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
                   HeadedRecordWriter* final_writer)
    : tempDir_(tmpDir), proteins_(proteins), modTable_(vmt),
      maxMods_(0), writer_(final_writer), totalWritten_(0) {
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
    // delete all temp file writers, since destructor writes end-of-records marker
    for (map< int, pair<string, RecordWriter*> >::iterator i = tempFiles_.begin();
         i != tempFiles_.end();
         i++) {
      delete i->second.second;
      i->second.second = NULL;
    }
    Merge();
    DeleteTempFiles();
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
      if (totalWritten_ % 10000 == 0) {
        carp(CARP_INFO, "Wrote %d peptides to temp files", totalWritten_);
      }
    }
  }

  // Return the total number of peptides written
  int64_t Total() const { return totalWritten_; }

 private:
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

    string file = GetTempName(tempDir_, tempFiles_.size());
    RecordWriter* writer = new RecordWriter(file, FLAGS_buf_size << 10);
    if (!writer->OK()) {
      DeleteTempFiles();
      CHECK(writer->OK());
    }
    tempFiles_[writerId] = make_pair(file, writer);
    return writer;
  }

  // Combine all temp files into the final file
  void Merge() {
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

  CHECK(reader->OK());
  HeadedRecordWriter writer(out_file, header, FLAGS_buf_size << 10);
  CHECK(writer.OK());

  ModsOutputter outputOrig(tmpDir, proteins, var_mod_table, &writer);
  ModsOutputterAlt outputAlt(tmpDir, proteins, var_mod_table, &writer);
  IModsOutputter* outputter;

  if (outputOrig.NumFiles() <= FLAGS_modsoutputter_file_threshold) {
    outputOrig.InitCountsMapper();
    outputter = &outputOrig;
  } else {
    // Switch to alternate ModsOutputter if the regular one would open too many files
    carp(CARP_DEBUG, "Using alternate ModsOutputter, original version would open %d files",
         outputOrig.NumFiles());
    outputter = &outputAlt;
  }

  pb::Peptide peptide;
  while (!reader->Done()) {
    CHECK(reader->Read(&peptide));
    outputter->Output(&peptide);
  }
  carp(CARP_INFO, "Created %d peptides.", outputter->Total());
  CHECK(reader->OK());
}

