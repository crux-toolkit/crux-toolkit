// Benjamin Diament
//

#include <stdio.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include <string>
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

using namespace std;

#define CHECK(x) GOOGLE_CHECK(x)

DEFINE_string(tmpfile_prefix, "/tmp/modified_peptides_partial_", "Temporary filename prefix.");
DEFINE_int32(buf_size, 1024, "Buffer size for files, in KBytes.");
DEFINE_int32(max_mods, 255, "Maximum number of modifications that can be applied "
                            "to a single peptide.");
DEFINE_int32(min_mods, 0, "Minimum number of modifications that can be applied "
                          "to a single peptide.");

static string GetTempName(int filenum) {
  char buf[10];
  sprintf(buf, "%d", filenum);
  return FLAGS_tmpfile_prefix + buf;
}

class ModsOutputter {
 public:
  unsigned long modpeptidecnt_;
  ModsOutputter(const vector<const pb::Protein*>& proteins,
		VariableModTable* var_mod_table,
		HeadedRecordWriter* final_writer)
    : proteins_(proteins),
      mod_table_(var_mod_table),
      max_counts_(*mod_table_->MaxCounts()),
      counts_mapper_vec_(max_counts_.size(), 0),
      final_writer_(final_writer),
      count_(0) {
    InitCountsMapper();
  }

  ~ModsOutputter() {
    for (int i = 0; i < writers_.size(); ++i)
      delete writers_[i];
    Merge();
  }

  void Output(pb::Peptide* peptide) {
    peptide_ = peptide;
    const pb::Location& loc = peptide->first_location();
    residues_ = proteins_[loc.protein_id()]->residues().data() + loc.pos();
    vector<int> counts(max_counts_.size(), 0);
    OutputNtermMods(0, counts);
//    OutputMods(0, counts);
  }

 private:
  void OutputMods(int pos, vector<int>& counts);
  void OutputNtermMods(int pos, vector<int>& counts);
  void OutputCtermMods(int pos, vector<int>& counts);
  void Merge();

  void InitCountsMapper() {
    int prod = 1;
    for (int i = 0; i < max_counts_.size(); ++i) {
      counts_mapper_vec_[i] = prod;
      if (max_counts_[i] == 0)
        prod *= (max_counts_[i]+2);
      else
        prod *= (max_counts_[i]+1);
    }

    writers_.resize(prod);
    for (int i = 0; i < prod; ++i) {
      writers_[i] = new RecordWriter(GetTempName(i), FLAGS_buf_size << 10);
      if (!writers_[i]->OK()) {
        // delete temporary files
        for (int j = 0; j < i; ++j)
          unlink(GetTempName(j).c_str());
        CHECK(writers_[i]->OK());
      }
    }

    const vector<double>& deltas = *mod_table_->OriginalDeltas();
    delta_by_file_.resize(prod);
    for (int i = 0; i < prod; ++i) {
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
    ++modpeptidecnt_;
    int index = DotProd(counts);
    double mass = peptide_->mass();
    peptide_->set_mass(delta_by_file_[index] + mass);
    writers_[index]->Write(peptide_);
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

//terminal modifications count as a modification and hence 
//it is taken into account the modification limit.
void ModsOutputter::OutputMods(int pos, vector<int>& counts) {
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

void ModsOutputter::OutputNtermMods(int pos, vector<int>& counts) {
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
    if (max_counts_[poss_max_ct] == 0) {
      ++counts[poss_max_ct];
      int delta_index = mod_table_->PossDeltIx(aa, i, NTPEP);
      peptide_->add_modifications(mod_table_->EncodeMod(pos, delta_index));
      OutputMods(1, counts);
      peptide_->mutable_modifications()->RemoveLast();
      --counts[poss_max_ct];
      any_term_modification = true;
    }
  }
  if (any_term_modification == false){
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

void ModsOutputter::OutputCtermMods(int pos, vector<int>& counts) {
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
    if (max_counts_[poss_max_ct] == 0) {
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
  if (any_term_modification == false){
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

struct greater_pepreader
  : public binary_function<PepReader*, PepReader*, bool> {
  bool operator()(PepReader* x, PepReader* y) {
    if (x == y)
      return false;
    return *y < *x;
  }
};

void ModsOutputter::Merge() {
  int num_files = writers_.size();
  vector<PepReader*> readers(num_files);
  for (int i = 0; i < num_files; ++i)
    readers[i] = new PepReader(GetTempName(i));

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
    unlink(GetTempName(i).c_str());
  }
}

void AddMods(HeadedRecordReader* reader, string out_file,
	     const pb::Header& header,
	     const vector<const pb::Protein*>& proteins, VariableModTable& var_mod_table) {
//  VariableModTable var_mod_table;
//  var_mod_table.Init(header.peptides_header().mods());
  CHECK(reader->OK());
  HeadedRecordWriter writer(out_file, header, FLAGS_buf_size << 10);
  CHECK(writer.OK());
  ModsOutputter outputter(proteins, &var_mod_table, &writer);
  outputter.modpeptidecnt_ = 0; 
  pb::Peptide peptide;
  while (!reader->Done()) {
    CHECK(reader->Read(&peptide));
    outputter.Output(&peptide);
  } 
//  cout << "no of modified peptides:\t" << outputter.modpeptidecnt_ << endl<<endl;
  CHECK(reader->OK());
}

void AddMods(HeadedRecordReader* reader, string out_file,
	     const pb::Header& header,
	     const vector<const pb::Protein*>& proteins) {
  VariableModTable var_mod_table;
  var_mod_table.Init(header.peptides_header().mods());
  CHECK(reader->OK());
  HeadedRecordWriter writer(out_file, header, FLAGS_buf_size << 10);
  CHECK(writer.OK());
  ModsOutputter outputter(proteins, &var_mod_table, &writer);

  pb::Peptide peptide;
  while (!reader->Done()) {
    CHECK(reader->Read(&peptide));
    outputter.Output(&peptide);
  }
  CHECK(reader->OK());
}
