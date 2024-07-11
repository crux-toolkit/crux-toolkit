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

// Original class to generate modified peptides. Writes to temporary files
// before merging them. As the number of possible modifications increases, the
// number of required temporary files can grow extremely large.
// This is a straigthforward, also naive implementation to combinatorially generate
// all the modified peptides. It recursively includes variable PTMs.
class ModsOutputter {//: public IModsOutputter {
 public:
  ModsOutputter(string tempDir,
                const vector<const pb::Protein*>& proteins,
                VariableModTable* var_mod_table,
                HeadedRecordWriter* final_writer,
                unsigned long long memory_limit)
    : tempDir_(tempDir),
      modPeptideCnt_(0),
      proteins_(proteins),
      mod_table_(var_mod_table),
      max_counts_(*mod_table_->MaxCounts()),
      counts_mapper_vec_(max_counts_.size(), 0),
      final_writer_(final_writer),
      count_(0),
      totalWritten_(0),
      memory_limit_(memory_limit) {
    numFiles_ = 1;
    for (int i = 0; i < max_counts_.size(); ++i) {
      counts_mapper_vec_[i] = numFiles_;
      if (max_counts_[i] == 0)
        numFiles_ *= (max_counts_[i]+2);
      else
        numFiles_ *= (max_counts_[i]+1);
    }
    InitCountsMapper();
  }

  ~ModsOutputter() {
  }
  // Return the total number of peptides written
  uint64_t Total() const { return totalWritten_; }

  void InitCountsMapper() {
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

  bool GetTempFileNames(vector<string>& filenames){
    DumpPeptides();
    for (vector<string>::iterator tf = temp_file_names_.begin(); tf != temp_file_names_.end(); ++tf) { 
      filenames.push_back(*tf);
    }
    return true;
  }

 private:
  string tempDir_;
  int numFiles_;
  int64_t modPeptideCnt_;
  
  unsigned long long memory_limit_;
  vector<pb::Peptide> pb_peptide_list_;
  unsigned long long temp_file_cnt_;
  vector<string> temp_file_names_;
  uint64_t totalWritten_;
  
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

  // Terminal modifications count as a modification and hence
  // it is taken into account the modification limit.
  void OutputMods(int pos, vector<int>& counts) {
    if (TotalMods(counts) > FLAGS_max_mods) {
      return;
    }
/*    if (pos == peptide_->length()) {
      OutputCtermMods(pos-1, counts);  // TODO: Why is this here? It seems never to be executed. 
    } else {
      if (pos == peptide_->length()-1){
        OutputCtermMods(pos, counts);
      } else {
*/        char aa = residues_[pos];
        int num_poss = mod_table_->NumPoss(aa);
        for (int i = 0; i < num_poss; ++i) {
          int poss_max_ct = mod_table_->PossMaxCt(aa, i);
          if (counts[poss_max_ct] < max_counts_[poss_max_ct]) {
            ++counts[poss_max_ct];
            int delta_index = mod_table_->PossDeltIx(aa, i);
            // We put regular variable modifications on the last amino acid
            // but then we additionall also can place terminal modifications            
            peptide_->add_modifications(mod_table_->EncodeMod(pos, delta_index));
            if (pos == peptide_->length()-1){ 
              OutputCtermMods(pos, counts);
            } else {
              OutputMods(pos+1, counts);
            }
            peptide_->mutable_modifications()->RemoveLast();
            --counts[poss_max_ct];
          }
        }
        // Having this call to OutputMods come last is, in fact, correct, but it's
        // tricky to see why. When modified peptides have equal mass, we want
        // modified positions toward the front of the peptide to appear before those
        // that come toward the end of the peptide. Having this call at the end
        // achieves that.
        
        // proceed without further mods
        if (pos == peptide_->length()-1){ 
          OutputCtermMods(pos, counts);
        } else {
          OutputMods(pos+1, counts); 
        }
  //     }
  //   }
  }

  void PlaceVariableNTermMod(int pos, char aa, mods_spec_type mod_spec, vector<int>& counts){
    int num_poss = mod_table_->NumPoss(aa, mod_spec);

    for (int i = 0; i < num_poss; ++i) {
      int poss_max_ct = mod_table_->PossMaxCt(aa, i, mod_spec);
      if (counts[poss_max_ct] < max_counts_[poss_max_ct]) {
        ++counts[poss_max_ct];
        int delta_index = mod_table_->PossDeltIx(aa, i, mod_spec);  
        peptide_->set_nterm_mod(mod_table_->EncodeMod(pos, delta_index));
        OutputMods(0, counts);
        peptide_->clear_nterm_mod();
        --counts[poss_max_ct];
      }
    }    
  }

  void OutputNtermMods(int pos, vector<int>& counts) {
    if (TotalMods(counts) > FLAGS_max_mods) {
      return;
    }
    bool any_term_modification = false;

    //add static N-terminal modifications
    //TODO: It looks like this code never runs. Static terminal mods are done 
    //through mass tables and they never make it into the possibles_ tables
    vector<char> aas{residues_[0], 'X'};
    for(auto aa = aas.begin(); aa < aas.end(); aa++){

      int num_poss = mod_table_->NumPoss(*aa, NTPEP);
      for (int i = 0; i < num_poss; ++i) {
        int poss_max_ct = mod_table_->PossMaxCt(*aa, i, NTPEP);
        if (max_counts_[poss_max_ct] == 0) {
          int delta_index = mod_table_->PossDeltIx(*aa, i, NTPEP);
          peptide_->set_nterm_mod(mod_table_->EncodeMod(pos, delta_index));
          OutputMods(0, counts);
          peptide_->clear_nterm_mod();
          any_term_modification = true;
        }
      }
    }

    if (!any_term_modification) {
      //if there were no static modificatinos add variable terminal modifications
      for(auto it = aas.begin(); it < aas.end(); it++){
        if(peptide_->first_location().pos() == 0){
          //this is protein N-terminal
            PlaceVariableNTermMod(pos, *it, NTPRO, counts);
          }
        PlaceVariableNTermMod(pos, *it, NTPEP, counts); //peptide terminal
      }
      OutputMods(0, counts);  //regular mods
    }
  }

  void PlaceCTermVariableMod(int pos, char aa, mods_spec_type mod_type, vector<int>& counts){
      int num_poss = mod_table_->NumPoss(aa, mod_type);
      for (int i = 0; i < num_poss; ++i) {
        int poss_max_ct = mod_table_->PossMaxCt(aa, i, mod_type);
        if (counts[poss_max_ct] < max_counts_[poss_max_ct]) {
          ++counts[poss_max_ct];
          int delta_index = mod_table_->PossDeltIx(aa, i, mod_type);
          peptide_->set_cterm_mod(mod_table_->EncodeMod(pos, delta_index));

          if (TotalMods(counts) >= FLAGS_min_mods) {
            peptide_->set_id(count_++);
            Write(counts);
          }
          peptide_->clear_cterm_mod();
          --counts[poss_max_ct];
        }
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
    vector<char> aas{residues_[pos], 'X'};        //we are looking up mod definitions for the current residue or X
    for(auto aa = aas.begin(); aa < aas.end(); aa++){
      int num_poss = mod_table_->NumPoss(*aa, CTPEP);
      for (int i = 0; i < num_poss; ++i) {
        int poss_max_ct = mod_table_->PossMaxCt(*aa, i, CTPEP);
        if (max_counts_[poss_max_ct] == 0) {                          //this condition is never satisfied because max_counts_ is initialized 
          int delta_index = mod_table_->PossDeltIx(*aa, i, CTPEP);    //with the var mods table and var mods max count is always >0  
          peptide_->set_cterm_mod(mod_table_->EncodeMod(pos, delta_index));

          if (TotalMods(counts) >= FLAGS_min_mods) {
            peptide_->set_id(count_++);
            Write(counts);
          }
          peptide_->clear_cterm_mod();
          any_term_modification = true;
        }
      }
    }

    if (!any_term_modification) {
      //if there were no static modifications add amino acid mods
      // char aa = residues_[pos];
      // //add any matching regular mods first
      // PlaceCTermVariableMod(pos, residues_[pos], MOD_SPEC, counts);

      //add variable c-terminal mods
      int prot_idx = peptide_->first_location().protein_id();
      auto prot_len = proteins_[prot_idx]->residues().length();

      for(auto aa = aas.begin(); aa < aas.end(); aa++) {
        //if this is protein's C_terminal
        if((peptide_->first_location().pos() + peptide_->length()) == prot_len)
          PlaceCTermVariableMod(pos, *aa, CTPRO, counts);
        PlaceCTermVariableMod(pos, *aa, CTPEP, counts);
      }

      if (TotalMods(counts) >= FLAGS_min_mods) {
        peptide_->set_id(count_++);
        Write(counts);
      }
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

  struct PbPeptideSort {
    PbPeptideSort() {}
    inline bool operator() (const pb::Peptide& x, const pb::Peptide& y) {
      return x.mass() < y.mass();
    }
  };

  void Write(const vector<int>& counts) {
    int index = DotProd(counts);
    double mass = peptide_->mass();
    peptide_->set_mass(delta_by_file_[index] + mass);
    pb_peptide_list_.push_back(*peptide_);
    peptide_->set_mass(mass);
    ++totalWritten_;
    if (totalWritten_ % 10000000 == 0) {
      carp(CARP_INFO, "Wrote %lu modified target peptides to temp files", totalWritten_);
    }    
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
   
  void DeleteTempFiles() {
    for (vector<string>::iterator tf = temp_file_names_.begin(); tf != temp_file_names_.end(); ++tf) { 
      unlink((*tf).c_str());
    }
  }

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
  ModsOutputter outputOrig(tmpDir, proteins, var_mod_table, &writer, memory_limit);

  pb::Peptide peptide;
  while (!reader->Done()) {
    CHECK(reader->Read(&peptide));
    outputOrig.Output(&peptide);
  }

  CHECK(reader->OK());
  unsigned long long peptide_num = outputOrig.Total();
  outputOrig.GetTempFileNames(temp_file_name);
  return peptide_num;
  
}

