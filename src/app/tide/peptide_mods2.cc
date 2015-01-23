// Benjamin Diament
//

#include <stdio.h>
#include <unistd.h>
#include <string>
#include <vector>
#include <algorithm>
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

DEFINE_string(mods, "", "Modification specification.");
DEFINE_string(mods_table, "", "Modification specification filename. May be "
	      "given instead of --mods.");
DEFINE_string(proteins, "", "File of proteins corresponding to peptides, as raw_proteins.proto");
DEFINE_string(in_peptides, "", "File of peptides, as peptides.proto");
DEFINE_string(out_peptides, "", "File to modified peptides");
DEFINE_string(tmpfile_prefix, "/tmp/modified_peptides_partial_", "Temporary filename prefix.");

class PepReader {
 public:
  PepReader(const string& filename, int index)
    : reader_(filename), index_(index) {
    CHECK(reader_.OK());
  }

  bool operator<(const PepReader& other) {
    double mass = current_.mass();
    double other_mass = other.current_.mass();
    if (mass < other_mass)
      return true;
    if (mass > other_mass)
      return false;
    return index_ < other.index_;
  }

  bool Advance() {
    if (reader_.Done())
      return false;
    reader_.Read(&current_);
    CHECK(reader_.OK());
    return true;
  }

  void ReadHeader(pb::Header* header) {
    CHECK(!reader_.Done());
    reader_.Read(header);
    CHECK(reader_.OK());
  }

  pb::Peptide* Current() { return &current_; }

 private:
  RecordReader reader_;
  int index_;
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

string& ReplaceAll(string& context, const string& from, const string& to) {		   
  size_t look_here = 0;
  size_t found_here;
  while ((found_here = context.find(from, look_here)) != string::npos) {
    context.replace(found_here, from.size(), to);
    look_here = found_here + to.size();
  }
  return context;
}

string BashEscape(const string& s) {
  string result = s;
  ReplaceAll(result, "'", "'\"'\"'");
  return result;
}

class Merger {
 public:
  Merger(const string& final_filename, const string& tmp_prefix,
	 const pb::Header& header)
    : num_files_(0), final_filename_(final_filename), tmp_prefix_(tmp_prefix) {
    current_writer_ = new RecordWriter(GetTempName(0));
    current_writer_->Write(&header);
  }

  RecordWriter* GetNew() {
    if (num_files_ > 0) {
      delete current_writer_;
      current_writer_ = new RecordWriter(GetTempName(num_files_));
    }
    ++num_files_;
    return current_writer_;
  }

  void Merge() {
    delete current_writer_;
    if (num_files_ == 1) {
      char buf[10000];
      sprintf(buf, "mv '%s' '%s'", BashEscape(GetTempName(0)).c_str(),
	      BashEscape(final_filename_).c_str());
      CHECK(0 == system(buf));
    } else {
      PepReader* readers[num_files_];
      for (int i = 0; i < num_files_; ++i)
	readers[i] = new PepReader(GetTempName(i), i);

      // Copy header record to output
      pb::Header header;
      readers[0]->ReadHeader(&header);
      HeadedRecordWriter writer(final_filename_, header);
      CHECK(writer.OK());

      // initialize heap
      PepReader** heap_end = readers + num_files_;
      for (PepReader** reader = readers; reader < heap_end; ++reader)
	if (!(*reader)->Advance())
	  swap(*reader--, *--heap_end);
      make_heap(readers, heap_end, greater_pepreader());

      // do heap merge
      int id = 0;
      while (heap_end > readers) {
	pop_heap(readers, heap_end, greater_pepreader());
	pb::Peptide* current = (*(heap_end-1))->Current();
	current->set_id(id++);
	writer.Write(current);
	CHECK(writer.OK());
	if ((*(heap_end-1))->Advance()) {
	  push_heap(readers, heap_end, greater_pepreader());
	} else {
	  --heap_end;
	}
      }

      // delete temporary files
      for (int i = 0; i < num_files_; ++i)
        unlink(GetTempName(i).c_str());
    }
  }

 private:
  string GetTempName(int filenum) {
    char buf[10];
    sprintf(buf, "%d", filenum);
    return tmp_prefix_ + buf;
  }

  int num_files_;
  string final_filename_;
  string tmp_prefix_;
  RecordWriter* current_writer_;
};


class ModifiedPeptide {
 public:
  ModifiedPeptide(pb::Peptide* pb_peptide, ModifiedPeptide* parent,
		  int mod_code)
    : pb_peptide_(pb_peptide), parent_(parent), mod_code_(mod_code) {
    // cerr << "Creating modified peptide. this: " << this << "  Peptide : " << pb_peptide_ << ";  ID: " << pb_peptide_->id() << ";  parent: " << parent_ << ";  mod_code: " << mod_code_ << endl;
  }

  void SetAddedMass(double added_mass) { added_mass_ = added_mass; }
  double Mass() { return pb_peptide_->mass() + added_mass_; }

  void Write(RecordWriter* writer) {
    double orig_mass = pb_peptide_->mass();
    pb_peptide_->mutable_modifications()->Clear();
    AddModCodes();
    pb_peptide_->set_mass(orig_mass + added_mass_);
    pb_peptide_->set_id(num_written_++);
    writer->Write(pb_peptide_);
    pb_peptide_->set_mass(orig_mass);
  }

 private:
  void AddModCodes() {
    if (mod_code_ == -1) {
      assert(parent_ == NULL);
      return;
    }
    pb_peptide_->add_modifications(mod_code_);
    if (parent_ != NULL)
      parent_->AddModCodes();
  }

  pb::Peptide* pb_peptide_;
  ModifiedPeptide* parent_;
  int mod_code_;
  double added_mass_;

  static int num_written_;
};

int ModifiedPeptide::num_written_ = 0;

struct less_modified_peptide
  : public binary_function<ModifiedPeptide*, ModifiedPeptide*, bool> {
  bool operator()(ModifiedPeptide* x, ModifiedPeptide* y) {
    if (x == y)
      return false;
    return x->Mass() < y->Mass();
  }
};

class ModifiedPeptideCollection {
 public:
  ModifiedPeptideCollection(Merger* merger, int max, double max_mass)
    : merger_(merger), max_(max), count_(0), max_mass_(max_mass) {
  }

  ~ModifiedPeptideCollection() {
    WriteAll();
    merger_->Merge();
  }

  void AddUnderlying(pb::Peptide* peptide) {
    if (++count_ > max_) {
      WriteAll();
      count_ = 1;
    }
    underlying_.push_back(peptide);
  }

  void Add(ModifiedPeptide* peptide) {
    // TODO: prevent overweight peptides from getting here in the first place.
    if (peptide->Mass() > max_mass_)
      return;
    coll_.push_back(peptide);
  }

  void WriteAll() {
    stable_sort(coll_.begin(), coll_.end(), less_modified_peptide());
    RecordWriter* writer = merger_->GetNew();
    vector<ModifiedPeptide*>::iterator i = coll_.begin();
    for (; i != coll_.end(); ++i)
      (*i)->Write(writer);

    for (i = coll_.begin(); i != coll_.end(); ++i)
      delete *i;
    coll_.clear();

    vector<pb::Peptide*>::iterator j = underlying_.begin();
    for (; j != underlying_.end(); ++j)
      delete *j;
    underlying_.clear();
  }

 private:
  vector<pb::Peptide*> underlying_;
  vector<ModifiedPeptide*> coll_;
  Merger* merger_;
  int max_;
  int count_;
  double max_mass_;
};


class ModsOutputter {
 public:
  ModsOutputter(const vector<const pb::Protein*>& proteins,
		ModifiedPeptideCollection* collection,
		VariableModTable* var_mod_table)
    : proteins_(proteins),
      collection_(collection),
      mod_table_(var_mod_table),
      max_counts_(*mod_table_->MaxCounts()),
      count_(0) {
  }

  void Output(pb::Peptide* peptide) {
    collection_->AddUnderlying(peptide);
    peptide_ = peptide;
    const pb::Location& loc = peptide->first_location();
    residues_ = proteins_[loc.protein_id()]->residues().data() + loc.pos();
    vector<int> counts(max_counts_.size(), 0);
    OutputMods(0 /*peptide->mass()*/, peptide->length() - 1, counts, NULL);
  }

 private:
  void OutputMods(double mass_so_far, int pos, vector<int>& counts,
		  ModifiedPeptide* parent);

  const vector<const pb::Protein*>& proteins_;

  ModifiedPeptideCollection* collection_;
  VariableModTable* mod_table_;
  const vector<int>& max_counts_;
  int count_;

  pb::Peptide* peptide_;
  const char* residues_;
};


void ModsOutputter::OutputMods(double mass_so_far, int pos,
			       vector<int>& counts,
			       ModifiedPeptide* parent) {
  if (pos == -1) {
    if (parent == NULL) {
      parent = new ModifiedPeptide(peptide_, NULL, -1);
      assert(mass_so_far == 0);
    }
    parent->SetAddedMass(mass_so_far);
    collection_->Add(parent);
  } else {
    OutputMods(mass_so_far, pos-1, counts, parent); // without further mods
    char aa = residues_[pos];
    int num_poss = mod_table_->NumPoss(aa);
    for (int i = 0; i < num_poss; ++i) {
      int poss_max_ct = mod_table_->PossMaxCt(aa, i);
      if (counts[poss_max_ct] < max_counts_[poss_max_ct]) {
	++counts[poss_max_ct];
	int delta_index = mod_table_->PossDeltIx(aa, i);
	int code = mod_table_->EncodeMod(pos, delta_index);
	ModifiedPeptide* next = new ModifiedPeptide(peptide_, parent, code);
	OutputMods(mass_so_far + mod_table_->PossDelta(aa, i), pos-1, counts, next);
	--counts[poss_max_ct];
      }
    }
  }
}


void Run(HeadedRecordReader* reader, const pb::Header& header,
	 const vector<const pb::Protein*>& proteins,
	 VariableModTable* var_mod_table) {
  Merger merger(FLAGS_out_peptides, FLAGS_tmpfile_prefix, header);
  ModifiedPeptideCollection collection(&merger, 1000000,
				       header.peptides_header().max_mass());
  ModsOutputter outputter(proteins, &collection, var_mod_table);

  pb::Peptide* peptide;
  while (!reader->Done()) {
    peptide = new pb::Peptide;
    CHECK(reader->Read(peptide));
    outputter.Output(peptide);
  }
}

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);

  vector<const pb::Protein*> proteins;
  CHECK(ReadRecordsToVector<pb::Protein>(&proteins, FLAGS_proteins));

  pb::Header orig_header, new_header, mods_header;
  // TODO: read mods_header from orig
  HeadedRecordReader reader(FLAGS_in_peptides, &orig_header);
  CHECK(orig_header.file_type() == pb::Header::PEPTIDES);
  CHECK(orig_header.has_peptides_header());

  // TODO: read mods_header from orig instead of this block (i.e. use commented line)
  // pb::ModTable* mod_table = orig_header.mutable_peptides_header()->mutable_mods();
  // Also ensure that serialized deltas make their way into new_header.
  VariableModTable var_mod_table;
  pb::ModTable mod_table;
  CHECK(FLAGS_mods.empty() != FLAGS_mods_table.empty());
  if (FLAGS_mods.empty()) {
    HeadedRecordReader mods_reader(FLAGS_mods, &mods_header);
    CHECK(mods_header.file_type() == pb::Header::MOD_TABLE);
    CHECK(!mods_reader.Done());
    CHECK(mods_reader.Read(&mod_table));
    CHECK(mods_reader.Done());
    
    CHECK(MassConstants::Init(&mod_table)); // do we really need this?
    CHECK(var_mod_table.Init(mod_table));
    CHECK(var_mod_table.SerializeUniqueDeltas(&mod_table));
  } else {
    CHECK(var_mod_table.Parse(FLAGS_mods.c_str()));
    mod_table.CopyFrom(*var_mod_table.ParsedModTable());
    CHECK(MassConstants::Init(&mod_table)); // do we really need this?
  }

  new_header.set_file_type(pb::Header::PEPTIDES);
  pb::Header_PeptidesHeader* subheader = new_header.mutable_peptides_header();
  subheader->CopyFrom(orig_header.peptides_header());
  // subheader->set_has_mods(true);
  pb::Header_Source* source = new_header.add_source();
  source->mutable_header()->CopyFrom(orig_header);
  source->set_filename(AbsPath(FLAGS_in_peptides));
  // TODO: This line not necessary (are we sure?) when mod_table is
  // sourced from orig_header.
  new_header.mutable_peptides_header()->mutable_mods()->CopyFrom(mod_table);
  CHECK(reader.OK());

  Run(&reader, new_header, proteins, &var_mod_table);

  CHECK(reader.OK());

  return 0;
}
