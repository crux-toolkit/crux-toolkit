// Benjamin Diament
//

#include <stdio.h>
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

DEFINE_string(mods, "", "Modification specification filename");
DEFINE_string(proteins, "", "File of proteins corresponding to peptides, as raw_proteins.proto");
DEFINE_string(in_peptides, "", "File of peptides, as peptides.proto");
DEFINE_string(out_peptides, "", "File to modified peptides");


class ModsOutputter {
 public:
  ModsOutputter(const vector<const pb::Protein*>& proteins,
		HeadedRecordWriter* writer,
		VariableModTable* var_mod_table)
    : proteins_(proteins),
      writer_(writer),
      mod_table_(var_mod_table),
      max_counts_(*mod_table_->MaxCounts()),
      count_(0) {
  }

  void Output(pb::Peptide* peptide) {
    peptide_ = peptide;
    const pb::Location& loc = peptide->first_location();
    residues_ = proteins_[loc.protein_id()]->residues().data() + loc.pos();
    vector<int> counts(max_counts_.size(), 0);
    OutputMods(peptide->mass(), peptide->length() - 1, counts);
  }

 private:
  void OutputMods(double mass_so_far, int pos, vector<int>& counts);

  const vector<const pb::Protein*>& proteins_;
  HeadedRecordWriter* writer_;
  VariableModTable* mod_table_;
  const vector<int>& max_counts_;
  int count_;

  pb::Peptide* peptide_;
  const char* residues_;
};


void ModsOutputter::OutputMods(double mass_so_far, int pos,
			       vector<int>& counts) {
  if (pos == -1) {
    peptide_->set_id(count_++);
    peptide_->set_mass(mass_so_far);
    writer_->Write(peptide_);
  } else {
    OutputMods(mass_so_far, pos-1, counts); // without further mods
    char aa = residues_[pos];
    int num_poss = mod_table_->NumPoss(aa);
    for (int i = 0; i < num_poss; ++i) {
      int poss_max_ct = mod_table_->PossMaxCt(aa, i);
      if (counts[poss_max_ct] < max_counts_[poss_max_ct]) {
	++counts[poss_max_ct];
	int delta_index = mod_table_->PossDeltIx(aa, i);
	peptide_->add_modifications(mod_table_->EncodeMod(pos, delta_index));
	OutputMods(mass_so_far + mod_table_->PossDelta(aa, i), pos-1, counts);
	peptide_->mutable_modifications()->RemoveLast();
	--counts[poss_max_ct];
      }
    }
  }
}

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);

  vector<const pb::Protein*> proteins;
  CHECK(ReadRecordsToVector<pb::Protein>(&proteins, FLAGS_proteins));

  pb::Header orig_header, new_header, mods_header;
  HeadedRecordReader reader(FLAGS_in_peptides, &orig_header);
  CHECK(orig_header.file_type() == pb::Header::PEPTIDES);
  CHECK(orig_header.has_peptides_header());

  HeadedRecordReader mods_reader(FLAGS_mods, &mods_header);
  CHECK(mods_header.file_type() == pb::Header::MOD_TABLE);
  pb::ModTable mod_table;
  CHECK(!mods_reader.Done());
  CHECK(mods_reader.Read(&mod_table));
  CHECK(mods_reader.Done());

  CHECK(MassConstants::Init(&mod_table));
  VariableModTable var_mod_table;
  var_mod_table.Init(mod_table);
  var_mod_table.SerializeUniqueDeltas(&mod_table);

  new_header.set_file_type(pb::Header::PEPTIDES);
  pb::Header_PeptidesHeader* subheader = new_header.mutable_peptides_header();
  subheader->CopyFrom(orig_header.peptides_header());
  // subheader->set_has_mods(true);
  pb::Header_Source* source = new_header.add_source();
  source->mutable_header()->CopyFrom(orig_header);
  source->set_filename(AbsPath(FLAGS_in_peptides));
  new_header.mutable_peptides_header()->mutable_mods()->CopyFrom(mod_table);
  HeadedRecordWriter writer(FLAGS_out_peptides, new_header);
  CHECK(reader.OK());
  CHECK(writer.OK());

  ModsOutputter outputter(proteins, &writer, &var_mod_table);

  pb::Peptide peptide;
  while (!reader.Done()) {
    CHECK(reader.Read(&peptide));
    outputter.Output(&peptide);
  }

  CHECK(reader.OK());
  CHECK(writer.OK());

  return 0;
}
