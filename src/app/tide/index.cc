// Benjamin Diament
//
// A crude cobbling together of read_fasta.cc,
// make_peptides_cmdline.cc, and peptide_peaks.cc to make a
// single-step indexer.

#include <sys/types.h>
#include <sys/stat.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include <errno.h>

#include <gflags/gflags.h>
#include "header.pb.h"
#include "records.h"
#include "records_to_vector-inl.h"
#include "peptide.h"
#include "theoretical_peak_set.h"
#include "abspath.h"

using namespace std;

#define CHECK(x) GOOGLE_CHECK(x)

DEFINE_string(fasta, "", "Input FASTA file");
DEFINE_string(proteins, "", "File of raw proteins to create, as raw_proteins.proto");
DEFINE_string(peptides, "", "File of peptides to create");
DEFINE_string(aux_locations, "", "File of peptide locations within proteins "
                                 "to create.");

DECLARE_int32(buf_size); // for peptide_mods temp files.

extern void TranslateFastaToPB(const string& fasta_filename,
			       const string& proteins_filename, string* command_line = NULL,
			       pb::Header* header = NULL);

extern bool MakePeptides(pb::Header* header, const string& peptides_file, 
                         const string& aux_locs_file);

extern void AddTheoreticalPeaks(const vector<const pb::Protein*>& proteins,
				const string& input_filename,
				const string& output_filename);

extern void SettingsFromFlags(pb::Header* header);

extern void AddMods(HeadedRecordReader* reader, string out_file,
		    const pb::Header& header,
		    const vector<const pb::Protein*>& proteins);

bool FileExists(const string& filename) {
#ifdef _MSC_VER
  // FIXME CEG Add code to check for file existance on Windows
#else
  struct stat dummy;
  if (stat(filename.c_str(), &dummy) == 0) {
    CHECK(S_ISREG(dummy.st_mode));
    fprintf(stderr, "INFO: File %s already exists.\n", filename.c_str());
    return true;
  }
  CHECK(errno == ENOENT);
#endif
  return false;
}

int main(int argc, char* argv[]) {
  string command_line;
  for (int i = 0; i < argc; i++) {
    command_line += argv[i];
    command_line += " ";
  }

  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);

  CHECK(!FLAGS_fasta.empty());
  if (FLAGS_proteins.empty())
    FLAGS_proteins = FLAGS_fasta + ".protix";
  if (FLAGS_peptides.empty())
    FLAGS_peptides = FLAGS_fasta + ".pepix";
  if (FLAGS_aux_locations.empty())
    FLAGS_aux_locations = FLAGS_fasta + ".auxlocs";

  pb::Header raw_proteins_header;
  if (FileExists(FLAGS_proteins)) {
    HeadedRecordReader reader(FLAGS_proteins, &raw_proteins_header);
  } else {
    cerr << "Reading " << FLAGS_fasta << endl;
    TranslateFastaToPB(FLAGS_fasta, FLAGS_proteins, &command_line, 
		       &raw_proteins_header);
  }

  pb::Header header_with_mods;
  SettingsFromFlags(&header_with_mods);
  header_with_mods.set_file_type(pb::Header::PEPTIDES);
  header_with_mods.set_command_line(command_line);
  pb::Header_Source* source = header_with_mods.add_source();
  source->mutable_header()->CopyFrom(raw_proteins_header);
  source->set_filename(AbsPath(FLAGS_proteins));

  pb::Header header_no_mods;
  header_no_mods.CopyFrom(header_with_mods);
  pb::ModTable* del = header_no_mods.mutable_peptides_header()->mutable_mods();
  del->mutable_variable_mod()->Clear();
  del->mutable_unique_deltas()->Clear();

  string modless_peptides = FLAGS_peptides + ".nomods.tmp"; 
  string peakless_peptides = FLAGS_peptides + ".nopeaks.tmp"; 
  bool need_mods
    = header_with_mods.peptides_header().mods().variable_mod_size() > 0;
  string basic_peptides = need_mods ? modless_peptides : peakless_peptides;

  if (!FileExists(basic_peptides) || !FileExists(FLAGS_aux_locations)) {
    cerr << "Computing unmodified peptides..." << endl;
    // TODO: pass proteins vector to MakePeptides
    CHECK(MakePeptides(&header_no_mods, basic_peptides, FLAGS_aux_locations));
  }

  vector<const pb::Protein*> proteins;
  CHECK(ReadRecordsToVector<pb::Protein>(&proteins, FLAGS_proteins));

  if (need_mods && !FileExists(peakless_peptides)) {
    cerr << "Computing modified peptides..." << endl;
    HeadedRecordReader reader(modless_peptides, NULL, FLAGS_buf_size << 10);			    
    AddMods(&reader, peakless_peptides, header_with_mods, proteins);
  }

  if (!FileExists(FLAGS_peptides)) {
    cerr << "Precomputing theoretical spectra..." << endl;
    AddTheoreticalPeaks(proteins, peakless_peptides, FLAGS_peptides);
  }

  return 0;
}
