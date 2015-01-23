#include <gflags/gflags.h>
#include <google/protobuf/stubs/common.h>

using namespace std;

#define CHECK(x) GOOGLE_CHECK(x)

DEFINE_string(fasta, "", "Input FASTA file");
DEFINE_string(proteins, "", "File of raw proteins to create, as raw_proteins.proto");

namespace pb { class Header; }
extern void TranslateFastaToPB(const string& fasta_filename,
			       const string& proteins_filename, string* command_line = NULL,
			       pb::Header* header = NULL);

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);

  CHECK(!FLAGS_fasta.empty());
  CHECK(!FLAGS_proteins.empty());

  TranslateFastaToPB(FLAGS_fasta, FLAGS_proteins);

  return 0;
}
