#include <string>
#include <gflags/gflags.h>
#include "header.pb.h"
#include "abspath.h"

using namespace std;

#define CHECK(x) GOOGLE_CHECK(x)

DEFINE_string(proteins, "", "File of raw proteins, raw_proteins.proto");
DEFINE_string(peptides, "", "File of peptides to create");
DEFINE_string(aux_locations, "", "File of auxiliary locations to create");

extern bool MakePeptides(pb::Header* header, const string& peptides_file, 
                         const string& aux_locs_file);

extern void SettingsFromFlags(pb::Header* header);

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);
 
  pb::Header header;
  pb::Header_Source* source = header.add_source();
  source->set_filename(AbsPath(FLAGS_proteins));
  SettingsFromFlags(&header);

  CHECK(!FLAGS_peptides.empty());
  CHECK(!FLAGS_aux_locations.empty());
  if (MakePeptides(&header, FLAGS_peptides, FLAGS_aux_locations))
    return 0;
  return 1;
}
