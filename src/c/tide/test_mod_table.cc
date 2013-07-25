// Benjamin Diament
//

#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>
#include <gflags/gflags.h>
#include "records.h"
#include "header.pb.h"
#include "mass_constants.h"
#include "modifications.h"

using namespace std;

#define CHECK(x) GOOGLE_CHECK(x)

DEFINE_string(mod_table, "", "Mod Table");

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);

  pb::Header orig_header, new_header;
  HeadedRecordReader reader(FLAGS_mod_table, &orig_header);
  CHECK(orig_header.file_type() == pb::Header::MOD_TABLE);
  pb::ModTable mod_table;
  CHECK(!reader.Done());
  CHECK(reader.Read(&mod_table));
  CHECK(reader.Done());

  CHECK(MassConstants::Init(&mod_table));
  VariableModTable var_mod_table;
  var_mod_table.Init(mod_table);
  CHECK(var_mod_table.SerializeUniqueDeltas(&mod_table));

  var_mod_table.Show();

  return 0;
}
