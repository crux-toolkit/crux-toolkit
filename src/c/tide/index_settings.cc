#include <string>
#include <gflags/gflags.h>
#include "header.pb.h"
#include "records.h"
#include "modifications.h"
#include "mass_constants.h"

using namespace std;

#define CHECK(x) GOOGLE_CHECK(x)

DEFINE_double(min_mass, 200, "Peptide mass inclusion threshold");
DEFINE_double(max_mass, 7200, "Peptide mass inclusion threshold");
DEFINE_int32(min_length, 6, "Peptide length inclusion threshold");
DEFINE_int32(max_length, 50, "Peptide length inclusion threshold");

DEFINE_string(enzyme, "none", "Digestion enzyme. May be none, trypsin, "
                      "chymotrypsin, elastase, clostripain, cyanogen-bromide, "
                      "idosobenzoate, proline-endopeptidase, staph-protease, "
                      "modified-chymotrypsin, elastase-trypsin-chymotrypsin, "
                      "aspn.");
DEFINE_string(digestion, "full-digest", "Digestion completeness. "
                         "May be full-digest or partial-digest.");
DEFINE_int32(max_missed_cleavages, 0, "Maximum number of missed cleavages to "
	                              "allow in enzymatic digestion.");
DEFINE_bool(monoisotopic_precursor, false, "Use monoisotopic precursor masses "
	    "rather than average masses for residues.");

DEFINE_string(mods_spec, "", "Expression for static modifications to include. "
	      "Specify a comma-separated list of the form --mods=C+57.0,...");
DEFINE_string(mods_table, "", "Modification specification filename. May be "
	      "given instead of --mods_spec.");

void SettingsFromFlags(pb::Header* header) {
  // Fill in header as PeptidesHeader, according to command-line flags

  pb::Header_PeptidesHeader& pep_header = *header->mutable_peptides_header();
  pep_header.Clear();
  pep_header.set_min_mass(FLAGS_min_mass);
  pep_header.set_max_mass(FLAGS_max_mass);
  pep_header.set_min_length(FLAGS_min_length);
  pep_header.set_max_length(FLAGS_max_length);
  pep_header.set_monoisotopic_precursor(FLAGS_monoisotopic_precursor);
  pep_header.set_enzyme(FLAGS_enzyme);
  CHECK(!FLAGS_enzyme.empty());
  if (FLAGS_enzyme != "none") {
    CHECK(FLAGS_digestion == "full-digest"
	  || FLAGS_digestion == "partial-digest");
    pep_header.set_full_digestion(FLAGS_digestion == "full-digest");
    pep_header.set_max_missed_cleavages(FLAGS_max_missed_cleavages);
  }

  CHECK(FLAGS_mods_table.empty() || FLAGS_mods_spec.empty())
    << "Can't specify both --mods_table and --mods_spec";

  VariableModTable var_mod_table;
  if (!FLAGS_mods_table.empty()) {
    pb::Header mods_header;
    pb::ModTable mod_table;
    HeadedRecordReader mods_reader(FLAGS_mods_table, &mods_header);
    CHECK(mods_header.file_type() == pb::Header::MOD_TABLE);
    CHECK(!mods_reader.Done());
    CHECK(mods_reader.Read(&mod_table));
    CHECK(mods_reader.Done());
    
    CHECK(var_mod_table.Init(mod_table));
    CHECK(var_mod_table.SerializeUniqueDeltas(&mod_table));
    pep_header.mutable_mods()->CopyFrom(mod_table);
//    CHECK(MassConstants::Init(&mod_table));
  } else {
    CHECK(var_mod_table.Parse(FLAGS_mods_spec.c_str()));
    pep_header.mutable_mods()->CopyFrom(*var_mod_table.ParsedModTable());
//    CHECK(MassConstants::Init(var_mod_table.ParsedModTable()));
  }
}
