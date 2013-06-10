#include "TideIndexApplication.h"

#include "tide/records_to_vector-inl.h"

#include "tide/index_settings.cc"

#include "carp.h"

extern void TranslateFastaToPB(const string& fasta_filename,
                  			       const string& proteins_filename,
                               string* command_line = NULL,
                  			       pb::Header* header = NULL);
extern bool MakePeptides(pb::Header* header,
                         const string& peptides_file,
                         const string& aux_locs_file);
extern void AddTheoreticalPeaks(const vector<const pb::Protein*>& proteins,
                        				const string& input_filename,
                        				const string& output_filename);
//extern void SettingsFromFlags(pb::Header* header);
extern void AddMods(HeadedRecordReader* reader,
                    string out_file,
            		    const pb::Header& header,
            		    const vector<const pb::Protein*>& proteins);

TideIndexApplication::TideIndexApplication() {
}

TideIndexApplication::~TideIndexApplication() {
}

int TideIndexApplication::main(int argc, char** argv) {

  const char* option_list[] = {
    "enzyme",
    "digestion",
    "max-missed-cleavages",
    "max-length",
    "max-mass",
    "min-length",
    "min-mass",
    "peptides",
    "proteins",
    "monoisotopic-precursor",
    "mods-spec"
  };
  int num_options = sizeof(option_list) / sizeof(char*);
  const char* arg_list[] = {
    "protein fasta file"
  };
  int num_args = sizeof(arg_list) / sizeof(char*);
  initialize(arg_list, num_args, option_list, num_options, argc, argv);

  carp(CARP_INFO, "Running tide-index...");

  string cmd_line = "crux tide-index";
  for (int i = 1; i < argc; ++i) {
    cmd_line += " ";
    cmd_line += argv[i];
  }

  string fasta = get_string_parameter_pointer("protein fasta file");

  string enzyme = enzyme_type_to_string(get_enzyme_type_parameter("enzyme"));
  if (enzyme == "no-enzyme")
    enzyme = "none";
  FLAGS_enzyme = enzyme;

  DIGEST_T digestion = get_digest_type_parameter("digestion");
  if (enzyme == "none" && digestion == NON_SPECIFIC_DIGEST) {
    digestion = FULL_DIGEST;
  } else if (digestion != FULL_DIGEST && digestion != PARTIAL_DIGEST) {
    carp(CARP_FATAL, "'digestion' must be 'full-digest' or 'partial-digest'");
  }
  FLAGS_digestion = digest_type_to_string(digestion);
  FLAGS_max_missed_cleavages = get_int_parameter("max-missed-cleavages");
  FLAGS_max_length = get_int_parameter("max-length");
  FLAGS_max_mass = get_double_parameter("max-mass");
  FLAGS_min_length = get_int_parameter("min-length");
  FLAGS_min_mass = get_double_parameter("min-mass");

  string out_proteins = get_string_parameter_pointer("proteins");
  if (out_proteins.empty())
    out_proteins = fasta + ".protix";
  string out_peptides = get_string_parameter_pointer("peptides");
  if (out_peptides.empty())
    out_peptides = fasta + ".pepix";
  string out_aux = fasta + ".auxlocs";

  FLAGS_monoisotopic_precursor = get_boolean_parameter("monoisotopic-precursor");
  FLAGS_mods_spec = get_string_parameter_pointer("mods-spec");

  pb::Header raw_proteins_header;
  if (file_exists(out_proteins)) {
    HeadedRecordReader reader(out_proteins, &raw_proteins_header);
  } else {
    carp(CARP_INFO, "Reading %s", fasta.c_str());
    TranslateFastaToPB(fasta, out_proteins, &cmd_line, &raw_proteins_header);
  }

  pb::Header header_with_mods;
  SettingsFromFlags(&header_with_mods);
  header_with_mods.set_file_type(pb::Header::PEPTIDES);
  header_with_mods.set_command_line(cmd_line);
  pb::Header_Source* source = header_with_mods.add_source();
  source->mutable_header()->CopyFrom(raw_proteins_header);
  source->set_filename(AbsPath(out_proteins));

  pb::Header header_no_mods;
  header_no_mods.CopyFrom(header_with_mods);
  pb::ModTable* del = header_no_mods.mutable_peptides_header()->mutable_mods();
  del->mutable_variable_mod()->Clear();
  del->mutable_unique_deltas()->Clear();

  string modless_peptides = out_peptides + ".nomods.tmp";
  string peakless_peptides = out_peptides + ".nopeaks.tmp";
  bool need_mods = header_with_mods.peptides_header().mods().variable_mod_size() > 0;
  string basic_peptides = need_mods ? modless_peptides : peakless_peptides;

  if (!file_exists(basic_peptides) || !file_exists(out_aux)) {
    carp(CARP_INFO, "Computing unmodified peptides...");
    MakePeptides(&header_no_mods, basic_peptides, out_aux);
  }

  vector<const pb::Protein*> proteins;
  if (!ReadRecordsToVector<pb::Protein>(&proteins, out_proteins)) {
    carp(CARP_FATAL, "Error reading proteins file");
  }

  if (need_mods && !file_exists(peakless_peptides)) {
    carp(CARP_INFO, "Computing modified peptides...");
    HeadedRecordReader reader(modless_peptides, NULL, 1024 << 10); // 1024kb,convert to bytes = default bufsize
    AddMods(&reader, peakless_peptides, header_with_mods, proteins);
  }

  if (!file_exists(out_peptides)) {
    carp(CARP_INFO, "Precomputing theoretical spectra...");
    AddTheoreticalPeaks(proteins, peakless_peptides, out_peptides);
  }

  return 0;
}

string TideIndexApplication::getName() {
  return "tide-index";
}

string TideIndexApplication::getDescription() {
  return "Runs tide-index";
}

bool TideIndexApplication::needsOutputDirectory() {
  return true;
}

COMMAND_T TideIndexApplication::getCommand() {
  return TIDE_INDEX_COMMAND;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
