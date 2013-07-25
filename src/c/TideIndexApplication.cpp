#include <cstdio>
#include "carp.h"
#include "CarpStreamBuf.h"
#include "TideIndexApplication.h"

#include "tide/modifications.h"
#include "tide/records_to_vector-inl.h"

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
    "missed-cleavages",
    "max-length",
    "max-mass",
    "min-length",
    "min-mass",
    "monoisotopic-precursor",
    "mods-spec",
    "output-dir",
    "overwrite",
    "parameter-file",
    "verbosity"
  };

  const string default_cysteine = "C+57.0214637206";

  // Crux command line parsing
  int num_options = sizeof(option_list) / sizeof(char*);
  const char* arg_list[] = {
    "protein fasta file",
    "index name"
  };
  int num_args = sizeof(arg_list) / sizeof(char*);
  initialize(arg_list, num_args, option_list, num_options, argc, argv);

  carp(CARP_INFO, "Running tide-index...");

  // Build command line string
  string cmd_line = "crux tide-index";
  for (int i = 1; i < argc; ++i) {
    cmd_line += " ";
    cmd_line += argv[i];
  }

  // Set up output paths
  string fasta = get_string_parameter_pointer("protein fasta file");
  string index = get_string_parameter_pointer("index name");
  bool overwrite = get_boolean_parameter("overwrite");

  if (!file_exists(fasta)) {
    carp(CARP_FATAL, "Fasta file %s does not exist", fasta.c_str());
  }

  string out_proteins = index + "/" + "protix";
  string out_peptides = index + "/" + "pepix";
  string out_aux = index + "/" + "auxlocs";
  string modless_peptides = out_peptides + ".nomods.tmp";
  string peakless_peptides = out_peptides + ".nopeaks.tmp";

  if (create_output_directory(index.c_str(), overwrite) != 0) {
    carp(CARP_FATAL, "Error creating index directory");
  } else if (file_exists(out_proteins) ||
             file_exists(out_peptides) ||
             file_exists(out_aux)) {
    if (overwrite) {
      carp(CARP_DEBUG, "Cleaning old index file(s)");
      remove(out_proteins.c_str());
      remove(out_peptides.c_str());
      remove(out_aux.c_str());
      remove(modless_peptides.c_str());
      remove(peakless_peptides.c_str());
    } else {
      carp(CARP_FATAL, "Index file(s) already exist, use --overwrite T or a "
                       "different index name");
    }
  }

  // Reroute stderr
  CarpStreamBuf buffer;
  streambuf* old = cerr.rdbuf();
  cerr.rdbuf(&buffer);

  // Start tide-index
  pb::Header raw_proteins_header;
  carp(CARP_INFO, "Reading %s", fasta.c_str());
  TranslateFastaToPB(fasta, out_proteins, &cmd_line, &raw_proteins_header);

  pb::Header header_with_mods;

  // Set up peptides header
  pb::Header_PeptidesHeader& pep_header = *(header_with_mods.mutable_peptides_header());
  pep_header.Clear();
  pep_header.set_min_mass(get_double_parameter("min-mass"));
  pep_header.set_max_mass(get_double_parameter("max-mass"));
  pep_header.set_min_length(get_int_parameter("min-length"));
  pep_header.set_max_length(get_int_parameter("max-length"));
  pep_header.set_monoisotopic_precursor(
    get_boolean_parameter("monoisotopic-precursor"));
  string enzyme = enzyme_type_to_string(get_enzyme_type_parameter("enzyme"));
  if (enzyme == "no-enzyme") {
    enzyme = "none";
  } else {
    DIGEST_T digestion = get_digest_type_parameter("digestion");
    if (digestion != FULL_DIGEST && digestion != PARTIAL_DIGEST) {
      carp(CARP_FATAL, "'digestion' must be 'full-digest' or 'partial-digest'");
    }
    pep_header.set_full_digestion(digestion == FULL_DIGEST);
    pep_header.set_max_missed_cleavages(get_int_parameter("missed-cleavages"));
  }
  pep_header.set_enzyme(enzyme);
  string mods_spec = get_string_parameter_pointer("mods-spec");
  if (mods_spec.find('C') == string::npos) {
    mods_spec = (mods_spec.empty()) ?
      default_cysteine : default_cysteine + ',' + mods_spec;
    carp(CARP_DEBUG, "Using default cysteine mod '%s' ('%s')",
         default_cysteine.c_str(), mods_spec.c_str());
  }

  VariableModTable var_mod_table;
  if (!var_mod_table.Parse(mods_spec.c_str())) {
    carp(CARP_FATAL, "Error parsing mods");
  }
  pep_header.mutable_mods()->CopyFrom(*(var_mod_table.ParsedModTable()));
  if (!MassConstants::Init(var_mod_table.ParsedModTable())) {
    carp(CARP_FATAL, "Error in MassConstants::Init");
  }

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

  bool need_mods = header_with_mods.peptides_header().mods().variable_mod_size() > 0;
  string basic_peptides = need_mods ? modless_peptides : peakless_peptides;
  carp(CARP_DETAILED_DEBUG, "basic_peptides=%s", basic_peptides.c_str());

  carp(CARP_INFO, "Computing unmodified peptides...");
  MakePeptides(&header_no_mods, basic_peptides, out_aux);

  vector<const pb::Protein*> proteins;
  if (!ReadRecordsToVector<pb::Protein>(&proteins, out_proteins)) {
    carp(CARP_FATAL, "Error reading proteins file");
  }

  if (need_mods) {
    carp(CARP_INFO, "Computing modified peptides...");
    HeadedRecordReader reader(modless_peptides, NULL, 1024 << 10); // 1024kb buffer
    AddMods(&reader, peakless_peptides, header_with_mods, proteins);
  }

  carp(CARP_INFO, "Precomputing theoretical spectra...");
  AddTheoreticalPeaks(proteins, peakless_peptides, out_peptides);

  // Recover stderr
  cerr.rdbuf(old);

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
