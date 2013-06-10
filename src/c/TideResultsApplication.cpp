#include "TideResultsApplication.h"

#include "tide/records_to_vector-inl.h"

#include "carp.h"
#include "parameter.h"

TideResultsApplication::TideResultsApplication() {
}

TideResultsApplication::~TideResultsApplication() {
}

DEFINE_string(spectrum_fields, "spectrum_num,mz,charge",
                               "A comma delimited set of fields to show for "
                               "an experimental spectrum in the order listed. "
                               "Available options are: spectrum_num,mz,charge");
DEFINE_string(match_fields, "xcorr,sequence",
                               "A comma delimited set of fields to show for "
                               "a matching peptide in the order listed. "
                               "Available options are: xcorr,sequence");
DEFINE_string(protein_fields, "protein_name,pos,aa_before,aa_after",
                               "A comma delimited set of fields to show for "
                               "a protein associated with a matching peptide. "
                               "Available options are: protein_name,pos, "
                               "aa_before,aa_after");
DEFINE_bool(show_mods, false, "Display modifications in the peptide "
                              "sequence.");
DEFINE_bool(show_all_proteins, false, "Display all the proteins the peptide "
                                      "was found in.");
int TideResultsApplication::main(int argc, char** argv) {

  const char* option_list[] = {
    "aux-locations",
    "out-filename",
    "show-all-proteins",
    "show-mods",
    "out-format",
    "spectrum-fields",
    "match-fields",
    "protein-fields"
  };
  int num_options = sizeof(option_list) / sizeof(char*);
  const char* arg_list[] = {
    "tide results file",
    "tide proteins file",
    "spectrum records file"
  };
  int num_args = sizeof(arg_list) / sizeof(char*);
  initialize(arg_list, num_args, option_list, num_options, argc, argv);

  carp(CARP_INFO, "Running tide-results...");

  string cmd_line = "crux tide-results";
  for (int i = 1; i < argc; ++i) {
    cmd_line += " ";
    cmd_line += argv[i];
  }

  string results_file = get_string_parameter_pointer("tide results file");
  string proteins_file = get_string_parameter_pointer("tide proteins file");
  string spectra_file = get_string_parameter_pointer("spectrum records file");

  string aux_locations = get_string_parameter_pointer("aux-locations");
  string out_file = get_string_parameter_pointer("out-filename");
  bool show_all_proteins = get_boolean_parameter("show-all-proteins");
  FLAGS_show_all_proteins = show_all_proteins;
  if (show_all_proteins && !file_exists(aux_locations)) {
    carp(CARP_FATAL, "show-all-proteins was true, but aux-locations did not "
                     "point to an existing file");
  }

  FLAGS_show_mods = get_boolean_parameter("show-mods");
  string out_format = get_string_parameter_pointer("out-format");
  if (out_format != "text" && out_format != "pep.xml" && out_format != "sqt") {
    carp(CARP_FATAL, "'out-format' must be 'text', 'pep.xml', or 'sqt'");
  }
  FLAGS_spectrum_fields = get_string_parameter_pointer("spectrum-fields");
  FLAGS_match_fields = get_string_parameter_pointer("match-fields");
  FLAGS_protein_fields = get_string_parameter_pointer("protein-fields");

  // Read the proteins into a vector
  ProteinVec proteins;
  if (!ReadRecordsToVector<pb::Protein, const pb::Protein>(&proteins, proteins_file)) {
    carp(CARP_FATAL, "Error reading proteins file");
  }

  // Read the auxiliary locations into a vector only if the user has specified
  // the "show_all_proteins" flag
  AuxLocVec aux_locs;
  if (show_all_proteins) {
    if (!ReadRecordsToVector<pb::AuxLocation, const pb::AuxLocation>(&aux_locs, aux_locations)) {
      carp(CARP_FATAL, "Error reading auxiliary locations file");
    }
  }

  // Create a reader for the results protocol buffer file
  pb::Header results_header;
  HeadedRecordReader results_reader(AbsPath(results_file), &results_header);
  if (!results_header.file_type() == pb::Header::RESULTS) {
    carp(CARP_FATAL, "Results file is of wrong type");
  } else if (!results_header.has_results_header()) {
    carp(CARP_FATAL, "Results file does not have results header");
  } else if (!results_header.results_header().has_peptides_header()) {
    carp(CARP_FATAL, "Results file does not have peptides header");
  }

  const pb::ModTable* mod_table = &results_header.results_header().peptides_header().mods();
  MassConstants::Init(mod_table);

  // The spectra referenced in the results protocol buffer don't have full peak
  // information. We can get the peak information by using a spectrum index
  // stored in the results protocol buffer to look into the spectrumrecords
  // protocol buffer.
  SpectrumCollection spectra;
  pb::Header spectrum_header;
  if (!spectra.ReadSpectrumRecords(spectra_file, &spectrum_header)) {
    carp(CARP_FATAL, "Error reading spectra file");
  }
  MaxMZ::SetGlobalMax(spectra.FindHighestMZ());

  // Conditionally create the writer based on FLAGS_out_format. If the user
  // doesn't give one, or an invalid format is specified, we'll default to the
  // text writer.
  ResultsWriter* res_writer = NULL;;
  ostream* out_stream = (!out_file.empty()) ?
    new ofstream(string(out_file + "." + out_format).c_str()) : &cout;

  if (out_format == "sqt") {
    res_writer = new SqtResultsWriter(proteins, aux_locs, spectra, mod_table,
                                      *out_stream, results_header, cmd_line);
  } else if (out_format == "pep.xml") {
    res_writer = new PepXMLResultsWriter(proteins, aux_locs, *out_stream);
  } else {
    res_writer = new TextResultsWriter(proteins, aux_locs, *out_stream);
  }

  // Go through each line in the results protocol buffer and write it out via
  // the results writer
  if (!results_reader.OK()) {
    carp(CARP_FATAL, "Error reading results file");
  }
  pb::Results current_pb_result;
  while (!results_reader.Done()) {
    results_reader.Read(&current_pb_result);
    res_writer->WriteResults(current_pb_result);
  }

  if (!results_reader.OK()) {
    carp(CARP_FATAL, "Error reading results file");
  }

  for (ProteinVec::const_iterator i = proteins.begin(); i != proteins.end(); ++i)
    delete const_cast<pb::Protein*>(*i);

  // res_writer keeps a reference to the out_stream, so delete it first before
  // we delete the out_stream
  delete res_writer;

  if (out_stream != &cout) {
    out_stream->flush();
    delete out_stream;
  }

  return 0;
}

string TideResultsApplication::getName() {
  return "tide-results";
}

string TideResultsApplication::getDescription() {
  return "Runs tide-results";
}

bool TideResultsApplication::needsOutputDirectory() {
  return true;
}

COMMAND_T TideResultsApplication::getCommand() {
  return TIDE_RESULTS_COMMAND;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
