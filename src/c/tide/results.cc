// Tahmina Baker
//
// Results() is a post-processer for the results protocol buffer file 
// generated during Search. The goal is to support multiple formats 
// (e.g. text, pepxml and sqt) and provide the user with flexible options
// for viewing the search results. 
  
#include <fstream>
#include "records_to_vector-inl.h"
#include "abspath.h"
#include "results_writer.h"

using namespace std;


DEFINE_string(results_file, "", "Results file generated via Search, as "
                                "results.proto");

DEFINE_string(proteins, "", "File of proteins corresponding to peptides, "
                            "as raw_proteins.proto");

DEFINE_string(aux_locations, "", "File of auxiliary locations corresponding "
                                 "to peptides in the results");

DEFINE_string(out_format, "text", "The output format to be generated. Can be "
                                  "text, pep.xml or sqt. Default is text.");

DEFINE_string(out_filename, "", "Name of the output file to generate. An "
                                "extension will be added based on the "
                                " out_format.");

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

DEFINE_string(spectra, "", "Spectrum input file");


string UsageMessage(string prog_name) {
  return prog_name + " [ settings ]\n" 
    "settings must include at least the following:\n"
    "  --results_file=<filename>\n"
    "  --proteins=<filename>\n"    
    "\n";
}

ostream* GetOutputStream() {
  // If the user specified an output filename, the output stream should point 
  // to a file stream; otherwise, direct the output to standard out
  if (!FLAGS_out_filename.empty()) {
    string out_file = FLAGS_out_filename + "." + FLAGS_out_format;
    return new ofstream(out_file.c_str());
  } else {
    return &cout;
  }
}

int main(int argc, char* argv[]) {
  string command_line;
  for (int i = 0; i < argc; i++) {
    command_line += argv[i];
    command_line += " ";
  }

  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::SetUsageMessage(UsageMessage(argv[0]));
  google::ParseCommandLineFlags(&argc, &argv, true);

  CHECK(!FLAGS_results_file.empty());
  CHECK(!FLAGS_proteins.empty());

  // Read the proteins into a vector
  ProteinVec proteins;
  CHECK((ReadRecordsToVector<pb::Protein, const pb::Protein>(&proteins, 
    FLAGS_proteins)));
  
  // Read the auxiliary locations into a vector only if the user has
  // specified the "show_all_proteins" flag.
  AuxLocVec aux_locs;
  if (FLAGS_show_all_proteins) {
    CHECK((ReadRecordsToVector<pb::AuxLocation, const pb::AuxLocation>(&aux_locs,
      FLAGS_aux_locations)));
  }

  // Create a reader for the results protocol buffer file and check the header
  pb::Header results_header;
  HeadedRecordReader results_reader(AbsPath(FLAGS_results_file), 
                                    &results_header);
  CHECK(results_header.file_type() == pb::Header::RESULTS);
  CHECK(results_header.has_results_header());
  CHECK(results_header.results_header().has_peptides_header());

  const pb::ModTable* mod_table = 
    &results_header.results_header().peptides_header().mods();
  MassConstants::Init(mod_table);
 
  // The spectra referenced in the results protocol buffer don't have full
  // peak information. We can get the peak information by using a spectrum
  // index stored in the results protocol buffer to looks into the
  // spectrumrecords protocol buffer.
  SpectrumCollection spectra;
  pb::Header spectrum_header;
  CHECK(!FLAGS_spectra.empty());
  CHECK(spectra.ReadSpectrumRecords(FLAGS_spectra, &spectrum_header));
  MaxMZ::SetGlobalMax(spectra.FindHighestMZ());

  // Conditionally create the writer based on FLAGS_out_format. If the user 
  // doesn't give one, or an invalid format is specified, we'll default to 
  // the text writer.
  ResultsWriter *res_writer = NULL;
  ostream* out_stream = GetOutputStream();
  if (FLAGS_out_format == "sqt") {
    res_writer = new SqtResultsWriter(proteins, aux_locs, spectra, 
                                      mod_table, *out_stream, 
                                      results_header, command_line);
  } else if (FLAGS_out_format == "pep.xml") {
    res_writer = new PepXMLResultsWriter(proteins, aux_locs, *out_stream);
  } else {
    res_writer = new TextResultsWriter(proteins, aux_locs, *out_stream);
  }

  // Go through each line in the results protocol buffer and write it out
  // via the results writer
  CHECK(results_reader.OK());
  pb::Results current_pb_result;
  while (!results_reader.Done()) {
    results_reader.Read(&current_pb_result);
    res_writer->WriteResults(current_pb_result);
  }

  CHECK(results_reader.OK());

  ProteinVec::const_iterator i = proteins.begin();
  for (; i != proteins.end(); ++i)
    delete const_cast<pb::Protein*>(*i);

  // res_writer keeps a reference to the out_stream, so delete it first
  // before we delete the out_stream
  delete res_writer;

  if (out_stream != &cout) {
    out_stream->flush();
    delete out_stream;
  }

  return 0;
}
