#include "TideSearchApplication.h"

#include "tide/records_to_vector-inl.h"
#include "tide/spectrum_preprocess.h"

#include "carp.h"
#include "parameter.h"

TideSearchApplication::TideSearchApplication() {
}

TideSearchApplication::~TideSearchApplication() {
}

int TideSearchApplication::main(int argc, char** argv) {

  const char* option_list[] = {
    "mass-window",
    "top-matches",
    "results"
  };
  int num_options = sizeof(option_list) / sizeof(char*);
  const char* arg_list[] = {
    "tide peptides file",
    "tide proteins file",
    "spectrum records file"
  };
  int num_args = sizeof(arg_list) / sizeof(char*);
  initialize(arg_list, num_args, option_list, num_options, argc, argv);

  carp(CARP_INFO, "Running tide-search...");

  string cmd_line = "crux tide-search";
  for (int i = 1; i < argc; ++i) {
    cmd_line += " ";
    cmd_line += argv[i];
  }

  string peptides_file = get_string_parameter_pointer("tide peptides file");
  string proteins_file = get_string_parameter_pointer("tide proteins file");
  string spectra_file = get_string_parameter_pointer("spectrum records file");

  double mass_window = get_double_parameter("mass-window");
  int top_matches = get_int_parameter("top-matches");
  string out_format = get_string_parameter_pointer("results");
  if (out_format != "protobuf" && out_format != "text") {
    carp(CARP_FATAL, "'results' must be 'protobuf' or 'text'");
  }

  // TODO add these as options
  string results_file = "results.tideres"; // only if writing results to protobuf
  double max_mz = 0.0;

  ProteinVec proteins;
  ActivePeptideQueue* active_peptide_queue = NULL;
  pb::Header protein_header;
  if (!ReadRecordsToVector<pb::Protein, const pb::Protein>(&proteins,
      proteins_file, &protein_header)) {
    carp(CARP_FATAL, "Error reading proteins file");
  }

  pb::Header peptides_header;
  HeadedRecordReader peptide_reader(peptides_file, &peptides_header);
  if (!peptides_header.file_type() == pb::Header::PEPTIDES) {
    carp(CARP_FATAL, "Peptides file is of wrong type");
  } else if (!peptides_header.has_peptides_header()) {
    carp(CARP_FATAL, "Peptides file does not have peptides header");
  }
  MassConstants::Init(&peptides_header.peptides_header().mods());
  active_peptide_queue = new ActivePeptideQueue(peptide_reader.Reader(), proteins);

  SpectrumCollection spectra;
  pb::Header spectrum_header;
  if (!spectra.ReadSpectrumRecords(spectra_file, &spectrum_header)) {
    carp(CARP_FATAL, "Error reading spectra file");
  }

  spectra.Sort();
  if (max_mz == 0) {
    MaxMZ::SetGlobalMax(spectra.FindHighestMZ());
  } else {
    MaxMZ::SetGlobalMaxFromFlag();
  }

  if (out_format == "protobuf") {
    pb::Header results_header;
    createResultsHeader(&protein_header, &peptides_header, &results_header,
                        proteins_file, peptides_file, spectra_file, mass_window,
                        top_matches, cmd_line);
    PBReporter pb_reporter(results_file, results_header);
    search(spectra.SpecCharges(), proteins, active_peptide_queue, &pb_reporter,
           mass_window, top_matches);
  } else {
    TextReporter text_reporter;
    search(spectra.SpecCharges(), proteins, active_peptide_queue, &text_reporter,
           mass_window, top_matches);
  }

  delete active_peptide_queue;

  for (ProteinVec::const_iterator i = proteins.begin(); i != proteins.end(); ++i)
    delete const_cast<pb::Protein*>(*i);

  return 0;
}

void TideSearchApplication::createResultsHeader(
  pb::Header* proteins_header,
  pb::Header* peptides_header,
  pb::Header* results_header,
  const string& proteins_file,
  const string& peptides_file,
  const string& spectra_file,
  double mass_window,
  int top_matches,
  const string& cmd_line
) {
  // Set the results file_type
  results_header->set_file_type(pb::Header::RESULTS);
  results_header->set_command_line(cmd_line);
  
  // Add the proteins source used to create these results
  pb::Header_Source* source = results_header->add_source();
  source->set_filename(AbsPath(proteins_file));
  source->mutable_header()->CopyFrom(*proteins_header);

  // Add the peptides source used to create these results
  source = results_header->add_source();
  source->set_filename(AbsPath(peptides_file));
  source->mutable_header()->CopyFrom(*peptides_header);

  // Add the spectrumrecords source used to create these results
  source = results_header->add_source();
  source->set_filename(AbsPath(spectra_file));
  source->set_filetype("spectrumrecords");
  
  // Add the results specific header info
  pb::Header_ResultsHeader& results_specific_header = *results_header->mutable_results_header();
  results_specific_header.set_mass_window(mass_window);
  results_specific_header.set_top_matches(top_matches); 
  results_specific_header.mutable_peptides_header()->CopyFrom(peptides_header->peptides_header());
}

void TideSearchApplication::search(
  const vector<SpectrumCollection::SpecCharge>* spec_charges,
  const ProteinVec& proteins,
  ActivePeptideQueue* active_peptide_queue,
  Reporter* reporter,
  double mass_window,
  int top_matches
) {
  // This is the main search loop.
  ObservedPeakSet observed;
  // cycle through spectrum-charge pairs, sorted by neutral mass
  for (vector<SpectrumCollection::SpecCharge>::const_iterator sc = spec_charges->begin();
       sc != spec_charges->end();
       ++sc) {
    const Spectrum* spectrum = sc->spectrum;
    int charge = sc->charge;
    double pre_mass = sc->neutral_mass;
    int spectrum_index = sc->spectrum_index;

    // Normalize the observed spectrum and compute the cache of
    // frequently-needed values for taking dot products with theoretical
    // spectra.
    observed.PreprocessSpectrum(*spectrum, charge);

    // The active peptide queue holds the candidate peptides for spectrum.
    // These reside in a window of 2 * FLAGS_mass_window Daltons around the
    // neutral precursor mass of the spectrum.
    int size = active_peptide_queue->SetActiveRange(pre_mass - mass_window, 
                                                    pre_mass + mass_window);
    MatchSet::Arr match_arr(size); // Scored peptides will go here.

    // Programs for taking the dot-product with the observed spectrum are laid
    // out in memory managed by the active_peptide_queue, one program for each
    // candidate peptide. The programs will store the results directly into
    // match_arr. We now pass control to those programs.
    //collectScoresCompiled(active_peptide_queue, spectrum, observed, &match_arr,
    //                      size, charge);
    MatchSet::Arr* match_arr_ptr = &match_arr;
    if (active_peptide_queue->HasNext()) {
      // prog gets the address of the dot-product program for the first peptide
      // in the active queue.
      const void* prog = active_peptide_queue->NextPeptide()->Prog(charge);
      const int* cache = observed.GetCache();
      // results will get (score, counter) pairs, where score is the dot product
      // of the observed peak set with a candidate peptide. The candidate
      // peptide is given by counter, which refers to the index within the
      // ActivePeptideQueue, counting from the back. This complication
      // simplifies the generated programs, which now simply dump the counter.
      pair<int, int>* results = (match_arr_ptr)->data();

      // See compiler.h for a description of the programs beginning at prog and
      // how they are generated. Here we initialize certain registers to the
      // values expected by the programs and call the first one (*prog).
      //
      // See gnu assembler format for more on this format. We tell the compiler
      // to set these registers:
      // edx/rdx points to the cache.
      // eax/rax points to the first program.
      // ecx/rcx is the counter and gets the size of the active queue.
      // edi/rdi points to the results buffer.
      //
      // The push and pop operations are a workaround for a compiler that
      // doesn't understand that %ecx and %edi (or %rcx and %rdi) get
      // clobbered. Since they're already input registers, they can't be
      // included in the clobber list.
      __asm__ __volatile__("cld\n" // stos operations increment edi
#ifdef __x86_64__
                           "push %%rcx\n"
                           "push %%rdi\n"
                           "call *%%rax\n"
                           "pop %%rdi\n"
                           "pop %%rcx\n"
#else
                           "push %%ecx\n"
                           "push %%edi\n"
                           "call *%%eax\n"
                           "pop %%edi\n"
                           "pop %%ecx\n"
#endif
                           : // no outputs
                           : "d" (cache),
                             "a" (prog),
                             "c" (size),
                             "D" (results)
      );

      // match_arr is filled by the compiled programs, not by calls to
      // push_back(). We have to set the final size explicitly.
      match_arr_ptr->set_size(size);
    }

    // matches will arrange the results in a heap by score, return the top
    // few, and recover the association between counter and peptide. We output
    // the top matches.
    MatchSet matches(&match_arr);
    matches.Report(reporter, top_matches, spectrum, charge, 
                   spectrum_index, active_peptide_queue, proteins);
  }
}

string TideSearchApplication::getName() {
  return "tide-search";
}

string TideSearchApplication::getDescription() {
  return "Runs tide-search";
}

bool TideSearchApplication::needsOutputDirectory() {
  return true;
}

COMMAND_T TideSearchApplication::getCommand() {
  return TIDE_SEARCH_COMMAND;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
