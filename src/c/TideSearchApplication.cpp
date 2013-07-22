#include <cstdio>
#include "tide/abspath.h"
#include "tide/records_to_vector-inl.h"

#include "carp.h"
#include "parameter.h"
#include "SpectrumRecordWriter.h"
#include "TideSearchApplication.h"

extern AA_MOD_T* list_of_mods[MAX_AA_MODS];
extern int num_mods;

TideSearchApplication::TideSearchApplication() {
}

TideSearchApplication::~TideSearchApplication() {
}

int TideSearchApplication::main(int argc, char** argv) {

  const char* option_list[] = {
    "precursor-window",
    "precursor-window-type",
    "top-match",
    "store-spectra",
    "compute-sp",
    "txt-output",
    "sqt-output",
    "pepxml-output",
    "fileroot",
    "output-dir",
    "overwrite",
    "parameter-file",
    "verbosity"
  };
  int num_options = sizeof(option_list) / sizeof(char*);
  const char* arg_list[] = {
    "tide spectra file",
    "tide database index"
  };
  int num_args = sizeof(arg_list) / sizeof(char*);
  initialize(arg_list, num_args, option_list, num_options, argc, argv);

  carp(CARP_INFO, "Running tide-search...");

  string cmd_line = "crux tide-search";
  for (int i = 1; i < argc; ++i) {
    cmd_line += " ";
    cmd_line += argv[i];
  }

  string index_dir = get_string_parameter_pointer("tide database index");
  string peptides_file = index_dir + "/pepix";
  string proteins_file = index_dir + "/protix";
  string spectra_file = get_string_parameter_pointer("tide spectra file");

  double window = get_double_parameter("precursor-window");
  WINDOW_TYPE_T window_type = get_window_type_parameter("precursor-window-type");
  int top_matches = get_int_parameter("top-match");

  bool compute_sp = get_boolean_parameter("compute-sp");
  if (get_boolean_parameter("sqt-output") && !compute_sp){
    compute_sp = true;
    carp(CARP_INFO, "Enabling parameter compute-sp since SQT output is enabled "
                    " (this will increase runtime).");
  }

  const double max_mz = 0.0;

  carp(CARP_INFO, "Reading index %s", index_dir.c_str());
  // Read proteins index file
  ProteinVec proteins;
  ActivePeptideQueue* active_peptide_queue = NULL;
  pb::Header protein_header;
  if (!ReadRecordsToVector<pb::Protein, const pb::Protein>(&proteins,
      proteins_file, &protein_header)) {
    carp(CARP_FATAL, "Error reading index (%s)", proteins_file.c_str());
  }
  carp(CARP_DEBUG, "Read %d proteins", proteins.size());

  // Read peptides index file
  pb::Header peptides_header;
  HeadedRecordReader peptide_reader(peptides_file, &peptides_header);
  if (!peptides_header.file_type() == pb::Header::PEPTIDES ||
      !peptides_header.has_peptides_header()) {
    carp(CARP_FATAL, "Error reading index (%s)", peptides_file.c_str());
  }
  MassConstants::Init(&peptides_header.peptides_header().mods());
  active_peptide_queue = new ActivePeptideQueue(peptide_reader.Reader(), proteins);

  carp(CARP_INFO, "Reading spectra file %s", spectra_file.c_str());
  // Try to read file as spectrumrecords file
  SpectrumCollection spectra;
  pb::Header spectrum_header;
  string delete_spectra_file = "";
  if (!spectra.ReadSpectrumRecords(spectra_file, &spectrum_header)) {
    // Failed, try converting to spectrumrecords file
    carp(CARP_INFO, "Converting %s to spectrumrecords format",
                    spectra_file.c_str());
    string converted_spectra_file = get_string_parameter_pointer("store-spectra");
    if (converted_spectra_file.empty()) {
      char tmpnam_buffer[L_tmpnam];
      tmpnam(tmpnam_buffer);
      delete_spectra_file = converted_spectra_file = tmpnam_buffer;
    }
    carp(CARP_DEBUG, "New spectrumrecords filename: %s",
                     converted_spectra_file.c_str());
    if (!SpectrumRecordWriter::convert(spectra_file, converted_spectra_file)) {
      carp(CARP_FATAL, "Error converting %s to spectrumrecords format",
                       spectra_file.c_str());
    }
    carp(CARP_DEBUG, "Reading converted spectra file %s",
                     spectra_file.c_str());
    // Re-read converted file as spectrumrecords file
    if (!spectra.ReadSpectrumRecords(converted_spectra_file, &spectrum_header)) {
      carp(CARP_DEBUG, "Deleting %s", converted_spectra_file.c_str());
      remove(converted_spectra_file.c_str());
      carp(CARP_FATAL, "Error reading spectra file %s",
                       converted_spectra_file.c_str());
    }
  }

  carp(CARP_INFO, "Sorting spectra");
  spectra.Sort();
  if (max_mz == 0) {
    double highest_mz = spectra.FindHighestMZ();
    carp(CARP_DEBUG, "Max m/z %f", highest_mz);
    MaxMZ::SetGlobalMax(highest_mz);
  } else {
    MaxMZ::SetGlobalMaxFromFlag();
  }

  // Do the search
  carp(CARP_INFO, "Running search");
  cleanMods();
  search(spectra.SpecCharges(), active_peptide_queue, proteins, window,
         window_type, top_matches, spectra.FindHighestMZ(), compute_sp);

  // Delete temporary spectrumrecords file
  if (!delete_spectra_file.empty()) {
    carp(CARP_DEBUG, "Deleting %s", delete_spectra_file.c_str());
    remove(delete_spectra_file.c_str());
  }

  // Clean up
  delete active_peptide_queue;
  for (ProteinVec::const_iterator i = proteins.begin(); i != proteins.end(); ++i) {
    delete const_cast<pb::Protein*>(*i);
  }

  return 0;
}

/**
 * Free all existing mods
 */
void TideSearchApplication::cleanMods() {
  for (int i = 0; i < MAX_AA_MODS; ++i) {
    free_aa_mod(list_of_mods[i]);
    list_of_mods[i] = NULL;
  }
  num_mods = 0;
}

void TideSearchApplication::search(
  const vector<SpectrumCollection::SpecCharge>* spec_charges,
  ActivePeptideQueue* active_peptide_queue,
  const ProteinVec& proteins,
  double precursor_window,
  WINDOW_TYPE_T window_type,
  int top_matches,
  double highest_mz,
  bool compute_sp
) {
  OutputFiles output_files(this);
  output_files.writeHeaders();
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
    // Calculate and set the window, depending on the window type.
    double min_mass, max_mass;
    switch (window_type) {
    case WINDOW_MASS:
      min_mass = pre_mass - precursor_window;
      max_mass = pre_mass + precursor_window;
      break;
    case WINDOW_MZ: {
      double mz_minus_proton = spectrum->PrecursorMZ() - MASS_PROTON;
      min_mass = (mz_minus_proton - precursor_window) * charge;
      max_mass = (mz_minus_proton + precursor_window) * charge;
      break;
      }
    case WINDOW_PPM: {
      double tiny_precursor = precursor_window * 1e-6;
      min_mass = pre_mass / (1.0 + tiny_precursor);
      max_mass = pre_mass / (1.0 - tiny_precursor);
      break;
      }
    default:
      carp(CARP_FATAL, "Invalid window type");
    }
    carp(CARP_DETAILED_DEBUG, "Scan %d (%f m/z, %f neutral mass, charge %d) "
         "mass window is [%f, %f]",
         spectrum->SpectrumNumber(), spectrum->PrecursorMZ(), pre_mass, charge,
         min_mass, max_mass);

    int size = active_peptide_queue->SetActiveRange(min_mass, max_mass);
    MatchSet::Arr match_arr(size); // Scored peptides will go here.

    // Programs for taking the dot-product with the observed spectrum are laid
    // out in memory managed by the active_peptide_queue, one program for each
    // candidate peptide. The programs will store the results directly into
    // match_arr. We now pass control to those programs.
    collectScoresCompiled(active_peptide_queue, spectrum, observed, &match_arr,
                          size, charge);

    // matches will arrange the results in a heap by score, return the top
    // few, and recover the association between counter and peptide. We output
    // the top matches.
    MatchSet matches(&match_arr, highest_mz);
    matches.report(&output_files, top_matches, spectrum, charge,
                   active_peptide_queue, proteins, compute_sp);
  }
  output_files.writeFooters();
}

void TideSearchApplication::collectScoresCompiled(
  ActivePeptideQueue* active_peptide_queue,
  const Spectrum* spectrum,
  const ObservedPeakSet& observed,
  MatchSet::Arr* match_arr,
  int queue_size,
  int charge
) {
  if (!active_peptide_queue->HasNext()) {
    return;
  }
  // prog gets the address of the dot-product program for the first peptide
  // in the active queue.
  const void* prog = active_peptide_queue->NextPeptide()->Prog(charge);
  const int* cache = observed.GetCache();
  // results will get (score, counter) pairs, where score is the dot product
  // of the observed peak set with a candidate peptide. The candidate
  // peptide is given by counter, which refers to the index within the
  // ActivePeptideQueue, counting from the back. This complication
  // simplifies the generated programs, which now simply dump the counter.
  pair<int, int>* results = match_arr->data();

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
                         "c" (queue_size),
                         "D" (results)
  );

  // match_arr is filled by the compiled programs, not by calls to
  // push_back(). We have to set the final size explicitly.
  match_arr->set_size(queue_size);
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
