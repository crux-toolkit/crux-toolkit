// Benjamin Diament

// Search() implements the main scoring loop to match each observed spectrum
// with the best few candidate theoretical peptides and to report the results.
  
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <gflags/gflags.h>
#include "raw_proteins.pb.h"
#include "peptides.pb.h"
#include "spectrum.pb.h"
#include "records.h"
#include "records_to_vector-inl.h"
#include "peptide.h"
#include "theoretical_peak_set.h"
#include "spectrum_preprocess.h"
#include "spectrum_collection.h"
#include "active_peptide_queue.h"
#include "max_mz.h"
#include "report.h"

using namespace std;
/*
typedef vector<const pb::Protein*> ProteinVec;

DEFINE_string(proteins, "", "File of proteins corresponding to peptides, as "
                            "raw_proteins.proto");
DEFINE_string(peptides, "", "File of unfragmented peptides, as "
                            "peptides.proto");
DEFINE_string(spectra, "", "Spectrum input file");

DEFINE_string(results_file, "results.tideres", "Results output file");

DEFINE_bool(stats, false, "Compute sample moments on match candidates");

// TODO: Once the results postprocessor is done, we no longer want to write
// out the text version, so this flag should probably be removed!
DEFINE_string(results, "text", "Results format. Can be text or protobuf.");

DEFINE_double(mass_window, 3.0, "Precursor mass tolerance in Daltons");
DEFINE_int32(top_matches, 5, "Number of matches to report for each "
                             "spectrum");

#ifndef NDEBUG
DEFINE_bool(debug_peaks, false, "Show all preprocessed peaks");
#endif

///////////////////////////////////////////////////////////////////////////////
// Set of peptide-spectrum matches. Supply a Reporter to Report() to see
// best few. 
class MatchSet {
 public:
  typedef pair<int, int> Pair;
  typedef FixedCapacityArray<Pair> Arr;

  // Matches will be an array of pairs, (score, counter), where counter refers
  // to the index within the ActivePeptideQueue, counting from the back.  This
  // slight complication is due to the way the generated machine code fills the
  // counter in the matches buffer by decrementing the counter.
  MatchSet(Arr* matches) : matches_(matches) {}

  void Report(Reporter* reporter, int top_n, const Spectrum* spectrum,
              int charge, int spectrum_index, 
              const ActivePeptideQueue* peptides, const ProteinVec& proteins) {
    if (top_n > matches_->size())
      top_n = matches_->size();
    GetTop(top_n);
    pb::Stats stats;
    stats.set_count(matches_->size());
    if (FLAGS_stats) {
      double sum = 0;
      double sum_squares = 0;
      for (Arr::iterator i = matches_->begin(); i != matches_->end(); ++i) {
	double score = i->first / 100000000.0;
	sum += score;
	sum_squares += score * score;
      }
      stats.set_sum(sum);
      stats.set_sum_squares(sum_squares);
    }
    reporter->ReportSpectrum(spectrum, charge, spectrum_index, &stats);
    for (Arr::iterator i = matches_->end() - 1; top_n; --i, --top_n) { 
      // i->second is counter as it was during scoring and corresponds to the
      // index of the peptide in the ActivePeptideQueue, counting from the
      // back. GetPeptide() retrieves the corresponding Peptide.
      const Peptide* peptide = peptides->GetPeptide(i->second);
      reporter->ReportMatch(i->first, *peptide);
    }
    reporter->WriteReport();
  }

 private:
  struct less_score : public binary_function<Pair, Pair, bool> {
    // Compare scores, ignore counters.
    bool operator()(Pair x, Pair y) { return x.first < y.first; }
  };

  void GetTop(int top_n) {
    assert(top_n <= matches_->size());
    // move top n elements to end of array, with largest element last
    make_heap(matches_->begin(), matches_->end(), less_score());
    for (int i = 0; i < top_n; ++i)
      pop_heap(matches_->begin(), matches_->end() - i, less_score());
  }

  Arr* matches_;
};

void CollectScoresCompiled(ActivePeptideQueue* active_peptide_queue,
                           const Spectrum* spectrum,
                           const ObservedPeakSet& observed,
                           MatchSet::Arr* match_arr,
                           int queue_size, int charge) {
  if (!active_peptide_queue->HasNext())
    return;
  // prog gets the address of the dot-product program for the first peptide in
  // the active queue.
  const void* prog = active_peptide_queue->NextPeptide()->Prog(charge);
  const int* cache = observed.GetCache();
  // results will get (score, counter) pairs, where score is the dot product of
  // the observed peak set with a candidate peptide. The candidate peptide is
  // given by counter, which refers to the index within the ActivePeptideQueue,
  // counting from the back. This complication simplifies the generated
  // programs, which now simply dump the counter.
  pair<int, int>* results = match_arr->data();

  // See compiler.h for a description of the programs beginning at prog and how
  // they are generated. Here we initialize certain registers to the values
  // expected by the programs and call the first one (*prog).
  //
  // See gnu assembler format for more on this format. We tell the compiler to
  // set these registers:
  // edx/rdx points to the cache.
  // eax/rax points to the first program.
  // ecx/rcx is the counter and gets the size of the active queue.
  // edi/rdi points to the results buffer.
  //
  // The push and pop operations are a workaround for a compiler that
  // doesn't understand that %ecx and %edi (or %rcx and %rdi) get
  // clobbered. Since they're already input registers, they can't be
  // included in the clobber list.
#ifdef _MSC_VER
#else
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
#endif

  // match_arr is filled by the compiled programs, not by calls to
  // push_back(). We have to set the final size explicitly.
  match_arr->set_size(queue_size);
}


void Search(const vector<SpectrumCollection::SpecCharge>* spec_charges,
            const ProteinVec& proteins,
            ActivePeptideQueue* active_peptide_queue, Reporter* reporter) {
  // This is the main search loop.
  ObservedPeakSet observed;
  // cycle through spectrum-charge pairs, sorted by neutral mass
  vector<SpectrumCollection::SpecCharge>::const_iterator sc;
  for (sc = spec_charges->begin(); sc != spec_charges->end(); ++sc) {
    const Spectrum* spectrum = sc->spectrum;
    int charge = sc->charge;
    double pre_mass = sc->neutral_mass;
    int spectrum_index = sc->spectrum_index;

    // Normalize the observed spectrum and compute the cache of
    // frequently-needed values for taking dot products with theoretical
    // spectra.
    observed.PreprocessSpectrum(*spectrum, charge);

#ifndef NDEBUG
    if (FLAGS_debug_peaks)
      observed.ShowPeaks();
#endif

    // The active peptide queue holds the candidate peptides for spectrum.
    // These reside in a window of 2 * FLAGS_mass_window Daltons around the
    // neutral precursor mass of the spectrum.
    int size = active_peptide_queue->SetActiveRange(pre_mass - FLAGS_mass_window, 
                                                    pre_mass + FLAGS_mass_window);
    MatchSet::Arr match_arr(size); // Scored peptides will go here.

    // Programs for taking the dot-product with the observed spectrum are laid
    // out in memory managed by the active_peptide_queue, one program for each
    // candidate peptide. The programs will store the results directly into
    // match_arr. We now pass control to those programs.
    CollectScoresCompiled(active_peptide_queue, spectrum, observed, &match_arr,
                          size, charge);

    // matches will arrange the results in a heap by score, return the top
    // few, and recover the association between counter and peptide. We output
    // the top matches.
    MatchSet matches(&match_arr);
    matches.Report(reporter, FLAGS_top_matches, spectrum, charge, 
                   spectrum_index, active_peptide_queue, proteins);
  }
}

void CreateResultsHeader(pb::Header* proteins_header,
                         pb::Header* peptides_header,
                         pb::Header* results_header,
                         string& command_line) {
  // Set the results file_type
  results_header->set_file_type(pb::Header::RESULTS);
  results_header->set_command_line(command_line);
  
  // Add the proteins source used to create these results
  pb::Header_Source* source = results_header->add_source();
  source->set_filename(AbsPath(FLAGS_proteins));
  source->mutable_header()->CopyFrom(*proteins_header);

  // Add the peptides source used to create these results
  source = results_header->add_source();
  source->set_filename(AbsPath(FLAGS_peptides));
  source->mutable_header()->CopyFrom(*peptides_header);

  // Add the spectrumrecords source used to create these results
  source = results_header->add_source();
  source->set_filename(AbsPath(FLAGS_spectra));
  source->set_filetype("spectrumrecords");
  
  // Add the results specific header info
  pb::Header_ResultsHeader& results_specific_header = *results_header->mutable_results_header();
  results_specific_header.set_mass_window(FLAGS_mass_window);
  results_specific_header.set_top_matches(FLAGS_top_matches); 
  results_specific_header.mutable_peptides_header()->CopyFrom(peptides_header->peptides_header());
}

string UsageMessage(string prog_name) {
  return prog_name + " [ settings ]\n" 
    "settings must include at least the following:\n"
    "  --proteins=<filename>\n"
    "  --peptides=<filename>\n"
    "  --spectra=<filename>\n"
    "\n";
}
*/
int main(int argc, char* argv[]) {
/*  string command_line;
  for (int i = 0; i < argc; i++) {
    command_line += argv[i];
    command_line += " ";
  }

  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::SetUsageMessage(UsageMessage(argv[0]));
  google::ParseCommandLineFlags(&argc, &argv, true);
  
  CHECK(!FLAGS_proteins.empty());
  CHECK(!FLAGS_peptides.empty());
  CHECK(!FLAGS_spectra.empty());

  ProteinVec proteins;
  ActivePeptideQueue* active_peptide_queue = NULL;
  pb::Header protein_header;
  bool ok = ReadRecordsToVector<pb::Protein, const pb::Protein>(&proteins, 
                                                                FLAGS_proteins, 
                                                                &protein_header);
  CHECK(ok);

  pb::Header peptides_header;
  HeadedRecordReader peptide_reader(FLAGS_peptides, &peptides_header);
  CHECK(peptides_header.file_type() == pb::Header::PEPTIDES);
  CHECK(peptides_header.has_peptides_header());
  MassConstants::Init(&peptides_header.peptides_header().mods());
  active_peptide_queue = new ActivePeptideQueue(peptide_reader.Reader(),
						proteins);

  SpectrumCollection spectra;
  pb::Header spectrum_header;
  CHECK(spectra.ReadSpectrumRecords(FLAGS_spectra, &spectrum_header));
  
  spectra.Sort();
  if (FLAGS_max_mz == 0) {
    MaxMZ::SetGlobalMax(spectra.FindHighestMZ());
  } else {
    MaxMZ::SetGlobalMaxFromFlag();
  }

  // TODO: Once the postprocessor is written, we should probably just remove
  // the text reporter and only write out PB results.
  if (FLAGS_results == "protobuf") {
    pb::Header results_header;
    CreateResultsHeader(&protein_header, &peptides_header, &results_header,
                        command_line);
    PBReporter pb_reporter(FLAGS_results_file, results_header);
    Search(spectra.SpecCharges(), proteins, active_peptide_queue, &pb_reporter);
  } else {
    TextReporter text_reporter;
    Search(spectra.SpecCharges(), proteins, active_peptide_queue, &text_reporter);
  }
  
  delete active_peptide_queue;
      
  ProteinVec::const_iterator i = proteins.begin();
  for (; i != proteins.end(); ++i)
    delete const_cast<pb::Protein*>(*i);

#ifdef MEM_STATS
  cerr << "Total Fifo Allocations: " << FifoAllocator::Total() << endl;
#endif
*/
  return 0;
}
