#include <cstddef>
#include <cstdio>
#include "app/tide/abspath.h"
#include "app/tide/records_to_vector-inl.h"

#include "io/carp.h"
#include "parameter.h"
#include "io/SpectrumRecordWriter.h"
#include "TideIndexApplication.h"
#include "SpectrumConvertApplication.h"
#include "ParamMedicApplication.h"
#include "PSMConvertApplication.h"
#include "tide/mass_constants.h"
#include "TideMatchSet.h"
#include "tide/spectrum_collection.h"
#include "util/Params.h"
#include "util/FileUtils.h"
#include "util/StringUtils.h"
#include <math.h> 
#include <map>
#include "tide/ActivePeptideQueue.h"
#include "residue_stats.pb.h"
#include "crux_version.h"

//here :O
//#include <iostream>
//#include <set>
/*
#include <vector>
#include <queue>
#include <sstream>
#include <cstdlib>
*/
/* This constant is the product of the original "magic number" (10000,
 * on line 4622 of search28.c) that was used to rescale the XCorr
 * score, and the integerization constant used by Benjamin Diament in
 * Tide.  In the Tide publication, that constant was reported as 10^7.
 * --WSN, 10 March 2015 */
const double SpectrumConvertApplication::XCORR_SCALING = 100000000.0;

/* This constant is used to put the refactored XCorr back into the
 * same range as the original XCorr score.  It is the XCorr "magic
 * number" (10000) divided by the EVIDENCE_SCALE_INT (defined in
 * tide/spectrum_preprocess2.cc). */
const double SpectrumConvertApplication::RESCALE_FACTOR = 20.0;

// Things done:
// -- handle terminal mod structure in the peptides internally. GH issue: #639
// -- missing charge states, 0 charege states,, override charge states handled. --> need to update doc. GH issue: #557, #607

#define CHECK(x) GOOGLE_CHECK(x)

SpectrumConvertApplication::SpectrumConvertApplication() {
  remove_index_ = "";
  spectrum_flag_ = NULL;
  decoy_num_ = 0;
  num_range_skipped_ = 0;
  num_precursors_skipped_ = 0;
  num_isotopes_skipped_ = 0;
  num_retained_ = 0;
  num_spectra_ = 0;
  num_spectra_searched_ = 0;
  total_candidate_peptides_ = 0;
  precursor_window_ = 0;
  spectrum_min_mz_ = 0; 
  spectrum_max_mz_ = 0;
  min_scan_ = 0;
  max_scan_ = 0;
  min_peaks_ = 0;
  min_precursor_charge_ = 0;
  max_precursor_charge_ = 0;
  out_tsv_target_  = NULL; // original tide-search output format in tab-delimited text files (txt)
  out_tsv_decoy_ = NULL;  // original tide-search output format in tab-delimited text files (txt) for the decoy psms only
  out_mztab_target_ = NULL;      // mzTAB output format
  out_mztab_decoy_ = NULL;      // mzTAB output format for the decoy psms only
  out_pin_target_ = NULL;        // pin output format for percolator
  out_pin_decoy_ = NULL;        // pin output format for percolator for the decoy psms only
  total_spectra_num_ = 0;       // The total number of spectra searched. This is counted during the spectrum conversion

  for (int i = 0; i < NUMBER_LOCK_TYPES; i++) {  // LOCK_TYPES are defined in model/objects.h
    locks_array_.push_back(new boost::mutex());
  }
}

SpectrumConvertApplication::~SpectrumConvertApplication() {
  for (int i = 0; i < NUMBER_LOCK_TYPES; i++) {
    delete locks_array_[i];
  }
}

int SpectrumConvertApplication::main(int argc, char** argv) {
  return main(Params::GetStrings("tide spectra file"));
}

int SpectrumConvertApplication::main(const vector<string>& input_files) {
  return main(input_files, Params::GetString("tide database"));
}

int SpectrumConvertApplication::main(const vector<string>& input_files, const string input_index) {

  carp(CARP_INFO, "Running spectrum-convert...");

  
  bin_width_  = Params::GetDouble("mz-bin-width");
  bin_offset_ = Params::GetDouble("mz-bin-offset");

  use_neutral_loss_peaks_ = Params::GetBool("use-neutral-loss-peaks");
  use_flanking_peaks_ = Params::GetBool("use-flanking-peaks");

  negative_isotope_errors_ = getNegativeIsotopeErrors();  

  window_type_= string_to_window_type(Params::GetString("precursor-window-type"));
  precursor_window_ = Params::GetDouble("precursor-window");

  spectrum_min_mz_ =  Params::GetDouble("spectrum-min-mz") ;
  spectrum_max_mz_ = Params::GetDouble("spectrum-max-mz") ;
  min_peaks_ = Params::GetInt("min-peaks");
  min_precursor_charge_ = Params::GetInt("min-precursor-charge");
  max_precursor_charge_ = Params::GetInt("max-precursor-charge");  

  fragTol_ = Params::GetDouble("fragment-tolerance");
  granularityScale_ = Params::GetInt("evidence-granularity");  

  top_matches_ = Params::GetInt("top-match");  
  num_threads_ = Params::GetInt("num-threads");
  if (num_threads_ < 1) {
    num_threads_ = boost::thread::hardware_concurrency(); // MINIMUM # = 1.
    // (Meaning just main thread) Do not make this value below 1.
  } else if (num_threads_ > 64) {
    // make sure that number of threads are reasonable, e.g. user did not specify millions of threads...
    carp(CARP_FATAL, "Requested more than 64 threads.");
  }
  carp(CARP_INFO, "Number of Threads: %d", num_threads_);


  // Check scan-number parameter
  string scan_range = Params::GetString("scan-number");
  if (scan_range.empty()) {
    min_scan_ = 0;
    max_scan_ = BILLION;
    carp(CARP_DEBUG, "Searching all scans");
  } else if (scan_range.find('-') == string::npos) {
    // Single scan
    min_scan_ = max_scan_ = atoi(scan_range.c_str());
    carp(CARP_INFO, "Searching single scan %d", min_scan_);
  } else {
    if (!get_range_from_string(scan_range.c_str(), min_scan_, max_scan_)) {
      carp(CARP_FATAL, "The scan number range '%s' is invalid. "
           "Must be of the form <first>-<last>.", scan_range.c_str());
    }
    if (min_scan_ > max_scan_) {
      carp(CARP_FATAL, "Invalid scan range: %d to %d.", min_scan_, max_scan_);
    }
    carp(CARP_INFO, "Searching scan range %d to %d.", min_scan_, max_scan_);
  }

  // Determine which score function to use for scoring PSMs and store in SCORE_FUNCTION enum. Change this in ./src/util/crux-utils.cpp: 68
  curScoreFunction_ = string_to_score_function_type(Params::GetString("score-function"));  

  if (curScoreFunction_ == PVALUES && bin_width_ < 1.0) {
    carp(CARP_FATAL, "Tide-convert does not support P-value calculation with bin-width less than 1.0 Da.");
  }
  if (curScoreFunction_ >= NUMBER_SCORE_FUNCTIONS) {
    carp(CARP_FATAL, "Invalid score function.");
  }
 
  // Get a peptide reader to the peptide index datasets along with proteins, auxlocs. 
  vector<const pb::Protein*> proteins;
  vector<const pb::AuxLocation*> locations;
  pb::Header peptides_header;
  string peptides_file = FileUtils::Join(input_index, "pepix");  
  HeadedRecordReader peptide_reader = HeadedRecordReader(peptides_file, &peptides_header);
  getPeptideIndexData(input_index, proteins, locations, peptides_header);
  tide_index_mzTab_file_path_ = FileUtils::Join(input_index, TideIndexApplication::tide_index_mzTab_filename_);

  TideMatchSet::curScoreFunction_ = curScoreFunction_;
  TideMatchSet::top_matches_ = top_matches_;
  TideMatchSet::decoy_num_ = decoy_num_;
  TideMatchSet::mass_precision_ =  Params::GetInt("mass-precision");
  TideMatchSet::score_precision_ = Params::GetInt("precision");
  TideMatchSet::mod_precision_ = Params::GetInt("mod-precision");
  TideMatchSet::concat_ = Params::GetBool("concat");  

  // Create the output files, print headers
  createOutputFiles(); 

  // Convert the original file names into spectrum records if needed 
  // Update the file names in the variable inputFiles_ locally.
  // Run spectrum file convertion in parallel.
  for (vector<string>::const_iterator original_file_name = input_files.begin(); original_file_name != input_files.end(); ++original_file_name) {
    inputFiles_.push_back(SpectrumConvertApplication::InputFile(*original_file_name, *original_file_name, false));
  }
  // Launch threads to convert files
  boost::thread_group threadgroup_input_files;
  for (int t = 1; t < num_threads_; ++t) {
    boost::thread * currthread = new boost::thread(boost::bind(&SpectrumConvertApplication::getInputFiles, this, t));
    threadgroup_input_files.add_thread(currthread);
  }
  getInputFiles(0);
  // Join threads
  threadgroup_input_files.join_all();

  if (total_spectra_num_ > 0) {
    carp(CARP_INFO, "There were a total of %d spectrum conversions from %d input spectrum files.",
         total_spectra_num_, inputFiles_.size());
  }
  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6);

  // Create the active_peptide_queues and peptide_readers for each threads
  vector<HeadedRecordReader*> peptide_reader_threads;
  vector<ActivePeptideQueue*> APQ;
  for (int i = 0; i < num_threads_; i++) {
    peptide_reader_threads.push_back(new HeadedRecordReader(peptides_file, &peptides_header));
    APQ.push_back(new ActivePeptideQueue(peptide_reader_threads.back()->Reader(), proteins, &locations));
  }

  // Delete stuffs
  if (out_tsv_target_ != NULL)
    delete out_tsv_target_;
  if (out_tsv_decoy_ != NULL)
    delete out_tsv_decoy_;
  if (out_mztab_target_ != NULL)
    delete out_mztab_target_;
  if (out_mztab_decoy_ != NULL)
    delete out_mztab_decoy_;
  if (out_pin_target_ != NULL)
    delete out_pin_target_;
  if (out_pin_decoy_ != NULL)
    delete out_pin_decoy_;

  return 0;
}


vector<int> SpectrumConvertApplication::getNegativeIsotopeErrors() {
  string isotope_errors_string = Params::GetString("isotope-error");
  if (isotope_errors_string[0] == ',') {
    carp(CARP_FATAL, "Error in isotope_error parameter formatting: (%s)",
         isotope_errors_string.c_str());
  }
  for (string::const_iterator i = isotope_errors_string.begin(); i != isotope_errors_string.end(); i++) {
    if (*i == ',' && (i+1 == isotope_errors_string.end() || *(i+1) == ',')) {
      carp(CARP_FATAL, "Error in isotope_error parameter formatting: (%s) ", isotope_errors_string.c_str());
    }
  }

  vector<int> negative_isotope_errors(1, 0);
  if (!isotope_errors_string.empty()) {
    vector<int> isotope_errors = StringUtils::Split<int>(isotope_errors_string, ',');
    for (vector<int>::iterator i = isotope_errors.begin(); i != isotope_errors.end(); i++) {
      if (*i < 0) {
        carp(CARP_FATAL, "Found a negative isotope error: %d.", *i);
      } else if (find(negative_isotope_errors.begin(), negative_isotope_errors.end(),
                      -(*i)) != negative_isotope_errors.end()) {
        carp(CARP_FATAL, "Found duplicate when parsing isotope_error parameter: %d", *i);
      }
      negative_isotope_errors.push_back(-(*i));
    }
  }
  sort(negative_isotope_errors.begin(), negative_isotope_errors.end());
  return negative_isotope_errors;
}


void SpectrumConvertApplication::getPeptideIndexData(const string input_index, ProteinVec& proteins, vector<const pb::AuxLocation*>& locations, pb::Header& peptides_header){

  string peptides_file = FileUtils::Join(input_index, "pepix");  
  string proteins_file = FileUtils::Join(input_index, "protix");
  string auxlocs_file = FileUtils::Join(input_index, "auxlocs");  
  string residue_stats_file = FileUtils::Join(input_index, "residue_stat");  

  // Read protein index file
  carp(CARP_INFO, "Reading index %s", input_index.c_str());

  pb::Header protein_header;

  if (!ReadRecordsToVector<pb::Protein, const pb::Protein>(&proteins, proteins_file, &protein_header)) {
    carp(CARP_FATAL, "Error reading index (%s)", proteins_file.c_str());
  }
  // There shouldn't be more than one header in the protein pb.
  pb::Header_Source headerSource = protein_header.source(0);  
  string decoy_prefix = "";
  if (headerSource.has_decoy_prefix()){
    decoy_prefix = headerSource.decoy_prefix();
  } else {
    carp(CARP_WARNING, "You are using an index generated by an old version of tide-index."
                       "This will not affect your results, but this index may need to be "
                       "re-created to work with future versions of tide-index. ");
  }
  TideMatchSet::decoy_prefix_ = decoy_prefix;  
  if (headerSource.has_filename()){
    TideMatchSet::fasta_file_name_ = headerSource.filename();
  }
  
  // Read auxlocs index file
  ReadRecordsToVector<pb::AuxLocation>(&locations, auxlocs_file);

  // Read peptides index file

  if ((peptides_header.file_type() != pb::Header::PEPTIDES) ||
      !peptides_header.has_peptides_header()) {
    carp(CARP_FATAL, "Error reading index (%s)", peptides_file.c_str());
  }

  const pb::Header::PeptidesHeader& pepHeader = peptides_header.peptides_header();
  decoy_num_ = pepHeader.has_decoys_per_target() ? pepHeader.decoys_per_target() : 0;

  // Initizalize the Mass Constants class
  MassConstants::Init(&pepHeader.mods(), &pepHeader.nterm_mods(), &pepHeader.cterm_mods(),
      &pepHeader.nprotterm_mods(), &pepHeader.cprotterm_mods(), bin_width_, bin_offset_);

  if (curScoreFunction_ == PVALUES) {

    // Get the terminal mass bins for RES-EV
    //TODO assumption is that there is one nterm static mod per peptide
    if (pepHeader.nterm_mods().static_mod_size() > 0) {
      nTermMassBin_ = MassConstants::mass2bin(
                        MassConstants::mono_h + pepHeader.nterm_mods().static_mod(0).delta());
      nTermMass_ = MassConstants::mono_h + pepHeader.nterm_mods().static_mod(0).delta();
    } else {
      nTermMassBin_ = MassConstants::mass2bin(MassConstants::mono_h);
      nTermMass_ = MassConstants::mono_h;
    }

    //TODO assumption is that there is one cterm mod per peptide
    if (pepHeader.cterm_mods().static_mod_size() > 0) {
      cTermMassBin_ = MassConstants::mass2bin(MassConstants::mono_oh + pepHeader.cterm_mods().static_mod(0).delta());
      cTermMass_ = MassConstants::mono_oh + pepHeader.cterm_mods().static_mod(0).delta();
    } else {
      cTermMassBin_ = MassConstants::mass2bin(MassConstants::mono_oh);
      cTermMass_ = MassConstants::mono_oh;
    }

    if (FileUtils::Exists(residue_stats_file)) {  // The amino acid frequency calculation is done in tide-index.
      RecordReader residue_stat_reader = RecordReader(residue_stats_file);
      bool done;
      while (!(done = residue_stat_reader.Done())) {
        pb::ResidueStats last_residue_stat;
        residue_stat_reader.Read(&last_residue_stat);    
        double aa_mass = last_residue_stat.aamass();
        dAAMass_.push_back(aa_mass);
        dAAFreqN_.push_back(last_residue_stat.aafreqn());
        dAAFreqI_.push_back(last_residue_stat.aafreqi());
        dAAFreqC_.push_back(last_residue_stat.aafreqc());
        iAAMass_.push_back(MassConstants::mass2bin(aa_mass));
        string aa_str = last_residue_stat.aa_str();
        mMass2AA_[aa_mass] = aa_str;
      }
    } else {
      carp(CARP_INFO, "You are using an old format of the peptide index data. Please recreate your peptide index data");
      // The following (this whole else branch) part should be removed later.
      // Calculate the Amino Acid Frequencies for the P-value calculation if it wasn't done with tide-index.
      unsigned int len;
      unsigned int i;
      unsigned int cntTerm = 0;
      unsigned int cntInside = 0;
      string peptide_seq;
      unsigned int residue_bin;
      string tempAA;
      int mod_precision  = Params::GetInt("mod-precision");
      const unsigned int MaxModifiedAAMassBin = MassConstants::ToFixPt(2000.0);;   //2000 is the maximum mass of a modified amino acid
      unsigned int* nvAAMassCounterN = new unsigned int[MaxModifiedAAMassBin];   //N-terminal amino acids
      unsigned int* nvAAMassCounterC = new unsigned int[MaxModifiedAAMassBin];   //C-terminal amino acids
      unsigned int* nvAAMassCounterI = new unsigned int[MaxModifiedAAMassBin];   //inner amino acids in the peptides
      memset(nvAAMassCounterN, 0, MaxModifiedAAMassBin * sizeof(unsigned int));
      memset(nvAAMassCounterC, 0, MaxModifiedAAMassBin * sizeof(unsigned int));
      memset(nvAAMassCounterI, 0, MaxModifiedAAMassBin * sizeof(unsigned int));

      pb::Peptide current_pb_peptide_;
      HeadedRecordReader peptide_reader = HeadedRecordReader(peptides_file, &peptides_header);
      while (!(peptide_reader.Done())) { //read all peptides form index
        peptide_reader.Read(&current_pb_peptide_);
        Peptide peptide(current_pb_peptide_, proteins);
        vector<double> residue_masses = peptide.getAAMasses(); //retrieves the amino acid masses, modifications included
        peptide_seq = peptide.Seq();
        len = current_pb_peptide_.length();

        vector<double> residue_mods(len, 0);  // Initialize a vecotr of peptide length  with zeros. 
        
        // Handle variable modifications
        if (current_pb_peptide_.has_nterm_mod()){ // Handle N-terminal modifications
          int index;
          double delta;
          MassConstants::DecodeMod(ModCoder::Mod(current_pb_peptide_.nterm_mod()), &index, &delta);
          residue_mods[index] = delta;
        }

        for (i = 0; i < current_pb_peptide_.modifications_size(); ++i) {
          int index;
          double delta;
          MassConstants::DecodeMod(current_pb_peptide_.modifications(i), &index, &delta);
          residue_mods[index] = delta;
        }
        if (current_pb_peptide_.has_cterm_mod()){  // Handle C-terminal modifications
          int index;
          double delta;
          MassConstants::DecodeMod(ModCoder::Mod(current_pb_peptide_.cterm_mod()), &index, &delta);
          residue_mods[index] = delta;
        }
        // count AA masses
        residue_bin = MassConstants::ToFixPt(residue_masses[0]);
        ++nvAAMassCounterN[residue_bin];  // N-temrianl
        if (nvAAMassCounterN[residue_bin] == 1){
          tempAA = peptide_seq[0];
          if (residue_mods[0] != 0) {
            tempAA += "[" + StringUtils::ToString(residue_mods[0], mod_precision) + ']';
          }
          mMass2AA_[MassConstants::ToDouble(residue_bin)] = tempAA;
        }
        for (i = 1; i < len-1; ++i) {
          residue_bin = MassConstants::ToFixPt(residue_masses[i]);
          ++nvAAMassCounterI[residue_bin];  // non-terminal
          if (nvAAMassCounterI[residue_bin] == 1){
            tempAA = peptide_seq[i];
            if (residue_mods[i] != 0) {
              tempAA += "[" + StringUtils::ToString(residue_mods[i], mod_precision) + ']';
            }
            mMass2AA_[MassConstants::ToDouble(residue_bin)] = tempAA;
          }		
          ++cntInside;
        }
        residue_bin = MassConstants::ToFixPt(residue_masses[len - 1]);
        ++nvAAMassCounterC[residue_bin];  // C-temrinal
        if (nvAAMassCounterC[residue_bin] == 1){
          tempAA = peptide_seq[len - 1];
          if (residue_mods[len - 1] != 0) {
            tempAA += "[" + StringUtils::ToString(residue_mods[len - 1], mod_precision) + ']';
          }
          mMass2AA_[MassConstants::ToDouble(residue_bin)] = tempAA;
        }
        ++cntTerm;
      }
      unsigned int uiUniqueMasses = 0;
      double aa_mass;
      for (i = 0; i < MaxModifiedAAMassBin; ++i) {
        if (nvAAMassCounterN[i] || nvAAMassCounterI[i] || nvAAMassCounterC[i]) {
          ++uiUniqueMasses;
          aa_mass = MassConstants::ToDouble(i);
          dAAMass_.push_back(aa_mass);
          dAAFreqN_.push_back((double)nvAAMassCounterN[i] / cntTerm);
          dAAFreqI_.push_back((double)nvAAMassCounterI[i] / cntInside);
          dAAFreqC_.push_back((double)nvAAMassCounterC[i] / cntTerm);
          iAAMass_.push_back(MassConstants::mass2bin(aa_mass));          
        }
      }
      // Write the statistics to a Protocol Buffer, so that next time it is at hand
      RecordWriter residue_stat_writer = RecordWriter(residue_stats_file);
      CHECK(residue_stat_writer.OK());  
      for (int i = 0; i < dAAMass_.size(); ++i){
        pb::ResidueStats last_residue_stat;
        last_residue_stat.set_aamass(dAAMass_[i]);
        last_residue_stat.set_aafreqn(dAAFreqN_[i]);
        last_residue_stat.set_aafreqi(dAAFreqI_[i]);
        last_residue_stat.set_aafreqc(dAAFreqC_[i]);
        string aa_str = mMass2AA_[dAAMass_[i]];
        last_residue_stat.set_aa_str(aa_str);
        CHECK(residue_stat_writer.Write(&last_residue_stat));
      }       
      delete[] nvAAMassCounterN;
      delete[] nvAAMassCounterI;
      delete[] nvAAMassCounterC;
    } // Finished calculating AA frequencies
  }  
}

string SpectrumConvertApplication::getOutputFileName() {
  return output_file_name_;
}

// In order to add more options, you need to add them to ./src/util/Params.cpp
vector<string> SpectrumConvertApplication::getOptions() const {
  string arr[] = {
    "auto-mz-bin-width",
    "auto-precursor-window",
    "concat",
    "deisotope",
    "elution-window-size",
    "fileroot",
    "fragment-tolerance",
    "isotope-error",
    "mass-precision",
    "max-precursor-charge",
    "min-precursor-charge",
    "min-peaks",
    "mod-precision",
    "mz-bin-offset",
    "mz-bin-width",
    "mzid-output",
    "mztab-output",
    "num-threads",
    "output-dir",
    "override-charges",
    "overwrite",
    "parameter-file",
    "pepxml-output",
    "pin-output",
    "pm-charges",
    "pm-max-frag-mz",
    "pm-max-precursor-delta-ppm",
    "pm-max-precursor-mz",
    "pm-max-scan-separation",
    "pm-min-common-frag-peaks",
    "pm-min-frag-mz",
    "pm-min-peak-pairs",
    "pm-min-precursor-mz",
    "pm-min-scan-frag-peaks",
    "pm-pair-top-n-frag-peaks",
    "pm-top-n-frag-peaks",
    "precision",
    "precursor-window",
    "precursor-window-type",
    "print-search-progress",
    "remove-precursor-peak",
    "remove-precursor-tolerance",
    "scan-number",
    "score-function",
    "skip-preprocessing",
    "spectrum-max-mz",
    "spectrum-min-mz",
    "spectrum-parser",
    "sqt-output",
    "store-index",
    "store-spectra",
    "top-match",
    "txt-output",
    "use-flanking-peaks",
    "use-neutral-loss-peaks",
    "use-z-line",
    "verbosity"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}
/*
Parameters to be removed from ./Params.cpp"
brief-output
peptide-centric-search
use-tailor-calibration
exact-p-value
charge-state
evidence-granularity
*/

string SpectrumConvertApplication::getName() const {
  return "tide-search";
}

string SpectrumConvertApplication::getDescription() const {
  return
    "[[nohtml:Search a collection of spectra against a sequence database, "
    "returning a collection of peptide-spectrum matches (PSMs). This is a "
    "fast search engine but requires that you first build an index with "
    "tide-index.]]"
    "[[html:<p>Tide is a tool for identifying peptides from tandem mass "
    "spectra. It is an independent reimplementation of the SEQUEST<sup>&reg;"
    "</sup> algorithm, which assigns peptides to spectra by comparing the "
    "observed spectra to a catalog of theoretical spectra derived from a "
    "database of known proteins. Tide's primary advantage is its speed. Our "
    "published paper provides more detail on how Tide works. If you use Tide "
    "in your research, please cite:</p><blockquote>Benjamin J. Diament and "
    "William Stafford Noble. <a href=\"http://dx.doi.org/10.1021/pr101196n\">"
    "&quot;Faster SEQUEST Searching for Peptide Identification from Tandem "
    "Mass Spectra&quot;</a>. <em>Journal of Proteome Research</em>. "
    "10(9):3871-9, 2011.</blockquote> "
    "<p>When <code>tide-search</code> runs, it performs "
    "several intermediate steps, as follows:</p><ol>"
    "<li>If a FASTA file was provided, convert it to an index using "
    "<code>tide-index</code>.</li>"
    "<li>Convert the given "
    "fragmentation spectra to a binary format.</li><li>Search the spectra "
    "against the database and store the results in binary format.</li><li>"
    "Convert the results to one or more requested output formats.</li></ol><p>"
    "By default, the intermediate binary files are stored in the output "
    "directory and deleted when Tide finishes execution. If you plan to search "
    "against given database more than once or search a given set of spectra "
    "more than once, then you can direct Tide to save the binary spectrum "
    "files using the <code>--store-index</code> and "
    "<code>--store-spectra</code> options. "
    "Subsequent runs of the program will go faster "
    "if provided with inputs in binary format.</p>]]";
}

vector<string> SpectrumConvertApplication::getArgs() const {
  string arr[] = {
    "tide spectra file+",
    "tide database"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}


vector< pair<string, string> > SpectrumConvertApplication::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("tide-search.target.txt",
    "a tab-delimited text file containing the target PSMs. See <a href=\""
    "../file-formats/txt-format.html\">txt file format</a> for a list of the fields."));
  outputs.push_back(make_pair("tide-search.decoy.txt",
    "a tab-delimited text file containing the decoy PSMs. This file will only "
    "be created if the index was created with decoys."));
  outputs.push_back(make_pair("tide-search.params.txt",
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other Crux programs."));
  outputs.push_back(make_pair("tide-search.log.txt",
    "a log file containing a copy of all messages that were printed to the "
    "screen during execution."));
  return outputs;
}
bool SpectrumConvertApplication::needsOutputDirectory() const {
  return true;
}

COMMAND_T SpectrumConvertApplication::getCommand() const {
  return TIDE_SEARCH_COMMAND;
}

void SpectrumConvertApplication::processParams() {
  const string index = Params::GetString("tide database");

  if (!FileUtils::Exists(index)) {
    carp(CARP_FATAL, "'%s' does not exist", index.c_str());
  
  } else if (FileUtils::IsRegularFile(index)) {
    // Index is FASTA file

    carp(CARP_INFO, "Creating index from '%s'", index.c_str());
    string targetIndexName = Params::GetString("store-index");
    if (targetIndexName.empty()) {
      targetIndexName = FileUtils::Join(Params::GetString("output-dir"),
                                        "tide-search.tempindex");
      remove_index_ = targetIndexName;
    }
    TideIndexApplication indexApp;
    indexApp.processParams();
    if (indexApp.main(index, targetIndexName) != 0) {
      carp(CARP_FATAL, "tide-index failed.");
    }
    Params::Set("tide database", targetIndexName);
  
  } else {
    // Index is Tide index directory

    pb::Header peptides_header;
    string peptides_file = FileUtils::Join(index, "pepix");
    HeadedRecordReader peptide_reader(peptides_file, &peptides_header);
    if ((peptides_header.file_type() != pb::Header::PEPTIDES) ||
        !peptides_header.has_peptides_header()) {
      carp(CARP_FATAL, "Error reading index (%s).", peptides_file.c_str());
    }

    const pb::Header::PeptidesHeader& pepHeader = peptides_header.peptides_header();

    Params::Set("enzyme", pepHeader.enzyme());
    const char* digestString =
      digest_type_to_string(pepHeader.full_digestion() ? FULL_DIGEST : PARTIAL_DIGEST);
    Params::Set("digestion", digestString);
    Params::Set("isotopic-mass", pepHeader.monoisotopic_precursor() ? "mono" : "average");
  }
  // run param-medic?
  const string autoPrecursor = Params::GetString("auto-precursor-window");
  const string autoFragment = Params::GetString("auto-mz-bin-width");
  if (autoPrecursor != "false" || autoFragment != "false") {
    if (autoPrecursor != "false" && Params::GetString("precursor-window-type") != "ppm") {
      carp(CARP_FATAL, "Automatic peptide mass tolerance detection is only supported with ppm "
                       "units. Please re-run with auto-precursor-window set to 'false' or "
                       "precursor-window-type set to 'ppm'.");
    }
    ParamMedic::RunAttributeResult errorCalcResult;
    ParamMedicApplication::processFiles(Params::GetStrings("tide spectra file"),
      true, false, &errorCalcResult, NULL);

    if (autoPrecursor != "false") {
      string fail = errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_PRECURSOR_FAILURE);
      if (fail.empty()) {
        double sigma = StringUtils::FromString<double>(
          errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_PRECURSOR_SIGMA));
        double prediction = StringUtils::FromString<double>(
          errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_PRECURSOR_PREDICTION));
        carp(CARP_INFO, "precursor ppm standard deviation: %f", sigma);
        carp(CARP_INFO, "precursor error estimate (ppm): %.2f", prediction);
        Params::Set("precursor-window", prediction);
      } else {
        carp(autoPrecursor == "fail" ? CARP_FATAL : CARP_ERROR,
             "failed to calculate precursor error: %s", fail.c_str());
      }
    }
    if (autoFragment != "false") {
      string fail = errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_FRAGMENT_FAILURE);
      if (fail.empty()) {
        double sigma = StringUtils::FromString<double>(
          errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_FRAGMENT_SIGMA));
        double prediction = StringUtils::FromString<double>(
          errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_FRAGMENT_PREDICTION));
        carp(CARP_INFO, "fragment ppm standard deviation: %f", sigma);
        carp(CARP_INFO, "Fragment bin size estimate (Th): %.4f", prediction);
        Params::Set("mz-bin-width", prediction);
      } else {
        carp(autoFragment == "fail" ? CARP_FATAL : CARP_ERROR,
             "failed to calculate fragment error: %s", fail.c_str());
      }
    }
  }  
  
}
void SpectrumConvertApplication::getInputFiles(int thread_id) {
  // Try to read all spectrum files as spectrumrecords, convert those that fail
  if (thread_id > inputFiles_.size())
    return;
  for (vector<SpectrumConvertApplication::InputFile>::iterator original_file_name = inputFiles_.begin()+thread_id; 
       original_file_name < inputFiles_.begin() + (inputFiles_.size()); 
       original_file_name = original_file_name + num_threads_) 
    {
    carp(CARP_DEBUG, "Start processing input files");
    bool keepSpectrumrecords = true;
    string original_name = (*original_file_name).OriginalName;
    string spectrumrecords = original_name;
    // Check if the input file is spectrum records of google protocol buffer
    pb::Header header;
    HeadedRecordReader reader(original_name, &header);
    if (header.file_type() != pb::Header::SPECTRA) {
      // converting to spectrumrecords file 

      carp(CARP_INFO, "Converting %s to spectrumrecords format", original_name.c_str());
      carp(CARP_DEBUG, "Elapsed time starting conversion: %.3g s", wall_clock() / 1e6);
      
      spectrumrecords = Params::GetString("store-spectra");
      keepSpectrumrecords = !spectrumrecords.empty();
      if (!keepSpectrumrecords) {
        spectrumrecords = make_file_path(FileUtils::BaseName( original_name) + ".spectrumrecords.tmp");
      } else if (inputFiles_.size() > 1) {
        carp(CARP_FATAL, "Cannot use store-spectra option with multiple input "
                         "spectrum files");
      }
      carp(CARP_DEBUG, "New spectrumrecords filename: %s", spectrumrecords.c_str());
      int spectra_num = 0;
      if (!SpectrumRecordWriter::convert(original_name, spectrumrecords, spectra_num)) {
        carp(CARP_FATAL, "Error converting %s to spectrumrecords format", original_name.c_str());
      }
      locks_array_[LOCK_SPECTRUM_READING]->lock();
      total_spectra_num_ += spectra_num;
      locks_array_[LOCK_SPECTRUM_READING]->unlock();

    }
    (*original_file_name).SpectrumRecords  = spectrumrecords;
    (*original_file_name).Keep = keepSpectrumrecords;
    carp(CARP_DEBUG, "Finish converting");
  }
}

void SpectrumConvertApplication::createOutputFiles() {
  
  // Create output files for the search results
   
  bool overwrite = Params::GetBool("overwrite");  
  bool concat = Params::GetBool("concat");

  string concat_file_name;
  string target_file_name;
  string decoy_file_name;

  // Get output files for tsv format
  concat_file_name = make_file_path("tide-search.txt");
  target_file_name = make_file_path("tide-search.target.txt");
  decoy_file_name  = make_file_path("tide-search.decoy.txt");

  if (overwrite) {
    remove(concat_file_name.c_str());  
    remove(target_file_name.c_str());  
    remove(decoy_file_name.c_str());  
  }

  if (Params::GetBool("txt-output") == true) {  // original tide-search output format in tab-delimited text files (txt)
    string header = TideMatchSet::getHeader(TIDE_SEARCH_TSV, tide_index_mzTab_file_path_);

    if (concat) {

      out_tsv_target_ = create_stream_in_path(concat_file_name.c_str(), NULL, overwrite);
      output_file_name_ = concat_file_name;
      *out_tsv_target_ << header; 

    } else {

      out_tsv_target_ = create_stream_in_path(target_file_name.c_str(), NULL, overwrite);
      output_file_name_ = target_file_name;
      *out_tsv_target_ << header; 
      if (decoy_num_ > 0) {
        out_tsv_decoy_ = create_stream_in_path(decoy_file_name.c_str(), NULL, overwrite);
        *out_tsv_decoy_ << header;
      }
    }  
  }

  // Get output files for mzTAB format
  concat_file_name = make_file_path("tide-search.mzTab");
  target_file_name = make_file_path("tide-search.target.mzTab");
  decoy_file_name  = make_file_path("tide-search.decoy.mzTab");

  if (overwrite) {
    remove(concat_file_name.c_str());  
    remove(target_file_name.c_str());  
    remove(decoy_file_name.c_str());  
  }

  if ( Params::GetBool("mztab-output") == true) {   // mzTAB output format

    // Get column headers
    string header = TideMatchSet::getHeader(TIDE_SEARCH_MZTAB_TSV, tide_index_mzTab_file_path_);  // Gets the column headers
    if (concat) {

      out_mztab_target_ = create_stream_in_path(concat_file_name.c_str(), NULL, overwrite);
      output_file_name_ = concat_file_name;
      *out_mztab_target_ << header; 

    } else {

      out_mztab_target_ = create_stream_in_path(target_file_name.c_str(), NULL, overwrite);
      output_file_name_ = target_file_name;
      *out_mztab_target_ << header; 

      if (decoy_num_ > 0) {
        out_mztab_decoy_ = create_stream_in_path(decoy_file_name.c_str(), NULL, overwrite);
        *out_mztab_decoy_ << header;
      }
    }  
  }

}
