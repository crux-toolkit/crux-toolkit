#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>

#include <cstdio>
#include "tide/abspath.h"
#include "tide/records_to_vector-inl.h"

#include "carp.h"
#include "parameter.h"
#include "SpectrumRecordWriter.h"
#include "TideSearchApplication.h"

extern AA_MOD_T* list_of_mods[MAX_AA_MODS];
extern int num_mods;

bool TideSearchApplication::HAS_DECOYS = false;

TideSearchApplication::TideSearchApplication() {
  exact_pval_search = false;
}

TideSearchApplication::~TideSearchApplication() {
}

int TideSearchApplication::main(int argc, char** argv) {

  const char* option_list[] = {
    "precursor-window",
    "precursor-window-type",
    "spectrum-min-mz",
    "spectrum-max-mz",
    "min-peaks",
    "spectrum-charge",
    "scan-number",
    "top-match",
    "store-spectra",
    "compute-sp",
    "txt-output",
    "sqt-output",
    "pepxml-output",
    "mzid-output",
    "pinxml-output",
    "fileroot",
    "output-dir",
    "overwrite",
    "parameter-file",
    "exact-p-value",
    "mz-bin-width",
    "mz-bin-offset",
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

  // Check spectrum-charge parameter
  string charge_string = get_string_parameter_pointer("spectrum-charge");
  int charge_to_search;
  if (charge_string == "all") {
    carp(CARP_DEBUG, "Searching all charge states");
    charge_to_search = 0;
  } else {
    charge_to_search = atoi(charge_string.c_str());
    if (charge_to_search < 1 || charge_to_search > 6) {
      carp(CARP_FATAL, "Invalid spectrum-charge value %s",
           charge_string.c_str());
    }
    carp(CARP_DEBUG, "Searching charge state %d", charge_to_search);
  }

  // Check scan-number parameter
  string scan_range = get_string_parameter_pointer("scan-number");
  int min_scan, max_scan;
  if (scan_range == "__NULL_STR") {
    min_scan = 0;
    max_scan = BILLION;
    carp(CARP_DEBUG, "Searching all scans");
  }
  else if (scan_range.find('-') == string::npos) {
    // Single scan
    min_scan = max_scan = atoi(scan_range.c_str());
    carp(CARP_DEBUG, "Searching single scan %d", min_scan);
  } else {
    if (!get_range_from_string(scan_range.c_str(), min_scan, max_scan)) {
      carp(CARP_FATAL, "The scan number range '%s' is invalid. "
           "Must be of the form <first>-<last>", scan_range.c_str());
    } else {
      if (min_scan > max_scan) {
        int tmp_scan = min_scan;
        min_scan = max_scan;
        max_scan = tmp_scan;
        carp(CARP_DEBUG, "Switched scan range min and max");
      }
      carp(CARP_DEBUG, "Searching scan range %d-%d", min_scan, max_scan);
    }
  }
  //check to compute exact P-value 
  exact_pval_search = get_boolean_parameter("exact-p-value");
  double binWidth   = get_double_parameter("mz-bin-width");
  double binOffset  = get_double_parameter("mz-bin-offset");

  // Check compute-sp parameter
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

  // Check for decoys and set static variable if found
  for (vector<const pb::Protein*>::const_iterator i = proteins.begin();
       i != proteins.end();
       ++i) {
    if (TideMatchSet::isDecoy((*i)->name())) {
      HAS_DECOYS = true;
      break;
    }
  }

  // Read peptides index file
  pb::Header peptides_header;
  HeadedRecordReader peptide_reader(peptides_file, &peptides_header);
  if (!peptides_header.file_type() == pb::Header::PEPTIDES ||
      !peptides_header.has_peptides_header()) {
    carp(CARP_FATAL, "Error reading index (%s)", peptides_file.c_str());
  }
  MassConstants::Init(&peptides_header.peptides_header().mods());
  active_peptide_queue = new ActivePeptideQueue(peptide_reader.Reader(), proteins);

  active_peptide_queue->SetBinSize(binWidth, binOffset);

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

  bool txt_only = !get_boolean_parameter("sqt-output") &&
                  !get_boolean_parameter("pepxml-output") &&
                  !get_boolean_parameter("mzid-output") &&
                  !get_boolean_parameter("pinxml-output");
  OutputFiles* output_files = NULL;
  ofstream* target_file = NULL;
  ofstream* decoy_file = NULL;
  if (!txt_only) {
    carp(CARP_DEBUG, "Using OutputFiles to write matches");
    output_files = new OutputFiles(this);
  } else {
    carp(CARP_DEBUG, "Using TideMatchSet to write matches");
    bool overwrite = get_boolean_parameter("overwrite");
    string target_file_name = make_file_path("tide-search.target.txt");
    target_file = create_stream_in_path(target_file_name.c_str(), NULL, overwrite);
    if (HAS_DECOYS) {
      string decoy_file_name = make_file_path("tide-search.decoy.txt");
      decoy_file = create_stream_in_path(decoy_file_name.c_str(), NULL, overwrite);
    }
  }

  // Do the search
  carp(CARP_INFO, "Running search");
  cleanMods();
  search(spectra.SpecCharges(), active_peptide_queue, proteins, window,
         window_type, get_double_parameter("spectrum-min-mz"),
         get_double_parameter("spectrum-max-mz"), min_scan, max_scan,
         get_int_parameter("min-peaks"), charge_to_search,
         get_int_parameter("top-match"), spectra.FindHighestMZ(),
         output_files, target_file, decoy_file, compute_sp);

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
  if (output_files) {
    delete output_files;
  }
  if (target_file) {
    target_file->close();
    delete target_file;
    if (decoy_file) {
      decoy_file->close();
      delete decoy_file;
    }
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
  double spectrum_min_mz,
  double spectrum_max_mz,
  int min_scan,
  int max_scan,
  int min_peaks,
  int search_charge,
  int top_matches,
  double highest_mz,
  OutputFiles* output_files,
  ofstream* target_file,
  ofstream* decoy_file,
  bool compute_sp
) {

  if (output_files) {
    output_files->exact_pval_search = exact_pval_search;
    output_files->writeHeaders();
  } else if (target_file) {
    TideMatchSet::writeHeaders(target_file, false, compute_sp, exact_pval_search);
    TideMatchSet::writeHeaders(decoy_file, true, compute_sp, exact_pval_search);
  }

  // This is the main search loop.
  ObservedPeakSet observed;
  // cycle through spectrum-charge pairs, sorted by neutral mass
  for (vector<SpectrumCollection::SpecCharge>::const_iterator sc = spec_charges->begin();
       sc != spec_charges->end();
       ++sc) {
    const Spectrum* spectrum = sc->spectrum;

    double precursor_mz = spectrum->PrecursorMZ();
    int charge = sc->charge;
    int scan_num = spectrum->SpectrumNumber();

    if (precursor_mz < spectrum_min_mz || precursor_mz > spectrum_max_mz ||
        scan_num < min_scan || scan_num > max_scan ||
        spectrum->Size() < min_peaks ||
        (search_charge != 0 && charge != search_charge)) {
      continue;
    }

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
    if (exact_pval_search == false) {  //execute original tide-search program

      int size = active_peptide_queue->SetActiveRange(min_mass, max_mass);
      TideMatchSet::Arr2 match_arr2(size); // Scored peptides will go here.

      // Programs for taking the dot-product with the observed spectrum are laid
      // out in memory managed by the active_peptide_queue, one program for each
      // candidate peptide. The programs will store the results directly into
      // match_arr. We now pass control to those programs.
      collectScoresCompiled(active_peptide_queue, spectrum, observed, &match_arr2,
                            size, charge);

      // matches will arrange the results in a heap by score, return the top
      // few, and recover the association between counter and peptide. We output
      // the top matches.

      TideMatchSet::Arr match_arr(size);
      for (TideMatchSet::Arr2::iterator it = match_arr2.begin(); it != match_arr2.end(); ++it){
          TideMatchSet::Pair pair;
          pair.first.first = (double)(it->first/100000000.0);
	  pair.first.second = 0.0;
	  pair.second = it->second;
	  match_arr.push_back(pair);
      }

      TideMatchSet matches(&match_arr, highest_mz);
      matches.exact_pval_search = exact_pval_search;
      if (output_files) {
        matches.report(output_files, top_matches, spectrum, charge,
                       active_peptide_queue, proteins, compute_sp);
      } else {
        matches.report(target_file, decoy_file, top_matches, spectrum, charge,
                       active_peptide_queue, proteins, compute_sp);
      }
    } else {  //execute exact-pval-search 
      int size = active_peptide_queue -> SetActiveRangeBIons( min_mass, max_mass );
      TideMatchSet::Arr match_arr( size ); // Scored peptides will go here.
      deque< TheoreticalPeakSetBIons >::const_iterator iter1_;
      vector< unsigned int >::const_iterator it;
      //here comes Jeff's exactPvalue calculation
	  int count = 0;
      for ( iter1_ = active_peptide_queue->b_ion_queue_.begin(); iter1_ != active_peptide_queue->b_ion_queue_.end(); ++iter1_ ) {
		cout << count << "   ***   ";
		count += 1;
 	    for ( it = iter1_ -> unordered_peak_list_.begin(); it != iter1_ -> unordered_peak_list_.end(); it++ ) {
		   cout << *it << "\t";
	    } 
        cout << endl;
      }

      //&& for test only
      for ( TideMatchSet::Arr::iterator it = match_arr.begin(); it != match_arr.end(); it++ ) {
        it -> first.first = 1.0;
	    it -> first.second = 2.0;
	    it -> second = 10;
      }
      for ( TideMatchSet::Arr::iterator it = match_arr.begin(); it != match_arr.end(); it++ ) {
        cout << it -> first.first << "\t";
	    cout << it -> first.second << "\t";
	    cout << it -> second << "\t";
		cout << endl;
      }
	  //&& for test only
	  
//************************************************************************************	  
/* For one observed spectrum, calculates:
 *  - vector of cleavage evidence
 *  - score count vectors for a range of integer masses
 *  - p-values of XCorr match scores between spectrum and all selected candidate target and decoy peptides
 * Written by Jeff Howbert, October, 2013.
 * Ported to and integrated with Tide by Jeff Howbert, November, 2013.
 */
      //&& need to decide which of these constants to move to parameter file
      const double evidenceIntScale = 500.0;   // make score bins 10 times larger than with usual value = 50.0;
      const double binWidth = 1.0005079;       // monoisotopic
      const double binOffset = 0.40;           // changed from previous standard value 0.68
      const int minDeltaMass = 57;
      const int maxDeltaMass = 186;
      //&& end need to decide

      // right hand side variables
      int maxPrecurMass;
      double precurMz;
      int precurCharge;
      double* ionSeriesMass;
      double* ionSeriesIntens;
      int* pepIdxTarget;
      int* pepIdxDecoy;
      double* pepMassMonoTarget;
      double* pepMassMonoDecoy;
      int* intensArrayTheorTarget;
      int* intensArrayTheorDecoy;
      int* nPeakIntensArrayTheorTarget;
      int* nPeakIntensArrayTheorDecoy;
      double* aaProb1;
      double* aaProb2;

      // left hand side variables
      double pValueSpectrumBestMatchTar;
      double pValueSpectrumBestMatchDec;
      int pepPValueSpectrumBestMatchTar;
      int pepPValueSpectrumBestMatchDec;

      // internal variables
      int pe;
      double pepMass;
      int pepMassInt;
      int ma;

	//&& working here
    // int nIon = ( int )mxGetM( prhs[ 3 ] );
    // int nPepIdxTarget = ( int )mxGetM( prhs[ 5 ] );
    // int nPepIdxDecoy = ( int )mxGetM( prhs[ 6 ] );
    // int nPepSeqTarget = ( int )mxGetM( prhs[ 9 ] );
    // int nPepSeqDecoy = ( int )mxGetM( prhs[ 10 ] );
    // int maxTheorPeak = ( int )mxGetN( prhs[ 9 ] );
    
    // double* pValueTarget = new double[ nPepIdxTarget ];
    // double* pValueDecoy = new double[ nPepIdxDecoy ];

    // std::vector< int > pepMassIntUnique;
    // pepMassIntUnique.reserve( nPepIdxTarget + nPepIdxDecoy );
    // for ( pe = 0; pe < nPepIdxTarget; pe++ ) {
        // pepMass = pepMassMonoTarget[ pepIdxTarget[ pe ] - 1 ];     //&& take account of MATLAB base-one indexing
        // pepMassInt = ( int )floor( pepMass / binWidth + 1.0 - binOffset );
        // pepMassIntUnique.push_back( pepMassInt );
    // }
    // for ( pe = 0; pe < nPepIdxDecoy; pe++ ) {
        // pepMass = pepMassMonoDecoy[ pepIdxDecoy[ pe ] - 1 ];     //&& take account of MATLAB base-one indexing
        // pepMassInt = ( int )floor( pepMass / binWidth + 1.0 - binOffset );
        // pepMassIntUnique.push_back( pepMassInt );
    // }
    // std::sort( pepMassIntUnique.begin(), pepMassIntUnique.end(), sortAscInt );
    // std::vector< int >::iterator last = std::unique( pepMassIntUnique.begin(), pepMassIntUnique.end() );
    // pepMassIntUnique.erase( last, pepMassIntUnique.end() );
    // int nPepMassIntUniq = ( int )pepMassIntUnique.size();

    // int** evidenceObs = new int* [ nPepMassIntUniq ];
    // int* scoreOffsetObs = new int [ nPepMassIntUniq ];
    // double** pValueScoreObs = new double* [ nPepMassIntUniq ];
    // int* intensArrayTheor = new int [ maxPrecurMass ];       // initialized later in loop
    // for ( pe = 0; pe < nPepMassIntUniq; pe++ ) {
        // evidenceObs[ pe ] = new int[ maxPrecurMass ];
        // for ( ma = 0; ma < maxPrecurMass; ma++ ) {
            // evidenceObs[ pe ][ ma ] = 0;
        // }
        // scoreOffsetObs[ pe ] = 0;
        // pepMassInt = pepMassIntUnique[ pe ];
        // // preprocess to create one integerized evidence vector for each cluster of masses among selected peptides
        // double pepMassMonoMean = ( pepMassInt - 1.0 + binOffset ) * binWidth + 0.5;
        // createEvidenceArrayObserved( nIon, ionSeriesMass, ionSeriesIntens, evidenceIntScale, binWidth, binOffset, precurMz, precurCharge, pepMassMonoMean, maxPrecurMass, evidenceObs[ pe ] );

        // // NOTE: will have to go back to separate dynamic programming for target and decoy if they have different probNI and probC
        // int maxEvidence = *std::max_element( evidenceObs[ pe ], evidenceObs[ pe ] + maxPrecurMass );
        // int minEvidence = *std::min_element( evidenceObs[ pe ], evidenceObs[ pe ] + maxPrecurMass );
        // // estimate maxScore and minScore
        // int maxNResidue = ( int )floor( ( double )pepMassInt / ( double )minDeltaMass );
        // std::vector< int > sortEvidenceObs ( evidenceObs[ pe ], evidenceObs[ pe ] + maxPrecurMass );
        // std::sort( sortEvidenceObs.begin(), sortEvidenceObs.end(), sortDescInt );
        // int maxScore = 0;
        // int minScore = 0;
        // for ( int sc = 0; sc < maxNResidue; sc++ ) {
            // maxScore += sortEvidenceObs[ sc ];
        // }
        // for ( int sc = maxPrecurMass - maxNResidue; sc < maxPrecurMass; sc++ ) {
            // minScore += sortEvidenceObs[ sc ];
        // }
        // int bottomRowBuffer = maxEvidence + 1;
        // int topRowBuffer = -minEvidence;
        // int nRowDynProg = bottomRowBuffer - minScore + 1 + maxScore + topRowBuffer;
        // pValueScoreObs[ pe ] = new double[ nRowDynProg ];
        // scoreOffsetObs[ pe ] = calcScoreCount( maxPrecurMass, evidenceObs[ pe ], binWidth, binOffset, pepMassInt, maxEvidence, minEvidence, maxScore, minScore, aaProb1, aaProb2, pValueScoreObs[ pe ] );
    // }

    // // ***** find best target match ***********************************
    // double pValueSpecBestMatchTar = 1.0;
    // int pepPValueSpecBestMatchTar = 0;
    // for ( pe = 0; pe < nPepIdxTarget; pe++ ) {
        // int pepIdxTar = pepIdxTarget[ pe ] - 1;
        // pepMassInt = ( int )floor( pepMassMonoTarget[ pepIdxTar ] / binWidth + 1.0 - binOffset );   // for indexing into preprocessed evidence vectors
        // int pepMassIntIdx = 0;
        // for ( ma = 0; ma < nPepMassIntUniq; ma++ ) {
            // if ( pepMassIntUnique[ ma ] == pepMassInt ) {
                // pepMassIntIdx = ma;
                // break;
            // }
        // }      
        // // score XCorr for target peptide with integerized evidenceObs array
        // for ( ma = 0; ma < maxPrecurMass; ma++ ) {
            // intensArrayTheor[ ma ] = 0;
        // }
        // int nPeakIntensTheor = nPeakIntensArrayTheorTarget[ pepIdxTar ];
        // for ( ma = 0; ma < nPeakIntensTheor; ma++ ) {
            // intensArrayTheor[ intensArrayTheorTarget[ ma * nPepSeqTarget + pepIdxTar ] ] = 1;
        // }
        // int scoreRefactInt = 0;
        // for ( ma = 0; ma < maxPrecurMass; ma++ ) {
            // scoreRefactInt += evidenceObs[ pepMassIntIdx ][ ma ] * intensArrayTheor[ ma ];
        // }
        // int scoreCountIdx = scoreRefactInt + scoreOffsetObs[ pepMassIntIdx ];
        // double pValueTar = pValueScoreObs[ pepMassIntIdx ][ scoreCountIdx ];
        // if ( pValueSpecBestMatchTar > pValueTar ) {
            // pValueSpecBestMatchTar = pValueTar;
            // pepPValueSpecBestMatchTar = pepIdxTar;
        // }
// //         nSpecTarget( pepIdxTar ) = nSpecTarget( pepIdxTar ) + 1;
// //         if ( pValueTar < pValuePeptideBestMatchTarget( pepIdxTar ) )
// //             pValuePeptideBestMatchTarget( pepIdxTar ) = pValueTar;
// //             specPValuePeptideBestMatchTarget( pepIdxTar ) = spChIdx;
// //         end
    // }

    // // ***** find best decoy match ***********************************
    // double pValueSpecBestMatchDec = 1.0;
    // int pepPValueSpecBestMatchDec = 0;
    // for ( pe = 0; pe < nPepIdxDecoy; pe++ ) {
        // int pepIdxDec = pepIdxDecoy[ pe ] - 1;
        // pepMassInt = ( int )floor( pepMassMonoDecoy[ pepIdxDec ] / binWidth + 1.0 - binOffset );   // for indexing into preprocessed evidence vectors
        // int pepMassIntIdx = 0;
        // for ( ma = 0; ma < nPepMassIntUniq; ma++ ) {
            // if ( pepMassIntUnique[ ma ] == pepMassInt ) {
                // pepMassIntIdx = ma;
                // break;
            // }
        // }      
        // // score XCorr for target peptide with integerized evidenceObs array
        // for ( ma = 0; ma < maxPrecurMass; ma++ ) {
            // intensArrayTheor[ ma ] = 0;
        // }
        // int nPeakIntensTheor = nPeakIntensArrayTheorDecoy[ pepIdxDec ];
        // for ( ma = 0; ma < nPeakIntensTheor; ma++ ) {
            // intensArrayTheor[ intensArrayTheorDecoy[ ma * nPepSeqDecoy + pepIdxDec ] ] = 1;
        // }
        // int scoreRefactInt = 0;
        // for ( ma = 0; ma < maxPrecurMass; ma++ ) {
            // scoreRefactInt += evidenceObs[ pepMassIntIdx ][ ma ] * intensArrayTheor[ ma ];
        // }
        // int scoreCountIdx = scoreRefactInt + scoreOffsetObs[ pepMassIntIdx ];
        // double pValueDec = pValueScoreObs[ pepMassIntIdx ][ scoreCountIdx ];
        // if ( pValueSpecBestMatchDec > pValueDec ) {
            // pValueSpecBestMatchDec = pValueDec;
            // pepPValueSpecBestMatchDec = pepIdxDec;
        // }
// //     nSpecDecoy( pepIdxDec ) = nSpecDecoy( pepIdxDec ) + 1;
// //     if ( pValueDec < pValuePeptideBestMatchDecoy( pepIdxDec ) )
// //         pValuePeptideBestMatchDecoy( pepIdxDec ) = pValueDec;
// //         specPValuePeptideBestMatchDecoy( pepIdxDec ) = spChIdx;
// //     end
// }

    // pValueSpectrumBestMatchTar[ 0 ] = pValueSpecBestMatchTar;
    // pValueSpectrumBestMatchDec[ 0 ] = pValueSpecBestMatchDec;   
    // pepPValueSpectrumBestMatchTar[ 0 ] = pepPValueSpecBestMatchTar + 1; //&& adding 1 just to match MATLAB base-one indexing; will probably get rid of after full port to C++
    // pepPValueSpectrumBestMatchDec[ 0 ] = pepPValueSpecBestMatchDec + 1;

    // // clean up
    // delete [] pValueTarget;
    // delete [] pValueDecoy;
    // delete [] scoreOffsetObs;
    // for( pe = 0; pe < nPepMassIntUniq; pe++ ) {
        // delete [] evidenceObs[ pe ];
        // delete [] pValueScoreObs[ pe ];
    // }
    // delete [] evidenceObs;
    // delete [] pValueScoreObs;
    // delete [] intensArrayTheor;
//************************************************************************************	  

      //report matches
      TideMatchSet matches(&match_arr, highest_mz);
      matches.exact_pval_search = exact_pval_search;
      if (output_files) {
        matches.report(output_files, top_matches, spectrum, charge,
	               active_peptide_queue, proteins, compute_sp);
      } else {
        matches.report(target_file, decoy_file, top_matches, spectrum, charge,
	               active_peptide_queue, proteins, compute_sp);
      } 
    }
  }
  if (output_files) {
    output_files->writeFooters();
  }
}

void TideSearchApplication::collectScoresCompiled(
  ActivePeptideQueue* active_peptide_queue,
  const Spectrum* spectrum,
  const ObservedPeakSet& observed,
  TideMatchSet::Arr2* match_arr,
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

bool TideSearchApplication::hasDecoys() {
  return HAS_DECOYS;
}

string TideSearchApplication::getName() {
  return "tide-search";
}

string TideSearchApplication::getDescription() {
  return
  "Search a collection of spectra against a sequence "
  "database, returning a collection of peptide-spectrum "
  "matches (PSMs) scored by XCorr.";
}

bool TideSearchApplication::needsOutputDirectory() {
  return true;
}

COMMAND_T TideSearchApplication::getCommand() {
  return TIDE_SEARCH_COMMAND;
}

int TideSearchApplication::calcScoreCount( int numelEvidenceObs, int* evidenceObs, double binWidth, double binOffset, int pepMassInt, int maxEvidence, int minEvidence, int maxScore, int minScore, double* aaProb1, double* aaProb2, double* pValueScoreObs ) {
/* Calculates counts of peptides with various XCorr scores, given a preprocessed
 * MS2 spectrum, using dynamic programming.
 * Written by Jeff Howbert, October, 2012 (as function calcScoreCount).
 * Ported to and integrated with Tide by Jeff Howbert, November, 2013.
 */

// int calcScoreCount( int numelEvidenceObs, int* evidenceObs, double binWidth, double binOffset,
//		int pepMassInt, int maxEvidence, int minEvidence, int maxScore, int minScore,
//		double* aaProb1, double* aaProb2, double* pValueScoreObs ) {

    const int nDeltaMass = 18;
    // expected delta masses in b/y ion ladders, single amino acid residues only
    int deltaMass[ nDeltaMass ] = {  
                       57, // G
                       71, // A
                       87, // S
                       97, // P
                       99, // V
                      101, // T
                      113, // I, L
                      114, // N
                      115, // D
                      128, // K, Q
                      129, // E
                      131, // M
                      137, // H
                      147, // F
                      156, // R
                      160, // C with acetamide mod
                      163, // Y
                      186  // W
                    };
    int minDeltaMass = deltaMass[ 0 ];
    int maxDeltaMass = deltaMass[ nDeltaMass - 1 ];

    // internal variables
    int row;
    int col;
    int ma;
    int evidence;
    int de;
    int evidenceRow;
    double sumScore;
    
    int bottomRowBuffer = maxEvidence + 1;
    int topRowBuffer = -minEvidence;
    int colBuffer = maxDeltaMass;
    int colStart = ( int )floor( 1.0 / binWidth + 1.0 - binOffset );
    int scoreOffsetObs = bottomRowBuffer - minScore;

    int nRow = bottomRowBuffer - minScore + 1 + maxScore + topRowBuffer;
    int nCol = colBuffer + pepMassInt;
    int rowFirst = bottomRowBuffer;
    int rowLast = rowFirst - minScore + maxScore;
    int colFirst = colStart + 1;
    int colLast = pepMassInt - 17;
    int initCountRow = bottomRowBuffer - minScore;
    int initCountCol = maxDeltaMass + colStart;

    double** dynProgArray = 0;
    dynProgArray = new double* [ nRow ];
    for ( row = 0; row < nRow; row++ ) {
        dynProgArray[ row ] = new double[ nCol ];
        for ( col = 0; col < nCol; col++ ) {
            dynProgArray[ row ][ col ] = 0.0;
        }
    }
    double* scoreCountBinAdjust = 0;
    scoreCountBinAdjust = new double [ nRow ];
    for ( row = 0; row < nRow; row++ ) {
        scoreCountBinAdjust[ row ] = 0.0;
    }
 
    dynProgArray[ initCountRow ][ initCountCol ] = 1.0;    // initial count of peptides with mass = 1
    int deltaMassCol[ nDeltaMass ];
    for ( ma = colFirst; ma < colLast; ma++ ) {
        col = maxDeltaMass + ma;
        evidence = evidenceObs[ ma ];
        for ( de = 0; de < nDeltaMass; de++ ) {
            deltaMassCol[ de ] = col - deltaMass[ de ];
        }
        for ( row = rowFirst; row <= rowLast; row++ ) {
            evidenceRow = row - evidence;
            sumScore = 0.0;
            for ( de = 0; de < nDeltaMass; de++ ) {
                sumScore += dynProgArray[ evidenceRow ][ deltaMassCol[ de ] ] * aaProb1[ de ];
            }
            dynProgArray[ row ][ col ] = sumScore;
        }
    }
    ma = colLast;
    col = maxDeltaMass + ma;
    evidence = evidenceObs[ ma ];
    for ( de = 0; de < nDeltaMass; de++ ) {
        deltaMassCol[ de ] = col - deltaMass[ de ];
    }
    for ( row = rowFirst; row <= rowLast; row++ ) {
        evidenceRow = row - evidence;
        sumScore = 0.0;
        for ( de = 0; de < nDeltaMass; de++ ) {
            sumScore += dynProgArray[ evidenceRow ][ deltaMassCol[ de ] ] * aaProb2[ de ];
        }
        dynProgArray[ row ][ col ] = sumScore;
    }

    int colScoreCount = maxDeltaMass + colLast;
    double totalCount = 0.0;
    for ( row = 0; row < nRow; row++ ) {
        // at this point pValueScoreObs just holds counts from last column of dynamic programming array
        pValueScoreObs[ row ] = dynProgArray[ row ][ colScoreCount ];
        totalCount += pValueScoreObs[ row ];
        scoreCountBinAdjust[ row ] = pValueScoreObs[ row ] / 2.0;
    }
    // convert from counts to cumulative sum of counts
    for ( row = nRow - 2; row >= 0; row-- ) {
        pValueScoreObs[ row ] += pValueScoreObs[ row + 1 ];
    }
    for ( row = 0; row < nRow; row++ ) {
        // adjust counts to reflect center of bin, not edge
        pValueScoreObs[ row ] -= scoreCountBinAdjust[ row ];
        // normalize distribution; use exp( log ) to avoid potential underflow
        double logTotalCount = log( totalCount );
        pValueScoreObs[ row ] = exp( log( pValueScoreObs[ row ] ) - logTotalCount );
    }
    
    // clean up
    for( int row = 0; row < nRow; row++ ) {
        delete [] dynProgArray[ row ];
    }
    delete [] dynProgArray;
    delete [] scoreCountBinAdjust;
    
    return scoreOffsetObs;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
