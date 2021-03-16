#include <cstdio>
#include "app/tide/abspath.h"
#include "app/tide/records_to_vector-inl.h"

#include "io/carp.h"
#include "parameter.h"
#include "io/SpectrumRecordWriter.h"
#include "TideIndexApplication.h"
#include "TideSearchApplication.h"
#include "DIAmeterApplication.h"
#include "ParamMedicApplication.h"
#include "PSMConvertApplication.h"
#include "tide/mass_constants.h"
#include "TideMatchSet.h"
#include "util/Params.h"
#include "util/FileUtils.h"
#include "util/StringUtils.h"


const double DIAmeterApplication::XCORR_SCALING = 100000000.0;
const double DIAmeterApplication::RESCALE_FACTOR = 20.0;

DIAmeterApplication::DIAmeterApplication() { /* do nothing */ }
DIAmeterApplication::~DIAmeterApplication() { /* do nothing */ }


int DIAmeterApplication::main(int argc, char** argv) {
  return main(Params::GetStrings("tide spectra file"));
}
int DIAmeterApplication::main(const vector<string>& input_files) {
  return main(input_files, Params::GetString("tide database"));
}

int DIAmeterApplication::main(const vector<string>& input_files, const string input_index) {
  carp(CARP_INFO, "Running tide-search...");

  const string index = input_index;
  string peptides_file = FileUtils::Join(index, "pepix");
  string proteins_file = FileUtils::Join(index, "protix");
  string auxlocs_file = FileUtils::Join(index, "auxlocs");

  // params
  double bin_width_  = Params::GetDouble("mz-bin-width");
  double bin_offset_ = Params::GetDouble("mz-bin-offset");
  vector<int> negative_isotope_errors = TideSearchApplication::getNegativeIsotopeErrors();
  // bool use_neutral_loss_peaks_ = Params::GetBool("use-neutral-loss-peaks");
  // bool use_flanking_peaks_ = Params::GetBool("use-flanking-peaks");

  // Read proteins index file
  ProteinVec proteins;
  pb::Header protein_header;
  if (!ReadRecordsToVector<pb::Protein, const pb::Protein>(&proteins, proteins_file, &protein_header)) { carp(CARP_FATAL, "Error reading index (%s)", proteins_file.c_str()); }

  // Read auxlocs index file
  vector<const pb::AuxLocation*> locations;
  if (!ReadRecordsToVector<pb::AuxLocation>(&locations, auxlocs_file)) { carp(CARP_FATAL, "Error reading index (%s)", auxlocs_file.c_str()); }

  // Read peptides index file
  pb::Header peptides_header;
  HeadedRecordReader* peptide_reader = new HeadedRecordReader(peptides_file, &peptides_header);
  if ((peptides_header.file_type() != pb::Header::PEPTIDES) || !peptides_header.has_peptides_header()) { carp(CARP_FATAL, "Error reading index (%s)", peptides_file.c_str()); }

  const pb::Header::PeptidesHeader& pepHeader = peptides_header.peptides_header();
  DECOY_TYPE_T headerDecoyType = (DECOY_TYPE_T)pepHeader.decoys();
  int decoysPerTarget = pepHeader.has_decoys_per_target() ? pepHeader.decoys_per_target() : 0;
  carp(CARP_DEBUG, "decoysPerTarget = %d \t headerDecoyType=%d.", decoysPerTarget, headerDecoyType);

  MassConstants::Init(&pepHeader.mods(), &pepHeader.nterm_mods(), &pepHeader.cterm_mods(), &pepHeader.nprotterm_mods(), &pepHeader.cprotterm_mods(), bin_width_, bin_offset_);
  ModificationDefinition::ClearAll();
  TideMatchSet::initModMap(pepHeader.mods(), ANY);
  TideMatchSet::initModMap(pepHeader.nterm_mods(), PEPTIDE_N);
  TideMatchSet::initModMap(pepHeader.cterm_mods(), PEPTIDE_C);
  TideMatchSet::initModMap(pepHeader.nprotterm_mods(), PROTEIN_N);
  TideMatchSet::initModMap(pepHeader.cprotterm_mods(), PROTEIN_C);


  string output_file_name_ = make_file_path("tide-search.txt");
  ofstream* output_file = create_stream_in_path(output_file_name_.c_str(), NULL, Params::GetBool("overwrite"));
  TideMatchSet::writeHeadersDIA(output_file, Params::GetBool("compute-sp"));

  map<string, double> peptide_predrt_map;
  getPeptidePredRTMapping(&peptide_predrt_map);


  vector<InputFile> ms1_spectra_files = getInputFiles(input_files, 1);
  vector<InputFile> ms2_spectra_files = getInputFiles(input_files, 2);


  // Loop through spectrum files
  for (pair<vector<InputFile>::const_iterator, vector<InputFile>::const_iterator> f(ms1_spectra_files.begin(), ms2_spectra_files.begin());
		  f.first != ms1_spectra_files.end() && f.second != ms2_spectra_files.end(); ++f.first, ++f.second) {

	  string ms1_spectra_file = (f.first)->SpectrumRecords;
	  string ms2_spectra_file = (f.second)->SpectrumRecords;
	  string origin_file = (f.second)->OriginalName;

	  map<int, pair<double*, double*>> ms1scan_intensity_rank_map;
	  avg_noise_intensity_logrank_ = loadMS1Spectra(ms1_spectra_file, &ms1scan_intensity_rank_map);
	  SpectrumCollection* spectra = loadSpectra(ms2_spectra_file);

	  // insert the search code here and will split into a new function later
	  double highest_ms2_mz = spectra->FindHighestMZ();
	  carp(CARP_DEBUG, "Maximum observed MS2 m/z = %f.", highest_ms2_mz);
	  MaxBin::SetGlobalMax(highest_ms2_mz);


	  // Active queue to process the indexed peptides
	  ActivePeptideQueue* active_peptide_queue = new ActivePeptideQueue(peptide_reader->Reader(), proteins);
	  active_peptide_queue->SetBinSize(bin_width_, bin_offset_);
	  active_peptide_queue->SetOutputs(NULL, &locations, Params::GetInt("top-match"), true, output_file, NULL, highest_ms2_mz);

	  // initialize fields required for output
	  // Not sure if it is necessary, will check later
	  const vector<SpectrumCollection::SpecCharge>* spec_charges = spectra->SpecCharges();
	  int* sc_index = new int(-1);
	  FLOAT_T sc_total = (FLOAT_T)spec_charges->size();
	  int print_interval = Params::GetInt("print-search-progress");
	  // Keep track of observed peaks that get filtered out in various ways.
	  long int num_range_skipped = 0;
	  long int num_precursors_skipped = 0;
	  long int num_isotopes_skipped = 0;
	  long int num_retained = 0;

	  // This is the main search loop.
	  // DIAmeter supports spectrum centric match report only
	  ObservedPeakSet observed(bin_width_, bin_offset_, Params::GetBool("use-neutral-loss-peaks"), Params::GetBool("use-flanking-peaks") );

	  for (vector<SpectrumCollection::SpecCharge>::const_iterator sc = spec_charges->begin();sc < spec_charges->begin() + (spec_charges->size()); sc++) {
		  ++(*sc_index);
		  if (print_interval > 0 && *sc_index > 0 && *sc_index % print_interval == 0) { carp(CARP_INFO, "%d spectrum-charge combinations searched, %.0f%% complete", *sc_index, *sc_index / sc_total * 100); }

		  Spectrum* spectrum = sc->spectrum;
		  double precursor_mz = spectrum->PrecursorMZ();
		  int scan_num = spectrum->SpectrumNumber();
		  int ms1_scan_num = spectrum->MS1SpectrumNum();
		  int charge = sc->charge;

		  // The active peptide queue holds the candidate peptides for spectrum.
		  // Calculate and set the window, depending on the window type.
		  vector<double>* min_mass = new vector<double>();
		  vector<double>* max_mass = new vector<double>();
		  vector<bool>* candidatePeptideStatus = new vector<bool>();
		  double min_range, max_range;

		  // TideSearchApplication::computeWindow(*sc, string_to_window_type(Params::GetString("precursor-window-type")), Params::GetDouble("precursor-window"), Params::GetInt("max-precursor-charge"), &negative_isotope_errors, min_mass, max_mass, &min_range, &max_range);
		  computeWindowDIA(*sc, Params::GetInt("max-precursor-charge"), &negative_isotope_errors, min_mass, max_mass, &min_range, &max_range);
		  carp(CARP_DETAILED_DEBUG, "MS1Scan:%d \t MS2Scan:%d \t precursor_mz:%f \t charge:%d \t mass_range:[%f,%f]", ms1_scan_num, scan_num, precursor_mz, charge, min_range, max_range);


		  // Normalize the observed spectrum and compute the cache of frequently-needed
		  // values for taking dot products with theoretical spectra.
		  observed.PreprocessSpectrum(*spectrum, charge, &num_range_skipped, &num_precursors_skipped, &num_isotopes_skipped, &num_retained);
		  int nCandPeptide = active_peptide_queue->SetActiveRange(min_mass, max_mass, min_range, max_range, candidatePeptideStatus);
		  int candidatePeptideStatusSize = candidatePeptideStatus->size();
		  carp(CARP_DETAILED_DEBUG, "nCandPeptide:%d \t candidatePeptideStatusSize:%d \t mass range:[%f,%f]", nCandPeptide, candidatePeptideStatusSize, min_range, max_range);
		  if (nCandPeptide == 0) { continue; }

		  TideMatchSet::Arr2 match_arr2(candidatePeptideStatusSize); // Scored peptides will go here.
		  // Programs for taking the dot-product with the observed spectrum are laid
		  // out in memory managed by the active_peptide_queue, one program for each
		  // candidate peptide. The programs will store the results directly into
		  // match_arr. We now pass control to those programs.
		  TideSearchApplication::collectScoresCompiled(active_peptide_queue, spectrum, observed, &match_arr2, candidatePeptideStatusSize, charge);

		  // The denominator used in the Tailor score calibration method
		  double quantile_score = getTailorQuantile(&match_arr2);
		  // carp(CARP_DETAILED_DEBUG, "Tailor quantile_score:%f", quantile_score);

		  TideMatchSet::Arr match_arr(nCandPeptide);
		  for (TideMatchSet::Arr2::iterator it = match_arr2.begin(); it != match_arr2.end(); ++it) {
			  /// Yang: I cannot understand here that the peptide index and the rank is complement.
			  /// More attention is needed.
			  int peptide_idx = candidatePeptideStatusSize - (it->second);
			  if ((*candidatePeptideStatus)[peptide_idx]) {
				  TideMatchSet::Scores curScore;
			      curScore.xcorr_score = (double)(it->first / XCORR_SCALING);
			      curScore.rank = it->second;
			      curScore.tailor = ((double)(it->first / XCORR_SCALING) + 5.0) / quantile_score;
			      match_arr.push_back(curScore);
			      // carp(CARP_DETAILED_DEBUG, "peptide_idx:%d \t xcorr_score:%f \t rank:%d", peptide_idx, curScore.xcorr_score, curScore.rank);
			  }
		  }
		  TideMatchSet matches(&match_arr, highest_ms2_mz);
		  if (!match_arr.empty()) { reportDIA(output_file, origin_file, *sc, active_peptide_queue, proteins, locations, &matches, &ms1scan_intensity_rank_map, &peptide_predrt_map); }

	  }

	  delete spectra;
	  delete sc_index;
	  for (map<int, pair<double*, double*>>::const_iterator i = ms1scan_intensity_rank_map.begin(); i != ms1scan_intensity_rank_map.end(); i++) {
		  delete[] (i->second).first;
		  delete[] (i->second).second;
	  }
	  /* */
  }


  return 0;
}

void DIAmeterApplication::reportDIA(
	ofstream* output_file,  ///< output file to write to
	const string& spectrum_filename, ///< name of spectrum file
	const SpectrumCollection::SpecCharge& sc, ///< spectrum and charge for matches
	const ActivePeptideQueue* peptides, ///< peptide queue
	const ProteinVec& proteins, ///< proteins corresponding with peptides
	const vector<const pb::AuxLocation*>& locations,  ///< auxiliary locations
	TideMatchSet* matches, ///< object to manage PSMs
	map<int, pair<double*, double*>>* ms1scan_intensity_rank_map,
	map<string, double>* peptide_predrt_map
) {
	Spectrum* spectrum = sc.spectrum;
	int charge = sc.charge;
	int ms1_scan_num = spectrum->MS1SpectrumNum();

	// get top-n targets and decoys by the heap
	vector<TideMatchSet::Arr::iterator> targets, decoys;
	matches->gatherTargetsAndDecoys(peptides, proteins, targets, decoys, Params::GetInt("top-match"), true);
	carp(CARP_DETAILED_DEBUG, "Gathered targets:%d \t decoy:%d", targets.size(), decoys.size());

	// calculate precursor intensity logrank
	double *intensity_arr = NULL, *intensity_rank_arr = NULL;
	map<int, pair<double*, double*>>::iterator intensityIter = ms1scan_intensity_rank_map->find(ms1_scan_num);
	if (intensityIter == ms1scan_intensity_rank_map->end()) {
		carp(CARP_DETAILED_DEBUG, "No intensity found in MS1 scan:%d !!!", ms1_scan_num);
	}
	else {
		intensity_arr = (intensityIter->second).first;
		intensity_rank_arr = (intensityIter->second).second;
	}
	map<TideMatchSet::Arr::iterator, boost::tuple<double, double, double>> intensity_map;
	computePrecIntRank(targets, peptides, intensity_rank_arr, &intensity_map, charge);
	computePrecIntRank(decoys, peptides, intensity_rank_arr, &intensity_map, charge);

	// calculate delta_cn and delta_lcn
	map<TideMatchSet::Arr::iterator, FLOAT_T> delta_cn_map;
	map<TideMatchSet::Arr::iterator, FLOAT_T> delta_lcn_map;
	TideMatchSet::computeDeltaCns(targets, &delta_cn_map, &delta_lcn_map);
	TideMatchSet::computeDeltaCns(decoys, &delta_cn_map, &delta_lcn_map);

	// calculate SpScore if necessary
	map<TideMatchSet::Arr::iterator, pair<const SpScorer::SpScoreData, int> > sp_map;
	if (Params::GetBool("compute-sp")) {
		SpScorer sp_scorer(proteins, *spectrum, charge, matches->max_mz_);
		TideMatchSet::computeSpData(targets, &sp_map, &sp_scorer, peptides);
		TideMatchSet::computeSpData(decoys, &sp_map, &sp_scorer, peptides);
	}

	matches->writeToFileDIA(output_file, Params::GetInt("top-match"), targets, spectrum_filename, spectrum, charge, peptides, proteins, locations,
			&delta_cn_map, &delta_lcn_map, Params::GetBool("compute-sp")? &sp_map : NULL, &intensity_map, peptide_predrt_map);

	matches->writeToFileDIA(output_file, Params::GetInt("top-match"), decoys, spectrum_filename, spectrum, charge, peptides, proteins, locations,
				&delta_cn_map, &delta_lcn_map, Params::GetBool("compute-sp")? &sp_map : NULL, &intensity_map, peptide_predrt_map);

}


void DIAmeterApplication::computePrecIntRank(
	const vector<TideMatchSet::Arr::iterator>& vec,
	const ActivePeptideQueue* peptides,
	const double* intensity_rank_arr,
	map<TideMatchSet::Arr::iterator, boost::tuple<double, double, double>>* intensity_map,
	int charge
) {
	for (vector<TideMatchSet::Arr::iterator>::const_iterator i = vec.begin(); i != vec.end(); ++i) {
		const Peptide& peptide = *(peptides->GetPeptide((*i)->rank));
		double peptide_mz_m0 = peptide.Mass() / (charge * 1.0) + MASS_PROTON;
		unsigned int peptide_mzbin_m0 = MassConstants::mz2bin(peptide_mz_m0);
		unsigned int peptide_mzbin_m1 = MassConstants::mz2bin(peptide_mz_m0 + 1.0/(charge * 1.0));
		unsigned int peptide_mzbin_m2 = MassConstants::mz2bin(peptide_mz_m0 + 2.0/(charge * 1.0));

		double intensity_rank_m0 = avg_noise_intensity_logrank_;
		double intensity_rank_m1 = avg_noise_intensity_logrank_;
		double intensity_rank_m2 = avg_noise_intensity_logrank_;

		if (intensity_rank_arr != NULL) {
			intensity_rank_m0 = intensity_rank_arr[peptide_mzbin_m0];
			intensity_rank_m1 = intensity_rank_arr[peptide_mzbin_m1];
			intensity_rank_m2 = intensity_rank_arr[peptide_mzbin_m2];
		}
		carp(CARP_DETAILED_DEBUG, "Peptide: %s \t mass:%f \t mz:%f \t intensity_rank:%f,%f,%f \t rank:%d \t xcorr:%f", peptide.Seq().c_str(), peptide.Mass(), peptide_mz_m0, intensity_rank_m0, intensity_rank_m1, intensity_rank_m2, (*i)->rank, (*i)->xcorr_score );
		intensity_map->insert(make_pair((*i), boost::make_tuple(intensity_rank_m0, intensity_rank_m1, intensity_rank_m2)));
	}

}


vector<InputFile> DIAmeterApplication::getInputFiles(const vector<string>& filepaths, int ms_level) const {
  vector<InputFile> input_sr;

  if (Params::GetString("spectrum-parser") != "pwiz") { carp(CARP_FATAL, "spectrum-parser must be pwiz instead of %s", Params::GetString("spectrum-parser").c_str() ); }

  for (vector<string>::const_iterator f = filepaths.begin(); f != filepaths.end(); f++) {
	  string spectrum_input_url = *f;
	  string spectrumrecords_url = make_file_path(FileUtils::BaseName(spectrum_input_url) + ".spectrumrecords.ms" + to_string(ms_level));
	  carp(CARP_INFO, "Converting %s to spectrumrecords %s", spectrum_input_url.c_str(), spectrumrecords_url.c_str());
      carp(CARP_DEBUG, "New MS%d spectrumrecords filename: %s", ms_level, spectrumrecords_url.c_str());

      if (!FileUtils::Exists(spectrumrecords_url)) {
    	  if (!SpectrumRecordWriter::convert(spectrum_input_url, spectrumrecords_url, ms_level, true)) {
    		  carp(CARP_FATAL, "Error converting MS2 spectrumrecords from %s", spectrumrecords_url.c_str());
    	  }
      }

      /*SpectrumCollection spectra;
	  pb::Header spectrum_header;
      if (!spectra.ReadSpectrumRecords(spectrumrecords_url, &spectrum_header)) {
    	  FileUtils::Remove(spectrumrecords_url);
    	  carp(CARP_FATAL, "Error reading spectrumrecords: %s", spectrumrecords_url.c_str());
    	  continue;
      }*/
      input_sr.push_back(InputFile(*f, spectrumrecords_url, true));
  }

  return input_sr;
}


void DIAmeterApplication::getPeptidePredRTMapping(map<string, double>* peptide_predrt_map) {
	carp(CARP_INFO, "predrt-files: %s ", Params::GetString("predrt-files").c_str());

	double max_predrt=-10000.0, min_predrt=10000.0;
	map<string, double> tmp_map;

	// it's possible that multiple mapping files are provided and concatenated by comma
	vector<string> mapping_paths = StringUtils::Split(Params::GetString("predrt-files"), ",");
	for(int file_idx = 0; file_idx<mapping_paths.size(); file_idx++) {
		if (!FileUtils::Exists(mapping_paths.at(file_idx))) { carp(CARP_FATAL, "The mapping file %s does not exist! \n", mapping_paths.at(file_idx).c_str()); }
		else { carp(CARP_DEBUG, "parsing the mapping file: %s", mapping_paths.at(file_idx).c_str()); }

		std::ifstream file_stream(mapping_paths.at(file_idx).c_str());
		string next_data_string;
		if (file_stream.is_open()) {
			unsigned int line_cnt = 0;
			while (getline(file_stream, next_data_string)) {
				vector<string> column_values = StringUtils::Split(StringUtils::Trim(next_data_string), "\t");
				if (column_values.size() < 2) { carp(CARP_FATAL, "Each row should contains two columns! (observed %d) \n", column_values.size()); }

				line_cnt++;
				// check if the first row is the header or the real mapping
				if (line_cnt <= 1 && !StringUtils::IsNumeric(column_values.at(1), true, true)) { continue; }

				double predrt = stod(column_values.at(1));
				tmp_map.insert(make_pair(column_values.at(0), predrt));
				// carp(CARP_DETAILED_DEBUG, "Peptide:%s \t predrt:%f", column_values.at(0).c_str(), predrt );

				if (predrt > max_predrt) { max_predrt=predrt; }
				if (predrt < min_predrt) { min_predrt=predrt; }
			}
			file_stream.close();
		}
		carp(CARP_DETAILED_DEBUG, "min_predrt:%f \t max_predrt:%f", min_predrt, max_predrt );

		for (map<string, double>::iterator it = tmp_map.begin(); it != tmp_map.end(); it++) {
			peptide_predrt_map->insert(make_pair(it->first, (it->second - min_predrt)/(max_predrt - min_predrt) ));
		}
		/*for (map<string, double>::iterator it = peptide_predrt_map->begin(); it != peptide_predrt_map->end(); it++) {
			carp(CARP_DETAILED_DEBUG, "Peptide:%s \t predrt:%f", it->first.c_str(), it->second );
		}*/
	}
	carp(CARP_DETAILED_DEBUG, "peptide_predrt_map size:%d", peptide_predrt_map->size());
}


double DIAmeterApplication::loadMS1Spectra(const std::string& file, map<int, pair<double*, double*>>* ms1scan_intensity_rank_map) {
	SpectrumCollection* spectra = loadSpectra(file);
	double highest_mz = spectra->FindHighestMZ();
	unsigned int highest_mzbin = MassConstants::mz2bin(highest_mz+1.0);

	double accumulated_intensity_logrank = 0.0;
	const vector<SpectrumCollection::SpecCharge>* spec_charges = spectra->SpecCharges();
	for (vector<SpectrumCollection::SpecCharge>::const_iterator sc = spec_charges->begin();sc < spec_charges->begin() + (spec_charges->size()); sc++) {
		Spectrum* spectrum = sc->spectrum;
		int ms1_scan_num = spectrum->MS1SpectrumNum();
		int peak_num = spectrum->Size();
		double noise_intensity_logrank = log(1.0+peak_num);
		carp(CARP_DETAILED_DEBUG, "MS1Scan:%d \t peak_size:%d \t noise_intensity_logrank:%f", ms1_scan_num, peak_num, noise_intensity_logrank);

		vector<double> sorted_intensity_vec = spectrum->DescendingSortedPeakIntensity();
		double* intensity_arr = new double[highest_mzbin];
		double* intensity_rank_arr = new double[highest_mzbin];

		fill_n(intensity_arr, highest_mzbin, 0);
		fill_n(intensity_rank_arr, highest_mzbin, noise_intensity_logrank);
		accumulated_intensity_logrank += noise_intensity_logrank;

		for (int peak_idx=0; peak_idx<peak_num; ++peak_idx) {
			double peak_mz = spectrum->M_Z(peak_idx);
			double peak_intensity = spectrum->Intensity(peak_idx);
			unsigned int peak_mzbin = MassConstants::mz2bin(peak_mz);
			double peak_intensity_logrank = log(1.0+std::count_if(sorted_intensity_vec.begin(), sorted_intensity_vec.end(),[&](int val){ return val >= peak_intensity; }));
			// carp(CARP_DETAILED_DEBUG, "peak_idx:%d \t peak_mz:%f \t peak_mzbin:%d \t peak_intensity:%f \t peak_intensity_logrank:%f", peak_idx, peak_mz, peak_mzbin, peak_intensity, peak_intensity_logrank);

			intensity_arr[peak_mzbin] = max(intensity_arr[peak_mzbin], peak_intensity);
			intensity_rank_arr[peak_mzbin] = min(intensity_rank_arr[peak_mzbin], peak_intensity_logrank);

			// analogous to XCorr by filling the flanking bin with half intensity
			if (Params::GetBool("use-flanking-peaks")) {
				double flanking_intensity_logrank = log(1.0+std::count_if(sorted_intensity_vec.begin(), sorted_intensity_vec.end(),[&](int val){ return val >= (0.5*peak_intensity); }));

				intensity_arr[peak_mzbin-1] = max(intensity_arr[peak_mzbin-1], 0.5*peak_intensity);
				intensity_arr[peak_mzbin+1] = max(intensity_arr[peak_mzbin+1], 0.5*peak_intensity);
				intensity_rank_arr[peak_mzbin-1] = min(intensity_rank_arr[peak_mzbin-1], flanking_intensity_logrank);
				intensity_rank_arr[peak_mzbin+1] = min(intensity_rank_arr[peak_mzbin+1], flanking_intensity_logrank);
			}
		}
		(*ms1scan_intensity_rank_map)[ms1_scan_num] = make_pair(intensity_arr, intensity_rank_arr);
	}
	delete spectra;

	// calculate the average noise intensity logrank, which is used as default value when MS1 scan is empty.
	return accumulated_intensity_logrank / max(1.0, 1.0*spec_charges->size());
}

SpectrumCollection* DIAmeterApplication::loadSpectra(const std::string& file) {
	SpectrumCollection* spectra = new SpectrumCollection();
	pb::Header spectrum_header;

	if (!spectra->ReadSpectrumRecords(file, &spectrum_header)) {
	    carp(CARP_FATAL, "Error reading spectrum file %s", file.c_str());
	}
	// if (string_to_window_type(Params::GetString("precursor-window-type")) != WINDOW_MZ) {
	//	carp(CARP_FATAL, "Precursor-window-type must be mz in DIAmeter!");
	// }
	// spectra->Sort();
	// spectra->Sort<ScSortByMz>(ScSortByMz(Params::GetDouble("precursor-window")));

	// Precursor-window-type must be mz in DIAmeter, based upon which spectra are sorted
	// Precursor-window is the half size of isolation window
	spectra->Sort<ScSortByMzDIA>(ScSortByMzDIA());
	spectra->SetNormalizedObvRTime();
	return spectra;
}

void DIAmeterApplication::computeWindowDIA(
  const SpectrumCollection::SpecCharge& sc,
  int max_charge,
  vector<int>* negative_isotope_errors,
  vector<double>* out_min,
  vector<double>* out_max,
  double* min_range,
  double* max_range
) {
	double unit_dalton = BIN_WIDTH;
	double mz_minus_proton = sc.spectrum->PrecursorMZ() - MASS_PROTON;
	double precursor_window = fabs(sc.spectrum->IsoWindowUpperMZ()-sc.spectrum->IsoWindowLowerMZ()) / 2;

	for (vector<int>::const_iterator ie = negative_isotope_errors->begin(); ie != negative_isotope_errors->end(); ++ie) {
		out_min->push_back((mz_minus_proton - precursor_window) * sc.charge + (*ie * unit_dalton));
		out_max->push_back((mz_minus_proton + precursor_window) * sc.charge + (*ie * unit_dalton));
	}
	*min_range = (mz_minus_proton*sc.charge + (negative_isotope_errors->front() * unit_dalton)) - precursor_window*max_charge;
	*max_range = (mz_minus_proton*sc.charge + (negative_isotope_errors->back() * unit_dalton)) + precursor_window*max_charge;
	// carp(CARP_DETAILED_DEBUG, "Scan=%d Charge=%d Mass window=[%f, %f]", sc.spectrum->SpectrumNumber(), sc.charge, (*out_min)[0], (*out_max)[0]);
}

double DIAmeterApplication::getTailorQuantile(TideMatchSet::Arr2* match_arr2) {
	//Implementation of the Tailor score calibration method
	double quantile_score = 1.0;
	vector<double> scores;
	double quantile_th = 0.01;
	// Collect the scores for the score tail distribution
	for (TideMatchSet::Arr2::iterator it = match_arr2->begin(); it != match_arr2->end(); ++it) {
		scores.push_back((double)(it->first / XCORR_SCALING));
	}
	sort(scores.begin(), scores.end(), greater<double>());  //sort in decreasing order
	int quantile_pos = (int)(quantile_th*(double)scores.size()+0.5);

	if (quantile_pos < 3) { quantile_pos = 3; }
	quantile_score = scores[quantile_pos]+5.0; // Make sure scores positive

	return quantile_score;
}


string DIAmeterApplication::getName() const {
  return "diameter";
}

string DIAmeterApplication::getDescription() const {
  return "DIAmeter description";
}

vector<string> DIAmeterApplication::getArgs() const {
  string arr[] = {
    "tide spectra file+",
    "tide database"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector<string> DIAmeterApplication::getOptions() const {
  string arr[] = {
    // "file-column",
    // "fileroot",
    "max-precursor-charge",
    "min-peaks",
    "mod-precision",
    "mz-bin-offset",
    "mz-bin-width",
    "output-dir",
    "overwrite",
    "fragment-tolerance",
    "precursor-window",
	// "precursor-window-type",
	// "spectrum-charge",
    // "spectrum-parser",
	"predrt-files",
    "top-match",
    "use-flanking-peaks",
    "use-neutral-loss-peaks",
    "verbosity"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector< pair<string, string> > DIAmeterApplication::getOutputs() const {
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
bool DIAmeterApplication::needsOutputDirectory() const {
  return true;
}

COMMAND_T DIAmeterApplication::getCommand() const {
  return DIAMETER_COMMAND;
}

void DIAmeterApplication::processParams() {

}

