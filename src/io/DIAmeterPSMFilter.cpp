
//#include "MatchColumns.h"
#include "DIAmeterPSMFilter.h"
#include "DelimitedFile.h"
#include "carp.h"

using namespace std;

DIAmeterPSMFilter::DIAmeterPSMFilter(const char* file_name) {
	fileReader_ = new DelimitedFileReader(file_name, true);
    parseHeader();
}

DIAmeterPSMFilter::~DIAmeterPSMFilter() {
	delete fileReader_;
}

bool DIAmeterPSMFilter::psm_sorter(const PSMByScanCharge & psm1, const PSMByScanCharge & psm2) {
	return psm1.xcorr_score_ > psm2.xcorr_score_;
}

void DIAmeterPSMFilter::loadAndFilter(const char* output_file_name, bool filter) {
	 fileReader_->reset();
	 std::vector<PSMByScanCharge> pooled_data_vec;

	 ofstream* output_file = create_stream_in_path(output_file_name, NULL, Params::GetBool("overwrite"));
	 *output_file << StringUtils::Join(fileReader_->getColumnNames(), '\t').c_str() << endl;

	 int prev_scan = -1, prev_charge = -1;
	 while (fileReader_->hasNext()) {
		 int scan = fileReader_->getInteger(scan_idx_);
		 int charge = fileReader_->getInteger(charge_idx_);
		 std::vector<std::string> data = fileReader_->getCurrentRowData();

		 double xcorr = fileReader_->getDouble(xcorr_idx_);

		 double ensemble = 0.0;
		 for (int idx = 0; idx < toagg_column_indices_.size(); idx++) {
			 int curr_column_idx = toagg_column_indices_.at(idx);
			 double column_val = fileReader_->getDouble(curr_column_idx);
			 double curr_column_coeff = toagg_column_coeffs_.at(idx);
			 ensemble += column_val * curr_column_coeff;
		 }
		 // carp(CARP_DEBUG, "scan:%d \t charge:%d \t xcorr:%f \t ensemble:%f \t data_str:%s", scan, charge, xcorr, ensemble, StringUtils::Join(data, '\t').c_str());


		 double old_ensemble = StringUtils::FromString<double>(data.at(agg_idx_));
		 data[agg_idx_] = StringUtils::ToString<double>(ensemble, 6);
		 // carp(CARP_DEBUG, "old_ensemble:%f \t new_ensemble:%s", old_ensemble, data[agg_idx_].c_str());

		 // deal with the PSMs corresponding to the same (scan, charge) tuple
		 if (prev_scan != scan || prev_charge != charge) {
			 if (scan < prev_scan) { carp(CARP_FATAL, "Error! The scan must be sorted ascendingly!"); }

			 if (pooled_data_vec.size() > 0) { processPSMGroup(&pooled_data_vec, output_file, filter); }
			 pooled_data_vec.clear();
		 }

		 pooled_data_vec.push_back(PSMByScanCharge(xcorr, ensemble, data));
		 prev_scan = scan;
		 prev_charge = charge;

		 fileReader_->next();
	 }
	 if (pooled_data_vec.size() > 0) { processPSMGroup(&pooled_data_vec, output_file, filter); }

	 if (output_file) { output_file->close(); delete output_file; }
}

void DIAmeterPSMFilter::processPSMGroup(std::vector<PSMByScanCharge>* pooled_data_vec, ofstream* output_file, bool filter) {
	// carp(CARP_DEBUG, "pooled_data_vec:%d", pooled_data_vec->size());
	sort(pooled_data_vec->begin(), pooled_data_vec->end(), psm_sorter);

	PSMByScanCharge psm_baseline = pooled_data_vec->at(0);
	double ensemble_baseline = psm_baseline.ensemble_score_ + 0.000001;
	*output_file << StringUtils::Join(psm_baseline.data_, '\t').c_str() << endl;

	for (int idx = 1; idx < pooled_data_vec->size(); idx++) {
		PSMByScanCharge psm = pooled_data_vec->at(idx);
		double ensemble = psm.ensemble_score_;
		if (filter && (ensemble >= ensemble_baseline)) {
			*output_file << StringUtils::Join(psm.data_, '\t').c_str() << endl;
			// carp(CARP_DEBUG, "idx:%d \t xcorr:%f \t ensemble:%f", idx, psm.xcorr_score_, psm.ensemble_score_);
		}
	}

}



void DIAmeterPSMFilter::parseHeader() {
  for (int idx = 0; idx < NUMBER_MATCH_COLUMNS; idx++) {
    match_indices_[idx] = fileReader_->findColumn(get_column_header(idx));
  }
  carp(CARP_DEBUG, "ColumnNames:%s", StringUtils::Join(fileReader_->getColumnNames(), ',').c_str() );

  toagg_column_ids_.clear();
  toagg_column_indices_.clear();
  toagg_column_coeffs_.clear();

  double coeff_precursor = Params::GetDouble("coeff-precursor");
  double coeff_fragment = Params::GetDouble("coeff-fragment");
  double coeff_rtdiff = Params::GetDouble("coeff-rtdiff");
  double coeff_elution = Params::GetDouble("coeff-elution");
  carp(CARP_DETAILED_DEBUG, "coeff_precursor:%f \t coeff_fragment:%f \t coeff_rtdiff:%f \t coeff_elution:%f", coeff_precursor, coeff_fragment, coeff_rtdiff, coeff_elution );

  MATCH_COLUMNS_T toagg_columns_[] = {
		  TAILOR_COL, RT_DIFF_COL,
		  PRECURSOR_INTENSITY_RANK_M0_COL, PRECURSOR_INTENSITY_RANK_M1_COL, PRECURSOR_INTENSITY_RANK_M2_COL,
		  DYN_FRAGMENT_PVALUE_COL, STA_FRAGMENT_PVALUE_COL,
		  COELUTE_MS1_COL, COELUTE_MS1_MS2_COL  };

  double toagg_coeffs_[] = {
		  1.0, -coeff_rtdiff,
		  -coeff_precursor/3.0, -coeff_precursor/3.0, -coeff_precursor/3.0,
		  coeff_fragment/2.0, coeff_fragment/2.0,
		  coeff_elution/2.0, coeff_elution/2.0  };

  for (int idx = 0; idx < sizeof(toagg_columns_)/sizeof(toagg_columns_[0]); idx++) {
	  MATCH_COLUMNS_T curr_column_id = toagg_columns_[idx];
	  const char* curr_column_name = get_column_header(curr_column_id);
	  int curr_column_idx = fileReader_->findColumn(curr_column_name);
	  double curr_column_coeff = toagg_coeffs_[idx];
	  carp(CARP_DEBUG, "ColumnID:%d \t ColumnIndex:%d \t ColumnName:%s \t ColumnCoeff:%f", curr_column_id, curr_column_idx, curr_column_name, curr_column_coeff );

	  if (curr_column_idx >= 0) {
		  toagg_column_ids_.push_back(curr_column_id);
		  toagg_column_indices_.push_back(curr_column_idx);
		  toagg_column_coeffs_.push_back(curr_column_coeff);
	  }
  }

  agg_idx_ = fileReader_->findColumn(get_column_header(ENSEMBLE_SCORE_COL));
  scan_idx_ = fileReader_->findColumn(get_column_header(SCAN_COL));
  charge_idx_ = fileReader_->findColumn(get_column_header(CHARGE_COL));
  xcorr_idx_ = fileReader_->findColumn(get_column_header(XCORR_SCORE_COL));
  carp(CARP_DEBUG, "ensemble_idx:%d \t scan_idx:%d \t charge_idx:%d \t xcorr_idx:%d", agg_idx_, scan_idx_, charge_idx_, xcorr_idx_ );

}
