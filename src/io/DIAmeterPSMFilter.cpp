
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

int DIAmeterPSMFilter::getKey(int scan, int charge) { return 10*scan+charge; }

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

  MATCH_COLUMNS_T toagg_columns_[] = {TAILOR_COL, RT_DIFF_COL, PRECURSOR_INTENSITY_RANK_M0_COL, DYN_FRAGMENT_PVALUE_COL, COELUTE_MS1_COL };
  double toagg_coeffs_[] = {1.0, -coeff_rtdiff, -coeff_precursor, coeff_fragment, coeff_elution  };

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
  carp(CARP_DETAILED_DEBUG, "ensemble_idx:%d \t scan_idx:%d \t charge_idx:%d \t xcorr_idx:%d", agg_idx_, scan_idx_, charge_idx_, xcorr_idx_ );

}

void DIAmeterPSMFilter::calcBaseline() {
  fileReader_->reset();
  scan_charge_scores_map.clear();

  while (fileReader_->hasNext()) {
    int scan = fileReader_->getInteger(scan_idx_);
    int charge = fileReader_->getInteger(charge_idx_);
    int key = getKey(scan, charge);

    std::vector<std::string> data = fileReader_->getCurrentRowData();

    double xcorr = fileReader_->getDouble(xcorr_idx_);
    double ensemble = 0.0;
    for (int idx = 0; idx < toagg_column_indices_.size(); idx++) {
      int curr_column_idx = toagg_column_indices_.at(idx);
      double column_val = fileReader_->getDouble(curr_column_idx);
      double curr_column_coeff = toagg_column_coeffs_.at(idx);
      ensemble += column_val * curr_column_coeff;
    }

    map<int, boost::tuple<double, double>>::iterator baselineIter = scan_charge_scores_map.find(key);
    if (baselineIter == scan_charge_scores_map.end()) { scan_charge_scores_map[key] = boost::make_tuple(xcorr, ensemble); }
    else {
      double xcorr_old = (baselineIter->second).get<0>();
      if (xcorr_old < xcorr) {
        scan_charge_scores_map[key] = boost::make_tuple(xcorr, ensemble);
      }
    }

    fileReader_->next();
  }
}

void DIAmeterPSMFilter::loadAndFilter(const char* output_file_name, bool filter) {
  calcBaseline();
  fileReader_->reset();

  ofstream* output_file = create_stream_in_path(output_file_name, NULL, Params::GetBool("overwrite"));
  *output_file << StringUtils::Join(fileReader_->getColumnNames(), '\t').c_str() << endl;

  while (fileReader_->hasNext()) {
    int scan = fileReader_->getInteger(scan_idx_);
    int charge = fileReader_->getInteger(charge_idx_);
    int key = getKey(scan, charge);

    std::vector<std::string> data = fileReader_->getCurrentRowData();

    double xcorr = fileReader_->getDouble(xcorr_idx_);
    double ensemble = 0.0;
    for (int idx = 0; idx < toagg_column_indices_.size(); idx++) {
      int curr_column_idx = toagg_column_indices_.at(idx);
      double column_val = fileReader_->getDouble(curr_column_idx);
      double curr_column_coeff = toagg_column_coeffs_.at(idx);
      ensemble += column_val * curr_column_coeff;
    }
    data[agg_idx_] = StringUtils::ToString<double>(ensemble, 6);

    map<int, boost::tuple<double, double>>::iterator baselineIter = scan_charge_scores_map.find(key);
    if (baselineIter == scan_charge_scores_map.end()) { carp(CARP_FATAL, "The key must exist in scan_charge_scores_map! %d", key); }

    double xcorr_baseline = (baselineIter->second).get<0>();
    double ensemble_baseline = (baselineIter->second).get<1>() - 0.000001;
    if ((!filter) || (ensemble >= ensemble_baseline)) {
      *output_file << StringUtils::Join(data, '\t').c_str() << endl;
    }

    fileReader_->next();
  }

  if (output_file) { output_file->close(); delete output_file; }
}

