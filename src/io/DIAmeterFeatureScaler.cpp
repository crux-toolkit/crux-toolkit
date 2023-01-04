
#include "DIAmeterFeatureScaler.h"
#include "DelimitedFile.h"
#include "carp.h"
#include "util/MathUtil.h"

using namespace std;

DIAmeterFeatureScaler::DIAmeterFeatureScaler(const char* file_name) {
  fileReader_ = new DelimitedFileReader(file_name, true);
  parseHeader();
}

DIAmeterFeatureScaler::~DIAmeterFeatureScaler() {
  delete fileReader_;
  for(std::vector<int>::iterator it = toscale_column_indices_.begin(); it != toscale_column_indices_.end(); ++it) {
    delete toscale_match_values_[*it];
  }
}

void DIAmeterFeatureScaler::writeScaledFile(const char* output_file_name) {
  fileReader_->reset();
  std::vector<std::string> output_vec;
  std::vector<std::string> column_names = fileReader_->getColumnNames();

  ofstream* output_file = create_stream_in_path(output_file_name, NULL, Params::GetBool("overwrite"));
  *output_file << StringUtils::Join(column_names, '\t').c_str() << endl;

  int idx = 0;
  while (fileReader_->hasNext()) {
    output_vec.clear();
    for (idx = 0; idx < column_names.size(); idx++) { output_vec.push_back(fileReader_->getString(idx)); }

    for (idx = 0; idx < toscale_column_indices_.size(); idx++) {
      int curr_column_idx = toscale_column_indices_.at(idx);
      double quantile_low_score = toscale_column_quantiles_.at(idx).first;
      double quantile_high_score = toscale_column_quantiles_.at(idx).second;
      double denominator = quantile_high_score - quantile_low_score;

      double old_score = 0.0;
      if (StringUtils::ToLower(output_vec.at(curr_column_idx)) != "nan") {
      	old_score = StringUtils::FromString<double>(output_vec.at(curr_column_idx));
      }
      double new_score = old_score;
      if (!MathUtil::AlmostEqual(denominator, 0.0, 4)) {
        new_score = (old_score - quantile_low_score) / denominator;
      }
      output_vec[curr_column_idx] = StringUtils::ToString<double>(new_score, 6);
    }

    *output_file << StringUtils::Join(output_vec, '\t').c_str() << endl;
    fileReader_->next();
  }

  if (output_file) { output_file->close(); delete output_file; }
}

void DIAmeterFeatureScaler::calcDataQuantile(double quantile_low, double quantile_high) {
  int record_cnt = 0;
  fileReader_->reset();
  while (fileReader_->hasNext()) {
    record_cnt ++;

    for (int idx = 0; idx < toscale_column_indices_.size(); idx++) {
      int curr_column_idx = toscale_column_indices_.at(idx);
      double column_val = fileReader_->getDouble(curr_column_idx);
      toscale_match_values_[curr_column_idx]->push_back(column_val);
    }
    fileReader_->next();
  }
  toscale_column_quantiles_.clear();

  carp(CARP_DETAILED_DEBUG, "Record:%d ", record_cnt );
  for (int idx = 0; idx < toscale_column_indices_.size(); idx++) {
    int curr_column_idx = toscale_column_indices_.at(idx);
    std::vector<double>* curr_match_values = toscale_match_values_[curr_column_idx];
    sort(curr_match_values->begin(), curr_match_values->end());

    if (curr_match_values->size() <= 0) { toscale_column_quantiles_.push_back(make_pair(0.0, 1.0)); }
    else {
        int quantile_low_pos = (int)(quantile_low*(double)curr_match_values->size()+0.5); // +0.5 is for rounding purpose
        int quantile_high_pos = (int)(quantile_high*(double)curr_match_values->size()+0.5);
        double quantile_low_score = curr_match_values->at(quantile_low_pos);
        double quantile_high_score = curr_match_values->at(quantile_high_pos);
        toscale_column_quantiles_.push_back(make_pair(quantile_low_score, quantile_high_score));
        carp(CARP_DETAILED_DEBUG, "ColumnIndex:%d \t ColumnName:%s \t Size:%d \t quantile_low:%f \t quantile_high:%f", curr_column_idx, get_column_header(toscale_column_ids_.at(idx)), toscale_match_values_[curr_column_idx]->size(), quantile_low_score, quantile_high_score );
    }

  }
}

void DIAmeterFeatureScaler::parseHeader() {
  for (int idx = 0; idx < NUMBER_MATCH_COLUMNS; idx++) {
    match_indices_[idx] = fileReader_->findColumn(get_column_header(idx));
    toscale_match_values_[idx] = NULL;
  }
  carp(CARP_DETAILED_DEBUG, "ColumnNames:%s", StringUtils::Join(fileReader_->getColumnNames(), ',').c_str() );

  toscale_column_ids_.clear();
  toscale_column_indices_.clear();

  MATCH_COLUMNS_T toscale_columns_[] = { PRECURSOR_INTENSITY_RANK_M0_COL, PRECURSOR_INTENSITY_RANK_M1_COL, PRECURSOR_INTENSITY_RANK_M2_COL, DYN_FRAGMENT_PVALUE_COL, STA_FRAGMENT_PVALUE_COL };
  for (int idx = 0; idx < sizeof(toscale_columns_)/sizeof(toscale_columns_[0]); idx++) {
    MATCH_COLUMNS_T curr_column_id = toscale_columns_[idx];
    const char* curr_column_name = get_column_header(curr_column_id);
    int curr_column_idx = fileReader_->findColumn(curr_column_name);
    carp(CARP_DETAILED_DEBUG, "ColumnID:%d \t ColumnIndex:%d \t ColumnName:%s", curr_column_id, curr_column_idx, curr_column_name );
    if (curr_column_idx >= 0) {
      toscale_column_ids_.push_back(curr_column_id);
      toscale_column_indices_.push_back(curr_column_idx);
      toscale_match_values_[curr_column_idx] = new std::vector<double>();
    }
  }
}

std::vector<bool> DIAmeterFeatureScaler::getMatchColumnsPresent() {
  std::vector<bool> col_is_present;
  col_is_present.assign(NUMBER_MATCH_COLUMNS, false);
  for(int col_idx = 0; col_idx < NUMBER_MATCH_COLUMNS; col_idx++) {
    col_is_present[col_idx] = (match_indices_[col_idx] > -1);
  }
  return col_is_present;
}
