#include "DIAmeterCVSelector.h"
#include "DelimitedFile.h"
#include "carp.h"

using namespace std;

DIAmeterCVSelector::DIAmeterCVSelector(const char* file_name) {
  fileReader_ = new DelimitedFileReader(file_name, true);
  parseHeader();
}

DIAmeterCVSelector::~DIAmeterCVSelector() {
  delete fileReader_;
}

int DIAmeterCVSelector::getKey(int scan, int charge) { return 10*scan+charge; }

void DIAmeterCVSelector::parseHeader() {
  for (int idx = 0; idx < NUMBER_MATCH_COLUMNS; idx++) {
    match_indices_[idx] = fileReader_->findColumn(get_column_header(idx));
  }
  carp(CARP_DETAILED_DEBUG, "ColumnNames:%s", StringUtils::Join(fileReader_->getColumnNames(), ',').c_str() );

  agg_idx_ = fileReader_->findColumn(get_column_header(ENSEMBLE_SCORE_COL));
  scan_idx_ = fileReader_->findColumn(get_column_header(SCAN_COL));
  charge_idx_ = fileReader_->findColumn(get_column_header(CHARGE_COL));
  peptide_idx_ = fileReader_->findColumn(get_column_header(SEQUENCE_COL));
  td_idx_ = fileReader_->findColumn(get_column_header(TARGET_DECOY_COL));

  tailor_idx_ = fileReader_->findColumn(get_column_header(TAILOR_COL));
  rtdiff_idx_ = fileReader_->findColumn(get_column_header(RT_DIFF_COL));
  precursor_idx_ = fileReader_->findColumn(get_column_header(PRECURSOR_INTENSITY_RANK_M0_COL));
  fragment_idx_ = fileReader_->findColumn(get_column_header(DYN_FRAGMENT_PVALUE_COL));
  elution_idx_ = fileReader_->findColumn(get_column_header(COELUTE_MS1_COL));
  carp(CARP_DETAILED_DEBUG, "agg_idx_:%d\tscan_idx:%d\tcharge_idx:%d\tpeptide_idx_:%d\ttd_idx_:%d", agg_idx_, scan_idx_, charge_idx_, peptide_idx_, td_idx_ );
  carp(CARP_DETAILED_DEBUG, "tailor_idx:%d\trtdiff_idx:%d\tprecursor_idx:%d\tfragment_idx:%d\telution_idx:%d", tailor_idx_, rtdiff_idx_, precursor_idx_, fragment_idx_, elution_idx_ );

}

void DIAmeterCVSelector::loadData(const char* output_file_name) {
  fileReader_->reset();
  scPSMList.clear();

  int curr_key = 0, curr_scan = 0, curr_charge = 0;
  vector<PSMFeatEnsemble> curr_psms;

  while (fileReader_->hasNext()) {
    int scan = fileReader_->getInteger(scan_idx_);
    int charge = fileReader_->getInteger(charge_idx_);
    int key = getKey(scan, charge);

    std::string peptide = fileReader_->getString(peptide_idx_);
    std::string td = fileReader_->getString(td_idx_);

    bool is_target = true;
    if (StringUtils::IEquals(td, "target")) { is_target = true; }
    else if (StringUtils::IEquals(td, "decoy")) { is_target = false; }
    else { carp(CARP_FATAL, "td is neither target nor decoy: %s", td.c_str()); }

    double tailor_val = fileReader_->getDouble(tailor_idx_);
    double rtdiff_val = fileReader_->getDouble(rtdiff_idx_);
    double precursor_val = fileReader_->getDouble(precursor_idx_);
    double fragment_val = fileReader_->getDouble(fragment_idx_);
    double elution_val = fileReader_->getDouble(elution_idx_);
    std::string row_val = fileReader_->getString();

    // carp(CARP_DETAILED_DEBUG, "scan:%d\t charge:%d\t peptide:%s\t td:%s\t is_target:%d", scan, charge, peptide.c_str(), td.c_str(), int(is_target) );
    // carp(CARP_DETAILED_DEBUG, "tailor:%f\t rtdiff:%f\t precursor:%f\t fragment:%f\t elution:%f", tailor_val, rtdiff_val, precursor_val, fragment_val, elution_val);

    if (curr_key != key) {
      if (curr_psms.size() > 0) {
        ScanChargePSM scPSM(curr_scan, curr_charge, curr_psms);
        scPSMList.push_back(scPSM);
      }

      curr_psms.clear();
      curr_key = key;
      curr_scan = scan;
      curr_charge = charge;
    }

    PSMFeatEnsemble featEnsemble(tailor_val, precursor_val, fragment_val, rtdiff_val, elution_val, is_target, peptide, row_val);
    curr_psms.push_back(featEnsemble);

    fileReader_->next();
  }

  if (curr_psms.size() > 0) {
    ScanChargePSM scPSM(curr_scan, curr_charge, curr_psms);
    scPSMList.push_back(scPSM);
  }
  carp(CARP_DETAILED_DEBUG, "scPSMList:%d", scPSMList.size());
}

void DIAmeterCVSelector::FoldFilter(const char* output_file_name, std::vector<double>* paramRangeList, int totalFold) {
  ofstream* output_file = create_stream_in_path(output_file_name, NULL, Params::GetBool("overwrite"));
  *output_file << StringUtils::Join(fileReader_->getColumnNames(), '\t').c_str() << endl;

  vector<int> train_indices, test_indices;
  for (int targetFold=0; targetFold<totalFold; ++targetFold) {
    train_indices.clear(); test_indices.clear();

    for (int psm_idx=0; psm_idx<scPSMList.size(); ++psm_idx) {
      ScanChargePSM scPSM = scPSMList.at(psm_idx);
      int psm_fold = scPSMList.at(psm_idx).getFold(totalFold);

      if (psm_fold == targetFold) { test_indices.push_back(psm_idx); }
      else { train_indices.push_back(psm_idx); }
    }
    carp(CARP_DETAILED_DEBUG, "targetFold:%d\t train_indices:%d\t test_indices:%d", targetFold, train_indices.size(), test_indices.size() );

    boost::tuple<double, double, double, double> opt_param = selectFoldParam(paramRangeList, &train_indices);
    double coeff_precursor = opt_param.get<0>();
    double coeff_frag = opt_param.get<1>();
    double coeff_rtdiff = opt_param.get<2>();
    double coeff_elution = opt_param.get<3>();

    bool filter = true;
    if (coeff_precursor < 0 && coeff_frag < 0 && coeff_rtdiff < 0 && coeff_elution < 0) { filter=false; }

    for (int idx=0; idx<test_indices.size(); ++idx) {
      int psm_idx = test_indices.at(idx);
      ScanChargePSM scPSM = scPSMList.at(psm_idx);
      vector<PSMFeatEnsemble> psms = scPSM.psms_;
      if (psms.size() <= 0) { continue; }

      double tailor_baseline = psms.at(0).tailor_;
      double ensemble_baseline = psms.at(0).getEnsembleScore(coeff_precursor, coeff_frag, coeff_rtdiff, coeff_elution) - 0.000001;

      for (int idx2=0; idx2<psms.size(); ++idx2) {
        if (psms.at(idx2).tailor_ > tailor_baseline) { carp(CARP_FATAL, "tailor %f shouldn't beat baseline %f!", psms.at(idx2).tailor_, tailor_baseline); }
        double ensemble = psms.at(idx2).getEnsembleScore(coeff_precursor, coeff_frag, coeff_rtdiff, coeff_elution);
        if ((!filter) || (ensemble >= ensemble_baseline)) { *output_file << psms.at(idx2).data_.c_str() << endl; }
      }
    }
  }

  if (output_file) { output_file->close(); delete output_file; }
}

boost::tuple<double, double, double, double> DIAmeterCVSelector::selectFoldParam(std::vector<double>* paramRangeList, vector<int>* train_indices) {
  vector<boost::tuple<double, double, double, double>> param_combos;
  for (int prec_idx=0; prec_idx<paramRangeList->size(); ++prec_idx) {
    double coeff_precursor = paramRangeList->at(prec_idx);
    for (int frag_idx=0; frag_idx<paramRangeList->size(); ++frag_idx) {
      double coeff_frag = paramRangeList->at(frag_idx);
      for (int rt_idx=0; rt_idx<paramRangeList->size(); ++rt_idx) {
        double coeff_rtdiff = paramRangeList->at(rt_idx);
        for (int elut_idx=0; elut_idx<paramRangeList->size(); ++elut_idx) {
          double coeff_elution = paramRangeList->at(elut_idx);
          param_combos.push_back(boost::make_tuple(coeff_precursor, coeff_frag, coeff_rtdiff, coeff_elution));
        }
      }
    }
  }
  carp(CARP_DETAILED_DEBUG, "param_combos:%d", param_combos.size() );

  int max_targetCnt = 0, max_paramIdx = 0;
  vector<PSMRecord> filtered_records;
  for (int param_idx=0; param_idx<param_combos.size(); ++param_idx) {
    boost::tuple<double, double, double, double> curr_param = param_combos.at(param_idx);
    double coeff_precursor = curr_param.get<0>();
    double coeff_frag = curr_param.get<1>();
    double coeff_rtdiff = curr_param.get<2>();
    double coeff_elution = curr_param.get<3>();
    carp(CARP_DETAILED_DEBUG, "coeff_precursor:%f\t coeff_frag:%f\t coeff_rtdiff:%f\t coeff_elution:%f", coeff_precursor, coeff_frag, coeff_rtdiff, coeff_elution);

    filtered_records.clear();
    for (int idx=0; idx<train_indices->size(); ++idx) {
      int psm_idx = train_indices->at(idx);
      ScanChargePSM scPSM = scPSMList.at(psm_idx);
      // carp(CARP_DETAILED_DEBUG, "psm_idx:%d\t ms2scan:%d\t charge:%d", psm_idx, scPSM.ms2scan_, scPSM.charge_);

      vector<PSMFeatEnsemble> psms = scPSM.psms_;
      if (psms.size() <= 0) { continue; }

      double tailor_baseline = psms.at(0).tailor_;
      double ensemble_baseline = psms.at(0).getEnsembleScore(coeff_precursor, coeff_frag, coeff_rtdiff, coeff_elution) - 0.000001;

      for (int idx2=0; idx2<psms.size(); ++idx2) {
        if (psms.at(idx2).tailor_ > tailor_baseline) { carp(CARP_FATAL, "tailor %f shouldn't beat baseline %f!", psms.at(idx2).tailor_, tailor_baseline); }
        double ensemble = psms.at(idx2).getEnsembleScore(coeff_precursor, coeff_frag, coeff_rtdiff, coeff_elution);
        if (ensemble >= ensemble_baseline) {
          PSMRecord record(ensemble, psms.at(idx2).is_target_, psms.at(idx2).peptide_);
          filtered_records.push_back(record);
        }
      }
    }
    sort(filtered_records.begin(), filtered_records.end());
    int target_cnt = getTargetFDR(&filtered_records);
    // carp(CARP_DETAILED_DEBUG, "filtered_records:%d\t target_cnt:%d", filtered_records.size(), target_cnt );

    if (target_cnt > max_targetCnt) {
      max_targetCnt = target_cnt;
      max_paramIdx = param_idx;
    }
  }
  carp(CARP_DETAILED_DEBUG, "max_targetCnt:%d\t max_paramIdx:%d", max_targetCnt, max_paramIdx );

  // the case when no psm filtering
  filtered_records.clear();
  for (int idx=0; idx<train_indices->size(); ++idx) {
    int psm_idx = train_indices->at(idx);
    ScanChargePSM scPSM = scPSMList.at(psm_idx);
    vector<PSMFeatEnsemble> psms = scPSM.psms_;
    for (int idx2=0; idx2<psms.size(); ++idx2) {
      PSMRecord record(psms.at(idx2).tailor_, psms.at(idx2).is_target_, psms.at(idx2).peptide_);
      filtered_records.push_back(record);
    }
  }
  sort(filtered_records.begin(), filtered_records.end());
  int unfiltered_target_cnt = getTargetFDR(&filtered_records);
  carp(CARP_DETAILED_DEBUG, "unfiltered_target_cnt:%d", unfiltered_target_cnt);

  if (max_targetCnt > unfiltered_target_cnt) { return param_combos.at(max_paramIdx); }
  else { return boost::make_tuple(-1, -1, -1, -1); }

}

int DIAmeterCVSelector::getTargetFDR(std::vector<PSMRecord>* records, double fdr_thres) {
  map<std::string, int> peptide_map;
  double prev_score = 1000000, target_cnt = 0, decoy_cnt = 0;

  for (int idx=0; idx<records->size(); ++idx) {
    PSMRecord record = records->at(idx);

    double score = record.score_;
    bool is_target = record.is_target_;
    std::string peptide = record.peptide_;

    if (score > prev_score) { carp(CARP_FATAL, "score %f should be decreasing! %f!", score, prev_score); }
    prev_score = score;

    map<std::string, int>::iterator pepIter = peptide_map.find(peptide);
    if (pepIter != peptide_map.end()) { continue; }
    peptide_map[peptide] = 1;

    if (is_target) {
      target_cnt++;
      double fdr = decoy_cnt * 1.0 / target_cnt;
      if (fdr > fdr_thres) { break; }
    }
    else { decoy_cnt++; }
  }

  return int(target_cnt);
}


