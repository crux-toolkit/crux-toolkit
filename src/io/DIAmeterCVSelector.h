/**
 * DIAmeterCVSelector.h
 * DATE: Sep 06, 2021
 * AUTHOR: Yang Lu
 * DESCRIPTION: Object for choosing the optimal hyperparam in DIAmeter.
 **************************************************************************/

#ifndef DIAMETERCVSELECTOR_H
#define DIAMETERCVSELECTOR_H

#include "DelimitedFileReader.h"
#include "MatchColumns.h"
#include "PSMReader.h"
#include "boost/tuple/tuple.hpp"

struct PSMRecord {
  bool is_target_;
  double score_;
  std::string peptide_;

  PSMRecord(double score, bool is_target, std::string peptide)
    : score_(score), is_target_(is_target), peptide_(peptide) { }

  bool operator<(const PSMRecord& other) const {
    return (score_ > other.score_);
  }
};

struct PSMFeatEnsemble {
  bool is_target_;
  double tailor_, precursor_, fragment_, rtdiff_, elution_;
  std::string peptide_, data_;

  PSMFeatEnsemble(double tailor, double precursor, double fragment, double rtdiff, double elution, bool is_target, std::string peptide, std::string data)
    : tailor_(tailor), precursor_(precursor), fragment_(fragment), rtdiff_(rtdiff), elution_(elution),
		  is_target_(is_target), peptide_(peptide), data_(data) { }

  PSMFeatEnsemble(const PSMFeatEnsemble& other) {
  	tailor_ = other.tailor_;
  	precursor_ = other.precursor_;
  	fragment_ = other.fragment_;
  	rtdiff_ = other.rtdiff_;
  	elution_ = other.elution_;
  	is_target_ = other.is_target_;
  	peptide_ = other.peptide_;
  	data_ = other.data_;
  }

  bool operator<(const PSMFeatEnsemble& other) const {
    return (tailor_ > other.tailor_);
  }

  double getEnsembleScore(double coeff_precursor, double coeff_frag, double coeff_rtdiff, double coeff_elution) {
  	double ensemble = tailor_;
  	ensemble += (-coeff_rtdiff * rtdiff_);
  	ensemble += (-coeff_precursor * precursor_);
  	ensemble += coeff_frag * fragment_;
  	ensemble += coeff_elution * elution_;
  	return ensemble;
  }
};

struct ScanChargePSM {
	int ms2scan_, charge_;
	std::vector<PSMFeatEnsemble> psms_;

	ScanChargePSM(int ms2scan, int charge, std::vector<PSMFeatEnsemble> psms)
	: ms2scan_(ms2scan), charge_(charge) {
		for (int idx=0; idx<psms.size(); ++idx) { psms_.push_back(psms.at(idx)); }
	}

	int getFold(int totalFold) { return (10 * ms2scan_ + charge_ + 3) % totalFold; }
};

class DIAmeterCVSelector {
  protected:
    int match_indices_[NUMBER_MATCH_COLUMNS];
    std::vector<MATCH_COLUMNS_T> toagg_column_ids_;
    int agg_idx_, scan_idx_, charge_idx_, peptide_idx_, td_idx_;
    int tailor_idx_, precursor_idx_, fragment_idx_, rtdiff_idx_, elution_idx_;
    DelimitedFileReader* fileReader_;

    std::vector<ScanChargePSM> scPSMList;

    void parseHeader();
    int getKey(int scan, int charge);

  public:
    DIAmeterCVSelector(const char* file_name);
    ~DIAmeterCVSelector();

    void loadData(const char* output_file_name);
    void FoldFilter(const char* output_file_name, std::vector<double>* paramRangeList, int totalFold=3);
    boost::tuple<double, double, double, double> selectFoldParam(std::vector<double>* paramRangeList, vector<int>* train_indices);
    int getTargetFDR(std::vector<PSMRecord>* records, double fdr_thres=0.01);
};

#endif //DIAMETERCVSELECTOR_H

