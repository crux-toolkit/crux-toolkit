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


struct PSMFeatEnsemble {
    bool is_target_;
    double tailor_, precursor_, fragment_, rtdiff_, elution_;
    std::string peptide_, record_;

    PSMFeatEnsemble(double tailor, double precursor, double fragment, double rtdiff, double elution, bool is_target, std::string peptide, std::string record)
        : tailor_(tailor), precursor_(precursor), fragment_(fragment), rtdiff_(rtdiff), elution_(elution),
		  is_target_(is_target), peptide_(peptide), record_(record) { }

    PSMFeatEnsemble(const PSMFeatEnsemble& other) {
    	tailor_ = other.tailor_;
    	precursor_ = other.precursor_;
    	fragment_ = other.fragment_;
    	rtdiff_ = other.rtdiff_;
    	elution_ = other.elution_;
    	is_target_ = other.is_target_;
    	peptide_ = other.peptide_;
    	record_ = other.record_;
    }

    bool operator<(const PSMFeatEnsemble& other) const {
        return (tailor_ > other.tailor_);
    }

    double getEnsembleScore(double coeff_precursor, double coeff_fragment, double coeff_rtdiff, double coeff_elution) {
    	double ensemble = tailor_;
    	ensemble += (-coeff_rtdiff * rtdiff_);
    	ensemble += (-coeff_precursor * precursor_);
    	ensemble += coeff_fragment * fragment_;
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

	int getFold(int totalFold) { return (10 * ms2scan_ + charge_) % totalFold; }
};


class DIAmeterCVSelector {
    protected:
        int match_indices_[NUMBER_MATCH_COLUMNS];
        std::vector<MATCH_COLUMNS_T> toagg_column_ids_;
        int agg_idx_, scan_idx_, charge_idx_, peptide_idx_, td_idx_;
        int tailor_idx_, precursor_idx_, fragment_idx_, rtdiff_idx_, elution_idx_;
        DelimitedFileReader* fileReader_;

        void parseHeader();
        int getKey(int scan, int charge);

    public:
        DIAmeterCVSelector(const char* file_name);
        ~DIAmeterCVSelector();

};




/*

class DIAmeterCVSelector {
 protected:
    int match_indices_[NUMBER_MATCH_COLUMNS];
    std::vector<MATCH_COLUMNS_T> toagg_column_ids_;
    std::vector<int> toagg_column_indices_;
    std::vector<double> toagg_column_coeffs_;
    int agg_idx_, scan_idx_, charge_idx_, xcorr_idx_;
    DelimitedFileReader* fileReader_;

    // we use (scan*10+charge) as the key
    std::map<int, boost::tuple<double, double>> scan_charge_scores_map;

    void parseHeader();

    static bool psm_sorter(const PSMByScanCharge & psm1, const PSMByScanCharge & psm2);

 public:
    DIAmeterPSMFilter(const char* file_name);
    ~DIAmeterPSMFilter();

    void calcBaseline();
    void loadAndFilter(const char* output_file_name, bool filter=true);
};*/


#endif //DIAMETERCVSELECTOR_H

