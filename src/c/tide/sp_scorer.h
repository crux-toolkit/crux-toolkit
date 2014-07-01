// Tahmina Baker
//
// This file contains classes for calculating the sp score.
//
// The sp score is calculated as follows:
// sp = (intensity_sum*ion_match) * (1+(repeat_count*beta))/num_ions
//
// intensity_sum:
// Convert each ion mz into a bin and use as index into an array of 
// intensities for the spectrum which has been preprocessed for sp 
// calculations (see SpSpectrum for details). This will give us a
// single_intensity = spectrum_intensity_array[ion.bin]. We then sum up
// these intensities to get an intensity_sum.
//
// ion_match:
// The total count of ions satisfying the condition (single_intensity > 0)
// 
// repeat_count:
// Counts the links between CONSECUTIVE ions matching the spectrum. For example, 
// let's say the following y ion series matched with the spectrum at ion 
// charge 1: y1y2y3y6y12y13. We have consecutive ion matches between y1y2, y2y3
// and y12y13. The repeat_count then equals 3.
// 
// beta:
// Jimmy said this was just some number they came up with. It is simply 0.075.
// He doesn't know why.
// 
// num_ions:
// The total number of ions we tried to match with the spectrum.


#ifndef SP_SCORER_H
#define SP_SCORER_H

#include <vector>
#include "peptide.h"
#include "crux_sp_spectrum.h"

typedef vector<const pb::Protein*> ProteinVec;
typedef vector<const pb::AuxLocation*> AuxLocVec;


class SpScorer {
 public:
  struct SpScoreData {
    double sp_score;
    int sp_rank;
    int matched_ions;
    int total_ions;
    double intensity_sum;
    int repeat_count;
    
    SpScoreData() : sp_score(0.0), sp_rank(0), matched_ions(0), total_ions(0),
      intensity_sum(0.0), repeat_count(0) { }

    SpScoreData(double sp_score_in, int sp_rank_in, 
                int matched_ions_in, int total_ions_in, 
                double intensity_sum_in, int repeat_count_in)
      : sp_score(sp_score_in), sp_rank(sp_rank_in), 
      matched_ions(matched_ions_in), total_ions(total_ions_in), 
      intensity_sum(intensity_sum_in), repeat_count(repeat_count_in) { }

    void CalculateSpScore(double beta) {
      if (total_ions != 0) {
        sp_score = (intensity_sum * matched_ions) * 
                   (1 + (repeat_count * beta)) / total_ions;
      }
    }
  };
  
  SpScorer(const ProteinVec& proteins, const Spectrum& spectrum, 
           int charge, double max_mz);

  void Score(const pb::Peptide& pb_peptide, SpScoreData& sp_score_data);
  void RankSpScores(vector<SpScoreData>& scores, 
                    double* smallest_score = NULL);
  double TotalIonIntensity() {return sp_spectrum_.TotalIonIntensity();}

 private:
  typedef pair<double, int> SpScoreMatchPair;
  
  static bool CompareBySpScore(SpScoreMatchPair sp_score_1, 
                               SpScoreMatchPair sp_score_2) {
    if (sp_score_1.first > sp_score_2.first)
      return true;
    return false;
  }

  int GetBin(double mass, int charge) {
//    double mz = (mass + (charge - 1)*MassConstants::proton)/charge;
    return MassConstants::mass2bin(mass, charge);
  }

  bool IonLookup(double mass, int charge, bool previous_ion_matched,
                 SpScoreData& sp_score_data);

  
  const ProteinVec& proteins_;
  const Spectrum& spectrum_;
  SpSpectrum sp_spectrum_;
  int charge_;
  double max_mz_;
};

#endif // SP_SCORER_H
