// Tahmina Baker
//
// This file contains implementations for classes defined in sp_scorer.h.
// Please see the header file for details.

#include <list>
#include "sp_scorer.h"

SpScorer::SpScorer(const ProteinVec& proteins, const Spectrum& spectrum, 
                   int charge, double max_mz)
  : proteins_(proteins), spectrum_(spectrum), charge_(charge), max_mz_(max_mz),
  sp_spectrum_(spectrum, charge, max_mz) {
}

bool SpScorer::IonLookup(double mass, int charge, bool previous_ion_matched,
                         SpScoreData& sp_score_data) {
  bool matched = false;
  sp_score_data.total_ions++;
  int bin = GetBin(mass, charge);
  double intensity = 0.0;
  if (bin < max_mz_)
    intensity = sp_spectrum_.Intensity(bin);  
  if (intensity > 0) {
    // We have a match
    matched = true;
    sp_score_data.matched_ions++;
    sp_score_data.intensity_sum += intensity;
    
    if (previous_ion_matched) {
      sp_score_data.repeat_count++;
    }
  }

  return matched;
}

void SpScorer::Score(const pb::Peptide& pb_peptide, SpScoreData& sp_score_data) {
  Peptide peptide(pb_peptide, proteins_);
  vector<double> m_z(peptide.Len());
  string sequence = peptide.Seq();

  // Collect m/z values for each residue
  for (int i = 0; i < sequence.length(); i++)
    m_z[i] = MassConstants::mono_table[sequence[i]];

  // Account for modifications
  const ModCoder::Mod* mods;
  int num_mods = peptide.Mods(&mods);
  for (int i = 0; i < num_mods; i++) {
    int index;
    double delta;
    MassConstants::DecodeMod(mods[i], &index, &delta);
    m_z[index] += delta;
  }

  int precursor_charge = (charge_ == 1) ? 2 : charge_;
  for (int ion_charge = 1; ion_charge < precursor_charge; ion_charge++) {
    // Needed for keeping track of repeat_count
    bool previous_b_ion_matched = false;
    bool previous_y_ion_matched = false;
    
    double b_ion = MassConstants::proton;
    double y_ion = peptide.Mass() + MassConstants::proton;

    for (int i = 0; i < sequence.length(); i++) {
      // Calculate and look up b-ions
      if (i < sequence.length()-1) {
        b_ion += m_z[i];
        previous_b_ion_matched = IonLookup(b_ion, ion_charge, 
                                           previous_b_ion_matched, 
                                           sp_score_data);
      }

      // Calculate and look up y-ions
      if (i > 0) {
        y_ion -= m_z[i-1];
        previous_y_ion_matched = IonLookup(y_ion, ion_charge, 
                                           previous_y_ion_matched, 
                                           sp_score_data);
      }
    }
  }

  sp_score_data.CalculateSpScore(sp_spectrum_.Beta());
}

void SpScorer::RankSpScores(vector<SpScoreData>& scores, 
                            double* smallest_score) {

  if (scores.size() > 0) {
    *smallest_score = scores[0].sp_score;
  } else {
    *smallest_score = 0.0;
  }

  // We use this list to sort the matches according to sp score, then
  // use the match id to assign the rankings to the sp scores passed in
  list<SpScoreMatchPair> sp_score_match_list;
  for (int match = 0; match < scores.size(); match++) {
    sp_score_match_list.push_back(make_pair(scores[match].sp_score, match));
    if (scores[match].sp_score < *smallest_score)
      *smallest_score = scores[match].sp_score;
  }
  sp_score_match_list.sort(CompareBySpScore);
  
  // Assign rankings to the sp scores passed in
  list<SpScoreMatchPair>::iterator i;  
  int rank_count = 0;
  for (i=sp_score_match_list.begin(); i != sp_score_match_list.end(); ++i) {
    scores[i->second].sp_rank = ++rank_count;
  }

  if (smallest_score) {
    if (sp_score_match_list.size() > 0) {
      // After we go through the loop that assigns the rankings, i is at 
      // list.end(). If we iterate back one, we'll have the lowest sp score.
      *smallest_score = scores[(--i)->second].sp_score;
    } else {
      *smallest_score = 0.0;
    }
  }
}
