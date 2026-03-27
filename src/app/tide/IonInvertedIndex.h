#ifndef ION_INVERTED_INDEX
#define ION_INVERTED_INDEX
#include <map>
#include <deque>
#include "theoretical_peak_set.h"
#include "util/crux-utils.h"
#include "model/objects.h"
#include "util/Params.h"

class Peptide;

class IonInvertedIndex {
  public:
    // maps a peak index ix to a list of peptides which contains that peak ix as a theoretical peak.
    using IonsStorage = std::vector<std::deque<Peptide*>>;

    IonInvertedIndex(){
      // IonStorage should be a vector with a fixed length, not a map. Each vector component is a dequeue. 
      // Each vector index is associated to a peak_mz;
      capacity_ = 0;
      curScoreFunction_ = string_to_score_function_type(Params::GetString("score-function"));  

      if (curScoreFunction_ == HYPERSCORE) {
        capacity_ =  MassConstants::mass2bin(MAX_THEORETICAL_PEAK_MZ, 1)*2; 
        ions_b_.resize(capacity_);
        ions_y_.resize(capacity_);
        // ions_b2_.resize(capacity_);
        // ions_y2_.resize(capacity_);
      }
      // Create a vector with a length of capacity.
      // size_t top_matches_ = Params::GetInt("top-match");   // +1 for delta cn 

    }

    void insert_peaks(Peptide* peptide);

    void pop_peaks(Peptide* peptide);
    
    // Scoring with hyperscore
    void score_peaks(unsigned int peak_mz, double peak_int, int charge);  

    IonsStorage ions_b_;
    IonsStorage ions_y_;
    // IonsStorage ions_b2_;
    // IonsStorage ions_y2_;
    unsigned int capacity_ = 1;
    SCORE_FUNCTION_T curScoreFunction_;



};
#endif