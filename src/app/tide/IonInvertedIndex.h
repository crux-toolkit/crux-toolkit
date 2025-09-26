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
    using IonsStorage = std::vector<std::deque<std::pair<double, Peptide*>>>;

    IonInvertedIndex(){
      // IonStorage should be a vector with a fixed length, not a map. Each vector component is a dequeue. 
      // Each vector index is associated to a peak_mz;
       capacity_ = 0;
      SCORE_FUNCTION_T curScoreFunction = string_to_score_function_type(Params::GetString("score-function"));  
      if (curScoreFunction == HYPERSCORE) {
         capacity_ =  MassConstants::mass2bin(MAX_THEORETICAL_PEAK_MZ, 1); 
      }
      ions_b_.resize(capacity_);
      ions_y_.resize(capacity_);
      // Create a vector with a length of capacity.
    }

    void insert_peaks(Peptide* peptide);

    void pop_peaks(Peptide* peptide);
    
    // Scoring with hyperscore
    void score_peaks(double min_precursor_mass, unsigned int peak_mz, double peak_int);  

    // Iterator lowerBound(unsigned int peak);
    // Iterator upperBound(unsigned int peak);

    // ConstIterator lowerBound(unsigned int peak) const;
    // ConstIterator upperBound(unsigned int peak) const;

    // Iterator begin() { return ions_.begin(); }
    // ConstIterator begin() const { return ions_.begin(); }
    // Iterator end() { return ions_.end(); }
    // ConstIterator end() const { return ions_.end(); }
  private:
    IonsStorage ions_b_;
    IonsStorage ions_y_;
    unsigned int capacity_ = 1;

};
#endif