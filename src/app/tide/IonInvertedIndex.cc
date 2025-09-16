#include "IonInvertedIndex.h"
#include "mass_constants.h"
#include "peptide.h"

void IonInvertedIndex::insert_peaks(Peptide* peptide) {
  // peak is coming from the theoretical peptide, so it cannot be bigger than the capacity_
  if (capacity_ < 10)  // If capacity is small, then we are not using hyperscoring with inverted index
    return;
  for (const auto peak : peptide->peaks_1b) {
    ions_b_[peak].push_back(std::make_pair(peptide->Mass(), peptide));
  }
  for (const auto peak : peptide->peaks_2b) {
    ions_b_[peak].push_back(std::make_pair(peptide->Mass(), peptide));
  }
  for (const auto peak : peptide->peaks_1y) {
    ions_y_[peak].push_back(std::make_pair(peptide->Mass(), peptide));
  }
  for (const auto peak : peptide->peaks_2y) {
    ions_y_[peak].push_back(std::make_pair(peptide->Mass(), peptide));
  }
}
void IonInvertedIndex::score_peaks(double min_precursor_mass, unsigned int peak_mz, double peak_int) {
  if (peak_mz >= capacity_)  // peak_mz is coming from the experimental spectrum
    return;

  // Remove peptides from the inverted index, whose precursor mass in smaller than the lower precursor mass tolerance
  std::deque<std::pair<double, Peptide*>>* queue_ = &(ions_b_.at(peak_mz));
  while (!queue_->empty()) {
    double pept_mass = queue_->front().first;
    if (pept_mass >= min_precursor_mass) {
      break;
    }
    queue_->pop_front();
  }
  // Do the scoring
  for (auto peptide_pair : ions_b_[peak_mz] ) {
    Peptide* pept = peptide_pair.second;
    ++pept->Nb_;
    pept->Ib_ += peak_int;  
  }
}

// void IonInvertedIndex::erase(Peptide* peptide) {
//     for (const auto peak : peptide->peaks_0) {
//         assert(ions_.at(peak).front() == peptide);
//         ions_[peak].pop_front();
//     }
//     for (const auto peak : peptide->peaks_1) {
//         assert(ions_.at(peak).front() == peptide);
//         ions_[peak].pop_front();
//     }
// }

// IonInvertedIndex::Iterator IonInvertedIndex::lowerBound(unsigned int peak) {
//     return ions_.lower_bound(peak);
// }

// IonInvertedIndex::Iterator IonInvertedIndex::upperBound(unsigned int peak) {
//     return ions_.upper_bound(peak);
// }

// IonInvertedIndex::ConstIterator IonInvertedIndex::lowerBound(unsigned int peak) const {
//     return ions_.lower_bound(peak);
// }

// IonInvertedIndex::ConstIterator IonInvertedIndex::upperBound(unsigned int peak) const {
//     return ions_.upper_bound(peak);
// }
