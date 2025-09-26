#include "IonInvertedIndex.h"
#include "mass_constants.h"
#include "peptide.h"

void IonInvertedIndex::insert_peaks(Peptide* peptide) {
  // peak is coming from the theoretical peptide, so it cannot be bigger than the capacity_
  if (capacity_ < 10)  // If capacity is small, then we are not using hyperscoring with inverted index
    return;

  for (const auto peak : peptide->peaks_1b) {
    ions_b_.at(peak).emplace_back(peptide->Mass(), peptide);
  }
  // for (const auto peak : peptide->peaks_2b) {
  //   ions_b_.at(peak).emplace_back(peptide->Mass(), peptide);
  // }
  for (const auto peak : peptide->peaks_1y) {
    ions_y_.at(peak).emplace_back(peptide->Mass(), peptide);
  }
  // for (const auto peak : peptide->peaks_2y) {
  //   ions_y_.at(peak).emplace_back(peptide->Mass(), peptide);
  // }
}
void IonInvertedIndex::score_peaks(double min_precursor_mass, unsigned int peak_mz, double peak_int) {
  if (peak_mz >= capacity_)  // peak_mz is coming from the experimental spectrum
    return;

  // Remove peptides from the inverted index, whose precursor mass in smaller than the lower precursor mass tolerance
  auto begin_b = ions_b_.at(peak_mz).begin();
  while (begin_b != ions_b_.at(peak_mz).end()) {
    double pept_mass = begin_b->first;
    if (pept_mass >= min_precursor_mass) {
      break;
    }
    begin_b++;
  }
  auto begin_y = ions_y_.at(peak_mz).begin();
  while (begin_y != ions_y_.at(peak_mz).end()) {
    double pept_mass = begin_y->first;
    if (pept_mass >= min_precursor_mass) {
      break;
    }
    begin_y++;
  }
  // Do the scoring
  for (auto it = begin_b; it != ions_b_.at(peak_mz).end(); ++it) {
    Peptide* pept = it->second;
    ++pept->Nb_;
    pept->Ib_ += peak_int; 
  }

  for (auto it = begin_y; it != ions_y_.at(peak_mz).end(); ++it) {
    Peptide* pept = it->second;
    ++pept->Ny_;
    pept->Iy_ += peak_int; 
  }
}

void IonInvertedIndex::pop_peaks(Peptide* peptide) {
    for (const auto peak : peptide->peaks_1b) {
        assert(ions_b_.at(peak).front().second == peptide);
        ions_b_[peak].pop_front();
    }
    for (const auto peak : peptide->peaks_1y) {
        assert(ions_y_.at(peak).front().second == peptide);
        ions_y_[peak].pop_front();
    }
}

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
