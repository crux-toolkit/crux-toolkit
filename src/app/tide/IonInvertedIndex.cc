#include "IonInvertedIndex.h"
#include "mass_constants.h"
#include "peptide.h"

void IonInvertedIndex::insert_peaks(Peptide* peptide) {

  // peak is coming from the theoretical peptide, so it cannot be bigger than the capacity_
  if (capacity_ < 10)  // If capacity is small, then we are not using hyperscoring with inverted index
    return;

  for (const auto peak : peptide->peaks_1b) {
    ions_b_.at(peak).emplace_back(peptide);
  }
  for (const auto peak : peptide->peaks_1y) {
    ions_y_.at(peak).emplace_back(peptide);
  }
}

void IonInvertedIndex::pop_peaks(Peptide* peptide) {

  if (capacity_ < 10)  // peak_mz is coming from the experimental spectrum
    return;  

  for (const auto peak : peptide->peaks_1b) {
    assert(ions_b_.at(peak).front() == peptide);
    ions_b_[peak].pop_front();
  }
  for (const auto peak : peptide->peaks_1y) {
    assert(ions_y_.at(peak).front() == peptide);
    ions_y_[peak].pop_front();
  }
}

void IonInvertedIndex::score_peaks(unsigned int peak_mz, double peak_int, int charge) {

  if (peak_mz >= capacity_)  // peak_mz is coming from the experimental spectrum
    return;

}