#include "IonInvertedIndex.h"
#include "mass_constants.h"
#include "peptide.h"
void IonInvertedIndex::insert(Peptide* peptide) {
    for (const auto peak : peptide->peaks_0) {
        ions_[peak].push_back(peptide);
    }
    for (const auto peak : peptide->peaks_1) {
        ions_[peak].push_back(peptide);
    }
}
void IonInvertedIndex::erase(Peptide* peptide) {
    for (const auto peak : peptide->peaks_0) {
        assert(ions_.at(peak).front() == peptide);
        ions_[peak].pop_front();
    }
    for (const auto peak : peptide->peaks_1) {
        assert(ions_.at(peak).front() == peptide);
        ions_[peak].pop_front();
    }
}

IonInvertedIndex::Iterator IonInvertedIndex::lowerBound(unsigned int peak) {
    return ions_.lower_bound(peak);
}

IonInvertedIndex::Iterator IonInvertedIndex::upperBound(unsigned int peak) {
    return ions_.upper_bound(peak);
}

IonInvertedIndex::ConstIterator IonInvertedIndex::lowerBound(unsigned int peak) const {
    return ions_.lower_bound(peak);
}

IonInvertedIndex::ConstIterator IonInvertedIndex::upperBound(unsigned int peak) const {
    return ions_.upper_bound(peak);
}
