#include "IonInvertedIndex.h"
#include "mass_constants.h"
#include "peptide.h"
void IonInvertedIndex::insert(Peptide* peptide) {
    for (const auto peak : peptide->peaks_0) {
        ions_.insert({peak, peptide});
    }
    for (const auto peak : peptide->peaks_1) {
        ions_.insert({peak, peptide});
    }
}
void IonInvertedIndex::erase(Peptide* peptide) {
    for (const auto peak : peptide->peaks_0) {
        ions_.erase(ions_.find({peak, peptide}));
    }
    for (const auto peak : peptide->peaks_1) {
        ions_.erase(ions_.find({peak, peptide}));
    }
}

IonInvertedIndex::Iterator IonInvertedIndex::lowerBound(unsigned int peak) {
    return ions_.lower_bound({peak, nullptr});
}

IonInvertedIndex::Iterator IonInvertedIndex::upperBound(unsigned int peak) {
    return ions_.upper_bound({peak, nullptr});
}

IonInvertedIndex::ConstIterator IonInvertedIndex::lowerBound(unsigned int peak) const {
    return ions_.lower_bound({peak, nullptr});
}

IonInvertedIndex::ConstIterator IonInvertedIndex::upperBound(unsigned int peak) const {
    return ions_.upper_bound({peak, nullptr});
}
