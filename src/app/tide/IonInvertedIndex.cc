#include "IonInvertedIndex.h"
#include "mass_constants.h"
void IonInvertedIndex::insert(Peptide* peptide) {
    for (const auto& b_ion_mz : peptide->BIonMzbins()) {
        ions_.insert({b_ion_mz, peptide});
    }
    for (const auto& y_ion_mz : peptide->YIonMzbins()) {
        ions_.insert({y_ion_mz, peptide});
    }
}
void IonInvertedIndex::erase(Peptide* peptide) {
    for (const auto& b_ion_mz : peptide->BIonMzbins()) {
        ions_.erase(ions_.find({b_ion_mz, peptide}));
    }
    for (const auto& y_ion_mz : peptide->YIonMzbins()) {
        ions_.erase(ions_.find({y_ion_mz, peptide}));
    }
}

IonInvertedIndex::Iterator IonInvertedIndex::lowerBound(int mzbin) {
    return ions_.lower_bound({mzbin, nullptr});
}

IonInvertedIndex::Iterator IonInvertedIndex::upperBound(int mzbin) {
    return ions_.upper_bound({mzbin, nullptr});
}

IonInvertedIndex::ConstIterator IonInvertedIndex::lowerBound(int mzbin) const {
    return ions_.lower_bound({mzbin, nullptr});
}

IonInvertedIndex::ConstIterator IonInvertedIndex::upperBound(int mzbin) const {
    return ions_.upper_bound({mzbin, nullptr});
}
