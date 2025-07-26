#ifndef TOP_N_PEPTIDES_H
#define TOP_N_PEPTIDES_H
#include "peptide.h"
#include <set>
class IonInvertedIndex {
public:
    struct Ion {
        int mzbin;
        Peptide* peptide;
        friend bool operator<(const Ion& a, const Ion& b) {
            return a.mzbin < b.mzbin;
        }
    };
    using IonsStorage = std::multiset<Ion>;
    using Iterator = IonsStorage::iterator;
    using ConstIterator = IonsStorage::const_iterator;
    IonInvertedIndex() = default;
    void erase(Peptide* peptide);
    void insert(Peptide* peptide);
    Iterator lowerBound(int mzbin);
    Iterator upperBound(int mzbin);

    ConstIterator lowerBound(int mzbin) const;
    ConstIterator upperBound(int mzbin) const;

    Iterator begin() { return ions_.begin(); }
    ConstIterator begin() const { return ions_.begin(); }
    Iterator end() { return ions_.end(); }
    ConstIterator end() const { return ions_.end(); }
private:
    IonsStorage ions_;
};
#endif