#ifndef ION_INVERTED_INDEX
#define ION_INVERTED_INDEX
#include <set>
class Peptide;
class IonInvertedIndex {
public:
    struct Ion {
        unsigned int mzbin;
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
    Iterator lowerBound(unsigned int peak);
    Iterator upperBound(unsigned int peak);

    ConstIterator lowerBound(unsigned int peak) const;
    ConstIterator upperBound(unsigned int peak) const;

    Iterator begin() { return ions_.begin(); }
    ConstIterator begin() const { return ions_.begin(); }
    Iterator end() { return ions_.end(); }
    ConstIterator end() const { return ions_.end(); }
private:
    IonsStorage ions_;
};
#endif