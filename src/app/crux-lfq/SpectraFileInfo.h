
#pragma once

#include <string>
#include <tuple>

using std::string;

namespace CruxLFQ {
struct SpectraFileInfo {
    int BiologicalReplicate;
    int Fraction;
    string Condition;
    string FullFilePathWithExtension;
    int TechnicalReplicate;

    bool operator<(const SpectraFileInfo &other) const {
        return std::tie(BiologicalReplicate, Fraction, Condition, FullFilePathWithExtension, TechnicalReplicate) < std::tie(other.BiologicalReplicate, other.Fraction, other.Condition, other.FullFilePathWithExtension, TechnicalReplicate);
    }
    bool operator==(const SpectraFileInfo &other) const {
        return FullFilePathWithExtension == other.FullFilePathWithExtension;
    }
};
}  // namespace CruxQuant

namespace std {
template <>
struct hash<CruxLFQ::SpectraFileInfo> {
    size_t operator()(const CruxLFQ::SpectraFileInfo &spectraFileInfo) const {
        return std::hash<string>()(spectraFileInfo.FullFilePathWithExtension);
    }
};
}  // namespace std