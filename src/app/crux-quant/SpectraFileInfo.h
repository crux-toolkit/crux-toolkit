
#pragma once

#include <string>
#include <tuple>

using std::string;

namespace CruxQuant {
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
struct hash<CruxQuant::SpectraFileInfo> {
    size_t operator()(const CruxQuant::SpectraFileInfo &spectraFileInfo) const {
        return std::hash<string>()(spectraFileInfo.FullFilePathWithExtension);
    }
};
}  // namespace std