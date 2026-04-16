#pragma once

#include <string>
#include <unordered_map>
#include "SpectraFileInfo.h"

using std::string;

namespace CruxLFQ {

class ProteinGroup {
   private:
    std::unordered_map<SpectraFileInfo, double> Intensities;

   public:
    string ProteinGroupName;
    ProteinGroup() = default;
    ProteinGroup(const std::string& groupName) : ProteinGroupName(groupName) {}

    int GetHashCode() const {
        return std::hash<string>()(ProteinGroupName);
    };

    bool operator==(const ProteinGroup& other) const {
        return ProteinGroupName == other.ProteinGroupName;
    };
    void SetIntensity(SpectraFileInfo fileInfo, double intensity) {
        auto it = Intensities.find(fileInfo);
        if (it != Intensities.end()) {
            // fileInfo exists in the map
            it->second = intensity;
        } else {
            // fileInfo doesn't exist in the map
            Intensities[fileInfo] = intensity;
        }
    }
};
}  // namespace CruxQuant

namespace std {
template <>
struct hash<CruxLFQ::ProteinGroup> {
    size_t operator()(const CruxLFQ::ProteinGroup& proteinGroup) const {
        return proteinGroup.GetHashCode();
    }
};
}  // namespace std