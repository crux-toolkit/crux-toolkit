#pragma once

#include <map>
#include <string>
#include <unordered_set>

#include "ProteinGroup.h"
#include "Utils.h"

using std::map;
using std::string;

namespace CruxQuant {

class Peptides {
   private:
    string sequence;
    string modified_sequence;
    bool useForProteinQuant;
    map<string, DetectionType> detectionTypeMap;
    std::unordered_set<ProteinGroup> proteinGroups;
    map<string, double> intensityMap;

   public:
    Peptides() = default;
    Peptides(const string sequence, const string modified_sequence, const bool useForProteinQuant, const std::unordered_set<ProteinGroup> ProteinGroups) {
        this->sequence = sequence;
        this->modified_sequence = modified_sequence;
        this->useForProteinQuant = useForProteinQuant;
        this->proteinGroups = ProteinGroups;
    }
    void setSequence(const string sequence) {
        this->sequence = sequence;
    }

    void setModifiedSequence(const string modified_sequence) {
        this->modified_sequence = modified_sequence;
    }

    void insertProteinGroup(const ProteinGroup proteinGroup) {
        this->proteinGroups.insert(proteinGroup);
    }

    void setDetectionType(const string& file_name, const DetectionType& detectionType) {
        auto it = this->detectionTypeMap.find(file_name);
        if (it != this->detectionTypeMap.end()) {
            it->second = detectionType;
        } else {
            this->detectionTypeMap.insert(std::make_pair(file_name, detectionType));
        }
    }

    void setIntensity(const string& file_name, const double intensity) {
        auto it = this->intensityMap.find(file_name);
        if (it != this->intensityMap.end()) {
            it->second = intensity;
        } else {
            this->intensityMap.insert(std::make_pair(file_name, intensity));
        }
    }

    double getIntensity(const string& file_name) {
        auto it = this->intensityMap.find(file_name);
        if (it != this->intensityMap.end()) {
            return it->second;
        } else {
            return 0;
        }
    }

    DetectionType getDetectionType(const string& file_name) {
        auto it = this->detectionTypeMap.find(file_name);
        if (it != this->detectionTypeMap.end()) {
            return it->second;
        } else {
            return DetectionType::NotDetected;
        }
    }

    string getSequence() {
        return this->sequence;
    }

    bool UnambiguousPeptideQuant() {
        // Assuming Intensities and DetectionTypes are member variables of YourClass

        // Check if there is any intensity greater than 0
        bool anyIntensityGreaterThanZero = std::any_of(intensityMap.begin(), intensityMap.end(),
                                                       [](const std::pair<const std::string, double>& pair) { return pair.second > 0; });

        // Check if there is any DetectionType not equal to MSMSAmbiguousPeakfinding
        bool anyNonAmbiguousDetectionType = std::any_of(detectionTypeMap.begin(), detectionTypeMap.end(),
                                                        [](const std::pair<const std::string, DetectionType>& pair) { return pair.second != DetectionType::MSMSAmbiguousPeakfinding; });

        // Return true if both conditions are met
        return anyIntensityGreaterThanZero && anyNonAmbiguousDetectionType;
    }

    bool getUseForProteinQuant() {
        return this->useForProteinQuant;
    }

    std::unordered_set<ProteinGroup> getProteinGroups() {
        return this->proteinGroups;
    }
};
}  // namespace CruxQuant