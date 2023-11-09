#pragma once

#include <string>
#include <map>
#include <unordered_set>
#include "Utils.h"
#include "ProteinGroup.h"

using std::string;
using std::map;

namespace CruxQuant{

class Peptides{
    private:
        string sequence;
        string modified_sequence;
        bool useForProteinQuant;
        map<string, DetectionType> detectionTypeMap;
        std::unordered_set<ProteinGroup> proteinGroups;
        map<string, double> intensityMap;

    public:
        Peptides() = default;
        Peptides(const string sequence, const string modified_sequence, const bool useForProteinQuant, const std::unordered_set<ProteinGroup> ProteinGroups){
            this->sequence = sequence;
            this->modified_sequence = modified_sequence;
            this->useForProteinQuant = useForProteinQuant;
            this->proteinGroups = ProteinGroups;

        }
        void setSequence(const string sequence){
            this->sequence = sequence;
        }

        void setModifiedSequence(const string modified_sequence){
            this->modified_sequence = modified_sequence;
        }

        void insertProteinGroup(const ProteinGroup proteinGroup){
            this->proteinGroups.insert(proteinGroup);
        }

        void setDetectionType(const string& file_name, const DetectionType& detectionType){
            auto it = this->detectionTypeMap.find(file_name);
            if(it != this->detectionTypeMap.end()){
                it->second = detectionType;
            }else{
                this->detectionTypeMap.insert(std::make_pair(file_name, detectionType));
            }
        }

        void setIntensity(const string& file_name, const double intensity){
            auto it = this->intensityMap.find(file_name);
            if(it != this->intensityMap.end()){
                it->second = intensity;
            }else{
                this->intensityMap.insert(std::make_pair(file_name, intensity));
            }
        }

        double getIntensity(const string& file_name){
            auto it = this->intensityMap.find(file_name);
            if(it != this->intensityMap.end()){
                return it->second;
            }else{
                return 0;
            }
        }


};
}