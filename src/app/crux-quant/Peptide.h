#pragma once

#include <string>
#include <map>
#include <unordered_set>
#include "Utils.h"
#include "ProteinGroup.h"

using std::string;

namespace CruxQuant{

class Peptides{
    private:
        string sequence;
        string modified_sequence;
        bool useForProteinQuant;
        std::map<string, DetectionType> detectionTypeMap;
        std::unordered_set<ProteinGroup> proteinGroups;
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


};
}