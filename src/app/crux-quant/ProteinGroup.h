#pragma once

#include <string>


using std::string;

namespace CruxQuant{
    class ProteinGroup{
        public:
            string ProteinGroupName;
            ProteinGroup(const std::string& groupName) : ProteinGroupName(groupName) {}

            struct Hash{
                size_t operator()(const ProteinGroup& proteinGroup) const{
                    return std::hash<string>()(proteinGroup.ProteinGroupName);
                }
            };
    }
}