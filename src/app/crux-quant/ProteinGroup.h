#pragma once

#include <string>


using std::string;

namespace CruxQuant{
    class ProteinGroup{
        public:
            string ProteinGroupName;
            ProteinGroup() = default;
            ProteinGroup(const std::string& groupName) : ProteinGroupName(groupName) {}

            int GetHashCode() const{
                return std::hash<string>()(ProteinGroupName);
            };

            bool operator==(const ProteinGroup& other) const{
                return ProteinGroupName == other.ProteinGroupName;
            };
    };
}

namespace std{
    template<>
    struct hash<CruxQuant::ProteinGroup>{
        size_t operator()(const CruxQuant::ProteinGroup& proteinGroup) const{
            return proteinGroup.GetHashCode();
        }
    };
}