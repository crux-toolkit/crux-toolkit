#pragma once

#include <vector>
#include <map>
#include <string>

#include "Utils.h"

namespace CruxQuant{
    class CruxLFQResults{
        public:
            std::map<std::string, std::vector<ChromatographicPeak>> Peaks;

            // Constructor that accepts a list of spectra files
            CruxLFQResults(const std::vector<std::string>& spectraFiles){
                for (const std::string& spectraFile : spectraFiles) {
                    Peaks[spectraFile] = std::vector<ChromatographicPeak>();
                }
            }
                
    };
}