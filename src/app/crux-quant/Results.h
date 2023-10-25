#pragma once

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "Utils.h"
#include "io/carp.h"

namespace CruxQuant {
class CruxLFQResults {
   public:
    std::map<std::string, std::vector<ChromatographicPeak>> Peaks;

    // Constructor that accepts a list of spectra files
    CruxLFQResults(const std::vector<std::string>& spectraFiles) {
        for (const std::string& spectraFile : spectraFiles) {
            Peaks[spectraFile] = std::vector<ChromatographicPeak>();
        }
    }

    void writeResults(const string& results_file) {
        
        string results_header = "Peptides\tSpectrum File\tQuantity";

        std::ofstream outFile(results_file);
        if (outFile) {
            // Create a custom stream buffer with a larger buffer size (e.g., 8192 bytes)
            const std::size_t bufferSize = 8192;
            char buffer[bufferSize];
            outFile.rdbuf()->pubsetbuf(buffer, bufferSize);
            outFile << results_header << std::endl;
            for (const auto& peak : Peaks) {
                for (const auto& c_p : peak.second) {
                    for (const auto& id : c_p.identifications) {
                        outFile << id.sequence << "\t" << id.spectralFile << "\t" << c_p.identifications.size() << std::endl;
                    }
                }
            }
            outFile.close();
        } else {
            carp(CARP_FATAL, "Failed to open results file for writing.");
        }
    }
};
}  // namespace CruxQuant