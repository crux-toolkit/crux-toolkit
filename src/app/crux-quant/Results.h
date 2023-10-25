#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "ChromatographicPeak.h"
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

    string tabSeperatedHeader() {
        std::ostringstream oss;

        // Append the header fields
        oss << "File Name"
            << "\t"
            << "Base Sequence"
            << "\t"
            << "Full Sequence"
            << "\t"
            << "Peptide Monoisotopic Mass"
            << "\t"
            << "MS2 Retention Time"
            << "\t"
            << "Precursor Charge"
            << "\t"
            << "Theoretical MZ"
            << "\t"
            << "Peak intensity"
            << "\t"
            << "Num Charge States Observed"
            << "\t"
            << "PSMs Mapped"
            << "\t"
            << "Peak Split Valley RT"
            << "\t"
            << "Peak Apex Mass Error (ppm)"
            << "\t";
        std::string header = oss.str();

        return header;
    }

    void writeResults(const string& results_file) {
        carp(CARP_INFO, "Writing output...");

        string results_header = tabSeperatedHeader();

        std::ofstream outFile(results_file);
        if (outFile) {
            // Create a custom stream buffer with a larger buffer size (e.g., 8192 bytes)
            const std::size_t bufferSize = 8192;
            char buffer[bufferSize];
            outFile.rdbuf()->pubsetbuf(buffer, bufferSize);
            outFile << results_header << std::endl;

            for (auto& peak : Peaks) {
                auto& peak_data = peak.second;
                std::sort(peak_data.begin(), peak_data.end(), [](const ChromatographicPeak& a, const ChromatographicPeak& b) {
                    // First, sort by SpectraFileInfo.FilenameWithoutExtension in ascending order
                    if (a.spectralFile != b.spectralFile) {
                        return a.spectralFile < b.spectralFile;
                    }
                    // If filenames are the same, sort by Intensity in descending order
                    return a.intensity > b.intensity;
                });

                for (auto& c_p : peak_data) {
                    outFile << c_p.ToString() << std::endl;
                }
            }

            outFile.close();
        } else {
            carp(CARP_FATAL, "Failed to open results file for writing.");
        }
    }
};
}  // namespace CruxQuant