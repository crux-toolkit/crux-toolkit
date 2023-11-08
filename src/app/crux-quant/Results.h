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
#include "Peptide.h"
#include "ProteinGroup.h"

using std::string;
using std::vector;
using std::map;

namespace CruxQuant {
class CruxLFQResults {
   public:
    map<string, std::vector<ChromatographicPeak>> Peaks;
    map<string, Peptides> PeptideModifiedSequences;
    map<string, ProteinGroup> ProteinGroups;
    vector<string> spectraFiles;

    // Constructor that accepts a list of spectra files
    CruxLFQResults(const vector<string>& spectraFiles) {
        for (const string& spectraFile : spectraFiles) {
            Peaks[spectraFile] = vector<ChromatographicPeak>();
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

    void setPeptideModifiedSequencesAndProteinGroups(const vector<Identification>& identifications) {
        for(const Identification& id : identifications){
            auto it = PeptideModifiedSequences.find(id.modifications);
            if(it == PeptideModifiedSequences.end()){
                Peptides peptide(id.sequence, id.modifications, id.useForProteinQuant, id.proteinGroups);
                PeptideModifiedSequences[id.modifications] = peptide;
            }else{
                Peptides& peptide = it->second;
                for(const ProteinGroup& proteinGroup : id.proteinGroups){
                    peptide.insertProteinGroup(proteinGroup);
                }
            }

            for(const ProteinGroup& proteinGroup : id.proteinGroups){
               ProteinGroups[proteinGroup.ProteinGroupName] = proteinGroup;
            }

        }
    }

    void calculatePeptideResults(bool quantifyAmbiguousPeptides){
        for(auto& sequence : PeptideModifiedSequences){
            for(string& file: spectraFiles){
                sequence.second.setDetectionType(file, DetectionType::NOT_DETECTED);
                sequence.second.setIntensity(file, 0);
            }
        }

        for(auto& filePeaks: Peaks){
            map<string, vector<ChromatographicPeak>> groupedPeaks;
            for(auto& peak: filePeaks.second){
                if(peak.NumIdentificationsByFullSeq == 1){
                    groupedPeaks[peak.identifications.front().modifications].push_back(peak);
                }
            }

            for(auto& sequenceWithPeaks: groupedPeaks){
                string sequence = sequenceWithPeaks.first;
                double intensity =  std::max_element(
                    sequenceWithPeaks.second.begin(), sequenceWithPeaks.second.end(),
                    [](const ChromatographicPeak& a, const ChromatographicPeak& b) {
                        return a.intensity < b.intensity;
                    })->intensity;
                auto bestPeak = std::find_if(
                    sequenceWithPeaks.second.begin(), sequenceWithPeaks.second.end(),
                    [&intensity](const ChromatographicPeak& p) {
                        return p.intensity == intensity;
                    });

                DetectionType detectionType;
                if(bestPeak->isMbrPeak && intensity > 0){
                    detectionType = DetectionType::MBR;
                }else if(!bestPeak->isMbrPeak && intensity > 0){
                    detectionType = DetectionType::MSMS;
                }else if(!bestPeak->isMbrPeak && intensity == 0){
                    detectionType = DetectionType::MSMSIdentifiedButNotQuantified;
                }else{
                    detectionType = DetectionType::NOT_DETECTED;
                }

                PeptideModifiedSequences[sequence].setDetectionType(filePeaks.first, detectionType);
                PeptideModifiedSequences[sequence].setIntensity(filePeaks.first, intensity);
            }

            // report ambiguous quantification
            vector<ChromatographicPeak> ambiguousPeaks;
            for(auto& peak: filePeaks.second){
                if(peak.NumIdentificationsByFullSeq > 1){
                    ambiguousPeaks.push_back(peak);
                }
            }

            for(auto& ambiguousPeak : ambiguousPeaks){
                for(auto& id : ambiguousPeak.identifications){
                    string sequence = id.modifications;
                    double alreadyRecordedIntensity = PeptideModifiedSequences[sequence].getIntensity(filePeaks.first);
                    double fractionAmbiguous = ambiguousPeak.intensity / (alreadyRecordedIntensity + ambiguousPeak.intensity);

                    if(quantifyAmbiguousPeptides){
                        // If the peptide intensity hasn't been recorded, the intensity is set equal to the intensity of the ambiguous peak
                        if(std::abs(alreadyRecordedIntensity) < 0.01){
                            PeptideModifiedSequences[sequence].setIntensity(filePeaks.first, ambiguousPeak.intensity);
                            PeptideModifiedSequences[sequence].setDetectionType(filePeaks.first, DetectionType::MSMSAmbiguousPeakfinding);
                        }else if(fractionAmbiguous > 0.3){
                            PeptideModifiedSequences[sequence].setDetectionType(filePeaks.first, DetectionType::MSMSAmbiguousPeakfinding);
                        }
                    }else if(fractionAmbiguous > 0.3){
                        PeptideModifiedSequences[sequence].setDetectionType(filePeaks.first, DetectionType::MSMSAmbiguousPeakfinding);
                        PeptideModifiedSequences[sequence].setIntensity(filePeaks.first, 0);
                    }
                }
            }
        }
        // if(!quantifyAmbiguousPeptides){
        //     handleAmbiguityInFractions();
        // }
    }
};
}  // namespace CruxQuant