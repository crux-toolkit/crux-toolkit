#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "SpectraFileInfo.h"
#include "ChromatographicPeak.h"
#include "Peptide.h"
#include "ProteinGroup.h"
#include "Utils.h"
#include "io/carp.h"

using std::map;
using std::string;
using std::vector;

namespace CruxQuant {

class CruxLFQResults {
   public:
    map<string, std::vector<ChromatographicPeak>> Peaks;
    map<string, Peptides> PeptideModifiedSequences;
    map<string, ProteinGroup> ProteinGroups;
    vector<SpectraFileInfo> spectraFiles;

    // Constructor that accepts a list of spectra files
    CruxLFQResults(const vector<string> &spectra_files) {
        for (const string &spectra_file : spectra_files) {
            Peaks[spectra_file] = vector<ChromatographicPeak>();
            SpectraFileInfo spectraFileInfo;
            spectraFileInfo.FullFilePathWithExtension = spectra_file;
            spectraFiles.push_back(spectraFileInfo);
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
            << "Peak Detection Type"
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

    void writeResults(const string &results_file) {
        carp(CARP_INFO, "Writing output...");

        string results_header = tabSeperatedHeader();

        std::ofstream outFile(results_file);
        if (outFile) {
            // Create a custom stream buffer with a larger buffer size (e.g., 8192 bytes)
            const std::size_t bufferSize = 8192;
            char buffer[bufferSize];
            outFile.rdbuf()->pubsetbuf(buffer, bufferSize);
            outFile << results_header << std::endl;

            for (auto &peak : Peaks) {
                auto &peak_data = peak.second;
                std::sort(peak_data.begin(), peak_data.end(), [](const ChromatographicPeak &a, const ChromatographicPeak &b) {
                    // First, sort by SpectraFileInfo.FilenameWithoutExtension in ascending order
                    if (a.spectralFile != b.spectralFile) {
                        return a.spectralFile < b.spectralFile;
                    }
                    // If filenames are the same, sort by Intensity in descending order
                    return a.intensity > b.intensity; });

                for (auto &c_p : peak_data) {
                    outFile << c_p.ToString() << std::endl;
                }
            }

            outFile.close();
        } else {
            carp(CARP_FATAL, "Failed to open results file for writing.");
        }
    }

    void setPeptideModifiedSequencesAndProteinGroups(const vector<Identification> &identifications) {
        for (const Identification &id : identifications) {
            auto it = PeptideModifiedSequences.find(id.modifications);
            if (it == PeptideModifiedSequences.end()) {
                Peptides peptide(id.sequence, id.modifications, id.useForProteinQuant, id.proteinGroups);
                PeptideModifiedSequences[id.modifications] = peptide;
            } else {
                Peptides &peptide = it->second;
                for (const ProteinGroup &proteinGroup : id.proteinGroups) {
                    peptide.insertProteinGroup(proteinGroup);
                }
            }

            for (const ProteinGroup &proteinGroup : id.proteinGroups) {
                ProteinGroups[proteinGroup.ProteinGroupName] = proteinGroup;
            }
        }
    }

    void calculatePeptideResults(bool quantifyAmbiguousPeptides) {
        for (auto &sequence : PeptideModifiedSequences) {
            for (auto &file : spectraFiles) {
                sequence.second.setDetectionType(file.FullFilePathWithExtension, DetectionType::NotDetected);
                sequence.second.setIntensity(file.FullFilePathWithExtension, 0);
            }
        }

        for (auto &filePeaks : Peaks) {
            map<string, vector<ChromatographicPeak>> groupedPeaks;
            for (auto &peak : filePeaks.second) {
                if (peak.NumIdentificationsByFullSeq == 1) {
                    groupedPeaks[peak.identifications.front().modifications].push_back(peak);
                }
            }

            for (auto &sequenceWithPeaks : groupedPeaks) {
                string sequence = sequenceWithPeaks.first;
                double intensity = std::max_element(
                                       sequenceWithPeaks.second.begin(), sequenceWithPeaks.second.end(),
                                       [](const ChromatographicPeak &a, const ChromatographicPeak &b) {
                                           return a.intensity < b.intensity;
                                       })
                                       ->intensity;
                auto bestPeak = std::find_if(
                    sequenceWithPeaks.second.begin(), sequenceWithPeaks.second.end(),
                    [&intensity](const ChromatographicPeak &p) {
                        return p.intensity == intensity;
                    });

                DetectionType detectionType;
                if (bestPeak->isMbrPeak && intensity > 0) {
                    detectionType = DetectionType::MBR;
                } else if (!bestPeak->isMbrPeak && intensity > 0) {
                    detectionType = DetectionType::MSMS;
                } else if (!bestPeak->isMbrPeak && intensity == 0) {
                    detectionType = DetectionType::MSMSIdentifiedButNotQuantified;
                } else {
                    detectionType = DetectionType::NotDetected;
                }

                PeptideModifiedSequences[sequence].setDetectionType(filePeaks.first, detectionType);
                PeptideModifiedSequences[sequence].setIntensity(filePeaks.first, intensity);
            }

            // report ambiguous quantification
            vector<ChromatographicPeak> ambiguousPeaks;
            for (auto &peak : filePeaks.second) {
                if (peak.NumIdentificationsByFullSeq > 1) {
                    ambiguousPeaks.push_back(peak);
                }
            }

            for (auto &ambiguousPeak : ambiguousPeaks) {
                for (auto &id : ambiguousPeak.identifications) {
                    string sequence = id.modifications;
                    double alreadyRecordedIntensity = PeptideModifiedSequences[sequence].getIntensity(filePeaks.first);
                    double fractionAmbiguous = ambiguousPeak.intensity / (alreadyRecordedIntensity + ambiguousPeak.intensity);

                    if (quantifyAmbiguousPeptides) {
                        // If the peptide intensity hasn't been recorded, the intensity is set equal to the intensity of the ambiguous peak
                        if (std::abs(alreadyRecordedIntensity) < 0.01) {
                            PeptideModifiedSequences[sequence].setIntensity(filePeaks.first, ambiguousPeak.intensity);
                            PeptideModifiedSequences[sequence].setDetectionType(filePeaks.first, DetectionType::MSMSAmbiguousPeakfinding);
                        } else if (fractionAmbiguous > 0.3) {
                            PeptideModifiedSequences[sequence].setDetectionType(filePeaks.first, DetectionType::MSMSAmbiguousPeakfinding);
                        }
                    } else if (fractionAmbiguous > 0.3) {
                        PeptideModifiedSequences[sequence].setDetectionType(filePeaks.first, DetectionType::MSMSAmbiguousPeakfinding);
                        PeptideModifiedSequences[sequence].setIntensity(filePeaks.first, 0);
                    }
                }
            }
        }
        if (!quantifyAmbiguousPeptides) {
            handleAmbiguityInFractions();
        }
    }

    // May be redundant with our current implementation
    void handleAmbiguityInFractions() {
        // handle ambiguous peaks in fractionated data
        // if the largest fraction intensity is ambiguous, zero out the other fractions for the sample
        map<string, vector<SpectraFileInfo>> sampleGroups;
        for (auto &file : spectraFiles) {
            sampleGroups[file.Condition].push_back(file);
        }

        for (auto &sampleGroup : sampleGroups) {
            map<int, vector<SpectraFileInfo>> samples;
            for (auto &file : sampleGroup.second) {
                samples[file.BiologicalReplicate].push_back(file);
                for (auto &sample : samples) {
                    // skip unfractionated samples
                    std::set<int> uniqueFractions;
                    for (auto &p : sample.second) {
                        uniqueFractions.insert(p.Fraction);
                    }
                    if (uniqueFractions.size() == 1) {
                        continue;
                    }

                    map<std::pair<SpectraFileInfo, std::string>, std::vector<ChromatographicPeak>> peaksForEachSequence;

                    for (auto &file : sample.second) {
                        for (auto &peak : Peaks[file.FullFilePathWithExtension]) {
                            for (auto &id : peak.identifications) {
                                auto key = std::make_pair(file, id.modifications);
                                auto it = peaksForEachSequence.find(key);
                                if (it != peaksForEachSequence.end()) {
                                    it->second.push_back(peak);
                                } else {
                                    peaksForEachSequence[key] = std::vector<ChromatographicPeak>{peak};
                                }
                            }
                        }
                    }

                    vector<Peptides> peptides;
                    for (auto &item : PeptideModifiedSequences) {
                        peptides.push_back(item.second);
                    }

                    vector<std::pair<double, DetectionType>> fractionIntensitiesWithDetectionTypes;

                    for (auto &peptide : peptides) {
                        fractionIntensitiesWithDetectionTypes.clear();
                        bool ambiguityObservedInSample = false;

                        for (auto &file : sample.second) {
                            double fractionIntensity = peptide.getIntensity(file.FullFilePathWithExtension);
                            DetectionType detectionType = peptide.getDetectionType(file.FullFilePathWithExtension);

                            if (detectionType == DetectionType::MSMSAmbiguousPeakfinding) {
                                ambiguityObservedInSample = true;

                                auto key = std::make_pair(file, peptide.getSequence());
                                double sum = 0.0;
                                for (auto &p : peaksForEachSequence[key]) {
                                    sum += p.intensity;
                                }
                                fractionIntensity = sum;
                            }

                            fractionIntensitiesWithDetectionTypes.push_back(std::make_pair(fractionIntensity, detectionType));
                        }
                        //
                        if (ambiguityObservedInSample) {
                            // Define a comparison function
                            auto compare = [](
                                               const std::pair<double, DetectionType> &a,
                                               const std::pair<double, DetectionType> &b) {
                                return a.first < b.first;
                            };

                            // Use std::max_element to find the highest intensity
                            auto highestIntensity = *std::max_element(
                                fractionIntensitiesWithDetectionTypes.begin(),
                                fractionIntensitiesWithDetectionTypes.end(),
                                compare);

                            DetectionType highestFractionIntensityDetectionType = highestIntensity.second;
                            if (highestFractionIntensityDetectionType == DetectionType::MSMSAmbiguousPeakfinding) {
                                // highest fraction intensity is ambiguous - zero out the other fractions
                                for (SpectraFileInfo &file : sample.second) {
                                    peptide.setIntensity(file.FullFilePathWithExtension, 0);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void calculateProteinResultsMedianPolish(bool useSharedPeptides) {
        // Reset protein intensities to 0
        for (auto &proteinGroup : ProteinGroups) {
            for (const auto &file : spectraFiles) {
                proteinGroup.second.SetIntensity(file, 0);
            }
        }

        vector<Peptides> peptides;
        for (auto &item : PeptideModifiedSequences) {
            if(item.second.UnambiguousPeptideQuant()){
                peptides.push_back(item.second);
            }
        }

        std::unordered_map<ProteinGroup, vector<Peptides>> proteinGroupToPeptides;

        for (Peptides& peptide : peptides) {
            if (!peptide.getUseForProteinQuant() || (peptide.getProteinGroups().size() > 1 && !useSharedPeptides)) {
                continue;
            }

            for (ProteinGroup pg : peptide.getProteinGroups()) {
                auto it = proteinGroupToPeptides.find(pg);
                if (it != proteinGroupToPeptides.end()) {
                    // ProteinGroup found, add the peptide to the existing list
                    it->second.push_back(peptide);
                } else {
                    // ProteinGroup not found, create a new entry
                    proteinGroupToPeptides[pg] = {peptide};
                }
            }
        }
    }
};
}  // namespace CruxQuant