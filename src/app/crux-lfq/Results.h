#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "ChromatographicPeak.h"
#include "Peptide.h"
#include "ProteinGroup.h"
#include "SpectraFileInfo.h"
#include "Utils.h"
#include "CQStatistics.h"
#include "io/carp.h"

using std::map;
using std::string;
using std::vector;

namespace CruxLFQ {

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
            if (item.second.UnambiguousPeptideQuant()) {
                peptides.push_back(item.second);
            }
        }

        std::unordered_map<ProteinGroup, vector<Peptides>> proteinGroupToPeptides;

        for (Peptides &peptide : peptides) {
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

        std::map<std::string, std::vector<SpectraFileInfo>> filesGroupedByCondition;
        for (const auto &file : spectraFiles) {
            filesGroupedByCondition[file.Condition].push_back(file);
        }
        std::map<int, std::vector<SpectraFileInfo>> groupedByBiologicalReplicate;
        for (const auto &group : filesGroupedByCondition) {
            for (const auto &file : group.second) {
                groupedByBiologicalReplicate[file.BiologicalReplicate].push_back(file);
            }
        }
        // quantify each protein
        for (auto &proteinGroupPair : ProteinGroups) {
            ProteinGroup &proteinGroup = proteinGroupPair.second;

            auto it = proteinGroupToPeptides.find(proteinGroup);
            if (it != proteinGroupToPeptides.end()) {
                // set up peptide intensity table
                // top row is the column effects, left column is the row effects
                // the other cells are peptide intensity measurements
                std::vector<std::string> conditions;
                std::transform(spectraFiles.begin(), spectraFiles.end(), std::back_inserter(conditions), [](const CruxLFQ::SpectraFileInfo &file) {
                    return file.Condition;
                });
                int numSamples = std::unordered_set<std::string>(conditions.begin(), conditions.end()).size();
                std::vector<std::vector<double>> peptideIntensityMatrix(it->second.size() + 1, std::vector<double>(numSamples + 1));

                // populate matrix w/ log2-transformed peptide intensities
                // if a value is missing, it will be filled with NaN
                int sampleN = 0;
                for (const auto &group : filesGroupedByCondition) {
                    for (const auto &sample : groupedByBiologicalReplicate) {
                        for (Peptides &peptide : it->second) {
                            double sampleIntensity = 0;
                            double highestFractionIntensity = 0;

                            // the fraction w/ the highest intensity is used as the sample intensity for this peptide.
                            // if there is more than one replicate of the fraction, then the replicate intensities are averaged
                            std::map<int, std::vector<SpectraFileInfo>> groupedByFraction;
                            for (const auto &file : sample.second) {
                                groupedByFraction[file.Fraction].push_back(file);
                            }
                            for (const auto &fraction : sample.second) {
                                double fractionIntensity = 0;
                                int replicatesWithValidValues = 0;

                                for (const SpectraFileInfo &replicate : groupedByFraction[fraction.Fraction]) {
                                    double replicateIntensity = peptide.getIntensity(replicate.FullFilePathWithExtension);

                                    if (replicateIntensity > 0) {
                                        fractionIntensity += replicateIntensity;
                                        replicatesWithValidValues++;
                                    }
                                }

                                if (replicatesWithValidValues > 0) {
                                    fractionIntensity /= replicatesWithValidValues;
                                }

                                if (fractionIntensity > highestFractionIntensity) {
                                    highestFractionIntensity = fractionIntensity;
                                    sampleIntensity = highestFractionIntensity;
                                }
                            }

                            int sampleNumber = sample.first;

                            if (sampleIntensity == 0) {
                                sampleIntensity = std::numeric_limits<double>::quiet_NaN();
                            } else {
                                sampleIntensity = std::log2(sampleIntensity);
                            }

                            auto peptideIndex = std::distance(it->second.begin(), std::find(it->second.begin(), it->second.end(), peptide));
                            peptideIntensityMatrix[peptideIndex + 1][sampleN + 1] = sampleIntensity;
                        }

                        sampleN++;
                    }
                }

                // if there are any peptides that have only one measurement, mark them as NaN
                // unless we have ONLY peptides with one measurement
                int peptidesWithMoreThanOneMmt = std::count_if(peptideIntensityMatrix.begin() + 1, peptideIntensityMatrix.end(),
                                                               [](const std::vector<double> &row) {
                                                                   return std::count_if(row.begin() + 1, row.end(), [](double cell) { return !std::isnan(cell); }) > 1;
                                                               });

                if (peptidesWithMoreThanOneMmt > 0) {
                    for (size_t i = 1; i < peptideIntensityMatrix.size(); i++) {
                        int validValueCount = std::count_if(peptideIntensityMatrix[i].begin(), peptideIntensityMatrix[i].end(),
                                                            [](double p) { return !std::isnan(p) && p != 0; });

                        if (validValueCount < 2 && numSamples >= 2) {
                            for (size_t j = 1; j < peptideIntensityMatrix[0].size(); j++) {
                                peptideIntensityMatrix[i][j] = std::numeric_limits<double>::quiet_NaN();
                            }
                        }
                    }
                }

                // do median polish protein quantification
                // row effects in a protein can be considered ~ relative ionization efficiency
                // column effects are differences between conditions
                MedianPolish(peptideIntensityMatrix);

                double overallEffect = peptideIntensityMatrix[0][0];
                std::vector<double> columnEffects(peptideIntensityMatrix[0].begin() + 1, peptideIntensityMatrix[0].end());
                auto peptidesForThisProtein = it->second;
                double referenceProteinIntensity = std::pow(2, overallEffect) * peptidesForThisProtein.size();

                // check for unquantifiable proteins; these are proteins w/ quantified peptides, but
                // the protein is still not quantifiable because there are not peptides to compare across runs

                std::vector<std::string> possibleUnquantifiableSample;
                sampleN = 0;
                for (auto &group : filesGroupedByCondition) {
                    for (auto &sample : groupedByBiologicalReplicate) {
                        bool isMissingValue = true;
                        for (SpectraFileInfo spectraFile : sample.second) {
                            if (std::any_of(
                                    peptidesForThisProtein.begin(),
                                    peptidesForThisProtein.end(),
                                    [&](Peptides p) { return p.getIntensity(spectraFile.FullFilePathWithExtension) != 0; })) {
                                isMissingValue = false;
                                break;
                            }
                        }

                        if (!isMissingValue && columnEffects[sampleN] == 0) {
                            possibleUnquantifiableSample.push_back(group.first + "_" + std::to_string(sample.first));
                        }

                        sampleN++;
                    }
                }
                // set the sample protein intensities

                sampleN = 0;
                for (auto &group : filesGroupedByCondition)  
                {
                    for (auto &sample : groupedByBiologicalReplicate)
                    {
                        // this step un-logs the protein "intensity". in reality this value is more like a fold-change
                        // than an intensity, but unlike a fold-change it's not relative to a particular sample.
                        // by multiplying this value by the reference protein intensity calculated earlier, then we get
                        // a protein intensity value
                        double columnEffect = columnEffects[sampleN];
                        double sampleProteinIntensity = std::pow(2, columnEffect) * referenceProteinIntensity;

                        // the column effect can be 0 in some cases. sometimes it's a valid value and sometimes it's not.
                        // so we need to check to see if it is actually a valid value
                        bool isMissingValue = true;

                        for (SpectraFileInfo spectraFile : sample.second) {
                            if (std::any_of(
                                peptidesForThisProtein.begin(), 
                                peptidesForThisProtein.end(), [&](Peptides p) { return p.getIntensity(spectraFile.FullFilePathWithExtension) != 0; })) {
                                isMissingValue = false;
                                break;
                            }
                        }

                        if (!isMissingValue) {
                            if (possibleUnquantifiableSample.size() > 1 && std::find(possibleUnquantifiableSample.begin(), possibleUnquantifiableSample.end(), group.first + "_" + std::to_string(sample.first)) != possibleUnquantifiableSample.end()) {
                                proteinGroup.SetIntensity(sample.second.front(), std::numeric_limits<double>::quiet_NaN());
                            } else {
                                proteinGroup.SetIntensity(sample.second.front(), sampleProteinIntensity);
                            }
                        }

                        sampleN++;
                    }
                }
            }
        }
    };

    void MedianPolish(std::vector<vector<double>> &table, int maxIterations = 10, double improvementCutoff = 0.0001) {
        // technically, this is weighted mean polish and not median polish.
        // but it should give similar results while being more robust to issues
        // arising from missing values.
        // the weights are inverse square difference to median.

        // subtract overall effect
        std::vector<double> allValues;
        for (const auto &row : table) {
            for (const auto &cell : row) {
                if (!std::isnan(cell) && cell != 0) {
                    allValues.push_back(cell);
                }
            }
        }

        if (!allValues.empty()) {
            double overallEffect = Median(allValues);
            table[0][0] += overallEffect;
            for (size_t r = 1; r < table.size(); r++) {
                for (size_t c = 1; c < table[0].size(); c++) {
                    table[r][c] -= overallEffect;
                }
            }
        }

        double sumAbsoluteResiduals = std::numeric_limits<double>::max();

        for (int i = 0; i < maxIterations; i++) {
            // subtract row effects
            for (size_t r = 0; r < table.size(); r++) {
                std::vector<double> rowValues(table[r].begin() + 1, table[r].end());
                rowValues.erase(std::remove_if(rowValues.begin(), rowValues.end(), [](double p) { return std::isnan(p); }), rowValues.end());

                if (!rowValues.empty()) {
                    double rowMedian = Median(rowValues);
                    std::vector<double> weights;
                    for (auto &p : rowValues) {
                        weights.push_back(1.0 / std::max(0.0001, std::pow(p - rowMedian, 2)));
                    }
                    double rowEffect = std::inner_product(rowValues.begin(), rowValues.end(), weights.begin(), 0.0) / std::accumulate(weights.begin(), weights.end(), 0.0);
                    table[r][0] += rowEffect;

                    for (size_t c = 1; c < table[0].size(); c++) {
                        table[r][c] -= rowEffect;
                    }
                }
            }

            // subtract column effects
            for (size_t c = 0; c < table[0].size(); c++) {
                std::vector<double> colValues;
                for (size_t r = 1; r < table.size(); r++) {
                    if (!std::isnan(table[r][c])) {
                        colValues.push_back(table[r][c]);
                    }
                }

                if (!colValues.empty()) {
                    double colMedian = Median(colValues);
                    std::vector<double> weights;
                    for (auto &p : colValues) {
                        weights.push_back(1.0 / std::max(0.0001, std::pow(p - colMedian, 2)));
                    }
                    double colEffect = std::inner_product(colValues.begin(), colValues.end(), weights.begin(), 0.0) / std::accumulate(weights.begin(), weights.end(), 0.0);
                    table[0][c] += colEffect;

                    for (size_t r = 1; r < table.size(); r++) {
                        table[r][c] -= colEffect;
                    }
                }
            }

            // calculate sum of absolute residuals and end the algorithm if it is not improving
            double iterationSumAbsoluteResiduals = 0.0;
            for (size_t r = 1; r < table.size(); r++) {
                for (size_t c = 1; c < table[r].size(); c++) {
                    if (!std::isnan(table[r][c])) {
                        iterationSumAbsoluteResiduals += std::abs(table[r][c]);
                    }
                }
            }

            if (std::abs((iterationSumAbsoluteResiduals - sumAbsoluteResiduals) / sumAbsoluteResiduals) < improvementCutoff) {
                break;
            }

            sumAbsoluteResiduals = iterationSumAbsoluteResiduals;
        }
    }
};

}  // namespace CruxQuant