#pragma once

#include <iterator>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "Utils.h"

using std::string;
using std::unordered_set;
using std::vector;

namespace CruxLFQ {

inline string ChromatographicPeakTabSeperatedHeader() {
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

struct ChromatographicPeak {
    vector<Identification> identifications;
    string spectralFile;
    vector<IsotopicEnvelope> isotopicEnvelopes;
    double intensity;
    double massError;
    int numChargeStatesObserved;
    IsotopicEnvelope apex;
    bool isMbrPeak;
    int _index = std::numeric_limits<int>::max();
    bool SplitRT = 0;
    int NumIdentificationsByBaseSeq = 1;
    int NumIdentificationsByFullSeq = 1;
    double MbrScore;  // Not really used.

    ChromatographicPeak(const Identification &id, bool _isMbrPeak, const string &_spectralFile)
        : isMbrPeak(_isMbrPeak), spectralFile(_spectralFile) {
        isotopicEnvelopes.push_back(IsotopicEnvelope());
        identifications.push_back(id);
        massError = std::numeric_limits<double>::quiet_NaN();
        numChargeStatesObserved = 0;
        _index = std::numeric_limits<int>::max();
        SplitRT = 0;
        NumIdentificationsByBaseSeq = 1;
        NumIdentificationsByFullSeq = 1;
    }

    void calculateIntensityForThisFeature(bool integrate) {
        if (!isotopicEnvelopes.empty()) {
            double maxIntensity = std::max_element(
                                      isotopicEnvelopes.begin(),
                                      isotopicEnvelopes.end(),
                                      [](const IsotopicEnvelope &a, const IsotopicEnvelope &b) {
                                          return a.intensity < b.intensity;
                                      })
                                      ->intensity;

            apex = *std::find_if(
                isotopicEnvelopes.begin(),
                isotopicEnvelopes.end(),
                [maxIntensity](const IsotopicEnvelope &p) {
                    return p.intensity == maxIntensity;
                });

            if (integrate) {
                intensity = std::accumulate(isotopicEnvelopes.begin(), isotopicEnvelopes.end(), 0.0, [](double sum, const IsotopicEnvelope &p) { return sum + p.intensity; });
            } else {
                intensity = apex.intensity;
            }

            massError = std::numeric_limits<double>::quiet_NaN();

            for (const auto &id : identifications) {
                double massErrorForId = ((toMass(apex.indexedPeak.mz, apex.chargeState) - id.peakFindingMass) / id.peakFindingMass) * 1e6;

                if (std::isnan(massError) || std::abs(massErrorForId) < std::abs(massError)) {
                    massError = massErrorForId;
                }
            }

            std::unordered_set<int> uniqueChargeStates;
            for (const IsotopicEnvelope &envelope : isotopicEnvelopes) {
                uniqueChargeStates.insert(envelope.chargeState);
            }
            numChargeStatesObserved = static_cast<int>(uniqueChargeStates.size());
        } else {
            intensity = 0;
            massError = std::numeric_limits<double>::quiet_NaN();
            numChargeStatesObserved = 0;
            apex = {};
        }
    }

    string ToString() {
        std::set<string> uniqueBaseSequences;
        std::set<string> uniqueModifiedSequences;

        for (const Identification &id : identifications) {
            uniqueBaseSequences.insert(id.sequence);
            uniqueModifiedSequences.insert(id.modifications);
        }

        std::ostringstream oss;

        oss << spectralFile << "\t";
        std::copy(uniqueBaseSequences.begin(), uniqueBaseSequences.end(),
                  std::ostream_iterator<std::string>(oss, "|"));
        oss << '\t';
        std::copy(uniqueModifiedSequences.begin(), uniqueModifiedSequences.end(),
                  std::ostream_iterator<std::string>(oss, "|"));
        oss << '\t';
        oss << identifications.front().monoIsotopicMass << '\t';
        if (!isMbrPeak) {
            // Append Ms2RetentionTimeInMinutes followed by a tab character
            oss << identifications.front().ms2RetentionTimeInMinutes << '\t';
        } else {
            // Append only a tab character
            oss << '\t';
        }
        oss << identifications.front().precursorCharge << '\t';
        oss << toMz(
                   identifications.front().monoIsotopicMass,
                   identifications.front().precursorCharge)
            << '\t';

        oss << intensity << '\t';

        oss << numChargeStatesObserved << '\t';
        if (isMbrPeak) {
            oss << "MBR\t";
        } else {
            oss << "MSMS\t";
        }
        oss << identifications.size() << '\t';
        oss << SplitRT << '\t';
        oss << massError << '\t';

        std::string result = oss.str();

        return result;
    }

    void resolveIdentifications() {
        std::unordered_set<std::string> distinctBaseSequences;
        std::unordered_set<std::string> distinctModifiedSequences;

        for (const Identification &id : identifications) {
            distinctBaseSequences.insert(id.sequence);
            distinctModifiedSequences.insert(id.modifications);
        }
        NumIdentificationsByBaseSeq = distinctBaseSequences.size();
        NumIdentificationsByFullSeq = distinctModifiedSequences.size();
    }

    void mergeFeatureWith(ChromatographicPeak &otherFeature, bool integrate) {
        if (&otherFeature != this) {
            // Create a set of IndexedMassSpectralPeak objects from the IsotopicEnvelopes of this feature
            unordered_set<IndexedMassSpectralPeak> thisFeaturesPeaks;
            for (const IsotopicEnvelope &envelope : isotopicEnvelopes) {
                thisFeaturesPeaks.insert(envelope.indexedPeak);
            }

            // Union and distinct Identifications
            std::vector<Identification> mergedIdentifications = identifications;
            for (const Identification &id : otherFeature.identifications) {
                if (std::find(mergedIdentifications.begin(), mergedIdentifications.end(), id) == mergedIdentifications.end()) {
                    mergedIdentifications.push_back(id);
                }
            }

            // Order the mergedIdentifications by PosteriorErrorProbability
            // This is redundant in our case since the PosteriorErrorProbability will be zero for all Identifications
            // TODO perhaps this block of code should be removed.
            std::sort(mergedIdentifications.begin(), mergedIdentifications.end(), [](const Identification &a, const Identification &b) { return a.posteriorErrorProbability < b.posteriorErrorProbability; });

            // Set the mergedIdentifications as the new Identifications for this feature
            identifications = mergedIdentifications;

            // Resolve Identifications (You should implement this method)
            resolveIdentifications();

            // Add IsotopicEnvelopes from otherFeature that are not already present in this feature
            for (const IsotopicEnvelope &envelope : otherFeature.isotopicEnvelopes) {
                if (thisFeaturesPeaks.find(envelope.indexedPeak) == thisFeaturesPeaks.end()) {
                    isotopicEnvelopes.push_back(envelope);
                }
            }

            // Recalculate intensity for this feature
            calculateIntensityForThisFeature(INTEGRATE);
        }
    }
};
}  // namespace CruxLFQ