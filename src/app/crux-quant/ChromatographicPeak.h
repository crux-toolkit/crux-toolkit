#pragma once

#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>

#include "Utils.h"

using std::string;
using std::vector;

namespace CruxQuant {

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

    ChromatographicPeak(const Identification& id, bool _isMbrPeak, const string& _spectralFile)
        : isMbrPeak(_isMbrPeak), spectralFile(_spectralFile) {
        isotopicEnvelopes.push_back(IsotopicEnvelope());
        identifications.push_back(id);
        massError = std::numeric_limits<double>::quiet_NaN();
        numChargeStatesObserved = 0;
        _index = std::numeric_limits<int>::max();
        SplitRT = 0;
    }

    void calculateIntensityForThisFeature(bool integrate) {
        if (!isotopicEnvelopes.empty()) {
            // Find the IsotopicEnvelope with the maximum intensity (apex)
            auto maxIntensityEnvelopes = std::max_element(
                isotopicEnvelopes.begin(), isotopicEnvelopes.end(),
                [](const IsotopicEnvelope& a, const IsotopicEnvelope& b) {
                    return a.intensity < b.intensity;
                });

            if (maxIntensityEnvelopes != isotopicEnvelopes.end()) {
                apex = *maxIntensityEnvelopes;

                if (integrate) {
                    // Calculate intensity by summing up all IsotopicEnvelope intensities
                    intensity = std::accumulate(
                        isotopicEnvelopes.begin(), isotopicEnvelopes.end(), 0.0,
                        [](double sum, const IsotopicEnvelope& envelope) {
                            return sum + envelope.intensity;
                        });
                } else {
                    // Set intensity to the apex intensity
                    intensity = apex.intensity;
                }

                // Calculate massError for each Identification
                massError = std::numeric_limits<double>::quiet_NaN();  // Set massError to NaN initially

                for (const Identification& id : identifications) {
                    double massErrorForId = ((toMass(apex.indexedPeak.mz, apex.chargeState) - id.peakFindingMass) / id.peakFindingMass) * 1e6;

                    if (std::isnan(massError) || std::abs(massErrorForId) < std::abs(massError)) {
                        massError = massErrorForId;
                    }
                }
                std::unordered_set<int> uniqueChargeStates;
                for (const IsotopicEnvelope& envelope : isotopicEnvelopes) {
                    uniqueChargeStates.insert(envelope.chargeState);
                }
                numChargeStatesObserved = static_cast<int>(uniqueChargeStates.size());
            }
        } else {
            // No isotopicEnvelopes, so set intensity to 0 and massError to NaN
            intensity = 0.0;
            massError = std::numeric_limits<double>::quiet_NaN();
            numChargeStatesObserved = 0;
            apex = {};  // Reset apex to default-constructed IsotopicEnvelope
        }
    }

    string ToString(){
        std::set<string> uniqueBaseSequences;
        std::set<string> uniqueModifiedSequences;

        for (const Identification& id : identifications) {
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
            identifications.front().precursorCharge
        ) << '\t';

        oss << intensity << '\t';

        oss << numChargeStatesObserved << '\t';
        oss << identifications.size() << '\t';
        oss << SplitRT << '\t';
        oss << massError << '\t';

        std::string result = oss.str();

        return result;
    }
};
}  // namespace CruxQuant