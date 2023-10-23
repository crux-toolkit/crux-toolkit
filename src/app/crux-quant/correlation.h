#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

namespace CruxQuant {
double Pearson(const std::vector<double>& dataA, const std::vector<double>& dataB) {
    int n = 0;
    double r = 0.0;

    double meanA = 0;
    double meanB = 0;
    double varA = 0;
    double varB = 0;

    if (dataA.size() != dataB.size()) {
        throw std::invalid_argument("The vectors must have the same length.");
    }

    for (size_t i = 0; i < dataA.size(); ++i) {
        double currentA = dataA[i];
        double currentB = dataB[i];

        double deltaA = currentA - meanA;
        double scaleDeltaA = deltaA / ++n;

        double deltaB = currentB - meanB;
        double scaleDeltaB = deltaB / n;

        meanA += scaleDeltaA;
        meanB += scaleDeltaB;

        varA += scaleDeltaA * deltaA * (n - 1);
        varB += scaleDeltaB * deltaB * (n - 1);
        r += (deltaA * deltaB * (n - 1)) / n;
    }

    return r / std::sqrt(varA * varB);
}
}  // namespace CruxQuant
