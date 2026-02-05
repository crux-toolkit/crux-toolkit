#include "PpmTolerance.h"

namespace CruxLFQ {
    PpmTolerance::PpmTolerance(double value) : Value(std::abs(value)) {}

    std::string PpmTolerance::ToString() const {
        return "±" + std::to_string(Value) + " PPM";
    }

    double PpmTolerance::GetMinimumValue(double mean) const {
        return mean * (1 - (Value / 1e6));
    }

    double PpmTolerance::GetMaximumValue(double mean) const {
        return mean * (1 + (Value / 1e6));
    }

    bool PpmTolerance::Within(double experimental, double theoretical) const {
        return std::abs((experimental - theoretical) / theoretical * 1e6) <= Value;
    }
}