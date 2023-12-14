#pragma once

#include <string>
#include <cmath>



namespace CruxQuant {

    class PpmTolerance{
        public:
            PpmTolerance(double value);
            ~PpmTolerance() = default;
            std::string ToString() const;
            double GetMinimumValue(double mean) const;
            double GetMaximumValue(double mean) const;
            bool Within(double experimental, double theoretical) const;
        private:
            double Value;
    };
}