#include "IndexedMassSpectralPeak.h"
#include <functional>

namespace CruxQuant
{
    class IsotopicEnvelope
    {

     private:
        template<typename T>
        class property {
        public:
            property(
                std::function<void(T)> setter,
                std::function<T()> getter) :
                setter_(setter),getter_(getter) { }
            operator T() const { return getter_(); }
            property<T>& operator= (const T &value) { setter_(value); return *this; }
            T& value() { return value_; }
        private:
            std::function<void(T)> setter_;
            std::function<T()> getter_;
            T value_;
        };

          IndexedMassSpectralPeak IndexedPeak;
          int ChargeState;

     public:
        IsotopicEnvelope(IndexedMassSpectralPeak monoisotopicPeak, int chargeState, double intensity)
        {
            IndexedPeak = monoisotopicPeak;
            ChargeState = chargeState;
            Intensity = intensity / chargeState;
        }

        property<double> Intensity;

        void Normalize(double normalizationFactor)
        {
            Intensity *= normalizationFactor;
        }

        string ToString()
        {
            return "+" + ChargeState + "|" + Intensity.ToString("F0") + "|" + IndexedPeak.RetentionTime.ToString("F3") + "|" + IndexedPeak.ZeroBasedMs1ScanIndex;
        }

        bool Equals(object obj)
        {
            auto otherEnv = (IsotopicEnvelope) obj;

            return otherEnv != null
                && otherEnv.ChargeState == this.ChargeState
                && otherEnv.IndexedPeak.Equals(this.IndexedPeak);
        }

        int GetHashCode()
        {
            return ChargeState.GetHashCode() + IndexedPeak.GetHashCode();
        }
    }
}
