#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

#include "CQStatistics.h"
#include "ChromatographicPeak.h"
#include "IndexedMassSpectralPeak.h"
#include "LFQPeptide.h"
#include "PpmTolerance.h"
#include "ProteinGroup.h"
#include "Results.h"
#include "Utils.h"

using namespace CruxLFQ;

// ============================================================
// PpmTolerance
// ============================================================

TEST(PpmToleranceTest, WithinReturnsTrueForExactMatch) {
    PpmTolerance tol(10.0);
    EXPECT_TRUE(tol.Within(500.0, 500.0));
}

TEST(PpmToleranceTest, WithinReturnsTrueInsideTolerance) {
    PpmTolerance tol(10.0);
    // 10 ppm of 500 Da = 0.005 Da
    EXPECT_TRUE(tol.Within(500.004, 500.0));
}

TEST(PpmToleranceTest, WithinReturnsFalseOutsideTolerance) {
    PpmTolerance tol(10.0);
    EXPECT_FALSE(tol.Within(500.01, 500.0));
}

TEST(PpmToleranceTest, NegativeValueTreatedAsAbsolute) {
    PpmTolerance tol(-10.0);
    EXPECT_TRUE(tol.Within(500.004, 500.0));
}

TEST(PpmToleranceTest, GetMinimumValue) {
    PpmTolerance tol(10.0);
    double min = tol.GetMinimumValue(500.0);
    EXPECT_NEAR(min, 500.0 * (1.0 - 10.0 / 1e6), 1e-9);
}

TEST(PpmToleranceTest, GetMaximumValue) {
    PpmTolerance tol(10.0);
    double max = tol.GetMaximumValue(500.0);
    EXPECT_NEAR(max, 500.0 * (1.0 + 10.0 / 1e6), 1e-9);
}

TEST(PpmToleranceTest, GetMinMaxBracketMean) {
    PpmTolerance tol(20.0);
    double mean = 1000.0;
    EXPECT_LT(tol.GetMinimumValue(mean), mean);
    EXPECT_GT(tol.GetMaximumValue(mean), mean);
}

// ============================================================
// IndexedMassSpectralPeak
// ============================================================

TEST(IndexedMassSpectralPeakTest, DefaultConstruction) {
    IndexedMassSpectralPeak peak;
    EXPECT_DOUBLE_EQ(peak.mz, 0.0);
    EXPECT_DOUBLE_EQ(peak.intensity, 0.0);
    EXPECT_EQ(peak.zeroBasedMs1ScanIndex, 0);
    EXPECT_EQ(peak.nativeScanNumber, 0);
    EXPECT_DOUBLE_EQ(peak.retentionTime, 0.0);
}

TEST(IndexedMassSpectralPeakTest, ParameterisedConstruction) {
    IndexedMassSpectralPeak peak(500.5, 12345.0, 3, 42, 15.7);
    EXPECT_DOUBLE_EQ(peak.mz, 500.5);
    EXPECT_DOUBLE_EQ(peak.intensity, 12345.0);
    EXPECT_EQ(peak.zeroBasedMs1ScanIndex, 3);
    EXPECT_EQ(peak.nativeScanNumber, 42);
    EXPECT_DOUBLE_EQ(peak.retentionTime, 15.7);
}

TEST(IndexedMassSpectralPeakTest, EqualityByMzAndScanIndex) {
    IndexedMassSpectralPeak a(500.5, 1000.0, 3, 10, 5.0);
    IndexedMassSpectralPeak b(500.5, 9999.0, 3, 99, 99.0);  // same mz + scan, different rest
    EXPECT_TRUE(a == b);
}

TEST(IndexedMassSpectralPeakTest, InequalityDifferentMz) {
    IndexedMassSpectralPeak a(500.5, 1000.0, 3, 10, 5.0);
    IndexedMassSpectralPeak b(500.6, 1000.0, 3, 10, 5.0);
    EXPECT_TRUE(a != b);
}

TEST(IndexedMassSpectralPeakTest, InequalityDifferentScanIndex) {
    IndexedMassSpectralPeak a(500.5, 1000.0, 3, 10, 5.0);
    IndexedMassSpectralPeak b(500.5, 1000.0, 4, 10, 5.0);
    EXPECT_TRUE(a != b);
}

TEST(IndexedMassSpectralPeakTest, CopyAssignment) {
    IndexedMassSpectralPeak a(300.0, 500.0, 1, 2, 3.0);
    IndexedMassSpectralPeak b;
    b = a;
    EXPECT_DOUBLE_EQ(b.mz, 300.0);
    EXPECT_DOUBLE_EQ(b.intensity, 500.0);
    EXPECT_EQ(b.zeroBasedMs1ScanIndex, 1);
    EXPECT_EQ(b.nativeScanNumber, 2);
    EXPECT_DOUBLE_EQ(b.retentionTime, 3.0);
}

// ============================================================
// CQStatistics — Median
// ============================================================

TEST(MedianTest, OddNumberOfElements) {
    std::vector<double> v = {3.0, 1.0, 2.0};
    EXPECT_DOUBLE_EQ(Median(v), 2.0);
}

TEST(MedianTest, EvenNumberOfElements) {
    std::vector<double> v = {1.0, 2.0, 3.0, 4.0};
    EXPECT_DOUBLE_EQ(Median(v), 2.5);
}

TEST(MedianTest, SingleElement) {
    std::vector<double> v = {7.0};
    EXPECT_DOUBLE_EQ(Median(v), 7.0);
}

TEST(MedianTest, AlreadySorted) {
    std::vector<double> v = {10.0, 20.0, 30.0, 40.0, 50.0};
    EXPECT_DOUBLE_EQ(Median(v), 30.0);
}

TEST(MedianTest, UnsortedLarger) {
    std::vector<double> v = {5.0, 1.0, 3.0, 9.0, 7.0};
    EXPECT_DOUBLE_EQ(Median(v), 5.0);
}

// ============================================================
// CQStatistics — Pearson
// ============================================================

TEST(PearsonTest, PerfectPositiveCorrelation) {
    std::vector<double> a = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> b = {2.0, 4.0, 6.0, 8.0};
    EXPECT_NEAR(Pearson(a, b), 1.0, 1e-9);
}

TEST(PearsonTest, PerfectNegativeCorrelation) {
    std::vector<double> a = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> b = {4.0, 3.0, 2.0, 1.0};
    EXPECT_NEAR(Pearson(a, b), -1.0, 1e-9);
}

TEST(PearsonTest, MismatchedSizesThrows) {
    std::vector<double> a = {1.0, 2.0};
    std::vector<double> b = {1.0};
    EXPECT_THROW(Pearson(a, b), std::invalid_argument);
}

// ============================================================
// IsotopicEnvelope
// ============================================================

TEST(IsotopicEnvelopeTest, IntensityDividedByCharge) {
    IndexedMassSpectralPeak peak(500.0, 1000.0, 0, 1, 1.0);
    IsotopicEnvelope env(peak, 2, 4000.0);
    EXPECT_DOUBLE_EQ(env.intensity, 2000.0);  // 4000 / 2
}

TEST(IsotopicEnvelopeTest, Normalize) {
    IndexedMassSpectralPeak peak(500.0, 1000.0, 0, 1, 1.0);
    IsotopicEnvelope env(peak, 1, 1000.0);
    env.Normalize(0.5);
    EXPECT_DOUBLE_EQ(env.intensity, 500.0);
}

TEST(IsotopicEnvelopeTest, Equality) {
    IndexedMassSpectralPeak peak(500.0, 1000.0, 0, 1, 1.0);
    IsotopicEnvelope a(peak, 2, 2000.0);
    IsotopicEnvelope b(peak, 2, 9999.0);  // same peak + charge
    EXPECT_TRUE(a == b);
}

TEST(IsotopicEnvelopeTest, Inequality_DifferentCharge) {
    IndexedMassSpectralPeak peak(500.0, 1000.0, 0, 1, 1.0);
    IsotopicEnvelope a(peak, 2, 2000.0);
    IsotopicEnvelope b(peak, 3, 2000.0);
    EXPECT_TRUE(a != b);
}

// ============================================================
// toMz / toMass round-trip
// ============================================================

TEST(MzMassConversionTest, ToMzPositiveCharge) {
    // (mass + charge * proton) / charge
    const double mass = 1000.0;
    const int charge = 2;
    double mz = toMz(mass, charge);
    EXPECT_NEAR(mz, (mass / charge) + PROTONMASS, 1e-9);
}

TEST(MzMassConversionTest, RoundTrip) {
    const double mass = 1500.0;
    const int charge = 3;
    double mz = toMz(mass, charge);
    double recovered = toMass(mz, charge);
    EXPECT_NEAR(recovered, mass, 1e-6);
}

TEST(MzMassConversionTest, ChargeOneRoundTrip) {
    const double mass = 800.0;
    double mz = toMz(mass, 1);
    EXPECT_NEAR(toMass(mz, 1), mass, 1e-9);
}

// ============================================================
// ChromatographicPeak — calculateIntensityForThisFeature
// ============================================================

namespace {
Identification makeId(const std::string& seq = "PEPTIDE",
                      const std::string& file = "file.mzML",
                      double mass = 800.0,
                      int charge = 2) {
    return Identification(seq, mass, mass, charge, file, 5.0, 1, seq);
}
}  // namespace

TEST(ChromatographicPeakTest, EmptyEnvelopesYieldsZeroIntensity) {
    auto id = makeId();
    ChromatographicPeak peak(id, false, "file.mzML");
    peak.calculateIntensityForThisFeature(false);
    EXPECT_DOUBLE_EQ(peak.intensity, 0.0);
    EXPECT_EQ(peak.numChargeStatesObserved, 0);
}

TEST(ChromatographicPeakTest, ApexIsHighestEnvelope) {
    auto id = makeId();
    ChromatographicPeak peak(id, false, "file.mzML");

    IndexedMassSpectralPeak p1(500.0, 100.0, 0, 1, 1.0);
    IndexedMassSpectralPeak p2(500.0, 500.0, 1, 2, 2.0);
    IndexedMassSpectralPeak p3(500.0, 200.0, 2, 3, 3.0);

    peak.isotopicEnvelopes.emplace_back(p1, 2, 200.0);  // intensity = 200/2 = 100
    peak.isotopicEnvelopes.emplace_back(p2, 2, 1000.0); // intensity = 1000/2 = 500
    peak.isotopicEnvelopes.emplace_back(p3, 2, 400.0);  // intensity = 400/2 = 200

    peak.calculateIntensityForThisFeature(false);

    EXPECT_DOUBLE_EQ(peak.intensity, 500.0);
    EXPECT_EQ(peak.apex.indexedPeak.zeroBasedMs1ScanIndex, 1);
}

TEST(ChromatographicPeakTest, IntegrateMode) {
    auto id = makeId();
    ChromatographicPeak peak(id, false, "file.mzML");

    IndexedMassSpectralPeak p1(500.0, 100.0, 0, 1, 1.0);
    IndexedMassSpectralPeak p2(500.0, 300.0, 1, 2, 2.0);

    peak.isotopicEnvelopes.emplace_back(p1, 1, 100.0);
    peak.isotopicEnvelopes.emplace_back(p2, 1, 300.0);

    peak.calculateIntensityForThisFeature(true);

    EXPECT_DOUBLE_EQ(peak.intensity, 400.0);
}

TEST(ChromatographicPeakTest, NumChargeStatesObserved) {
    auto id = makeId();
    ChromatographicPeak peak(id, false, "file.mzML");

    IndexedMassSpectralPeak p1(500.0, 100.0, 0, 1, 1.0);
    IndexedMassSpectralPeak p2(500.5, 200.0, 1, 2, 2.0);

    peak.isotopicEnvelopes.emplace_back(p1, 2, 200.0);
    peak.isotopicEnvelopes.emplace_back(p2, 3, 400.0);

    peak.calculateIntensityForThisFeature(false);

    EXPECT_EQ(peak.numChargeStatesObserved, 2);
}

TEST(ChromatographicPeakTest, ResolveIdentifications) {
    auto id1 = makeId("PEPTIDEA");
    auto id2 = makeId("PEPTIDEB");
    ChromatographicPeak peak(id1, false, "file.mzML");
    peak.identifications.push_back(id2);

    peak.resolveIdentifications();

    EXPECT_EQ(peak.NumIdentificationsByBaseSeq, 2);
}

TEST(ChromatographicPeakTest, MergeFeatureWithAddsEnvelopes) {
    auto id = makeId();
    ChromatographicPeak peak1(id, false, "file.mzML");
    ChromatographicPeak peak2(id, false, "file.mzML");

    IndexedMassSpectralPeak p1(500.0, 100.0, 0, 1, 1.0);
    IndexedMassSpectralPeak p2(501.0, 200.0, 1, 2, 2.0);

    peak1.isotopicEnvelopes.emplace_back(p1, 2, 200.0);
    peak2.isotopicEnvelopes.emplace_back(p2, 2, 400.0);

    peak1.mergeFeatureWith(peak2, false);

    EXPECT_EQ(static_cast<int>(peak1.isotopicEnvelopes.size()), 2);
    EXPECT_DOUBLE_EQ(peak1.intensity, 200.0);  // apex of merged = 200
}

// ============================================================
// LFQPeptides
// ============================================================

TEST(LFQPeptidesTest, SetAndGetIntensity) {
    LFQPeptides pep("PEPTIDE", "PEPTIDE", true, {});
    pep.setIntensity("file1.mzML", 1234.5);
    EXPECT_DOUBLE_EQ(pep.getIntensity("file1.mzML"), 1234.5);
}

TEST(LFQPeptidesTest, UnknownFileIntensityIsZero) {
    LFQPeptides pep("PEPTIDE", "PEPTIDE", true, {});
    EXPECT_DOUBLE_EQ(pep.getIntensity("missing.mzML"), 0.0);
}

TEST(LFQPeptidesTest, UpdateExistingIntensity) {
    LFQPeptides pep("PEPTIDE", "PEPTIDE", true, {});
    pep.setIntensity("file1.mzML", 100.0);
    pep.setIntensity("file1.mzML", 999.0);
    EXPECT_DOUBLE_EQ(pep.getIntensity("file1.mzML"), 999.0);
}

TEST(LFQPeptidesTest, SetAndGetDetectionType) {
    LFQPeptides pep("PEPTIDE", "PEPTIDE", true, {});
    pep.setDetectionType("file1.mzML", DetectionType::MSMS);
    EXPECT_EQ(pep.getDetectionType("file1.mzML"), DetectionType::MSMS);
}

TEST(LFQPeptidesTest, UnknownFileDetectionTypeIsNotDetected) {
    LFQPeptides pep("PEPTIDE", "PEPTIDE", true, {});
    EXPECT_EQ(pep.getDetectionType("missing.mzML"), DetectionType::NotDetected);
}

TEST(LFQPeptidesTest, UnambiguousQuantReturnsTrueWhenIntensityAndNonAmbiguous) {
    LFQPeptides pep("PEPTIDE", "PEPTIDE", true, {});
    pep.setIntensity("file1.mzML", 500.0);
    pep.setDetectionType("file1.mzML", DetectionType::MSMS);
    EXPECT_TRUE(pep.UnambiguousPeptideQuant());
}

TEST(LFQPeptidesTest, UnambiguousQuantReturnsFalseWhenZeroIntensity) {
    LFQPeptides pep("PEPTIDE", "PEPTIDE", true, {});
    pep.setIntensity("file1.mzML", 0.0);
    pep.setDetectionType("file1.mzML", DetectionType::MSMS);
    EXPECT_FALSE(pep.UnambiguousPeptideQuant());
}

TEST(LFQPeptidesTest, UnambiguousQuantReturnsFalseWhenAllAmbiguous) {
    LFQPeptides pep("PEPTIDE", "PEPTIDE", true, {});
    pep.setIntensity("file1.mzML", 500.0);
    pep.setDetectionType("file1.mzML", DetectionType::MSMSAmbiguousPeakfinding);
    EXPECT_FALSE(pep.UnambiguousPeptideQuant());
}

// ============================================================
// Identification — construction and equality
// ============================================================

TEST(IdentificationTest, ConstructionSetsFields) {
    Identification id("PEPTIDE", 800.4, 801.4, 2, "sample.mzML", 12.3, 42, "PEPTIDE", "protA");
    EXPECT_EQ(id.sequence, "PEPTIDE");
    EXPECT_DOUBLE_EQ(id.monoIsotopicMass, 800.4);
    EXPECT_DOUBLE_EQ(id.peptideMass, 801.4);
    EXPECT_EQ(id.precursorCharge, 2);
    EXPECT_EQ(id.spectralFile, "sample.mzML");
    EXPECT_DOUBLE_EQ(id.ms2RetentionTimeInMinutes, 12.3);
    EXPECT_EQ(id.scanId, 42);
    EXPECT_EQ(id.protein_id, "protA");
}

TEST(IdentificationTest, EqualityAndInequality) {
    Identification a("PEPTIDE", 800.4, 801.4, 2, "s.mzML", 12.3, 42, "PEPTIDE");
    Identification b("PEPTIDE", 800.4, 801.4, 2, "s.mzML", 12.3, 42, "PEPTIDE");
    Identification c("OTHER",   800.4, 801.4, 2, "s.mzML", 12.3, 42, "OTHER");
    EXPECT_TRUE(a == b);
    EXPECT_FALSE(a == c);
}

TEST(IdentificationTest, CopyConstructor) {
    Identification orig("PEPTIDE", 800.4, 801.4, 2, "s.mzML", 12.3, 42, "PEPTIDE", "protA");
    Identification copy(orig);
    EXPECT_TRUE(orig == copy);
}

// ============================================================
// CruxLFQResults — setPeptideModifiedSequencesAndProteinGroups
// ============================================================

TEST(CruxLFQResultsTest, AddsSingleIdentification) {
    CruxLFQResults results(std::vector<std::string>{"file.mzML"});
    Identification id("PEPTIDE", 800.0, 800.0, 2, "file.mzML", 5.0, 1, "PEPTIDE", "protA");
    results.setPeptideModifiedSequencesAndProteinGroups({id});

    EXPECT_EQ(results.PeptideModifiedSequences.count("PEPTIDE"), 1u);
}

TEST(CruxLFQResultsTest, DeduplicatesSameSequence) {
    CruxLFQResults results(std::vector<std::string>{"file.mzML"});
    Identification id1("PEPTIDE", 800.0, 800.0, 2, "file.mzML", 5.0, 1, "PEPTIDE");
    Identification id2("PEPTIDE", 800.0, 800.0, 2, "file.mzML", 6.0, 2, "PEPTIDE");
    results.setPeptideModifiedSequencesAndProteinGroups({id1, id2});

    EXPECT_EQ(results.PeptideModifiedSequences.size(), 1u);
}

TEST(CruxLFQResultsTest, TracksMultipleSequences) {
    CruxLFQResults results(std::vector<std::string>{"file.mzML"});
    Identification id1("PEPTIDEA", 800.0, 800.0, 2, "file.mzML", 5.0, 1, "PEPTIDEA");
    Identification id2("PEPTIDEB", 900.0, 900.0, 2, "file.mzML", 6.0, 2, "PEPTIDEB");
    results.setPeptideModifiedSequencesAndProteinGroups({id1, id2});

    EXPECT_EQ(results.PeptideModifiedSequences.size(), 2u);
}

// ============================================================
// CruxLFQResults — calculatePeptideResults
// ============================================================

TEST(CruxLFQResultsTest, UndetectedPeptideHasZeroIntensity) {
    std::string file = "file.mzML";
    CruxLFQResults results(std::vector<std::string>{file});

    Identification id("PEPTIDE", 800.0, 800.0, 2, file, 5.0, 1, "PEPTIDE");
    results.setPeptideModifiedSequencesAndProteinGroups({id});

    // No peaks added — peptide should remain at zero intensity / NotDetected
    results.calculatePeptideResults(false);

    EXPECT_DOUBLE_EQ(results.PeptideModifiedSequences["PEPTIDE"].getIntensity(file), 0.0);
    EXPECT_EQ(results.PeptideModifiedSequences["PEPTIDE"].getDetectionType(file),
              DetectionType::NotDetected);
}

TEST(CruxLFQResultsTest, DetectedPeptideGetsIntensityFromPeak) {
    std::string file = "file.mzML";
    CruxLFQResults results(std::vector<std::string>{file});

    Identification id("PEPTIDE", 800.0, 800.0, 2, file, 5.0, 1, "PEPTIDE");
    results.setPeptideModifiedSequencesAndProteinGroups({id});

    ChromatographicPeak peak(id, false, file);
    IndexedMassSpectralPeak mp(500.0, 1000.0, 0, 1, 5.0);
    peak.isotopicEnvelopes.emplace_back(mp, 2, 2000.0);
    peak.calculateIntensityForThisFeature(false);
    results.Peaks[file].push_back(peak);

    results.calculatePeptideResults(false);

    EXPECT_GT(results.PeptideModifiedSequences["PEPTIDE"].getIntensity(file), 0.0);
    EXPECT_EQ(results.PeptideModifiedSequences["PEPTIDE"].getDetectionType(file),
              DetectionType::MSMS);
}

// ============================================================
// CruxLFQResults — MedianPolish
// ============================================================

TEST(MedianPolishTest, ConvergesOnUniformMatrix) {
    // All peptides identical across samples: column effects should all be equal.
    // 2 peptides (rows 1-2), 2 samples (cols 1-2)
    std::vector<std::vector<double>> table = {
        {0.0, 0.0, 0.0},
        {0.0, 4.0, 4.0},
        {0.0, 4.0, 4.0}
    };

    CruxLFQResults results(std::vector<std::string>{});
    results.MedianPolish(table);

    // Column effects (table[0][1] and table[0][2]) should be equal
    EXPECT_NEAR(table[0][1], table[0][2], 1e-6);
}

TEST(MedianPolishTest, HandlesNaNValues) {
    // One missing value — should not crash and residuals should converge
    std::vector<std::vector<double>> table = {
        {0.0, 0.0, 0.0},
        {0.0, 4.0, std::numeric_limits<double>::quiet_NaN()},
        {0.0, 4.0, 4.0}
    };

    CruxLFQResults results(std::vector<std::string>{});
    EXPECT_NO_THROW(results.MedianPolish(table));
}

TEST(MedianPolishTest, OverallEffectExtracted) {
    // Constant matrix: overall effect should be the common value, residuals ~0
    std::vector<std::vector<double>> table = {
        {0.0, 0.0, 0.0},
        {0.0, 8.0, 8.0},
        {0.0, 8.0, 8.0}
    };

    CruxLFQResults results(std::vector<std::string>{});
    results.MedianPolish(table);

    // After polish, residuals (cells [1..][1..]) should be near zero
    for (size_t r = 1; r < table.size(); ++r)
        for (size_t c = 1; c < table[0].size(); ++c)
            if (!std::isnan(table[r][c]))
                EXPECT_NEAR(table[r][c], 0.0, 1e-4);
}

// ============================================================
// binarySearchForIndexedPeak
// ============================================================

TEST(BinarySearchTest, FindsCorrectStartIndex) {
    // Build a sorted list of peaks by zeroBasedMs1ScanIndex
    std::vector<IndexedMassSpectralPeak> peaks;
    for (int i = 0; i < 10; ++i) {
        peaks.emplace_back(500.0 + i * 0.01, 100.0, i, i + 1, static_cast<double>(i));
    }

    int idx = binarySearchForIndexedPeak(&peaks, 5);
    EXPECT_GE(idx, 0);
    EXPECT_LT(idx, static_cast<int>(peaks.size()));
}

TEST(BinarySearchTest, EmptyVectorReturnsZero) {
    std::vector<IndexedMassSpectralPeak> peaks;
    int idx = binarySearchForIndexedPeak(&peaks, 3);
    EXPECT_EQ(idx, 0);
}
