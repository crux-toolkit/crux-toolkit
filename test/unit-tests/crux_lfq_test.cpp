#include <gtest/gtest.h>

#include <vector>

#include "CruxLFQApplication.h"
#include "crux-lfq/Results.h"
#include "crux-lfq/Utils.h"

using std::string;
using std::vector;

string spectrum_file = "C:\\Users\\acqua\\Code\\data\\PXD005590\\B02_12_161103_D3_HCD_OT_4ul.raw.mzXML";
string psm_file_tide_search = TEST_DATA_PATH "/test-data/tide-search.txt";
string psm_file_assign_confidence = TEST_DATA_PATH "/test-data/assign-confidence.target.txt";

// A basic assertion which ensures gtest is setup properly.
TEST(HelloTest, BasicAssertions) {
    // Expect two strings not to be equal.
    EXPECT_STRNE("hello", "world");
    // Expect equality
    EXPECT_EQ(7 * 6, 42);
}

// This shows a fundamental flaw in the code that must be fixed
TEST(CreatePsmMapTest, TideSearchFormat) {
    // Test case for psm_file_format = "tide-search"
    vector<CruxLFQ::PSM> psm = CruxLFQ::create_psm(psm_file_tide_search, "tide-search");
    EXPECT_EQ(psm.size(), 15);
}

TEST(CreatePsmMapTest, AssignConfidenceFormat) {
    // Test case for psm_file_format = "assign-confidence"
    vector<CruxLFQ::PSM> psm = CruxLFQ::create_psm(psm_file_assign_confidence, "assign-confidence");
    EXPECT_EQ(psm.size(), 15);

    psm = CruxLFQ::create_psm(psm_file_assign_confidence, "assign-confidence", true);
    EXPECT_EQ(psm.size(), 30);
}

TEST(CruxLFQApplicationTest, CreateIdentifications) {
    vector<CruxLFQ::PSM> psm = CruxLFQ::create_psm(psm_file_assign_confidence, "assign-confidence");
    vector<CruxLFQ::Identification> allIdentifications = CruxLFQApplication::createIdentifications(psm, spectrum_file);
    ASSERT_EQ(15, allIdentifications.size());
    ASSERT_EQ("RPQYSNPPVQGEVMEGADNQGAGEQGRPVR", allIdentifications[0].sequence);
}

TEST(CruxLFQApplicationTest, SetPeptideModifiedSequencesAndProteinGroups) {
    vector<CruxLFQ::PSM> psm = CruxLFQ::create_psm(psm_file_assign_confidence, "assign-confidence");
    vector<CruxLFQ::Identification> allIdentifications = CruxLFQApplication::createIdentifications(psm, spectrum_file);
    vector<CruxLFQ::Identification> single_identification = vector<CruxLFQ::Identification>(allIdentifications.begin(), allIdentifications.begin() + 1);
    vector<string> spec_files = {spectrum_file};
    CruxLFQ::CruxLFQResults lfqResults(spec_files);
    vector<CruxLFQ::Identification> empty_identifications;
    lfqResults.setPeptideModifiedSequencesAndProteinGroups(empty_identifications);
    ASSERT_EQ(lfqResults.PeptideModifiedSequences.size(), 0);
    lfqResults.setPeptideModifiedSequencesAndProteinGroups(single_identification);
    ASSERT_EQ(lfqResults.PeptideModifiedSequences.size(), 1);
    lfqResults.setPeptideModifiedSequencesAndProteinGroups(allIdentifications);
    ASSERT_EQ(lfqResults.PeptideModifiedSequences.size(), 10);
}

TEST(CalculateTheoreticalIsotopeDistributionsTest, ActualIdentifications) {
    vector<CruxLFQ::PSM> psm = CruxLFQ::create_psm(psm_file_assign_confidence, "assign-confidence");
    vector<CruxLFQ::Identification> allIdentifications = CruxLFQApplication::createIdentifications(psm, spectrum_file);

    // Seems redundant but checked to ensure the right data is passed to calculateTheoreticalIsotopeDistributions
    ASSERT_EQ(15, allIdentifications.size());
    ASSERT_EQ("RPQYSNPPVQGEVMEGADNQGAGEQGRPVR", allIdentifications[0].sequence);

    unordered_map<string, vector<pair<double, double>>> modifiedSequenceToIsotopicDistribution = CruxLFQ::calculateTheoreticalIsotopeDistributions(allIdentifications);
    ASSERT_EQ(modifiedSequenceToIsotopicDistribution.size(), 10);
}

TEST(CalculateTheoreticalIsotopeDistributionsTest, EmptyInput) {
    vector<Identification> allIdentifications;
    unordered_map<string, vector<pair<double, double>>> result = calculateTheoreticalIsotopeDistributions(allIdentifications);
    EXPECT_TRUE(result.empty());
}

TEST(CalculateTheoreticalIsotopeDistributionsTest, SingleIdentification) {
    vector<Identification> allIdentifications;
    Identification id;
    id.sequence = "ACDEFG";
    id.precursorCharge = 2;
    allIdentifications.push_back(id);
    unordered_map<string, vector<pair<double, double>>> result = calculateTheoreticalIsotopeDistributions(allIdentifications);
    ASSERT_EQ(result.size(), 1);
    // TODO - Figure out why this code passes
    EXPECT_EQ(result["ACDEFG"].size(), 3);  // Assuming NUM_ISOTOPES_REQUIRED is 3 
}