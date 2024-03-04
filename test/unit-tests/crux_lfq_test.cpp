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
