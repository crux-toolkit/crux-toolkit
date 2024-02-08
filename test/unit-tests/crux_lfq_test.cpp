#include <gtest/gtest.h>

#include <vector>

#include "crux-lfq/Utils.h"

using std::string;
using std::vector;

string spectrum_files = "C:\\Users\\acqua\\Code\\data\\PXD005590";
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