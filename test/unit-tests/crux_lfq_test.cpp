#include <gtest/gtest.h>

#include "crux-lfq/Utils.h"

std::string spectrum_file = TEST_DATA_PATH "/test-data/test_filtered.mzML";
std::string psm_file = TEST_DATA_PATH "/test-data/assign-confidence.target.txt";

// A basic assertion which ensures gtest is setup properly.
TEST(HelloTest, BasicAssertions) {
    // Expect two strings not to be equal.
    EXPECT_STRNE("hello", "world");
    // Expect equality
    EXPECT_EQ(7 * 6, 42);
}
