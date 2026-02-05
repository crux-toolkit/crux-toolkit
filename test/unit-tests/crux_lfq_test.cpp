#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "CruxLFQApplication.h"
#include "crux-lfq/Results.h"
#include "crux-lfq/Utils.h"

using std::string;
using std::vector;

string spectrum_file;
string psm_file_percolator = TEST_DATA_PATH "/test-data/percolator.target.psms.txt";

// A basic assertion which ensures gtest is setup properly.
TEST(HelloTest, BasicAssertions) {
    if (spectrum_file.empty()) {
        GTEST_SKIP() << "No spectrum file provided";
    }
    // Expect two strings not to be equal.
    EXPECT_STRNE("hello", "world");
    // Expect equality
    EXPECT_EQ(7 * 6, 42);
}

TEST(CruxLFQTest, CreatePSMs) {
    vector<PSM> psms = CruxLFQApplication::create_percolator_psm(psm_file_percolator);
    // Test that PSMs are created successfully
    EXPECT_FALSE(psms.empty());
}

TEST(CruxLFQTest, CreateIdentifications) {
    vector<PSM> psms = CruxLFQApplication::create_percolator_psm(psm_file_percolator);
    vector<CruxLFQ::Identification> allIdentifications = CruxLFQApplication::createIdentifications(psms, spectrum_file);
    ASSERT_EQ(4, allIdentifications.size());
    ASSERT_EQ("LLALNSLYSPK", allIdentifications[0].sequence);
}

TEST(CruxLFQTest, calculateTheoreticalIsotopeDistributions) {
    vector<PSM> psms = CruxLFQApplication::create_percolator_psm(psm_file_percolator);
    vector<CruxLFQ::Identification> allIdentifications = CruxLFQApplication::createIdentifications(psms, spectrum_file);
    auto distributions = CruxLFQ::calculateTheoreticalIsotopeDistributions(allIdentifications);
    ASSERT_FALSE(distributions.empty());

    // Write distributions to file
    std::ofstream outFile("theoretical_isotope_distributions.tsv");
    outFile << "Sequence\tMass\tIntensity\n";

    for (const auto& entry : distributions) {
        const std::string& sequence = entry.first;
        const std::vector<std::pair<double, double>>& isotopes = entry.second;

        for (const auto& isotope : isotopes) {
            double mass = isotope.first;
            double intensity = isotope.second;
            outFile << sequence << "\t" << mass << "\t" << intensity << "\n";
        }
    }

    outFile.close();
    std::cout << "Theoretical isotope distributions written to theoretical_isotope_distributions.tsv\n";
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    if (argc > 1) {
        spectrum_file = argv[1];
    } else {
        // Set a default for test discovery, but skip tests that need the file
        spectrum_file = "";
        std::cerr << "Warning: No spectrum file provided. Some tests will be skipped." << std::endl;
    }

    return RUN_ALL_TESTS();
}