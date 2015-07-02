#ifndef TIDEREGRESSIONTEST_H
#define TIDEREGRESSIONTEST_H

#include <map>
#include <string>

struct TideRegressionSettings;

static std::map<std::string, TideRegressionSettings> settings_;
static std::string tideIndex_;
static std::string tideSearch_;
static std::string crux_;
static std::string fasta_;
static std::string spectrumRecords_;

static void initSettings();
static void addTest(const std::string& name, const TideRegressionSettings& settings);
int main(int argc, char** argv);
static void setPaths(const char* cruxPath);
static void runTest(std::map<std::string, TideRegressionSettings>::const_iterator i);
static void listTests();
static void printUsage(const char* argv0);

#endif

