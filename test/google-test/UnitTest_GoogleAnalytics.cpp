#include <gtest/gtest.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <regex>
#include "util/FileUtils.h"
#include "util/crux-utils.h"

using namespace std;

#define PATH_SEPARATOR "/"

class GoogleAnalyticsTest : public ::testing::Test {
  protected:
    const string cruxHomeDir{".crux"};
    const string cruxClientIdFileName{"client_id"};
    const string testUUID{"123e4567-e89b-12d3-a456-426614174000"};
    const regex uuid_regex{"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}"};
    stringstream filePath;

    virtual void SetUp(){
        setenv("HOME", ".", 1);
        filePath << std::getenv("HOME") << PATH_SEPARATOR << cruxHomeDir;
        if(FileUtils::Exists(filePath.str())){
            FileUtils::Remove(filePath.str());
        }
        filePath  << PATH_SEPARATOR << cruxClientIdFileName;
    };
    virtual void TearDown(){};
};

TEST_F(GoogleAnalyticsTest, TestUUID){
    auto uuid = ensureClientId();
    ASSERT_TRUE(FileUtils::Exists(filePath.str()));
    EXPECT_TRUE(regex_match(uuid, uuid_regex));
    FileUtils::Remove(filePath.str());
    auto fileStream = FileUtils::GetWriteStream(filePath.str(), true);
    *fileStream  << testUUID;
    fileStream->close();
    uuid = ensureClientId();
    EXPECT_EQ(uuid, testUUID);
}
