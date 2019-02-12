#include <aws/s3/S3Client.h>
#include <aws/s3/model/PutObjectRequest.h>
#include <aws/s3/model/GetObjectRequest.h>
#include <aws/core/Aws.h>
#include <aws/core/utils/memory/stl/AWSStringStream.h>

#include <regex>
#include <boost/algorithm/string.hpp>

#include "io/carp.h"

#include "AwsS3System.h"

using namespace Aws::S3;
using namespace Aws::S3::Model;
using namespace std;

const regex AwsS3System::uri_pattern{"s3:/(/[^/]+)((/[^/]+)*(([^/]+)(\\.[^/.]+)?))/?$"};


AwsS3System::AwsS3System(){
    Aws::InitAPI(m_options);
    m_client = new S3Client();
}

bool AwsS3System::Exists(const string &path){
    AwsS3System::S3ObjectInfo object_info = parseUrl(path);
    return false;
}

bool AwsS3System::IsRegularFile(const string &path){
    //TODO_RC: Implement method
}

bool AwsS3System::IsDir(const string &path){
    //TODO_RC: Implement method
}

bool AwsS3System::Mkdir(const string &path){
    //TODO_RC: Implement method
}

void AwsS3System::Rename(const string &from, const string &to){
    //TODO_RC: Implement method
}

void AwsS3System::Remove(const string &path){
    //TODO_RC: Implement method
}

string AwsS3System::Join(const string &path1, const string &path2){
    //TODO_RC: Implement method
}

string AwsS3System::Read(const string &path){
    //TODO_RC: Implement method
}

ostream* AwsS3System::GetWriteStream(const string &path, bool overwrite){
    //TODO_RC: Implement method. We have to accumulate the stream data in memory before sending it to S3
}

istream* AwsS3System::GetReadStream(const string &path){
    Aws::String BUCKET{"proteo"};
    Aws::String KEY{"folder/file.gz"};
    GetObjectRequest getObjectRequest;
    getObjectRequest.WithBucket(BUCKET)
                .WithKey(KEY);

    auto getObjectOutcome = m_client->GetObject(getObjectRequest);

    if(getObjectOutcome.IsSuccess())
    {
        carp(CARP_DEBUG, "Successfully retrieved object %s from s3.", path);
        auto& result = (istream&)(getObjectOutcome.GetResult().GetBody());
        return &result;
    }
    else
    {
            carp(CARP_FATAL, "Error while getting object %s. Error message: %s, \n%s", path, 
                    getObjectOutcome.GetError().GetExceptionName(),
                    getObjectOutcome.GetError().GetMessage());
    }
}

string AwsS3System::BaseName(const string &path){
    
    //TODO_RC: Implement method
}

string AwsS3System::DirName(const string &path){
    //TODO_RC: Implement method
}

string AwsS3System::Stem(const string &path){
    //TODO_RC: Implement method
}

string AwsS3System::Extension(const string &path){
    //TODO_RC: Implement method
}

void AwsS3System::CopyLocal(const string &orig, const string &dest){
    //TODO_RC: Implement method
}


AwsS3System::~AwsS3System(){
    if(m_client)
        delete m_client;
    Aws::ShutdownAPI(m_options);
}


AwsS3System::S3ObjectInfo AwsS3System::parseUrl(const string& url){
    string lower_path = boost::algorithm::to_lower_copy(url);
    smatch matches;
    AwsS3System::S3ObjectInfo result;

    //extract protocal name from the URI path
    if(regex_search(lower_path, matches, AwsS3System::uri_pattern)){
        result.bucketName = matches[1];
        result.key = matches[2];
        result.fileName = matches[3];
        vector<string> str_test;
        for(auto& m : matches){
           str_test.push_back(string{m});
        }
        return result;
    }    
    else
    {
        carp(CARP_FATAL, "%s is an invalid S3 URL", &url);
    }

}