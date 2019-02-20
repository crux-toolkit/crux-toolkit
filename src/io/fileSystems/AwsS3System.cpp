#include <aws/s3/S3Client.h>
#include <aws/s3/model/PutObjectRequest.h>
#include <aws/s3/model/GetObjectRequest.h>
#include <aws/s3/model/HeadBucketRequest.h>
#include <aws/s3/model/HeadObjectRequest.h>
#include <aws/s3/model/ListObjectsRequest.h>
#include <aws/core/Aws.h>
#include <aws/core/utils/memory/stl/AWSStringStream.h>

#include <regex>
#include <boost/algorithm/string.hpp>

#include "io/carp.h"

#include "AwsS3System.h"

using namespace Aws::S3;
using namespace Aws::S3::Model;
using namespace std;

const regex AwsS3System::uri_pattern{"s3://([^/]+)(((/[^/]+)*)(/([^/.]+)(\\.([^/.]+))?))(/)?$"  
        , regex::icase};
const regex AwsS3System::bucket_pattern{"s3://([^/]+)$", regex::icase};


AwsS3System::AwsS3System(){
    Aws::InitAPI(m_options);
    m_client = new S3Client();
}

bool AwsS3System::Exists(const string &path){
    return IsRegularFile(path) || IsDir(path);
}

bool AwsS3System::IsRegularFile(const string &path){
    AwsS3System::S3ObjectInfo object_info = parseUrl(path);
    if(!object_info.isFile)
        return false;
    //at this point we know that this is not a plain bucket and it does not have last /
    HeadObjectRequest request;
    request.WithBucket(object_info.bucketName)
        .WithKey(object_info.key);
    auto result = m_client->HeadObject(request);
    return result.IsSuccess();
}

/**
 * isDir in S3 means either a bucket or a prefix
 * */
bool AwsS3System::IsDir(const string &path){
    AwsS3System::S3ObjectInfo object_info = parseUrl(path);
    if(object_info.isBucket){
        HeadBucketRequest request;
        request.WithBucket(object_info.bucketName);
        auto response = m_client->HeadBucket(request);
        return response.IsSuccess();    //return true if such bucket exists
    }
    else{
        ListObjectsRequest request;
        request.WithBucket(object_info.bucketName)
            .WithPrefix(object_info.key);
        auto response = m_client->ListObjects(request);
        return (response.GetResult().GetContents().size() > 0);
    }
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
    carp(CARP_FATAL, "Writing to S3 is not supported.");
}

istream* AwsS3System::GetReadStream(const string &path){
    AwsS3System::S3ObjectInfo object_info = parseUrl(path);

    if(!IsRegularFile(path))
        carp(CARP_FATAL, "%s is not a valid S3 file URL.", path.c_str);

    GetObjectRequest getObjectRequest;
    getObjectRequest.WithBucket(object_info.bucketName)
                .WithKey(object_info.key);

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
    AwsS3System::S3ObjectInfo object_info = parseUrl(path);
    return object_info.fileName;
}

string AwsS3System::DirName(const string &path){
    AwsS3System::S3ObjectInfo object_info = parseUrl(path);
    return object_info.bucketName + object_info.prefix;
}

string AwsS3System::Stem(const string &path){
    AwsS3System::S3ObjectInfo object_info = parseUrl(path);
    return object_info.bucketName;
}

string AwsS3System::Extension(const string &path){
    AwsS3System::S3ObjectInfo object_info = parseUrl(path);
    return object_info.fileExtension;
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
    smatch matches;

    //check if this is a bucket
    if(regex_search(url, matches, AwsS3System::bucket_pattern)){
        return AwsS3System::S3ObjectInfo{string{matches[uriMatches::BUCKET]}};
    }

    //create complete path info
    if(regex_search(url, matches, AwsS3System::uri_pattern)){
        vector<string> str_test;
        for(auto& m : matches){
           str_test.push_back(string{m});
        }
        return  AwsS3System::S3ObjectInfo{matches};
    }    
    else
    {
        carp(CARP_FATAL, "%s is an invalid S3 URL", url.c_str());
    }
}

