#include <aws/s3/S3Client.h>
#include <aws/s3/model/PutObjectRequest.h>
#include <aws/s3/model/GetObjectRequest.h>
#include <aws/s3/model/HeadBucketRequest.h>
#include <aws/s3/model/HeadObjectRequest.h>
#include <aws/s3/model/ListObjectsRequest.h>
#include <aws/core/Aws.h>
#include <aws/core/utils/memory/stl/AWSStringStream.h>
#include <aws/core/config/AWSProfileConfigLoader.h>

#include <regex>
#include <boost/algorithm/string.hpp>
#include <stdio.h>
#include <map>

#include "io/carp.h"
#include "util/Params.h"

#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_streambuf.hpp>
//#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>

#include "AwsS3System.h"

using namespace Aws::S3;
using namespace Aws::S3::Model;
using namespace std;

const regex AwsS3System::uri_pattern{"s3://([^/]+)(((/[^/]+)*)(/([^/.]+)(\\.([^/.]*))*))(/)?$"  
        , regex::icase};
const regex AwsS3System::bucket_pattern{"s3://([^/]+)/?$", regex::icase};
const regex AwsS3System::path_pattern{"/?([^/]+(/[^/]+)*)", regex::icase};

static int const gz_magic[2] = {0x1f, 0x8b}; /* gzip magic header */


AwsS3System::AwsS3System(){
    m_systemId = SystemIdEnum::AWS_S3;

    Aws::InitAPI(m_options);
    Aws::Client::ClientConfiguration clientConfig;

    char* aws_region = getenv("AWS_DEFAULT_REGION");
    if(aws_region != NULL){
        clientConfig.region = aws_region;
        carp(CARP_INFO, "Using AWS region %s from the environment variable.", aws_region);
    }
    else{
        stringstream configFileStream("");
#if defined(_WIN32) || defined(_WIN64)
        configFileStream << getenv("HOMEDRIVE") << getenv("HOMEPATH") << "\\.aws\\config";
#endif
#ifdef __linux__
        configFileStream << getenv("HOME") << "/.aws/config";
#endif
        if(boost::filesystem::exists(configFileStream.str())){
            Aws::Config::AWSConfigFileProfileConfigLoader configLoader(configFileStream.str());
            configLoader.Load();
            auto profiles = configLoader.GetProfiles();
            string aws_profile = Params::GetString("aws-profile");
            if(profiles.find(aws_profile) != profiles.end()){
                carp(CARP_DEBUG, "Using %s AWS profile", aws_profile.c_str());
                //if region is not specified in this profile the system will 
                //quietly use default.
                //TODO: Test this assumption
                clientConfig.region = profiles[aws_profile].GetRegion();
            }
            else{
                carp(CARP_INFO, "Profile %s is not found in the AWS config file. "
                "Using AWS region %s from the default AWS profile.", 
                        profiles["default"].GetRegion().c_str());
                clientConfig.region = profiles["default"].GetRegion();
            }
        }
        else
            carp(CARP_INFO, "AWS config file not found. Using default AWS region %s.", clientConfig.region.c_str());
    }
    m_client = new S3Client(clientConfig);
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
    if(!result.IsSuccess()){
        carp(CARP_ERROR, "Cannot get AWS file info. Error type: %s\n Error message: %s.\n", 
            result.GetError().GetExceptionName().c_str(), result.GetError().GetMessage().c_str());
    }
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
        if(!response.IsSuccess()){
            carp(CARP_ERROR, "Cannot get AWS file info. Error type: %s\n Error message: %s.\n", 
                response.GetError().GetExceptionName().c_str(), response.GetError().GetMessage().c_str());
        }
        return response.IsSuccess();    //return true if such bucket exists
    } 
    else{
        ListObjectsRequest request;
        request.WithBucket(object_info.bucketName)
            .WithPrefix(object_info.key);
        auto response = m_client->ListObjects(request);
        if(!response.IsSuccess()){
            carp(CARP_ERROR, "Cannot get AWS file info. Error type: %s\n Error message: %s.\n", 
                response.GetError().GetExceptionName().c_str(), response.GetError().GetMessage().c_str());
        }
        auto results = response.GetResult().GetContents();
        switch(results.size()){
            case 0:
                //if the key is not found - it's not a dir, does not exist.
                return false;
            case 1:
                //if a single result, if the key size is the same, then it's a file
                // if it's a dir, it will have at least a trailing / or more
                return (results[0].GetKey().size() > object_info.key.size());
            default:
                return true;    //if more then one result - this is a dir.
        }
    }
}

bool AwsS3System::Mkdir(const string &path){
    //TODO_RC: Implement method
    carp(CARP_FATAL, "Writing to S3 is not supported.");
}

void AwsS3System::Rename(const string &from, const string &to){
    //TODO_RC: Implement method
    carp(CARP_FATAL, "Writing to S3 is not supported.");
}

void AwsS3System::Remove(const string &path){
    //TODO_RC: Implement method
    carp(CARP_FATAL, "Writing to S3 is not supported.");
}
/**
 * The first path must be a valid S3 path (start with s3://), the second path 
 * should be a regular path string to append to the first one.
 */ 
string AwsS3System::Join(const string &path1, const string &path2){
    if(path2.size() == 0)
        return path1;
    //validate the first path
    AwsS3System::S3ObjectInfo object_info = parseUrl(path1);
    std::stringstream path_result;
    smatch matches;
    if(regex_search(path2, matches, path_pattern)){
        string s{matches[1]}; 
        path_result << "s3://" << object_info.bucketName;
        if(!object_info.isBucket) 
            path_result << "/" << object_info.key;
        path_result << "/" << matches[1];
        return path_result.str(); 
    }
    else
        carp(CARP_FATAL, "%s is an invalid path fragment, can be added to an existing path.", path2.c_str());
    
}

std::string AwsS3System::AbsPath(const string& path){
    return path;
}

string AwsS3System::Read(const string &path){
    stringstream res_stream;
    res_stream << GetReadStream(path).rdbuf();
    return res_stream.str();
}

string AwsS3System::Read(const string &path, int byteCount){
    stringstream res_stream;
    streambuf * buf = GetReadStream(path).rdbuf();
    int i = 0;
    do {
        i++;
        res_stream << (char)buf -> sgetc();
    }while(buf-> snextc() != EOF && i < byteCount);

    return res_stream.str();
}


ostream* AwsS3System::GetWriteStream(const string &path, bool overwrite){
    //TODO_RC: Implement method. We have to accumulate the stream data in memory before sending it to S3
    carp(CARP_FATAL, "Writing to S3 is not supported.");
}

istream* AwsS3System::getReadStreamImpl(const string &path, StreamRecord &rec){
    AwsS3System::S3ObjectInfo object_info = parseUrl(path);

    if(!IsRegularFile(path))
        carp(CARP_FATAL, "%s is not a valid S3 file URL.", path.c_str());

    GetObjectRequest getObjectRequest;
    getObjectRequest.WithBucket(object_info.bucketName)
                .WithKey(object_info.key);

    //moving the return object into the free storage
    GetObjectOutcome* getObjectOutcome = new GetObjectOutcome(m_client->GetObject(getObjectRequest));

    if(getObjectOutcome->IsSuccess())
    {
        carp(CARP_DEBUG, "Successfully retrieved object %s from s3.", path.c_str());
        auto& result = (istream&)(getObjectOutcome->GetResult().GetBody());
        rec.object_pointer = (void*)getObjectOutcome;
         if(result.good()){
             //checking if this is a gzipped file
             if(result.rdbuf()->sbumpc() == gz_magic[0] && result.rdbuf()->sbumpc() == gz_magic[1]){
                result.rdbuf()->pubseekpos(0);  //rewind back to the starting position.
                //Read from the first command line argument, assume it's gzipped
                boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
                inbuf.push(boost::iostreams::gzip_decompressor());
                inbuf.push(result);
                //Convert streambuf to istream
                std::istream* instream = new std::istream(&inbuf);
                return instream;
             }
             else  { //not a gzipped stream, returning as is.
                result.rdbuf()->pubseekpos(0);  //rewind back to the starting position.
                 return &result;
             }
         }
         else{
             carp(CARP_FATAL, "Error while handling AWS data stream for the object %s", path.c_str());
         }
    }
    else
    {
        carp(CARP_FATAL, "Error while getting object %s. Error message: %s, \n%s", path.c_str(), 
                getObjectOutcome->GetError().GetExceptionName(),
                getObjectOutcome->GetError().GetMessage());
        //delete getObjectOutcome; //aborting the process, no need to clean up
    }
}

void AwsS3System::_CloseStream(const StreamRecord& p_streamRec){
    if(p_streamRec.is_input){
        delete (GetObjectOutcome*)p_streamRec.object_pointer;
    }
    else
        delete (PutObjectOutcome*)p_streamRec.object_pointer;
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
        //TODO_RC: Debug code
        vector<string> str_test;
        for(int i = 0; i < matches.size(); i++){
           str_test.push_back(matches.str(i));
           int p = matches.position(i);
        }
        // end debug code
        return  AwsS3System::S3ObjectInfo{matches};
    }    
    else
    {
        carp(CARP_FATAL, "%s is an invalid S3 URL", url.c_str());
    }
}

