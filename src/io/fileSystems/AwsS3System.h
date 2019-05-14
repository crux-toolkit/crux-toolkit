#ifndef AWSS3SYSTEM_H
#define AWSS3SYSTEM_H
#include <aws/s3/S3Client.h>
#include <aws/core/Aws.h>

#include <regex>

#include "GenericStorageSystem.h"

using namespace Aws::S3;
using namespace std;

class AwsS3System : public GenericStorageSystem{

private:
  Aws::SDKOptions m_options;
  Aws::S3::S3Client* m_client;
  std::vector<Aws::S3::Model::GetObjectOutcome*> m_inputStreams;
  vector<Aws::S3::Model::PutObjectOutcome*> m_outputStreams;

  static const regex uri_pattern;
  static const regex bucket_pattern;
  static const regex path_pattern;

public:

  bool Exists(const string &path);
  bool IsRegularFile(const string &path);
  bool IsDir(const string &path);
  bool Mkdir(const string &path);
  void Rename(const string &from, const string &to);
  void Remove(const string &path);
  string Join(const string &path1, const string &path2);
  string AbsPath(const string& path);
  string Read(const string &path);
  string Read(const string &path, int byteCount);
  ostream *GetWriteStream(const string &path, bool overwrite);
  istream& GetReadStream(const string &path);
  string BaseName(const string &path);
  string DirName(const string &path);
  string Stem(const string &path);
  string Extension(const string &path);
  void CopyLocal(const string &orig, const string &dest);

  AwsS3System();
  ~AwsS3System();

//these positions are defined by the structure of uri_pattern regexp
enum uriMatches{BUCKET=1, FULL_KEY=2, PREFIX=3, FULL_FILE_NAME=5, BASE_NAME=6, EXTENSION=8, LAST_SLASH=9};

private:
  struct S3ObjectInfo{
    string bucketName;
    string key;
    string prefix;
    string fileName;
    string fileExtension;
    string fullFileName;    //filename.ext
    bool isFile;
    bool isBucket;

    S3ObjectInfo(const string& p_bucket): bucketName{p_bucket}, isBucket{true}{};
    S3ObjectInfo(const smatch& p_matches){
      bucketName = p_matches[uriMatches::BUCKET];
      string key_str{p_matches[uriMatches::FULL_KEY]};
      if(key_str.size() > 0)
        key = key_str.substr(1);
      string prefix_str{p_matches[uriMatches::PREFIX]};
      if(prefix_str.size() > 0)
        prefix = prefix_str.substr(1);
      int size_ = p_matches.position(uriMatches::EXTENSION) - p_matches.position(uriMatches::BASE_NAME) - 1;
      fileName = p_matches.str(0).substr( p_matches.position(uriMatches::BASE_NAME), size_);
      
      fileExtension = p_matches[uriMatches::EXTENSION];
      fullFileName = string{p_matches[FULL_FILE_NAME]}.substr(1);
      isFile = (p_matches[uriMatches::LAST_SLASH] == "");
      isBucket  = false;
    }
  };

  S3ObjectInfo parseUrl(const string& url);
  void _CloseStream(const StreamRecord& p_streamRec);  

};


#endif