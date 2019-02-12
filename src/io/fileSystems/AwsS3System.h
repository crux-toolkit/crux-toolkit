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

  static const regex uri_pattern;

public:

  bool Exists(const string &path);
  bool IsRegularFile(const string &path);
  bool IsDir(const string &path);
  bool Mkdir(const string &path);
  void Rename(const string &from, const string &to);
  void Remove(const string &path);
  string Join(const string &path1, const string &path2);
  string Read(const string &path);
  ostream *GetWriteStream(const string &path, bool overwrite);
  istream *GetReadStream(const string &path);
  string BaseName(const string &path);
  string DirName(const string &path);
  string Stem(const string &path);
  string Extension(const string &path);
  void CopyLocal(const string &orig, const string &dest);

  AwsS3System();
  ~AwsS3System();

private:
  struct S3ObjectInfo{
    string bucketName;
    string key;
    string fileName;
    string fileExtension;
    bool isFile;
    bool isBucket;
  };

  S3ObjectInfo parseUrl(const string& url);
};


#endif