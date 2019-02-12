#ifndef GENERICSTORAGESYSTEM_H
#define GENERICSTORAGESYSTEM_H

#include <fstream>
#include <string>

using namespace std;

/*
* Represents a generic storage system interface. Implementations include 
* local file system, S3, HTTP or any other URI-based resource access.
*/
class GenericStorageSystem
{
public:
  static GenericStorageSystem *getStorage( string p_path);
  static void cleanup();

  virtual bool Exists(const string &path) = 0;
  virtual bool IsRegularFile(const string &path) = 0;
  virtual bool IsDir(const string &path) = 0;
  virtual bool Mkdir(const string &path) = 0;
  virtual void Rename(const string &from, const string &to) = 0;
  virtual void Remove(const string &path) = 0;
  virtual string Join(const string &path1, const string &path2) = 0;
  virtual string Read(const string &path) = 0;
  virtual ostream *GetWriteStream(const string &path, bool overwrite) = 0;
  virtual istream *GetReadStream(const string &path) = 0;
  virtual string BaseName(const string &path) = 0;
  virtual string DirName(const string &path) = 0;
  virtual string Stem(const string &path) = 0;
  virtual string Extension(const string &path) = 0;
  virtual void CopyLocal(const string &orig, const string &dest) = 0;

  void Copy(const string &orig, const string &dest);

  static constexpr int MAX_STORAGES = 2;
  enum
  {
    FILE_SYSTEM = 0,
    AWS_S3 = 1
  };

protected:
  static GenericStorageSystem *m_system[MAX_STORAGES];

  GenericStorageSystem(){};
  
  virtual ~GenericStorageSystem(){};

private:
  static int getStorageIndex(const string &path);
};

#endif