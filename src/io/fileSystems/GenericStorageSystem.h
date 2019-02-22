#ifndef GENERICSTORAGESYSTEM_H
#define GENERICSTORAGESYSTEM_H

#include <fstream>
#include <string>
#include <vector>

using namespace std;

/*
* Represents a generic storage system interface. Implementations include 
* local file system, S3, HTTP or any other URI-based resource access.
*/
class GenericStorageSystem
{
public:
  static GenericStorageSystem *getStorage( string p_path);
  static void CloseStream(ios_base& stream);

  static void cleanup();

  virtual bool Exists(const string &path) = 0;
  virtual bool IsRegularFile(const string &path) = 0;
  virtual bool IsDir(const string &path) = 0;
  virtual bool Mkdir(const string &path) = 0;
  virtual void Rename(const string &from, const string &to) = 0;
  virtual void Remove(const string &path) = 0;
  virtual string Join(const string &path1, const string &path2) = 0;
  virtual string AbsPath(const string& path) = 0;
  virtual string Read(const string &path) = 0;
  virtual ostream *GetWriteStream(const string &path, bool overwrite) = 0;
  virtual istream& GetReadStream(const string &path) = 0;
  virtual string BaseName(const string &path) = 0;
  virtual string DirName(const string &path) = 0;
  virtual string Stem(const string &path) = 0;
  virtual string Extension(const string &path) = 0;
  virtual void CopyLocal(const string &orig, const string &dest) = 0;

  void Copy(const string &orig, const string &dest);

  static constexpr int MAX_STORAGES = 2;
  enum SystemIdEnum
  {
    FILE_SYSTEM = 0,
    AWS_S3 = 1
  };

  struct StreamRecord{
      SystemIdEnum system_id;
      bool is_input;       //true if input stream, false is output.
      ios_base* stream_pointer;
      void* object_pointer;    //[RC]this is brute force, but AWS library should not be imported into this scope
      StreamRecord(SystemIdEnum p_system_id, bool p_is_input, ios_base* p_stream_pointer, void* p_object_pointer){
        system_id = p_system_id;
        is_input = p_is_input;
        stream_pointer = p_stream_pointer;
        object_pointer = p_object_pointer;
      };
  };

protected:
  void _RegisterStream(const StreamRecord& p_streamRec);

  static GenericStorageSystem *m_system[MAX_STORAGES];

  static std::vector<StreamRecord> m_openStreams; 



  GenericStorageSystem(){};
  
  virtual ~GenericStorageSystem(){};

private:
  static int getStorageIndex(const string &path);
  virtual void _CloseStream(const StreamRecord& p_streamRec) = 0;


};

#endif