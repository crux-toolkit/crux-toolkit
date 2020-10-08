#ifndef GENERICSTORAGESYSTEM_H
#define GENERICSTORAGESYSTEM_H

#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <mutex>
#include <iostream>


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
  virtual string Read(const string &path, int byteCount) = 0;
  virtual ostream *GetWriteStream(const string &path, bool overwrite) = 0;
  istream& GetReadStream(const string &path);
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

  struct Identifiable{
    int streamId;
    static int globalStreamCounter;


     Identifiable() : streamId{++globalStreamCounter} {};
  };

  class IdentifiableInputStream: public istream, public Identifiable{
    public:
      IdentifiableInputStream(istream& pStream) : istream{pStream.rdbuf()} {}
  };

  class IdentifiableOutputStream : public ostream, public Identifiable{
    public:
      IdentifiableOutputStream(ostream& pStream) : ostream{pStream.rdbuf()}  {};
  };


  struct StreamRecord{
      SystemIdEnum system_id;
      bool is_input;           //true if input stream, false is output.
      int stream_id;           //unique integer identifying this stream
      void* object_pointer;    //[RC]this is brute force, but AWS library should not be imported into this scope
      string path;

      StreamRecord(SystemIdEnum p_system_id, bool p_is_input, int p_stream_id, void* p_object_pointer, const string& p_path): system_id{p_system_id}, path{p_path} {
        is_input = p_is_input;
        stream_id = p_stream_id;
        object_pointer = p_object_pointer;
      };
      bool operator==(IdentifiableInputStream &stream) const 
      { 
        return this->stream_id == stream.streamId;
      }
      bool operator==(const StreamRecord &rec) const
      { 
        return (this->stream_id == rec.stream_id);
      }

      string toString() const {
        stringstream strStreamInfo{"Stream: "};
        strStreamInfo << "storage system: " <<  system_id;
        strStreamInfo << "; is input: " << is_input;
        strStreamInfo << "; stream id: " << stream_id;
        strStreamInfo << "; file_path: " << path << "\n";
        return strStreamInfo.str();
      }
  };

protected:
  void _RegisterStream(const StreamRecord& p_streamRec);
  virtual istream* getReadStreamImpl(const string &path, StreamRecord& rec) = 0; 

  static GenericStorageSystem *m_system[MAX_STORAGES];

  static std::map<int, StreamRecord> m_openStreams; 
  static std::mutex m_registerMutex;

  SystemIdEnum m_systemId;



  GenericStorageSystem(){};
  
  virtual ~GenericStorageSystem(){};

private:
  static int getStorageIndex(const string &path);
  virtual void _CloseStream(const StreamRecord& p_streamRec) = 0;


};

#endif