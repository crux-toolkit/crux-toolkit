#include "BoostFileSystem.h"
#include "io/carp.h"

#include <fstream>

using namespace std;

string BoostFileSystem::Read(const string& path) {
  try {
    ifstream stream(path.c_str());
    stream.seekg(0, ios::end);
    size_t size = stream.tellg();
    string buffer(size, '\0');
    stream.seekg(0);
    stream.read(&buffer[0], size);
    return buffer;
  } catch (...) {
    throw runtime_error("Error reading file '" + path + "'");
  }
}

ostream* BoostFileSystem::GetWriteStream(const string& path, bool overwrite) {
  if (Exists(path) && !overwrite) {
    return NULL;
  }
  ofstream* stream = new ofstream(path.c_str());
  if (!stream->good()) {
    delete stream;
    return NULL;
  }
  return stream;
}
istream& BoostFileSystem::GetReadStream(const string &path){
  ifstream* stream = new ifstream(path);
  if(!stream->good()){
    delete stream;
    carp(CARP_FATAL, "Cannot open input stream for the file %s", path.c_str());
  }
  StreamRecord sr{GenericStorageSystem::SystemIdEnum::FILE_SYSTEM, true, 
      (ios_base*)stream, (void*)stream};
  _RegisterStream(sr);

  return *stream;
}
