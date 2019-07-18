#include "BoostFileSystem.h"
#include "io/carp.h"
#include "pwiz/utility/misc/random_access_compressed_ifstream.hpp"

#include <ios>
#include <fstream>

using namespace std;

string BoostFileSystem::Read(const string& path) {
  return BoostFileSystem::Read(path, -1);
}

string BoostFileSystem::Read(const std::string &path, int byteCount){
  try {
    ifstream stream(path.c_str());
    stream.seekg(0, ios::end);
    size_t size = stream.tellg();
    int bytesToRead = byteCount > size ? size : byteCount;
    string buffer(bytesToRead, '\0');
    stream.seekg(0);
    stream.read(&buffer[0], bytesToRead);
    return buffer;
  } catch (...) {
    throw runtime_error("Error reading file '" + path + "'");
  }
}
  
ostream* BoostFileSystem::GetWriteStream(const string& path, bool overwrite) {
  if (Exists(path) && !overwrite) {
    return NULL;
  }
  ostream* stream = (std::ostream*) new ofstream(path.c_str());
  if (! stream -> good()) {
    delete stream;
    return NULL;
  }
  return stream;
}

istream* BoostFileSystem::getReadStreamImpl(const string &path, StreamRecord &rec){
  istream* stream = new pwiz::util::random_access_compressed_ifstream(path.c_str());
  if(!stream->good()){
    delete stream;
    carp(CARP_FATAL, "Cannot open input stream for the file %s", path.c_str());
  }
  rec.object_pointer = (void*)stream;

  return stream;
}
