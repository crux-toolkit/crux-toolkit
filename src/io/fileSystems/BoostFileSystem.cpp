#include "BoostFileSystem.h"

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
istream* BoostFileSystem::GetReadStream(const string &path){
  ifstream* stream = new ifstream(path);
  if(!stream->good()){
    delete stream;
    return NULL;
  }
  return stream;
}
