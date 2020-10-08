#include "FileUtils.h"
#include "boost/filesystem.hpp"
#include <fstream>
#include <stdexcept>
#include "io/fileSystems/GenericStorageSystem.h"
#include "io/carp.h"

using namespace std;

bool FileUtils::Exists(const string& path) {
  return GenericStorageSystem::getStorage(path)->Exists(path);
  //return boost::filesystem::exists(path);
} 

bool FileUtils::IsRegularFile(const string& path) {
  return  GenericStorageSystem::getStorage(path)->IsRegularFile(path);
}

bool FileUtils::IsDir(const string& path) {
  return GenericStorageSystem::getStorage(path)->IsDir(path);
}

// returns true if a new directory was created, otherwise false
bool FileUtils::Mkdir(const string& path) {
  return GenericStorageSystem::getStorage(path)->Mkdir(path);;
}

void FileUtils::Rename(const string& from, const string& to) {
  GenericStorageSystem* s_from = GenericStorageSystem::getStorage(from);
  GenericStorageSystem* s_to = GenericStorageSystem::getStorage(to);
  if(s_from != s_to)
    carp(CARP_FATAL, "Rename operation must have both URIs from the same storage system. Got from: %s and to: %s.", from, to);
  if (s_from->Exists(from)) 
    s_from->Rename(from, to);
}

void FileUtils::Remove(const string& path) {
  GenericStorageSystem* s = GenericStorageSystem::getStorage(path);
  if (s->Exists(path)) {
    s->Remove(path);
  }
}

string FileUtils::Join(const string& path1, const string& path2) {
  GenericStorageSystem* s = GenericStorageSystem::getStorage(path1);
  return s->Join(path1, path2);
}

string FileUtils::AbsPath(const string& path){
  GenericStorageSystem* s = GenericStorageSystem::getStorage(path);
  return s->AbsPath(path);
}

string FileUtils::Read(const string& path) {
  GenericStorageSystem* s = GenericStorageSystem::getStorage(path);
  return s->Read(path);
}

string FileUtils::Read(const string& path, int byteCount) {
  GenericStorageSystem* s = GenericStorageSystem::getStorage(path);
  return s->Read(path);
}

ostream* FileUtils::GetWriteStream(const string& path, bool overwrite) {
  GenericStorageSystem* s = GenericStorageSystem::getStorage(path);
  ostream* res = s->GetWriteStream(path, overwrite);
  if(!res->good())
    carp(CARP_FATAL, "Cannot open a write stream with path %s", path);
  return res;
}

istream& FileUtils::GetReadStream(const string& path){
  GenericStorageSystem* s = GenericStorageSystem::getStorage(path);
  return s->GetReadStream(path);
}

void FileUtils::CloseStream(ios_base& stream){
  GenericStorageSystem::CloseStream(stream);
}

string FileUtils::BaseName(const string& path) {
  GenericStorageSystem* s = GenericStorageSystem::getStorage(path);

  return s->BaseName(path);
}

string FileUtils::DirName(const string& path) {
  GenericStorageSystem* s = GenericStorageSystem::getStorage(path);

  return s->DirName(path);
}

string FileUtils::Stem(const string& path) {
  GenericStorageSystem* s = GenericStorageSystem::getStorage(path);

  return s->Stem(path);
}

string FileUtils::Extension(const string& path) {
  GenericStorageSystem* s = GenericStorageSystem::getStorage(path);

  return s->Extension(path);
}

void FileUtils::Copy(const std::string& orig, const std::string& dest) {
  GenericStorageSystem* s_from = GenericStorageSystem::getStorage(orig);
  GenericStorageSystem* s_to = GenericStorageSystem::getStorage(dest);

  if(s_from == s_to)
    s_from->CopyLocal(orig, dest);
  else
  {
    //TODO_RC: Verify that this will work efficiently.
    istream& in_str = s_from->GetReadStream(orig);
    ostream* out_str = s_to->GetWriteStream(dest, false);
    *out_str << (in_str.rdbuf());
    out_str->flush();
    delete out_str;
    //delete in_str;
  }
  
}

