#include "FileUtils.h"
#include "boost/filesystem.hpp"
#include <fstream>
#include <stdexcept>

using namespace std;

bool FileUtils::Exists(const string& path) {
  return boost::filesystem::exists(path);
}

bool FileUtils::IsRegularFile(const string& path) {
  return boost::filesystem::is_regular_file(path);
}

bool FileUtils::IsDir(const string& path) {
  return boost::filesystem::is_directory(path);
}

// returns true if a new directory was created, otherwise false
bool FileUtils::Mkdir(const string& path) {
  return boost::filesystem::create_directory(path);
}

void FileUtils::Rename(const string& from, const string& to) {
  if (Exists(from)) {
    boost::filesystem::rename(from, to);
  }
}

void FileUtils::Remove(const string& path) {
  if (Exists(path)) {
    boost::filesystem::remove_all(path);
  }
}

string FileUtils::Join(const string& path1, const string& path2) {
  return (boost::filesystem::path(path1) / boost::filesystem::path(path2)).string();
}

string FileUtils::Read(const string& path) {
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

ofstream* FileUtils::GetWriteStream(const string& path, bool overwrite) {
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

string FileUtils::BaseName(const string& path) {
  boost::filesystem::path p(path);
  return p.has_filename() ? p.filename().string() : "";
}

string FileUtils::DirName(const string& path) {
  boost::filesystem::path p(path);
  return p.has_parent_path() ? p.parent_path().string() : "";
}

string FileUtils::Stem(const string& path) {
  boost::filesystem::path p(path);
  return p.has_stem() ? p.stem().string() : "";
}

string FileUtils::Extension(const string& path) {
  boost::filesystem::path p(path);
  return p.has_extension() ? p.extension().string() : "";
}
void FileUtils::Copy(const std::string& orig, const std::string& dest) {
  if (Exists(orig)) {
    boost::filesystem::copy(orig, dest);
  }
}

