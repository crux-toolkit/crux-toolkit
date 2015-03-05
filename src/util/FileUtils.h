#ifndef FILEUTILS_H
#define FILEUTILS_H

#include <fstream>
#include <string>

class FileUtils {
 public:
  static bool Exists(const std::string& path);
  static std::string Read(const std::string& path);
  static std::ofstream* GetWriteStream(const std::string& path, bool overwrite);
  static std::string BaseName(const std::string& path);
  static std::string DirName(const std::string& path);
  static std::string Stem(const std::string& path);
  static std::string Extension(const std::string& path);
 private:
  FileUtils();
  ~FileUtils();
};

#endif

