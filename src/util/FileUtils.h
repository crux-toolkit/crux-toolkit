#ifndef FILEUTILS_H
#define FILEUTILS_H

#include <fstream>
#include <string>

class FileUtils {
 public:
  static bool Exists(const std::string& path);
  static bool IsRegularFile(const std::string& path);
  static bool IsDir(const std::string& path);
  static bool Mkdir(const std::string& path);
  static void Rename(const std::string& from, const std::string& to);
  static void Remove(const std::string& path);
  static std::string Join(const std::string& path1, const std::string& path2);
  static std::string Read(const std::string& path);
  static std::ofstream* GetWriteStream(const std::string& path, bool overwrite);
  static std::string BaseName(const std::string& path);
  static std::string DirName(const std::string& path);
  static std::string Stem(const std::string& path);
  static std::string Extension(const std::string& path);
  static void Copy(const std::string& orig, const std::string& dest);
 private:
  FileUtils();
  ~FileUtils();
};

#endif

