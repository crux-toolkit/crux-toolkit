#ifndef BOOSTFILESYSTEM_H
#define BOOSTFILESYSTEM_H

#include "GenericStorageSystem.h"
#include "boost/filesystem.hpp"

class BoostFileSystem : public GenericStorageSystem
{

public:
  bool Exists(const std::string &path) { return boost::filesystem::exists(path); }
  bool IsRegularFile(const std::string &path) { return boost::filesystem::is_regular_file(path); }
  bool IsDir(const std::string &path) { return boost::filesystem::is_directory(path); }
  bool Mkdir(const std::string &path) { return boost::filesystem::create_directory(path); }
  void Rename(const std::string &from, const std::string &to)
  {
    if (Exists(from))
    {
      boost::filesystem::rename(from, to);
    }
  }
  void Remove(const std::string &path)
  {
    if (Exists(path))
    {
      boost::filesystem::remove_all(path);
    }
  }
  std::string Join(const std::string &path1, const std::string &path2)
  {
    return (boost::filesystem::path(path1) / boost::filesystem::path(path2)).string();
  }

  // Compute an absolute path from the parameter by prepending the current
  // working directory when the given path is a relative path.
  // Although the result is correct (according to 'man path_resolution')
  // no effort is made to decode symlinks or normalize /. and /..
  std::string AbsPath(const string &path)
  {
    return boost::filesystem::absolute(path).string();
  }

  std::string Read(const std::string &path);
  std::ostream *GetWriteStream(const std::string &path, bool overwrite);
  std::istream &GetReadStream(const string &path);

  std::string BaseName(const std::string &path)
  {
    boost::filesystem::path p(path);
    return p.has_filename() ? p.filename().string() : "";
  }
  std::string DirName(const std::string &path)
  {
    boost::filesystem::path p(path);
    return p.has_parent_path() ? p.parent_path().string() : "";
  }
  std::string Stem(const std::string &path)
  {
    boost::filesystem::path p(path);
    return p.has_stem() ? p.stem().string() : "";
  }
  std::string Extension(const std::string &path)
  {
    boost::filesystem::path p(path);
    return p.has_extension() ? p.extension().string() : "";
  }

  void CopyLocal(const std::string &orig, const std::string &dest)
  {
    if (Exists(orig))
    {
      boost::filesystem::copy(orig, dest);
    }
  }

  BoostFileSystem(){};

  ~BoostFileSystem(){};

private:
void _CloseStream(const StreamRecord& p_streamRec){
    if(p_streamRec.is_input)
        delete (ifstream*)p_streamRec.object_pointer;
    else
        delete (ofstream*)p_streamRec.object_pointer;
}


};

#endif