#include "abspath.h"
#include <boost/filesystem.hpp>

using namespace std;

// Compute an absolute path from the parameter by prepending the current
// working directory when the given path is a relative path.
// Although the result is correct (according to 'man path_resolution')
// no effort is made to decode symlinks or normalize /. and /..
string AbsPath(const string& path) {
  return boost::filesystem::absolute(path).string();
}
