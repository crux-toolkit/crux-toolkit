#ifdef _MSC_VER
#include <direct.h>
#define getcwd _getcwd
#else
#include <unistd.h>
#endif
#include <stdlib.h>
#include <string>
#include "abspath.h"

using namespace std;

// Compute an absolute path from the parameter by prepending the current
// working directory when the given path is a relative path.
// Although the result is correct (according to 'man path_resolution')
// no effort is made to decode symlinks or normalize /. and /..
string AbsPath(const string& path) {
  if (path.empty())
    return "";
  if (path[0] == '/')
    return path;
  char* cwd = getcwd(NULL, 0);
  string result(cwd);
  free(cwd);
  result += string("/") + path;
  return result;
}
