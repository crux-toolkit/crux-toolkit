#ifndef WINCRUX_H
#define WINCRUX_H

#ifdef _MSC_VER
#define WIN32_LEAN_AND_MEAN
#include <direct.h>
#include <fcntl.h>
#include <io.h>
#include <string.h>
#include <windows.h>
#include "utils.h"
#include "windirent.h"

// The Windows C runtime has a hard limit on the 
// number of simultaneously opened files
#define WIN_MAX_OPEN_FILES 2048

// Rename some functions to the windows version
#define access _access
#define isfinite _finite
#define isnan _isnan
#define pclose _pclose
#define popen _popen
#define chdir _chdir
#define getcwd _getcwd
#define mkdir(a, b) _mkdir(a)
#define mkstemp _mktemp_s
#define sleep(x) Sleep(1000 * (x))
#define snprintf _snprintf

#undef NO_ERROR
#undef max
#undef min

#define R_OK 04
#define W_OK 02
#define F_OK 00

#define S_ISDIR(mode)  (((mode) & S_IFMT) == S_IFDIR)

int gettimeofday(struct timeval *tv, struct timezone *tz);
char *realpath(const char * file_name, char * resolved_name);

typedef struct 
{ 
   HANDLE f; 
   HANDLE m; 
   void *p; 
} SIMPLE_UNMMAP; 

// map 'filename' and return a pointer to it.
void *stub_mmap(const char *filename, SIMPLE_UNMMAP *un);

// Release memory mapped file
void stub_unmmap(SIMPLE_UNMMAP *un);

int scandir(
  const char *dirname, 
  struct dirent ***namelist,
  int (*select)(struct dirent *),
  int (*compar)(const void *, const void *)
 );

int alphasort(const void *d1, const void *d2);

int isinf(FLOAT_T x);
float log2(float x);

#endif
#endif

