//
// sysdepend.h
//
// System specific dependencies, syntactic sugar hidden here.
//
#ifndef SYSDEPEND_INC
#define SYSDEPEND_INC

#ifdef _MSC_VER               // MS Visual Studio quirks
#pragma warning(disable:4305) // don't bark about double to float conversion
#pragma warning(disable:4244) // don't bark about double to float conversion
#pragma warning(disable:4786) // don't bark about "identifier was truncated to '255' characters in the browser information"
#pragma warning(disable:4996) // don't bark about "unsafe" functions
#endif

// platform independent
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/stat.h>
#include <stdio.h>	// for FILE declaration

#if defined(_MSC_VER) || defined(__MINGW32__)	// MinGW or MS Visual Studio?
#define WINDOWS_NATIVE
#endif

#ifndef TPPLIB
#define TPPLIB // so RAMP, for one, knows it can call on other tpplib functions
#endif


#ifdef WINDOWS_NATIVE	// --- Windows definitions ----------------------------

// DEPRECATED or IRRELEVANT
// #define USING_RELATIVE_WEBSERVER_PATH 1

// include this here so the Array class doesn't mess up SDK includes later
#include <winsock2.h> 
#include <windows.h> 

// rand isn't as random as rand48, but hobble by
#define drand48() ((double)rand()/(double)RAND_MAX))
#define lrand48() rand()
#define srand48(seed) srand(seed)
#define srandom(seed) srand(seed)
#define random() rand()

#include <direct.h>	// win32 filesystem directory manipulation functions
#include <io.h>		// win32 low level I/O functions

#define sleep(sec) Sleep(sec*1000)
typedef int gid_t;
typedef int uid_t;

#ifndef GNUPLOT_BINARY
#  define GNUPLOT_BINARY "wgnuplot.exe"
#endif

#ifndef _MSC_VER
#include <sys/types.h>
#include <unistd.h>
#endif


// we have to do a bit of preprocessing to deal with all the unixy system and
// pipe calls in TPP
//
#define system(a)  win32_system(a)         /* function is in win32_system.c */
#undef  popen
#define popen(a,b) win32_popen(a,b)        /* function is in win32_system.c */

// this may seem redundant, but not when "using ::system;" is in play as in some
// gcc hdr files
#define tpplib_system(a) win32_system(a)   /* function is in win32_system.c */
#define tpplib_popen(a,b) win32_popen(a,b) /* function is in win32_system.c */
#undef  pclose
#define pclose(a) _pclose(a)

int  win32_system(const char *cmd);
FILE *win32_popen(const char *cmd, const char *mode);

#endif			// --- Windows definitions ----------------------------


#ifdef _MSC_VER		// --- MS Visual Studio -------------------------------

// use Win32 threads
#define TPP_WIN32THREADS

typedef int mode_t;

#define lrint(num) (long)floor((num) + 0.5)
#define round(num) (long)floor((num) + 0.5)
#include <float.h>
#define isnan _isnan
#define isfinite _finite
#define isinf !_finite
#define tempnam _tempnam
#define mktemp(s) _mktemp(s)

#define S_ISREG(mode) ((mode)&_S_IFREG)
#define S_ISDIR(mode) ((mode)&_S_IFDIR)

#ifndef F_OK     // for use with access()
#  define F_OK 0 // exists
#  define W_OK 2 // write OK
#  define R_OK 4 // read OK
#endif

#define snprintf _snprintf
#ifndef strcasecmp
#  define strcasecmp stricmp
#  define strncasecmp strnicmp
#endif

inline unsigned _int64 strtoull(const char *str, char **end,int base) {
	return _strtoui64(str,end,base);
}

#endif			// --- MS Visual Studio -------------------------------


#ifndef WINDOWS_NATIVE	// -- Linux/Other -------------------------------------

#ifndef GNUPLOT_BINARY
#define GNUPLOT_BINARY "gnuplot"
#endif

#define tpplib_system(a)  system(a)
#define tpplib_popen(a,b) popen(a,b)

#endif			// -- Linux/Other -------------------------------------


#ifndef WEXITSTATUS
#define WEXITSTATUS(x) x
#endif

#ifndef Stat_t
#define Stat_t struct stat
#endif

#ifndef Boolean 
typedef unsigned int Boolean;
#endif

#ifndef u_char
typedef unsigned char u_char;
#endif

#ifndef u_short
typedef unsigned short u_short;
#endif

#ifndef u_int
typedef unsigned int u_int;
#endif

#ifndef u_long
typedef unsigned long u_long;
#endif

#ifndef True
#define True 1
#endif

#ifndef False
#define False 0
#endif

#define bDebug False

#define myfabs(a) ((a) < 0.0 ? (a) * -1.0 : (a))

#define DABS(a) ((a) < 0.0 ? (a) * -1.0 : (a))

//#define MIN(a,b) (a < b ? a : b)
//#define MAX(a,b) (a > b ? a : b)

// class to instatiate at program start to handle install dir issues etc
// MOSTLY DEPRECATED
#include "hooks_tpp.h"

/*
  like strdup(), but uses new instead of malloc
*/
inline char* strCopy(const char* orig) {
   if (!orig) {
      return NULL;
   }
   size_t len;
   char* output = new char[(len=strlen(orig)+1)];
   memcpy(output, orig, len);
   return output;
}

// faster and safer than #define style macro, which may evaluate twice
template<typename T> inline const T& Min(const T &a, const T &b) {
   return (a<b)?a:b; 
}
template<typename T> inline const T& Max(const T &a, const T &b) {
   return (a>b)?a:b; 
}

inline char* strlwr(char* orig) {
   if (!orig) {
       return NULL;
   }
   for (char *o=orig;*o;o++) {
     *o = tolower(*o);
   }
   return orig;
}

inline void cnvtUpper(char* o) {
   while (*o) {
     *o = toupper(*o);
	 o++;
   }
}

inline bool isEmptyFile(const char *file) {
  Stat_t statbuf;
  bool result=(0==stat(file,&statbuf));
  if (result) {
    if (statbuf.st_size <= 0) {
      result = true;
    }
    else {
      result = false;
    }
  }
  else {
    result = true;
  }
  return result;
}

inline double sqr(double x) {
   return x*x;
}

inline float sqr(float x) {
   return x*x;
}

inline const char *getCmdlineQuoteChar() {
#ifdef WINDOWS_NATIVE
    return "\"";
#else 	
    return "'";
#endif
}

#endif		// SYSDEPEND_INC
