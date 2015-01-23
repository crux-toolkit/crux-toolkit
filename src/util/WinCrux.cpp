/**
 * \file WinCrux.cpp
 * \brief Support functions need to compile Crux under native windows.
 */

#include <stdlib.h>
#include <time.h>
#include <Winsock2.h>
#include "carp.h"
#include "utils.h"
#include "WinCrux.h"
#include <iostream>
#include <vector>


using namespace std;

// Windows GetTimeOfDay() code from http://www.suacommunity.com/dictionary/index.php

#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
 
struct timezone
{
  int  tz_minuteswest; /* minutes W of Greenwich */
  int  tz_dsttime;     /* type of dst correction */
};
 
// Definition of a gettimeofday function
 
int gettimeofday(struct timeval *tv, struct timezone *tz)
{
// Define a structure to receive the current Windows filetime
  FILETIME ft;
 
// Initialize the present time to 0 and the timezone to UTC
  unsigned __int64 tmpres = 0;
  static int tzflag = 0;
 
  if (NULL != tv)
  {
    GetSystemTimeAsFileTime(&ft);
 
// The GetSystemTimeAsFileTime returns the number of 100 nanosecond 
// intervals since Jan 1, 1601 in a structure. Copy the high bits to 
// the 64 bit tmpres, shift it left by 32 then or in the low 32 bits.
    tmpres |= ft.dwHighDateTime;
    tmpres <<= 32;
    tmpres |= ft.dwLowDateTime;
 
// Convert to microseconds by dividing by 10
    tmpres /= 10;
 
// The Unix epoch starts on Jan 1 1970.  Need to subtract the difference 
// in seconds from Jan 1 1601.
    tmpres -= DELTA_EPOCH_IN_MICROSECS;
 
// Finally change microseconds to seconds and place in the seconds value. 
// The modulus picks up the microseconds.
    tv->tv_sec = (long)(tmpres / 1000000UL);
    tv->tv_usec = (long)(tmpres % 1000000UL);
  }
 
  if (NULL != tz)
  {
    if (!tzflag)
    {
      _tzset();
      tzflag++;
    }
  
// Adjust for the timezone west of Greenwich
      tz->tz_minuteswest = _timezone / 60;
    tz->tz_dsttime = _daylight;
  }
 
  return 0;
}

char *realpath(const char * file_name, char * resolved_name) {

  char * full_path_buffer = (char *) mymalloc(MAX_PATH * sizeof(char));
  size_t needed_buff_size 
    = GetFullPathName(file_name, MAX_PATH, full_path_buffer, NULL);

  if (needed_buff_size == 0) {
    // An error occurred.
    LPSTR lpMsgBuf;
    FormatMessage(
      FORMAT_MESSAGE_ALLOCATE_BUFFER 
      | FORMAT_MESSAGE_FROM_SYSTEM 
      | FORMAT_MESSAGE_IGNORE_INSERTS,
      NULL,
      GetLastError(),
      MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
      (LPTSTR) &lpMsgBuf,
      0,
      NULL 
    );
    carp(CARP_FATAL, lpMsgBuf);

  }

  if (needed_buff_size > MAX_PATH) {
    full_path_buffer = (char *) myrealloc(full_path_buffer, needed_buff_size);
    needed_buff_size 
    = GetFullPathName(file_name, MAX_PATH, full_path_buffer, NULL);
  }

  return full_path_buffer;
}

/* 
This code was placed in the public domain by the author, 
Sean Barrett, in November 2007. Do with it as you will. 
(Seee the page for stb_vorbis or the mollyrocket source 
page for a longer description of the public domain non-license). 
*/ 

// Public domain code from https://mollyrocket.com/forums/viewtopic.php?p=2529
// map 'filename' and return a pointer to it. fill out *length and *un if not-NULL 
void *stub_mmap(const char *filename, SIMPLE_UNMMAP *un) 
{ 
   HANDLE f = CreateFile(filename, GENERIC_READ, FILE_SHARE_READ,  NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL); 
   HANDLE m; 
   void *p; 
   if (!f) return NULL; 
   m = CreateFileMapping(f, NULL, PAGE_READONLY, 0,0, NULL); 
   if (!m) { CloseHandle(f); return NULL; } 
   p = MapViewOfFile(m, FILE_MAP_READ, 0,0,0); 
   if (!p) { CloseHandle(m); CloseHandle(f); return NULL; } 
   if (un) { 
      un->f = f; 
      un->m = m; 
      un->p = p; 
   } 
   return p; 
} 

void stub_unmmap(SIMPLE_UNMMAP *un) 
{ 
   UnmapViewOfFile(un->p); 
   CloseHandle(un->m); 
   CloseHandle(un->f); 
} 

int scandir(
  const char *dirname, 
  struct dirent ***namelist,
  int (*select)(struct dirent *),
  int (*compar)(const void *, const void *)
) {
  errno = 0;
  
  DIR *directory = opendir(dirname);
  if (!directory) {
    return -1;
  }

  struct dirent *entry;
  vector<struct dirent *> *entries = new vector<struct dirent *>;
  while(entry = readdir(directory)) {
    if (!select || select(entry)) {
      struct dirent *valid_entry = new struct dirent;
      copy_dirent(entry, valid_entry);
      entries->push_back(valid_entry);
    }
  }
  if (errno) {
    return -1;
  }

  *namelist = static_cast<struct dirent **>(malloc(sizeof(struct dirent *) * entries->size()));
  copy(entries->begin(), entries->end(), *namelist);

  return entries->size();
 }

int alphasort(const void *d1, const void *d2) {

  // FIXME this is only a stub
  return 0;
}

char *mkdtemp(char *temp) {
  
  // Create temporary file name from template
  char *dirname = _mktemp(temp);
  
  if (dirname) {
    // Create the directory
    int result = _mkdir(dirname);
    if (result < 0) {
      free(dirname);
      dirname = NULL;
    }
  }

  return dirname;

}

int isinf(FLOAT_T x) {
  return !_finite(x);
}

float log2(FLOAT_T x) {
  return log(x) / log(2.0);
}
