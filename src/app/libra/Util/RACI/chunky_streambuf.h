// RACI: short for Random Access Compressed ifstream
// This is a standalone version of Brian Pratt's implementation in ProteoWizard
// It has also been un-boosted for size and portability
// February 10 2020, Michael Hoopmann, ISB. All credit and licensing maintained below.

// $Id: random_access_compressed_ifstream.hpp 4140 2012-11-22 01:07:16Z pcbrefugee $
//
// This is just an istream which chooses its streambuf implementation
// based on whether the file is gzipped or not.  Could be extended
// to work with other compression formats, too.
//
// What makes it interesting compared to the classic gzstream implementation
// is the ability to perform seeks in a reasonably efficient manner.  In the
// event that a seek is requested (other than a rewind, or tellg()) the file 
// is decompressed once and snapshots of the compression state are made 
// every 1MB or so.  Further seeks are then quite efficient since they don't
// have to begin at the head of the file.
//
// It also features threaded readahead with adaptive buffering - it will read 
// increasingly larger chunks of the raw file as it perceives a sequential read 
// in progress, and will launch a thread to grab the next probable chunk of the 
// file while the previous chunk is being parsed.  This is especially helpful for
// files being read across a slow network connection.
// If lots of seeking is going on, and the buffer size proves excessive, the read 
// buffer size decreases.  This is also true for the non-compressed case, so this 
// is generally useful for file read speed optimization.
//
// Copyright (C) Insilicos LLC 2008,2011 All Rights Reserved.
//
// draws heavily on example code from the zlib distro, so
// for conditions of distribution and use, see copyright notice in zlib.h
//
// based on:
/* gzio.c -- IO on .gz files
* Copyright (C) 1995-2005 Jean-loup Gailly.
* For conditions of distribution and use, see copyright notice in zlib.h
  Main idea there is as follows:

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

*/

// efficient random access stuff based on
/* zran.c -- example of zlib/gzip stream indexing and random access
* Copyright (C) 2005 Mark Adler
* For conditions of distribution and use, see copyright notice in zlib.h
Version 1.0  29 May 2005  Mark Adler */



#ifndef _CHUNKY_STREAMBUF_H
#define _CHUNKY_STREAMBUF_H

#define DECOMPRESS_BUFSIZE 16384
#define DEFAULT_CHUNKY_BUFSIZE 32768
#define MAX_CHUNKY_BUFSIZE DEFAULT_CHUNKY_BUFSIZE*1024

#ifdef _WIN32
#define strdup _strdup
#endif

#include <fstream>
#include <string.h>
#include "zlib.h"

/****************************************
  chunky_streambuf
  *****************************************/
class chunky_streambuf : public std::streambuf {
#define N_INBUFS 3 // lookback, current, readahead
public:

  //constructors and destructors
  chunky_streambuf();
  virtual ~chunky_streambuf();

  //functions
  void close();
  int find_readbuf_for_gptr();
  int find_readbuf_for_pos(std::streamoff pos);
  std::streamoff get_next_read_pos(); 
  bool is_open() const;
  std::streampos my_seekg(std::streamoff offset, std::ios_base::seekdir whence, std::ios_base::openmode Mode); // streampos not the same as streamoff, esp in win32
  bool open(const char* path);
  bool physical_seek(const std::streampos& pos);
  bool physical_seek(const std::streamoff& off);
  bool pos_is_in_readbuf(std::streamoff pos);
  std::streamsize readBuffer(int inbuf, std::streamoff headpos, std::streamsize readlen);
  void set_inbuf(int bufnum, std::streamoff readpos); // for reusing an already loaded buffer
  void set_inbuf(int bufnum, std::streamoff headpos, std::streamsize readlen, std::streamoff readpos, std::streamsize requested_readlen);
  void update_istream_ptrs(std::streamoff new_headpos, std::streamsize new_buflen, std::streamsize new_posoffset = 0);

  virtual pos_type seekoff(off_type off, std::ios_base::seekdir way, std::ios_base::openmode which = std::ios_base::in); // we don't do out
  virtual pos_type seekpos(pos_type pos, std::ios_base::openmode which = std::ios_base::in); // we don't do out
  virtual int_type underflow(); // repopulate the input buffer
  

  //inline functions
  inline std::streamoff& chars_used() { return inbuf[current_inbuf].chars_used;  }/* number of chars actually read out of buffer */
  inline Byte*& filereadbuf() { return inbuf[current_inbuf].filereadbuf; } /* file read buffer */
  inline Byte*& filereadbuf(int n) { return inbuf[n].filereadbuf; } /* file read buffer */
  inline size_t& maxbufsize() { return inbuf[current_inbuf].maxbufsize; } /* size of read buffer */
  inline size_t& maxbufsize(int n) { return inbuf[n].maxbufsize; }/* size of read buffer */
  inline std::streamoff& readbuf_head_filepos() { return inbuf[current_inbuf].readbuf_head_filepos; }
  inline std::streamoff& readbuf_head_filepos(int n) { return inbuf[n].readbuf_head_filepos; }
  inline std::streamsize& readbuf_len() { return inbuf[current_inbuf].readbuf_len; } /* length of readbuf last time we populated it */
  inline std::streamsize& readbuf_len(int n) { return inbuf[n].readbuf_len; }/* length of readbuf last time we populated it */
  inline std::streamsize& readbuf_requested_len() { return inbuf[current_inbuf].readbuf_requested_len; } /* bytes we tried to read last time we populated it */
  inline std::streamsize& readbuf_requested_len(int n) { return inbuf[n].readbuf_requested_len; }
  inline void resize_readbufs(size_t newsize);
  inline bool shortread(int n) const { return ((inbuf[n].readbuf_requested_len > 0) && (inbuf[n].readbuf_requested_len != inbuf[n].readbuf_len));  }

  //data members
  std::filebuf* handle;   /* handle of file we're reading */
  char* path;   /* path name for debugging only */
  size_t   desired_readbuf_len; /* dynamic sizing of disk reads */
  std::streamoff last_seek_pos; /* last requested seek position */
  int current_inbuf;
  std::streamoff flen; /* unknown until first SEEK_END */

  //data member structures
  struct {
    Byte* filereadbuf; /* file read buffer */
    size_t   maxbufsize; /* size of read buffer */
    std::streamoff readbuf_head_filepos; /* filepos for head of readbuf */
    std::streamsize readbuf_len; /* length of readbuf last time we populated it */
    std::streamsize readbuf_requested_len; /* number of bytes we actually tried to load */
    std::streamoff chars_used; /* number of chars actually read out of buffer */
  } inbuf[N_INBUFS]; // read out of one while the other populates in another thread
  
};

#endif
