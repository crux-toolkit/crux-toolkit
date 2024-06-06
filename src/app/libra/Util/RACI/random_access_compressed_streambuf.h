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

#ifndef _RANDOM_ACCESS_COMPRESSED_STREAMBUF_H
#define _RANDOM_ACCESS_COMPRESSED_STREAMBUF_H

#include <cassert>
#include <fstream>
#include <vector>
#include "chunky_streambuf.h"

#ifndef Z_BUFSIZE
#define Z_BUFSIZE 16384
#endif

#define SPAN 1048576L       /* desired distance between access points */
#define WINSIZE 32768U      /* sliding window size */
#define RACI_CHUNK 16384         /* file input buffer size */

#define ALLOC(size) malloc(size)
#define gzio_raw_readerror(s) (s->infile->bad())

/* gzip flag byte */
#define ASCII_FLAG   0x01 /* bit 0 set: file probably ascii text */
#define HEAD_CRC     0x02 /* bit 1 set: header CRC present */
#define EXTRA_FIELD  0x04 /* bit 2 set: extra field present */
#define ORIG_NAME    0x08 /* bit 3 set: original file name present */
#define COMMENT      0x10 /* bit 4 set: file comment present */
#define RESERVED     0xE0 /* bits 5..7: reserved */

static int const gz_magic[2] = { 0x1f, 0x8b }; /* gzip magic header */

/****************************************
  access point entry
  *****************************************/
class synchpoint {
public:
  std::streamoff out;          /* corresponding offset in uncompressed data */
  std::streamoff in;           /* offset in input file of first full byte */
  z_stream* state; // stream state
};


class random_access_compressed_streambuf : public std::streambuf {
  friend class random_access_compressed_ifstream; // so we can modify some behaviors
public:

  //constructors and destructors
  random_access_compressed_streambuf(chunky_streambuf* rdbuf); // ctor
  virtual ~random_access_compressed_streambuf();

  //functions
  bool is_open() const;
  chunky_streambuf* close(); // close file and hand back readbuf

protected:

  //functions
  virtual pos_type seekoff(off_type off, std::ios_base::seekdir way, std::ios_base::openmode which = std::ios_base::in); // we don't do out
  virtual pos_type seekpos(pos_type pos, std::ios_base::openmode which = std::ios_base::in); // we don't do out
  virtual int_type underflow(); // repopulate the input buffer


private:

  //functions
  synchpoint* addIndexEntry(std::streamoff in, std::streamoff out); //Add an entry to the access point list.
  std::streampos get_next_read_pos();
  std::streampos my_seekg(std::streampos offset, std::ios_base::seekdir whence, std::ios_base::openmode Mode); // streampos not the same as streamoff, esp in win32
  bool pos_is_in_outbuf(std::streampos pos);
  void update_istream_ptrs(std::streampos new_headpos, std::streamoff new_buflen, std::streamoff new_posoffset = 0);

  // gzip stuff
  int build_index();
  void check_header();
  int destroy();
  int do_flush(int flush);
  int get_byte();
  int get_buf(int len);
  uLong getLong();
  
  //data members
  z_stream stream;
  int      z_err;   /* error code for last stream operation */
  int      z_eof;   /* set if end of input file */
  std::istream* infile;   /* raw .gz file we're reading */
  Byte* inbuf;  /* input buffer */
  Byte* outbuf; /* output buffer */
  uLong    crc;     /* crc32 of uncompressed data */
  std::streamoff  start;   /* start of compressed data in file (header skipped) */
  std::streamoff  uncompressedLength; /* total length of uncompressed file */
  std::streampos  last_seek_pos; /* last requested seek position (uncompressed) */
  std::streampos  outbuf_headpos; /* filepos for head of outbuf */
  std::streamoff	outbuf_len; /* length of outbuf last time we populated it */
  std::vector<synchpoint*> index; // index for random access

  
};

#endif
