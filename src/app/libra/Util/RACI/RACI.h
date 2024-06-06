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

#ifndef _RACI_H
#define _RACI_H

#include <iostream>
#include "chunky_streambuf.h"
#include "random_access_compressed_streambuf.h"

/****************************************
  Main object class - what the user interacts with
  *****************************************/
class RACI : public std::istream {
public:

  //constructors & destructor
  RACI();
  RACI(const char* fname);
  ~RACI();

  //functions
  void open(const char* fname); // for ease of use as ifstream replacement
  bool is_open() const; // for ease of use as ifstream replacement
  void close(); // for ease of use as ifstream replacement
  
  enum eCompressionType { NONE, GZIP }; // maybe add bz2 etc later?
  eCompressionType getCompressionType() const {return compressionType;}

private:
  
  //data members
  eCompressionType compressionType;
};

#endif