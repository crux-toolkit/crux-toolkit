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

#include "chunky_streambuf.h"

chunky_streambuf::chunky_streambuf() {
  handle = NULL;
  flen = 0;
  path = NULL;
  desired_readbuf_len = DEFAULT_CHUNKY_BUFSIZE;
  last_seek_pos = -1;
  for (current_inbuf = N_INBUFS; current_inbuf--;) {
    inbuf[current_inbuf].chars_used = 0;
    inbuf[current_inbuf].filereadbuf = NULL;
  }
  current_inbuf = 0;
}

chunky_streambuf::~chunky_streambuf() {
  close();
  for (current_inbuf = N_INBUFS; current_inbuf--;) {
    if(filereadbuf()!=NULL) free(filereadbuf());
  }
  if(path!=NULL) free(path);
}

void chunky_streambuf::close() { // for ifstream-ish-ness
  if (is_open()) {
    handle->close();
    delete handle;
    handle = NULL;
  }
}

int chunky_streambuf::find_readbuf_for_gptr() {
  for (int n = N_INBUFS; n--;) {
    if (inbuf[n].filereadbuf &&
      ((char*)inbuf[n].filereadbuf <= gptr()) &&
      (gptr() <= ((char*)inbuf[n].filereadbuf + inbuf[n].readbuf_len))) {
      return n;
    }
  }
  return -1;
}

int chunky_streambuf::find_readbuf_for_pos(std::streamoff pos) {
  for (int n = N_INBUFS; n--;) {
    if ( readbuf_head_filepos(n) <= pos &&  pos < (readbuf_head_filepos(n) + readbuf_len(n))) {
      return n;
    }
  }
  return -1;
}

std::streamoff chunky_streambuf::get_next_read_pos() {
  return readbuf_head_filepos() + (gptr() - (char*)filereadbuf());
}

bool chunky_streambuf::is_open() const { // for ifstream-ish-ness
  return (handle && handle->is_open());
}

std::streampos chunky_streambuf::my_seekg(std::streamoff offset, std::ios_base::seekdir whence, std::ios_base::openmode Mode) {
  if (!is_open()) {
    return -1;
  }
  // watch out for a simple rewind or tellg
  if (0 == offset) {
    if (whence == std::ios_base::cur) {
      if (last_seek_pos >= 0) {  // unsatisfied seek request
        return last_seek_pos; // that's where we'd be if we'd actually seeked
      } else {
        return get_next_read_pos(); // nothing to do
      }
    }
    if (whence == std::ios_base::beg) { // just a rewind
        // do we already have this loaded?
      int n = find_readbuf_for_pos(0);
      if (n >= 0) {
        set_inbuf(n, 0);
      } else { // rewind 
        physical_seek(0);
        set_inbuf(current_inbuf + 1, 0, 0, 0, 0); // blow the cache
      }
      last_seek_pos = -1; // no need to seek
      return 0;
    } // end just a rewind
  } // end if offset==0

   /* find where in stream to start */
   /* compute absolute position */
  std::streamoff pos = offset;
  if (whence == std::ios_base::cur) {
    pos += get_next_read_pos();
  } else if (whence == std::ios_base::end) {
    if (!flen) { // length is unknown
      flen = std::streamoff(handle->pubseekoff(0, std::ios_base::end));
    }
    pos = flen + pos;
  }

  // do we already have this loaded?
  int n = find_readbuf_for_pos(pos); // will wait for readhead if any
  if (n >= 0) {
    set_inbuf(n, pos);
    last_seek_pos = -1; // no need to actually seek
  } else {
    // just note the request for now, actually seek at read time
    last_seek_pos = pos;
    set_inbuf(current_inbuf + 1, 0, 0, 0, 0); // blow the cache
  }
  return pos;
}

bool chunky_streambuf::open(const char* path) {
  if (!path) {
    return false;
  }

  if(handle!=NULL){
    printf("chunky_streambuf::open file already open\n");
    return false;
  }

  handle = new std::filebuf;
  handle->open(path,std::ios::in|std::ios::binary);
  flen = 0;
  desired_readbuf_len = 0;
  // dynamic read buffer sizing - start small
  for (current_inbuf = N_INBUFS; current_inbuf--;) {
    chars_used() = 0;
    readbuf_len() = 0;
    readbuf_head_filepos() = -1;
    filereadbuf() = NULL;
    maxbufsize() = 0;
  }
  // allocate a modest read buffer (we'll grow it if that makes
  // sense for how the file is being read)
  resize_readbufs(4 * DEFAULT_CHUNKY_BUFSIZE); // alloc moderately
  resize_readbufs(DEFAULT_CHUNKY_BUFSIZE); // but start small
  for (current_inbuf = N_INBUFS; current_inbuf--;) {
    if (!filereadbuf()) {
      return false; // no memory
    }
  }
  this->path = strdup(path);
  if (this->path == NULL || !is_open()) {
    return false;
  }

  set_inbuf(0, 0, 0, 0, 0); // we're pointed at head of file, but no read yet
  return true;
}

bool chunky_streambuf::physical_seek(const std::streampos& pos) {
  return physical_seek(std::streamoff(pos));
}

bool chunky_streambuf::physical_seek(const std::streamoff& off) {
  if (off >= 0) {
    handle->pubseekoff(off, std::ios_base::beg);
    return true;
  }
  return false;
}

bool chunky_streambuf::pos_is_in_readbuf(std::streamoff pos) {
  return ((readbuf_head_filepos() <= pos) &&
    (pos < (readbuf_head_filepos() + readbuf_len())));
}

void chunky_streambuf::resize_readbufs(size_t newsize) {
  newsize = (newsize < DEFAULT_CHUNKY_BUFSIZE) ? DEFAULT_CHUNKY_BUFSIZE : newsize;
  newsize = (newsize > MAX_CHUNKY_BUFSIZE) ? MAX_CHUNKY_BUFSIZE : newsize;
  desired_readbuf_len = newsize;
  int hotbuf = find_readbuf_for_gptr();
  for (int n = N_INBUFS; n--;) {
    if (newsize > maxbufsize(n)) {
      Byte* newbuf = (Byte*)realloc(filereadbuf(n), newsize);
      if (newbuf) {
        filereadbuf(n) = newbuf;
        maxbufsize(n) = newsize;
        if (n == hotbuf) { // let parent know we moved the read buffer
            // declare first, next, last for istream use
          setg((char*)newbuf, (char*)newbuf + (gptr() - eback()), (char*)newbuf + readbuf_len(n));
        }
      } else {
        // Failed to reallocate. So, stick with the old size and stop trying.
        desired_readbuf_len = maxbufsize(n);
        break;
      }
    }
  }
}

std::streamsize chunky_streambuf::readBuffer(int inbuf, std::streamoff headpos, std::streamsize readlen) {
  inbuf = inbuf % N_INBUFS;
  std::streamsize nread = handle->sgetn((char*)filereadbuf(inbuf), readlen);
  // std::cout << headpos << " read " << nread << " ask " << readlen << "\n";
  if (nread > 0) { // in case of eof, leave current buffer state alone
    readbuf_len(inbuf) = nread;
    readbuf_head_filepos(inbuf) = headpos;
    readbuf_requested_len(inbuf) = readlen;
    return nread;
  } else {
    return 0;
  }
}

std::streampos chunky_streambuf::seekoff(std::streamoff offset, std::ios_base::seekdir whence, std::ios_base::openmode Mode) {
  return my_seekg(offset, whence, Mode);
}

std::streampos chunky_streambuf::seekpos(std::streampos pos, std::ios_base::openmode Mode) {
  return my_seekg(std::streamoff(pos), std::ios_base::beg, Mode);
}

// for reusing an already loaded buffer
void chunky_streambuf::set_inbuf(int bufnum, std::streamoff readpos) { 
  bufnum = bufnum % N_INBUFS;
  set_inbuf(bufnum, readbuf_head_filepos(bufnum), readbuf_len(bufnum), readpos, readbuf_requested_len(bufnum));
}

void chunky_streambuf::set_inbuf(int bufnum, std::streamoff headpos, std::streamsize readlen, std::streamoff readpos, std::streamsize requested_readlen) {
  current_inbuf = bufnum % N_INBUFS;
  update_istream_ptrs(headpos, readlen, readpos - headpos);
  inbuf[current_inbuf].readbuf_requested_len = requested_readlen;
}

void chunky_streambuf::update_istream_ptrs(std::streamoff new_headpos, std::streamsize new_buflen, std::streamsize new_posoffset) {
  readbuf_head_filepos() = new_headpos; // note the filepos corresponding to buf head
  readbuf_len() = new_buflen; // how many bytes in the buffer?
  if (new_buflen) { // not just a reset to provoke underflow
    int n = find_readbuf_for_gptr();
    if (n >= 0)
      inbuf[n].chars_used = gptr() - eback(); // note usage of old buffer
  }
  // declare first, next, last for istream use
  setg((char*)filereadbuf(), (char*)filereadbuf() + new_posoffset, (char*)filereadbuf() + new_buflen);
}


// this gets called each time ifstream uses up its input buffer
chunky_streambuf::int_type chunky_streambuf::underflow() {
  std::streamsize nread = 0;
  std::streamoff next_read_pos;
  if (last_seek_pos >= 0) { // we had a seek request
    int bufnum = find_readbuf_for_pos(last_seek_pos); // will wait for readahead thread
    if (bufnum >= 0) { // found it already loaded
      set_inbuf(bufnum, last_seek_pos);
      std::streamsize offset = last_seek_pos - readbuf_head_filepos();
      last_seek_pos = -1; // satisfied it
      return filereadbuf()[offset];
    } else { // actually need to read
        // did we leave a lot of chars unread in current buffer?
      if (chars_used() && (int)desired_readbuf_len > 2 * chars_used()) {
        // reduce the read size a bit
        size_t newsize = DEFAULT_CHUNKY_BUFSIZE * (1 + (chars_used() / DEFAULT_CHUNKY_BUFSIZE)); // ugh - why isn't min() portable
        resize_readbufs(newsize);
      }
      next_read_pos = last_seek_pos; // will be new outbuf head pos
      last_seek_pos = -1; // satisfied it
      if (physical_seek(next_read_pos)) {
        nread = readBuffer(current_inbuf + 1, next_read_pos, desired_readbuf_len);
        if (nread) {
          set_inbuf(current_inbuf + 1, next_read_pos);
        }
      }
    }
  } else {
    // we hit end of buffer - perhaps we are streaming?
    next_read_pos = get_next_read_pos();
    int bufnum = find_readbuf_for_pos(next_read_pos); // will wait for readahead thread
    bool readahead_success = false;
    if (bufnum >= 0) { // already loaded
      set_inbuf(bufnum, next_read_pos);
      readahead_success = true;
      if (shortread(bufnum) || // eof after this buffer
        find_readbuf_for_pos(readbuf_head_filepos() + readbuf_len()) > -1) {
        // fully cached - just get on with it
        return filereadbuf()[get_next_read_pos() - readbuf_head_filepos()];
      }
    }
    // are we making contiguous reads?
    int last_inbuf = find_readbuf_for_pos(readbuf_head_filepos() - 1);
    bool streaming = ((last_inbuf > -1) && (find_readbuf_for_pos(readbuf_head_filepos(last_inbuf) - 1) > -1));
    // we read all of last buffer, try a bigger read now
    resize_readbufs(desired_readbuf_len * 4);

    if (!readahead_success) { // need an immediate blocking read
      if (physical_seek(next_read_pos)) {
        nread = readBuffer(current_inbuf + 1, next_read_pos, desired_readbuf_len);
        if (nread > 0) { // eof?
          set_inbuf(current_inbuf + 1, next_read_pos);
        }
      }
    } else {
      nread = readbuf_len();
    }
    if (nread && streaming && !shortread(current_inbuf)) { // not eof, not already reading, in sequential read mode
        // read the next chunk asynchronously in hopes we'll want that next
      int readerThread_inbuf = (current_inbuf + 1) % N_INBUFS;
      if (physical_seek(readbuf_head_filepos() + nread)) {
        readBuffer(readerThread_inbuf, readbuf_head_filepos() + nread, desired_readbuf_len);
      }
    }
  }
  if (!nread) {
    return traits_type::eof();
  } else {
    return filereadbuf()[get_next_read_pos() - readbuf_head_filepos()];
  }
}

