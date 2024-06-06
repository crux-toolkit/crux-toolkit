#include "random_access_compressed_streambuf.h"

random_access_compressed_streambuf::random_access_compressed_streambuf(chunky_streambuf* rawbuf) {
  int err;

  stream.zalloc = (alloc_func)0;
  stream.zfree = (free_func)0;
  stream.opaque = (voidpf)0;
  stream.next_in = inbuf = Z_NULL;
  stream.next_out = outbuf = Z_NULL;
  stream.avail_in = stream.avail_out = 0;
  last_seek_pos = -1;
  infile = new std::istream(rawbuf); // dynamic disk buffer size
  z_err = Z_OK;
  z_eof = 0;
  crc = crc32(0L, Z_NULL, 0);

  if (infile->fail()) {
    destroy();
    return;
  }

  stream.next_in = inbuf = (Byte*)ALLOC(Z_BUFSIZE);
  outbuf = (Byte*)ALLOC(DECOMPRESS_BUFSIZE);

  err = inflateInit2(&stream, -MAX_WBITS);
  /* windowBits is passed < 0 to tell that there is no zlib header.
  * Note that in this case inflate *requires* an extra "dummy" byte
  * after the compressed stream in order to complete decompression and
  * return Z_STREAM_END. Here the gzip CRC32 ensures that 4 bytes are
  * present after the compressed stream.
  */
  if (err != Z_OK || inbuf == Z_NULL) {
    destroy();
    return;
  }
  stream.avail_out = DECOMPRESS_BUFSIZE;

  errno = 0;

  check_header(); /* skip the .gz header */
  infile->clear(); // clear eof flag for short files
  start = infile->tellg() - (std::streampos)stream.avail_in;
  update_istream_ptrs(0, 0); // we're pointed at head of file, but no read yet
}

random_access_compressed_streambuf::~random_access_compressed_streambuf() {
  destroy();
}

// Add an entry to the access point list. stealing from zran.c
synchpoint* random_access_compressed_streambuf::addIndexEntry(std::streamoff in, std::streamoff out) {
  synchpoint* next = new synchpoint(); /* fill in entry and increment how many we have */
  if (next) {
    next->in = in;
    next->out = out;
    next->state = new z_stream;
    inflateCopy(next->state, &stream);
    index.push_back(next);
  }
  return next;
}


int random_access_compressed_streambuf::build_index() {

/* Make one entire pass through the compressed stream and build an index, with
    access points about every span bytes of uncompressed output -- span is
    chosen to balance the speed of random access against the memory requirements
    of the list, about 32K bytes per access point.  Note that data after the end
    of the first zlib or gzip stream in the file is ignored.  build_index()
    returns 0 on success, Z_MEM_ERROR for out of memory, Z_DATA_ERROR for an error in
    the input file, or Z_ERRNO for a file read error.
    On success, *built points to the resulting index. */

  int ret;
  std::streamoff span = SPAN;
  std::streamoff totin, totout;        /* our own total counters to avoid 4GB limit */
  std::streamoff last;                 /* totout value of last access point */
  unsigned char* input = new unsigned char[RACI_CHUNK];
  unsigned char* window = new unsigned char[WINSIZE];
  z_stream& strm = stream;

  /* initialize inflate */
  stream.avail_in = 0;
  stream.avail_out = 0;
  stream.total_in = 0;
  stream.total_out = 0;
  stream.next_in = inbuf;
  crc = crc32(0L, Z_NULL, 0);
  ret = inflateReset(&strm);
  infile->clear(); // clear stale eof bit if any
  infile->seekg((std::streamoff)start); // rewind
  if (ret != Z_OK) {
    goto build_index_error;
  }

  /* inflate the input, maintain a sliding window, and build an index -- this
  also validates the integrity of the compressed data using the check
  information at the end of the gzip or zlib stream */
  totout = last = 0;
  totin = start;
  addIndexEntry(totin, totout); // note head of file

  do {
    /* get some compressed data from input file */
    infile->read((char*)input, RACI_CHUNK);
    strm.avail_in = (uInt)infile->gcount();
    if (gzio_raw_readerror(this)) {
      ret = Z_ERRNO;
      goto build_index_error;
    }
    if (strm.avail_in == 0) {
      ret = Z_DATA_ERROR;
      goto build_index_error;
    }
    strm.next_in = input;

    /* process all of that, or until end of stream */
    do {
      /* reset sliding window if necessary */
      if (strm.avail_out == 0) {
        strm.avail_out = WINSIZE;
        strm.next_out = window;
      }

      /* inflate until out of input, output, or at end of block --
      update the total input and output counters */
      totin += strm.avail_in; // note input filepos at start of inflate
      totout += strm.avail_out; // note uncompressed filepos prior to inflate
      ret = inflate(&strm, Z_BLOCK);      /* return at end of block */
      totin -= strm.avail_in; // if we use this as a synchpoint we'll have to repopulate the input buffer
      totout -= strm.avail_out;
      if (ret == Z_NEED_DICT) ret = Z_DATA_ERROR;
      if (ret == Z_MEM_ERROR || ret == Z_DATA_ERROR) goto build_index_error;
      if (ret == Z_STREAM_END) {
        // reached end successfully
        infile->clear(); // clear the fail bit
        break;
      }

      /* add an index entry every 'span' bytes   */
      if ((totout - last) > span) {
        if (!addIndexEntry(totin, totout)) {
          ret = Z_MEM_ERROR;
          goto build_index_error;
        }
        last = totout;
      }
    } while (strm.avail_in != 0);
  } while (ret != Z_STREAM_END);

  /* return index (release unused entries in list) */
  uncompressedLength = totout;

  /* return error */
build_index_error:
  delete[] window;
  delete[] input;
  return ret;
}

/* Check the gzip header of a random_access_gzstream opened for reading. set err
to Z_DATA_ERROR if the magic header is not present, or present but the rest of the header is incorrect.
IN assertion: the stream s has already been created sucessfully; stream.avail_in is zero for the first time, 
but may be non-zero for concatenated .gz files. */
void random_access_compressed_streambuf::check_header() {
  int method; /* method byte */
  int flags;  /* flags byte */
  uInt len;
  int c;

  /* Assure two bytes in the buffer so we can peek ahead -- handle case
  where first byte of header is at the end of the buffer after the last
  gzip segment */
  len = stream.avail_in;
  if (len < 2) {
    if (len) inbuf[0] = stream.next_in[0];
    errno = 0;
    infile->read((char*)(inbuf + len), Z_BUFSIZE >> len);
    len = (uInt)infile->gcount();
    if (len <= 0 && gzio_raw_readerror(this)) z_err = Z_ERRNO;
    stream.avail_in += len;
    stream.next_in = inbuf;
    if (stream.avail_in < 2) {
      if (stream.avail_in) z_err = Z_DATA_ERROR;
      return;
    }
  }

  /* Peek ahead to check the gzip magic header */
  if (stream.next_in[0] != gz_magic[0] || stream.next_in[1] != gz_magic[1]) {
    z_err = Z_DATA_ERROR;
    return;
  }
  stream.avail_in -= 2;
  stream.next_in += 2;

  /* Check the rest of the gzip header */
  method = get_byte();
  flags = get_byte();
  if (method != Z_DEFLATED || (flags & RESERVED) != 0) {
    z_err = Z_DATA_ERROR;
    return;
  }

  /* Discard time, xflags and OS code: */
  for (len = 0; len < 6; len++) {
    get_byte();
  }

  if ((flags & EXTRA_FIELD) != 0) { /* skip the extra field */
    len = (uInt)get_byte();
    len += ((uInt)get_byte()) << 8;
    /* len is garbage if EOF but the loop below will quit anyway */
    while (len-- != 0 && get_byte() != EOF);
  }
  if ((flags & ORIG_NAME) != 0) { /* skip the original file name */
    while ((c = get_byte()) != 0 && c != EOF);
  }
  if ((flags & COMMENT) != 0) {   /* skip the .gz file comment */
    while ((c = get_byte()) != 0 && c != EOF);
  }
  if ((flags & HEAD_CRC) != 0) {  /* skip the header crc */
    for (len = 0; len < 2; len++) {
      get_byte();
    }
  }
  z_err = z_eof ? Z_DATA_ERROR : Z_OK;
}

chunky_streambuf* random_access_compressed_streambuf::close() { // for ifstream-ish-ness
  chunky_streambuf* rawbuf = (chunky_streambuf*)infile->rdbuf(); // preserve
  infile->rdbuf(NULL);
  destroy();
  return rawbuf; // hand it back to parent
}

int random_access_compressed_streambuf::destroy() {
  int err = Z_OK;

  bool bClosedOK = true;
  if (stream.state != NULL) {
    err = inflateEnd(&stream);
    if(infile) delete infile;
    infile = NULL;
  }

  // clean up the seek index list if any
  for (int i = (int)index.size(); i--;) {
    inflateEnd(index[i]->state);
    delete index[i]->state;
    delete index[i];
  }
  index.clear(); // set length 0
  if (!bClosedOK) {

#ifdef ESPIPE
    if (errno != ESPIPE) /* fclose is broken for pipes in HP/UX */
#endif
      err = Z_ERRNO;
  }
  if (z_err < 0) {
    err = z_err;
  }

  if(inbuf) {
    free(inbuf);
    inbuf=NULL;
  }
  if(outbuf){
    free(outbuf);
    outbuf=NULL;
  }
  return err;
}

/* ===========================================================================
Read a byte from a random_access_gzstream; update next_in and avail_in. Return EOF
for end of file.
IN assertion: the stream s has been sucessfully opened for reading. */
int random_access_compressed_streambuf::get_byte() {
  if (z_eof) return EOF;
  if (stream.avail_in == 0) {
    errno = 0;
    infile->read((char*)inbuf, Z_BUFSIZE);
    stream.avail_in = (uInt)infile->gcount(); // how many did we read?
    if (stream.avail_in <= 0) {
      z_eof = 1;
      if (gzio_raw_readerror(this)) z_err = Z_ERRNO;
      else if (infile->eof()) infile->clear(); // clear eof flag
      return EOF;
    }
    stream.next_in = inbuf;
  }
  stream.avail_in--;
  return *(stream.next_in)++;
}

//Reads a long in LSB order from the given random_access_gzstream. Sets z_err in case of error.
uLong random_access_compressed_streambuf::getLong() {
  uLong x = (uLong)get_byte();
  int c;

  x += ((uLong)get_byte()) << 8;
  x += ((uLong)get_byte()) << 16;
  c = get_byte();
  if (c == EOF) {
    z_err = Z_DATA_ERROR;
  }
  x += ((uLong)c) << 24;
  return x;
}

std::streampos random_access_compressed_streambuf::get_next_read_pos() {
  return (unsigned long)outbuf_headpos + (gptr() - (char*)outbuf);
}

bool random_access_compressed_streambuf::is_open() const { // for ifstream-ish-ness
  return true; // only ever exist when file is open
}

std::streampos random_access_compressed_streambuf::my_seekg(std::streampos offset, std::ios_base::seekdir whence, std::ios_base::openmode Mode) {
  if (z_err == Z_ERRNO || z_err == Z_DATA_ERROR) {
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
    if ((whence == std::ios_base::beg) && // just a rewind
      !index.size()) { // no index yet
      // do we already have this decompressed?
      if (pos_is_in_outbuf(offset) || // already in buffer
        (offset == outbuf_headpos)) // or buffer is not yet loaded, let underflow() do it
      {
        update_istream_ptrs(outbuf_headpos, outbuf_len);
      } else { // rewind without provoking index build
        stream.avail_in = 0;
        stream.avail_out = 0;
        stream.total_in = 0;
        stream.total_out = 0;
        stream.next_in = inbuf;
        crc = crc32(0L, Z_NULL, 0);
        (void)inflateReset(&stream);
        last_seek_pos = -1; // no need to seek
        infile->seekg(start);
        update_istream_ptrs(0, 0); // blow the cache
      }
      return 0;
    } // end just a rewind and no index yet
  } // end if offset==0

  z_err = Z_OK;
  z_eof = 0;
  /* find where in stream to start */
  /* compute absolute position */
  std::streampos pos = offset;
  if (whence == std::ios_base::cur) {
    pos += get_next_read_pos();
  } else if (whence == std::ios_base::end) {
    if (!index.size()) { // first seek - build index
      build_index(); // sets uncompressedLength
    }
    pos += uncompressedLength;
  }

  // do we already have this decompressed?
  if (pos_is_in_outbuf(pos)) {
    update_istream_ptrs(outbuf_headpos, outbuf_len, pos - outbuf_headpos);
    last_seek_pos = -1; // no need to actually seek
  } else {
    // just note the request for now, actually seek at read time
    last_seek_pos = pos;
    if (!index.size()) { // first seek - build index
      build_index();
    }
    update_istream_ptrs(outbuf_headpos, 0); // blow the cache
  }
  return pos;
}

bool random_access_compressed_streambuf::pos_is_in_outbuf(std::streampos pos) {
  return ((outbuf_headpos <= pos) && (pos < (outbuf_headpos + outbuf_len)));
}

std::streampos random_access_compressed_streambuf::seekpos(std::streampos pos, std::ios_base::openmode Mode) {
  return my_seekg(pos, std::ios_base::beg, Mode);
}

std::streampos random_access_compressed_streambuf::seekoff(std::streamoff offset, std::ios_base::seekdir whence, std::ios_base::openmode Mode) {
  return my_seekg(std::streampos(offset), whence, Mode);
}

// this gets called each time ifstream uses up its input buffer
random_access_compressed_streambuf::int_type random_access_compressed_streambuf::underflow() {
  int nread = 0;
  int len = DECOMPRESS_BUFSIZE; // we'll try to read the next full chunk
  if (last_seek_pos >= 0) { // we had a seek request

      //here's where we tear out the inefficient default seek behavior and replace
      //with the stuff from zran.c - bpratt

      //Use the index to read len bytes from offset into buf, return bytes read or
      //negative for error (Z_DATA_ERROR or Z_MEM_ERROR).  If data is requested past
      //the end of the uncompressed data, then extract() will return a value less
      //than len, indicating how much as actually read into buf.  This function
      //should not return a data error unless the file was modified since the index
      //was generated.  extract() may also return Z_ERRNO if there is an error on
      //reading or seeking the input file. 

    int skip, ret;
    unsigned char* discard = new unsigned char[WINSIZE];
    std::streamoff offset = last_seek_pos; // last requested absolute uncompress filepos
    std::streampos next_read_pos = last_seek_pos; // will be new outbuf head pos
    last_seek_pos = -1; // satisfied it

    // first locate the index entry which will get us at or just before target
    size_t ind = index.size();
    while (--ind && index[ind]->out > offset);
    // and prepare to decompress
    synchpoint* synch = index[ind];
    inflateEnd(&stream); // tidy up old ptrs
    inflateCopy(&stream, synch->state);
    z_stream& strm = stream;
    assert(strm.total_in == synch->in - start);
    infile->clear(); // clear eof flag if any
    infile->seekg(synch->in);

    // skip uncompressed bytes until offset reached 
    offset -= synch->out;  // now offset is the number of uncompressed bytes we need to skip
    strm.avail_in = 0; // inflateCopy doesn't retain the input buffer
    skip = 1;                               // while skipping to offset 
    do {
      // define where to put uncompressed data, and how much 
      if (offset == 0 && skip) {          // at offset now
        strm.avail_out = len;
        strm.next_out = (Bytef*)outbuf;
        skip = 0;                       // only do this once
      }
      if (offset > WINSIZE) {             // skip WINSIZE bytes
        strm.avail_out = WINSIZE;
        strm.next_out = discard;
        offset -= WINSIZE;
      } else if (offset != 0) {             // last skip 
        strm.avail_out = (unsigned)offset;
        strm.next_out = discard;
        offset = 0;
      }

      // uncompress until avail_out filled, or end of stream 
      do {
        if (strm.avail_in == 0) {
          infile->read((char*)inbuf, RACI_CHUNK);
          strm.avail_in = infile->gcount();
          if (gzio_raw_readerror(this)) {
            ret = Z_ERRNO;
            goto perform_seek_ret;
          } else if (infile->eof()) {
            infile->clear(); // clear eof flag
          }
          if (strm.avail_in == 0) {
            ret = Z_DATA_ERROR;
            goto perform_seek_ret;
          }
          strm.next_in = inbuf;
        }
        ret = inflate(&strm, Z_NO_FLUSH);       // normal inflate
        if (ret == Z_NEED_DICT) ret = Z_DATA_ERROR;
        if (ret == Z_MEM_ERROR || ret == Z_DATA_ERROR) goto perform_seek_ret;
        if (ret == Z_STREAM_END) break;
      } while (strm.avail_out != 0);

      // if reach end of stream, then don't keep trying to get more 
      if (ret == Z_STREAM_END) {
        // reached EOF, clear fail bit
        infile->clear();
        break;
      }

      // do until offset reached and requested data read, or stream ends 
    } while (skip);
    nread = skip ? 0 : len - strm.avail_out;
    update_istream_ptrs(next_read_pos, nread);
    // return error if any
  perform_seek_ret:
    delete[]discard;

  } else {
    std::streampos buftailpos = get_next_read_pos(); // buf tail will be buf head unless we seek

    Bytef* start = (Bytef*)outbuf; // starting point for crc computation 
    Byte* next_out; // == stream.next_out but not forced far (for MSDOS) 

    next_out = (Byte*)start;
    stream.next_out = next_out;
    stream.avail_out = len;

    while (stream.avail_out != 0) {

      if (stream.avail_in == 0 && !z_eof) {

        errno = 0;
        infile->read((char*)inbuf, Z_BUFSIZE);
        stream.avail_in = infile->gcount();
        if (stream.avail_in <= 0) {
          z_eof = 1;
          if (gzio_raw_readerror(this)) {
            z_err = Z_ERRNO;
            break;
          } else if (infile->eof()) {
            infile->clear(); // clear eof flag
          }
        }
        stream.next_in = inbuf;
      }
      z_err = inflate(&(stream), Z_NO_FLUSH);

      if (z_err != Z_OK || z_eof) {
        break;
      }
    }
    crc = crc32(crc, start, (uInt)(stream.next_out - start));

    if (len == (int)stream.avail_out &&
      (z_err == Z_DATA_ERROR || z_err == Z_ERRNO)) {
      return traits_type::eof();;
    }
    nread = (int)(len - stream.avail_out); // how many bytes actually read?
    update_istream_ptrs(buftailpos, nread); // update outbuf head position, length
  }
  return (nread > 0) ? outbuf[0] : traits_type::eof();
}

void random_access_compressed_streambuf::update_istream_ptrs(std::streampos new_headpos, std::streamoff new_buflen, std::streamoff new_posoffset) {
  outbuf_headpos = new_headpos; // note the decompressed filepos corresponding to buf head
  outbuf_len = new_buflen; // how many bytes in the buffer are consumable?
  // declare first, next, last for istream use
  setg((char*)outbuf, (char*)outbuf + new_posoffset, (char*)outbuf + new_buflen);
}

