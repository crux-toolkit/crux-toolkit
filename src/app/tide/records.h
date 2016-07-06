// Benjamin Diament
//
// Utility routines for reading and writing a file of protocol buffers as
// records. Each protocol buffer message is encoded into a record which is
// prepended by its length in the varint format. See protocol buffer
// documentation: 
//     http://code.google.com/apis/protocolbuffers/docs/overview.html.
// A '\0' is appended to the end of the file.
// According to the documentation it is inappropriate (also ineffective) to
// encode a series of records as a protocol buffer containing a single repeated
// field. Records provide a very basic way to store vastly many messages in
// a file.
//
// The corresponding program records.py can be used to read and write in
// the same format in Python, provided that all messages are the same type.
//
// See peptide_downselect.cc for a very simple example of the use of records.
//
// Each file of records is prepended with a 4-byte arbitrarily-chosen string,
// called a magic number, to provide a check that a given file is indeed a
// file of records.
//
// HeadedRecordReader and HeadedRecordWriter provide an interface for 
// prepending a file of records with a single record of type Header.proto,
// containing client-provided meta-information about the file.
// All tide generated files are expected to be HeadedRecords.
//
// RecordReader and RecordWriter make use of the classes:
//     google::protobuf::io::FileInputStream
//     google::protobuf::io::FileOutputStream
//     google::protobuf::io::CodedInputStream
//     google::protobuf::io::CodedOutputStream,
// which are described in the protocol buffer documentation.
//
// TODO 259: This class probably should be rewritten, paying more
// attention to the various idioms for use of ZeroCopyStream and
// CodedStream! Things could be greatly simplified.
//
// Note that CodedInputStream isn't built to handle large streams of
// input, so it should be reconstructed at each record. Perhaps the
// underlying ZeroCopyStream should handle EOF determination


#ifndef RECORDS_H
#define RECORDS_H

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <fcntl.h>
#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#endif
#include <iostream>
#include <string>
#include <google/protobuf/message.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/coded_stream.h>
#include "header.pb.h"
#include "io/carp.h"

using namespace std;

#ifndef UINT32_MAX
#define UINT32_MAX 0xfffffffful
#endif
#define MAGIC_NUMBER  0xfead1234ul

class RecordWriter {
 public:
  explicit RecordWriter(const string& filename, int buf_size = -1)
    : raw_output_(NULL), coded_output_(NULL) {
    if ((fd_ = open(filename.c_str(), O_CREAT | O_TRUNC | O_WRONLY, 0644)) < 0) {
      carp(CARP_FATAL, "Couldn't open file %s for write (errno %d: %s).",
	   filename.c_str(), errno, strerror(errno));
      return;
    }
    raw_output_ = new google::protobuf::io::FileOutputStream(fd_, buf_size);
    Init();
  }
  
  explicit RecordWriter(google::protobuf::io::ZeroCopyOutputStream* raw_output)
    : fd_(-1), raw_output_(raw_output) {
    Init();
    raw_output_ = NULL; // we do not own (and will not delete) raw_output
  }

  ~RecordWriter() {
    coded_output_->WriteVarint32(0); // end-of-records marker
    delete coded_output_;
    delete raw_output_;
    if (fd_ > -1)
      close(fd_);
  }

  // client should check once after construction
  bool OK() const { return NULL != coded_output_; }

  bool Write(const google::protobuf::Message* message) {
    coded_output_->WriteVarint32(message->ByteSize());
    if (coded_output_->HadError()) {
      delete coded_output_;
      coded_output_ = NULL;
      return false;
    }
    message->SerializeWithCachedSizes(coded_output_);
    return !coded_output_->HadError();
  }

 private:
  void Init() {
    coded_output_ = new google::protobuf::io::CodedOutputStream(raw_output_);
    coded_output_->WriteLittleEndian32(MAGIC_NUMBER);
    if (coded_output_->HadError()) {
      delete coded_output_;
      coded_output_ = NULL;
    }
  }

  int fd_;
  google::protobuf::io::ZeroCopyOutputStream* raw_output_;
  google::protobuf::io::CodedOutputStream* coded_output_;
};


class RecordReader {
 public:
  explicit RecordReader(const string& filename, int buf_size = -1)
    : raw_input_(NULL), coded_input_(NULL), size_(UINT32_MAX), valid_(false) {
    fd_ = open(filename.c_str(), O_RDONLY);
    if (fd_ < 0)
      return;
    raw_input_ = new google::protobuf::io::FileInputStream(fd_, buf_size);
    google::protobuf::io::CodedInputStream coded_input(raw_input_);
    google::protobuf::uint32 magic_number;
    if (coded_input.ReadLittleEndian32(&magic_number) 
        && magic_number == MAGIC_NUMBER)
      valid_ = true;
  }

  ~RecordReader() {
    if (coded_input_)
      delete coded_input_;
    delete raw_input_;
    if (fd_ >= 0)
      close(fd_);
  }

  bool OK() const { return valid_; }

  bool Done() {
    if (!valid_)
      return true;
    coded_input_ = new google::protobuf::io::CodedInputStream(raw_input_);
    if (!coded_input_->ReadVarint32(&size_))
      return valid_ = false;
    return (size_ == 0);
    // TODO 260: there should be an easy way to tell if we're at the end
    // of the file, but I can't seem to find a reliable way to do that
    // easily...
  }

  bool Read(google::protobuf::Message* message) {
    if (!valid_)
      return false;
    assert(size_ != UINT32_MAX);
    google::protobuf::io::CodedInputStream::Limit limit
      = coded_input_->PushLimit(size_);
    if (!message->ParseFromCodedStream(coded_input_))
      return valid_ = false;
    if (!coded_input_->ConsumedEntireMessage() ||
        coded_input_->BytesUntilLimit() != 0)
      return valid_ = false;
    coded_input_->PopLimit(limit); // for completeness; perhaps remove
    delete coded_input_;
    coded_input_ = NULL;
    size_ = UINT32_MAX;
    return true;
  }

 private:
  int fd_;
  google::protobuf::io::ZeroCopyInputStream* raw_input_;
  google::protobuf::io::CodedInputStream* coded_input_;
  google::protobuf::uint32 size_;
  bool valid_;
};

class HeadedRecordWriter {
 public:
  HeadedRecordWriter(const string& filename, const pb::Header& header,
                     int buf_size = -1) 
    : writer_(filename, buf_size) {
    if (!writer_.OK())
      carp(CARP_FATAL, "Cannot create the file %s\n", filename.c_str());
    Write(&header);
  }

  HeadedRecordWriter(google::protobuf::io::ZeroCopyOutputStream* raw_output,
                     const pb::Header& header) 
    : writer_(raw_output) {
    Write(&header);
  }

  RecordWriter* Writer() { return &writer_; }

  // client should check once after construction
  bool OK() const { return writer_.OK(); }

  bool Write(const google::protobuf::Message* message) {
    return writer_.Write(message);
  }

 private:
  RecordWriter writer_;
};

class HeadedRecordReader {
 public:
  HeadedRecordReader(const string& filename, pb::Header* header = NULL,
                     int buf_size = -1)
    : reader_(filename, buf_size), header_(header),
    del_header_(header == NULL) {
    if (header == NULL)
      header_ = new pb::Header;
    if (!Done())
      Read(header_);
  }
  ~HeadedRecordReader() { if (del_header_) delete header_; }

  RecordReader* Reader() { return &reader_; }

  bool OK() const { return reader_.OK(); }
  bool Done() { return reader_.Done(); }
  bool Read(google::protobuf::Message* message) { 
    return reader_.Read(message);
  }
  const pb::Header* GetHeader() const { return header_; }

 private:
  void ReadMagicNumber();
  RecordReader reader_;
  pb::Header* header_;
  bool del_header_;
};

#endif // RECORDS_H
