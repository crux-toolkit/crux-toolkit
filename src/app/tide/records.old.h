#ifndef RECORDS_H
#define RECORDS_H

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string>
#include <google/protobuf/message.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/coded_stream.h>

using namespace std;

#define UINT32_MAX 0xfffffffful

class RecordWriter {
 public:
 RecordWriter(const string& filename)
   : raw_output_(NULL), coded_output_(NULL) {
    fd_ = open(filename.c_str(), O_CREAT | O_TRUNC | O_WRONLY, 0644);
    if (fd_ < 0)
      return;
    raw_output_ = new google::protobuf::io::FileOutputStream(fd_);
    coded_output_ = new google::protobuf::io::CodedOutputStream(raw_output_);
  }

  ~RecordWriter() {
    coded_output_->WriteVarint32(0); // end-of-records marker
    delete coded_output_;
    delete raw_output_;
    close(fd_);
  }

  // client should check once after construction
  bool OK() const { return NULL != coded_output_; }

  bool Write(const google::protobuf::Message* message) {
    if (!coded_output_->WriteVarint32(message->ByteSize()))
      return false;
    return message->SerializeWithCachedSizes(coded_output_);
  }

 private:
  int fd_;
  google::protobuf::io::ZeroCopyOutputStream* raw_output_;
  google::protobuf::io::CodedOutputStream* coded_output_;
};

class RecordReader {
 public:
  RecordReader(const string& filename)
    : raw_input_(NULL), coded_input_(NULL), size_(UINT32_MAX), valid_(false) {
    fd_ = open(filename.c_str(), O_RDONLY);
    if (fd_ < 0)
      return;
    raw_input_ = new google::protobuf::io::FileInputStream(fd_);
    coded_input_ = new google::protobuf::io::CodedInputStream(raw_input_);
    valid_ = true;
  }

  ~RecordReader() {
    delete coded_input_;
    delete raw_input_;
    close(fd_);
  }

  bool OK() const { return valid_; }

  bool Done() {
    if (!valid_)
      return true;
    assert(coded_input_->BytesUntilLimit() == -1);
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
    if (!coded_input_->ExpectAtEnd())
      return valid_ = false;
    coded_input_->PopLimit(limit);
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

#endif // RECORDS_H
