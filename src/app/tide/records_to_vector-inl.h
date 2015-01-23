//  Benjamin Diament
//
//  Templatized function to read a file of records (as in records.{h,cc}) to
//  store in a vector<>.
//
//  The template parameters are:
//  Protobuf:   any protocol buffer message type.
//  VecElement: the type of the vector elements (vector<> actually contains
//              pointers). Usually same as Protobuf, but may be const. Note that
//              VecElement is implicit in the first function parameter.
//
//  Example:
//
//    vector<const pb::Protein*> proteins;
//    ReadRecordsToVector<pb::Protein>(&proteins, filename);
//
#include "records.h"
template<class Protobuf, class VecElement>
bool ReadRecordsToVector(vector<VecElement*>* vec, const string& filename,
			 pb::Header* header = NULL) {
  HeadedRecordReader reader(filename, header);
  while (!reader.Done()) {
    Protobuf* protobuf = new Protobuf;
    reader.Read(protobuf);
    vec->push_back(protobuf);
  }
  if (!reader.OK()) { // Discard everything if we fail.
    for (int i = 0; i < vec->size(); ++i)
      delete (*vec)[i];
    vec->clear();
    return false;
  }
  return true;
}

template<class Protobuf>
int CountRecords(const string& filename) {
  int count = 0;
  HeadedRecordReader reader(filename);
  while (!reader.Done()) {
    Protobuf protobuf;
    reader.Read(&protobuf);
    ++count;
  }
  if (!reader.OK())
    return -1;

  return count;
}
