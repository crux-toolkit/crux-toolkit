#! /usr/bin/python
# Benjamin Diament

# Utility routines for reading and writing a file of protocol buffers as
# records. Each protocol buffer message is encoded into a record which is
# prepended by its length in the varint format. See protocol buffer
# documentation: 
#     http://code.google.com/apis/protocolbuffers/docs/overview.html.
# A '\0' is appended to the end of the file.
# According to the documentation it is inappropriate (also ineffective) to
# encode a series of records as a protocol buffer containing a single repeated
# field. Records provide a very basic way to store vastly many messages in
# a file.

# Any protocol buffer message type may be used, but is the same for the whole
# file. The corresponding program records.h can be used to read and write in
# the same format in C++.

# Example: assume a protocol buffer message "Protein" is defined as follows:
# 
#     message Protein {
#       optional string name = 1;
#       optional int32 id = 2;
#       optional string residues = 3;
#     }
#
# Then we could write some records of Protein to a file thus:
#
#     import records
#     from protoobj import raw_proteins_pb2
#     writer = records.RecordWriter(sys.stdout)
#     protein = raw_proteins_pb2.Protein()
#     protein.name = "My favorite protein"
#     protein.id = 23
#     protein.residues = "EYTRRP"
#     writer.write(protein)
#     ... (write some more records)
#
# Such a file could be read back thus:
#
#    import records
#    from protoobj import raw_proteins_pb2
#
#    for protein in records.Records(sys.stdin, raw_proteins_pb2.Protein):
#      print "id: %s name: %s" % (protein.id, protein.name)
#      print protein.residues
#      print
#
# The additional routines HeadedRecords and HeadedRecordsWriter are
# provided for reading and writing files that begin with a Header.proto
# entry.

import sys, types
from protoobj import header_pb2

# All HeadedRecord files begin with this arbitrarily-chosen length 4
# identifying string.
MAGIC_NUMBER = '\x34\x12\xad\xfe'

def GetSize(file):
  # Read a varint from the file representing the size of the
  # next record, then decode a record of type message_type.
  digit = file.read(1)
  if not digit:
    return
  multiplier = 1
  size = 0
  digit = ord(digit)
  while digit >= 128:
    size += multiplier*(digit - 128)
    multiplier *= 128
    digit = ord(file.read(1))
  size += multiplier*digit
  return size  


def ReadRecord(file, message_type):
  size = GetSize(file)
  if size == 0:
    return None
  message = message_type()
  message.ParseFromString(file.read(size))
  return message


def SkipRecord(file):
  size = GetSize(file)
  file.seek(size, 1)


def SkipManyRecords(file, num):
  for i in range(num):
    SkipRecord(file)


def CheckMagicNumber(file):
  assert(file.read(4) == MAGIC_NUMBER)


def WriteMagicNumber(file):
  file.write(MAGIC_NUMBER)


def Records(file, message_type):
  while True:
    message = ReadRecord(file, message_type)
    if type(message) is types.NoneType:
      return
    yield message


def HeadedRecords(file, message_type):
  # returns a pair: (header, records_iterator)
  CheckMagicNumber(file)
  return (ReadRecord(file, header_pb2.Header), Records(file, message_type))


class RecordWriter:
  def __init__(self, file):
    self.file = file

  def __del__(self):
    if self.file:
      self.file.write('\0')
      self.file.close()

  def TakeFile(self):
    f = self.file
    self.file = None
    return f

  def write(self, record):
    data = record.SerializeToString()
    length = len(data)
    while length > 0:
      digit = length & 127
      length >>= 7
      if length > 0:
        digit += 128
      self.file.write(chr(digit))
    self.file.write(data)

  def WriteMagicNumber(self):
    WriteMagicNumber(self.file)


def HeadedRecordWriter(file, header):
  WriteMagicNumber(file)
  writer = RecordWriter(file)
  writer.write(header)
  return writer
