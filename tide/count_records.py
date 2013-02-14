#! /usr/bin/python
# Benjamin Diament

import sys, types
from protoobj import header_pb2

# All HeadedRecord files begin with this arbitrarily-chosen length 4
# identifying string.
MAGIC_NUMBER = '\x34\x12\xad\xfe'

def ReadRecordLen(file):
  # Read a varint from the file representing the size of the next record
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

def CheckMagicNumber(file):
  assert(file.read(4) == MAGIC_NUMBER)

file = open(sys.argv[1], "rb")
CheckMagicNumber(file)
count = 0
while True:
  size = ReadRecordLen(file)
  if (size == 0):
    break
  file.seek(size, 1)
  count += 1
  if (count % 1e6 == 0):
    sys.stderr.write("Count = %d million. " % (count / 1e6))
    sys.stderr.write("Position = %s\n" % file.tell())
print "Final file position: %s. %s Records." % (file.tell(), count - 1)
assert(file.read(1) == '')
