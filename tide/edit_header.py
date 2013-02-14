#!/usr/bin/env python

import sys, shutil, code
import records
from protoobj import header_pb2

infile = open(sys.argv[1])
records.CheckMagicNumber(infile)
header = records.ReadRecord(infile, header_pb2.Header)

print "header =", header

def Done(filename = sys.argv[1] + ".newheader"):
  out = open(filename, "wb")
  writer = records.HeadedRecordWriter(out, header)
  shutil.copyfileobj(infile, writer.TakeFile(), 2**25)

banner = """
Edit value of header and call Done(filename).
Default filename appends ".newheader" to old filename.
"""

code.interact(banner, local=globals())
