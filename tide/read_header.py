#!/usr/bin/env python

import sys
import records
from protoobj import header_pb2

file = open(sys.argv[1])
records.CheckMagicNumber(file)
header = records.ReadRecord(file, header_pb2.Header)
file.close()

print header
