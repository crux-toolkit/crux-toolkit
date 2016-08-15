#!/usr/bin/env python


import operator
import socket
from sys import argv
from collections import defaultdict

host_counts = defaultdict(int)
app_counts = defaultdict(int)

for filename in argv[1:]:
    input_file = open(filename)
    for line in input_file:
        if "zip" in line or "gz" in line:
            values = line.split(",")
            date_time = values[0]
            ip_address = values[1]
            file = values[2]
            platform = None
            if "Windows.i386" in line: platform = "Windows 32-bit"
            if "Windows.ix86" in line: platform = "Windows 32-bit"
            if "Windows.AMD64" in line: platform = "Windows 64-bit"
            if "Linux.i686" in line: platform = "Linux 32-bit"
            if "Linux.x86_64" in line: platform = "Linux 64-bit"
            if "Darwin" in line: platform = "Mac OS"
            if "Source" in line: platform = "Source"
            if platform == None:
                print line
            hostname = "Unknown"
            try:
                hostname = socket.gethostbyaddr(ip_address)[0]
            except:
                pass
            host_counts[hostname] += 1
            app_counts[platform] += 1
    input_file.close()

sorted_host_counts = sorted(host_counts.items(), key=operator.itemgetter(1))
for key, count in sorted_host_counts[-20:]:
    print key, count

sorted_app_counts = sorted(app_counts.items(), key=operator.itemgetter(1))
sum = 0
for key, count in sorted_app_counts:
    sum += count
    print key, count
print "Total Downloads", sum
