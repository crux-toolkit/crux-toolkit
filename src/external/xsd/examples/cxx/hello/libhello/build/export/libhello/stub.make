# file      : examples/cxx/hello/libhello/build/export/libhello/stub.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(src_root)/libhello/makefile)

$(call export,\
  l: $(out_root)/libhello/hello.l,\
  cpp-options: $(out_root)/libhello/hello.l.cpp-options)
