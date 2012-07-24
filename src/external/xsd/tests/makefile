# file      : tests/makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
# license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))../build/bootstrap.make

default   := $(out_base)/
test      := $(out_base)/.test
clean     := $(out_base)/.clean

$(default): $(out_base)/cxx/
$(test): $(out_base)/cxx/.test
$(clean): $(out_base)/cxx/.clean

$(call import,$(src_base)/cxx/makefile)
