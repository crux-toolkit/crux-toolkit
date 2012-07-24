# file      : tests/cxx/makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
# license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))../../build/bootstrap.make

default   := $(out_base)/
test      := $(out_base)/.test
clean     := $(out_base)/.clean

$(default): $(out_base)/parser/ $(out_base)/tree/
$(test): $(out_base)/parser/.test $(out_base)/tree/.test
$(clean): $(out_base)/parser/.clean $(out_base)/tree/.clean

$(call import,$(src_base)/parser/makefile)
$(call import,$(src_base)/tree/makefile)
