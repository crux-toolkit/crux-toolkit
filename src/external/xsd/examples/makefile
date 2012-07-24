# file      : examples/makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
# license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))../build/bootstrap.make

default  := $(out_base)/
install  := $(out_base)/.install
dist     := $(out_base)/.dist
dist-win := $(out_base)/.dist-win
clean    := $(out_base)/.clean

$(default): $(out_base)/cxx/parser/ $(out_base)/cxx/tree/
$(install): $(out_base)/cxx/parser/.install $(out_base)/cxx/tree/.install
$(dist): $(out_base)/cxx/parser/.dist $(out_base)/cxx/tree/.dist
$(dist-win): $(out_base)/cxx/parser/.dist-win $(out_base)/cxx/tree/.dist-win
$(clean): $(out_base)/cxx/parser/.clean $(out_base)/cxx/tree/.clean

$(call import,$(src_base)/cxx/parser/makefile)
$(call import,$(src_base)/cxx/tree/makefile)
