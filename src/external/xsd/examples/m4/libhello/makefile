# file      : examples/cxx/hello/libhello/makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))build/bootstrap.make

default   := $(out_base)/
install   := $(out_base)/.install
clean     := $(out_base)/.clean
test      := $(out_base)/.test

.PHONY: $(default) $(install) $(test) $(clean)

$(default): $(out_base)/hello/hello.l $(out_base)/test/driver

$(install): $(out_base)/hello/.install

$(clean): $(out_base)/hello/.clean $(out_base)/test/.clean

$(test): $(out_base)/test/.test


# Aliases.
#
ifdef %interactive%

.PHONY: install test clean

install: $(install)
clean: $(clean)
test: $(test)

endif

# Import build rules.
#
$(call import,$(src_base)/hello/makefile)
$(call import,$(src_base)/test/makefile)
