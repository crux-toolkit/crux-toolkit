# file      : examples/cxx/hello/libhello/test/makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))../build/bootstrap.make

cxx_tun := driver.cxx
cxx_obj := $(addprefix $(out_base)/,$(cxx_tun:.cxx=.o))
cxx_od  := $(cxx_obj:.o=.o.d)

hello.l             := $(out_root)/libhello/hello.l
hello.l.cpp-options := $(out_root)/libhello/hello.l.cpp-options

driver  := $(out_base)/driver
clean   := $(out_base)/.clean
test    := $(out_base)/.test


# Build.
#
$(driver): $(cxx_obj) $(hello.l)

$(cxx_obj) $(cxx_od): $(hello.l.cpp-options)

$(call -include,$(cxx_od))


# Test.
#
.PHONY: $(test)

$(test): $(driver)
	$<


# Clean.
#
.PHONY: $(clean)

$(clean): $(driver).o.clean          \
  $(addsuffix .cxx.clean,$(cxx_obj)) \
  $(addsuffix .cxx.clean,$(cxx_od))


# Aliases.
#
ifdef %interactive%

.PHONY: clean test

test: $(test)
clean: $(clean)

endif


# How to.
#
$(call include,$(bld_root)/cxx/o-e.make)
$(call include,$(bld_root)/cxx/cxx-o.make)
$(call include,$(bld_root)/cxx/cxx-d.make)


# Import build rules.
#
$(call import,$(src_root)/libhello/makefile)
