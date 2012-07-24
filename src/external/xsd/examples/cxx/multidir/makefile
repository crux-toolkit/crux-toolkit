# file      : examples/cxx/multidir/makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))build/bootstrap.make

cxx_tun  := driver.cxx foo/foo.cxx bar/bar.cxx
cxx_obj  := $(addprefix $(out_base)/,$(cxx_tun:.cxx=.o))
cxx_od   := $(cxx_obj:.o=.o.d)

driver   := $(out_base)/driver
clean    := $(out_base)/.clean


# Build.
#
$(driver): $(cxx_obj)

$(call -include,$(cxx_od))


# Clean
#
.PHONY: $(clean)

$(clean): $(driver).o.clean           \
  $(addsuffix .cxx.clean,$(cxx_obj))  \
  $(addsuffix .cxx.clean,$(cxx_od))


# Aliases.
#
ifdef %interactive%

.PHONY: clean

clean: $(clean)

endif


# How to.
#
$(call include,$(bld_root)/cxx/o-e.make)
$(call include,$(bld_root)/cxx/cxx-o.make)
$(call include,$(bld_root)/cxx/cxx-d.make)
