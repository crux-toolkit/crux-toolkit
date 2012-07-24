# file      : examples/m4/hello/libhello/hello/makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))../build/bootstrap.make

cxx_tun   := hello.cxx
cxx_obj   := $(addprefix $(out_base)/,$(cxx_tun:.cxx=.o))
cxx_od    := $(cxx_obj:.o=.o.d)

hello.l             := $(out_base)/hello.l
hello.l.cpp-options := $(out_base)/hello.l.cpp-options

install   := $(out_base)/.install
clean     := $(out_base)/.clean


# Build.
#
$(hello.l): $(cxx_obj)

$(cxx_obj) $(cxx_od): $(hello.l.cpp-options)

#@@ not necessary
#$(out_base)/hello.hxx: $(src_base)/hello.hxx.m4

$(hello.l.cpp-options): prefix := hello/ $(out_root)/
$(hello.l.cpp-options): value := -I$(src_root) -I$(out_root)

$(call -include,$(cxx_od))


# Clean.
#
.PHONY: $(clean)

$(clean): $(hello.l).o.clean       \
  $(addsuffix .cxx.clean,$(cxx_obj)) \
  $(addsuffix .cxx.clean,$(cxx_od))  \
  $(out_base)/hello.hxx.m4.clean \
  $(hello.l.cpp-options).clean


# !! broken: need to install hello.hxx from out_root Install.
#
.PHONY: $(install)

$(install): $(hello.l)
	$(call install-lib,$<,$(install_lib_dir)/$(ld_lib_prefix)hello$(ld_lib_suffix))
	$(call install-dir,$(src_base),$(install_inc_dir)/hello,-name '*.hxx')

# Aliases.
#
ifdef %interactive%

.PHONY: clean install

clean: $(clean)
install: $(install)

endif


# How to.
#
$(call include,$(bld_root)/cxx/o-l.make)
$(call include,$(bld_root)/cxx/cxx-o.make)
$(call include,$(bld_root)/cxx/cxx-d.make)
$(call include,$(bld_root)/m4/m4.make)
$(call include,$(bld_root)/install.make)
