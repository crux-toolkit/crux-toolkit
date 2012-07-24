# file      : examples/cxx/hello/hello/makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))build/bootstrap.make

cxx_tun := hello.cxx
cxx_obj := $(addprefix $(out_base)/,$(cxx_tun:.cxx=.o))
cxx_od  := $(cxx_obj:.o=.o.d)

hello   := $(out_base)/hello
clean   := $(out_base)/.clean
install := $(out_base)/.install


# Secure default target.
#
#$(hello):


# Import information about libhello.
#
$(call import,\
  $(scf_root)/import/libhello/stub.make,\
  l: hello.l,cpp-options: hello.l.cpp-options)


# Build.
#
$(hello): $(cxx_obj) $(hello.l)

$(cxx_obj) $(cxx_od): $(hello.l.cpp-options)

$(call -include,$(cxx_od))


# Clean.
#
.PHONY: $(clean)

$(clean): $(hello).o.clean           \
  $(addsuffix .cxx.clean,$(cxx_obj)) \
  $(addsuffix .cxx.clean,$(cxx_od))


# Install.
#
.PHONY: $(install)

$(install): $(hello)
	$(call install-exec,$<,$(install_bin_dir)/hello)

# Aliases.
#
ifdef %interactive%

.PHONY: clean install

clean: $(clean)
install: $(install)

endif


# How to.
#
$(call include,$(bld_root)/cxx/o-e.make)
$(call include,$(bld_root)/cxx/cxx-o.make)
$(call include,$(bld_root)/cxx/cxx-d.make)
$(call include,$(bld_root)/install.make)
