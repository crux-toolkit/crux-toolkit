# file      : examples/cxx/parser/makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
# license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))../../../build/bootstrap.make

examples := generated hello library mixin mixed multiroot performance \
polymorphism polyroot wildcard

default  := $(out_base)/
install  := $(out_base)/.install
dist     := $(out_base)/.dist
dist-win := $(out_base)/.dist-win
clean    := $(out_base)/.clean

# Build
#
$(default): $(addprefix $(out_base)/,$(addsuffix /,$(examples)))


# Install & Dist.
#
$(install) $(dist) $(dist-win): path := $(subst $(src_root)/,,$(src_base))

$(install): $(addprefix $(out_base)/,$(addsuffix /.install,$(examples)))
	$(call install-data,$(src_base)/README,$(install_doc_dir)/xsd/$(path)/README)

$(dist): $(addprefix $(out_base)/,$(addsuffix /.dist,$(examples)))
	$(call install-data,$(src_base)/README,$(dist_prefix)/$(path)/README)

$(dist-win): $(addprefix $(out_base)/,$(addsuffix /.dist-win,$(examples)))
	$(call install-data,$(src_base)/README,$(dist_prefix)/$(path)/README.txt)
	$(call message,,unix2dos $(dist_prefix)/$(path)/README.txt)

# Clean
#
$(clean): $(addprefix $(out_base)/,$(addsuffix /.clean,$(examples)))

$(call include,$(bld_root)/install.make)

$(foreach e,$(examples),$(call import,$(src_base)/$e/makefile))
