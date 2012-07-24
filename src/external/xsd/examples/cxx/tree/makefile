# file      : examples/cxx/tree/makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
# license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))../../../build/bootstrap.make

all_examples := binary caching embedded custom hello library messaging \
mixed multiroot performance polymorphism streaming wildcard compression \
dbxml xpath


build_examples := binary caching embedded custom hello library messaging \
mixed multiroot performance polymorphism streaming wildcard

ifeq ($(xsd_with_zlib),y)
build_examples += compression
endif

ifeq ($(xsd_with_dbxml),y)
build_examples += dbxml
endif

ifeq ($(xsd_with_xqilla),y)
build_examples += xpath
endif

default  := $(out_base)/
install  := $(out_base)/.install
dist     := $(out_base)/.dist
dist-win := $(out_base)/.dist-win
clean    := $(out_base)/.clean

# Build.
#
$(default): $(addprefix $(out_base)/,$(addsuffix /,$(build_examples)))

# Install & Dist.
#
$(install) $(dist) $(dist-win): path := $(subst $(src_root)/,,$(src_base))

$(install): $(addprefix $(out_base)/,$(addsuffix /.install,$(all_examples)))
	$(call install-data,$(src_base)/README,$(install_doc_dir)/xsd/$(path)/README)

$(dist): $(addprefix $(out_base)/,$(addsuffix /.dist,$(all_examples)))
	$(call install-data,$(src_base)/README,$(dist_prefix)/$(path)/README)

$(dist-win): $(addprefix $(out_base)/,$(addsuffix /.dist-win,$(all_examples)))
	$(call install-data,$(src_base)/README,$(dist_prefix)/$(path)/README.txt)
	$(call message,,unix2dos $(dist_prefix)/$(path)/README.txt)

# Clean.
#
$(clean): $(addprefix $(out_base)/,$(addsuffix /.clean,$(build_examples)))

$(call include,$(bld_root)/install.make)

ifneq ($(filter $(MAKECMDGOALS),dist dist-win install),)
$(foreach e,$(all_examples),$(call import,$(src_base)/$e/makefile))
else
$(foreach e,$(build_examples),$(call import,$(src_base)/$e/makefile))
endif
