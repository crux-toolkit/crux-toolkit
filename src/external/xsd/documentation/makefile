# file      : documentation/makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
# license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))../build/bootstrap.make

default  := $(out_base)/
install  := $(out_base)/.install
dist     := $(out_base)/.dist
dist-win := $(out_base)/.dist-win
cleandoc := $(out_base)/.cleandoc

$(default): $(out_base)/cxx/

# Install.
#
$(install): $(out_base)/cxx/.install	
	$(call install-data,$(src_base)/default.css,$(install_doc_dir)/xsd/default.css)
	$(call install-data,$(src_base)/future.xhtml,$(install_doc_dir)/xsd/future.xhtml)
	$(call install-data,$(src_base)/schema-authoring-guide.xhtml,$(install_doc_dir)/xsd/schema-authoring-guide.xhtml)
	$(call install-data,$(src_base)/xsd.xhtml,$(install_doc_dir)/xsd/xsd.xhtml)
	$(call install-data,$(src_base)/custom-literals.xsd,$(install_doc_dir)/xsd/custom-literals.xsd)
	$(call install-data,$(src_base)/xsd.1,$(install_man_dir)/man1/xsd.1)

# Dist.
#
dist-common := $(out_base)/.dist-common

$(dist-common):
	$(call install-data,$(src_base)/default.css,$(dist_prefix)/documentation/default.css)
	$(call install-data,$(src_base)/xsd.xhtml,$(dist_prefix)/documentation/xsd.xhtml)
	$(call install-data,$(src_base)/xsd.1,$(dist_prefix)/documentation/xsd.1)
	$(call install-data,$(src_base)/future.xhtml,$(dist_prefix)/documentation/future.xhtml)
	$(call install-data,$(src_base)/schema-authoring-guide.xhtml,$(dist_prefix)/documentation/schema-authoring-guide.xhtml)
	$(call install-data,$(src_base)/custom-literals.xsd,$(dist_prefix)/documentation/custom-literals.xsd)

$(dist): $(dist-common) $(out_base)/cxx/.dist
$(dist-win): $(dist-common) $(out_base)/cxx/.dist-win

# Clean.
#
$(cleandoc): $(src_base)/cxx/.cleandoc

$(call include,$(bld_root)/install.make)

$(call import,$(src_base)/cxx/makefile)
