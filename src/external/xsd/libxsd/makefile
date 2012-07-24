# file      : libxsd/makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
# license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))../build/bootstrap.make

install  := $(out_base)/.install
dist     := $(out_base)/.dist
dist-win := $(out_base)/.dist-win


# Install.
#
$(install):
	$(call install-dir,$(src_base)/xsd,$(install_inc_dir)/xsd)
	$(call install-data,$(src_base)/FLOSSE,$(install_doc_dir)/libxsd/FLOSSE)
	$(call install-data,$(src_base)/GPLv2,$(install_doc_dir)/libxsd/GPLv2)
	$(call install-data,$(src_base)/LICENSE,$(install_doc_dir)/libxsd/LICENSE)
	$(call install-data,$(src_base)/README,$(install_doc_dir)/libxsd/README)

# Dist.
#
dist-common := $(out_base)/.dist-common

$(dist-common):
	$(call install-dir,$(src_base)/xsd,$(dist_prefix)/libxsd/xsd)

$(dist): $(dist-common)
	$(call install-data,$(src_base)/GPLv2,$(dist_prefix)/libxsd/GPLv2)
	$(call install-data,$(src_base)/FLOSSE,$(dist_prefix)/libxsd/FLOSSE)
	$(call install-data,$(src_base)/LICENSE,$(dist_prefix)/libxsd/LICENSE)
	$(call install-data,$(src_base)/README,$(dist_prefix)/libxsd/README)

$(dist-win): $(dist-common)
	$(call install-data,$(src_base)/GPLv2,$(dist_prefix)/libxsd/GPLv2.txt)
	$(call message,,unix2dos $(dist_prefix)/libxsd/GPLv2.txt)
	$(call install-data,$(src_base)/FLOSSE,$(dist_prefix)/libxsd/FLOSSE.txt)
	$(call message,,unix2dos $(dist_prefix)/libxsd/FLOSSE.txt)
	$(call install-data,$(src_base)/LICENSE,$(dist_prefix)/libxsd/LICENSE.txt)
	$(call message,,unix2dos $(dist_prefix)/libxsd/LICENSE.txt)
	$(call install-data,$(src_base)/README,$(dist_prefix)/libxsd/README.txt)
	$(call message,,unix2dos $(dist_prefix)/libxsd/README.txt)

$(call include,$(bld_root)/install.make)
