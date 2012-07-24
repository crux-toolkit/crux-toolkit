# file      : makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
# license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))build/bootstrap.make

default  := $(out_base)/
test     := $(out_base)/.test
install  := $(out_base)/.install
dist     := $(out_base)/.dist
dist-win := $(out_base)/.dist-win
clean    := $(out_base)/.clean
cleandoc := $(out_base)/.cleandoc

$(default): $(out_base)/xsd/      \
            $(out_base)/tests/    \
            $(out_base)/examples/ \
            $(out_base)/documentation/

# Test.
#
$(test): $(out_base)/tests/.test


# Install.
#
$(install): $(out_base)/xsd/.install           \
            $(out_base)/libxsd/.install        \
	    $(out_base)/examples/.install      \
            $(out_base)/documentation/.install
	$(call install-dir,$(src_base)/dist/examples/build,$(install_doc_dir)/xsd/examples/build)
	$(call install-dir,$(src_base)/dist/examples/cxx,$(install_doc_dir)/xsd/examples/cxx,-name makefile)
	$(call install-data,$(src_base)/dist/examples/makefile,$(install_doc_dir)/xsd/examples/makefile)
	$(call install-data,$(src_base)/FLOSSE,$(install_doc_dir)/xsd/FLOSSE)
	$(call install-data,$(src_base)/GPLv2,$(install_doc_dir)/xsd/GPLv2)
	$(call install-data,$(src_base)/LICENSE,$(install_doc_dir)/xsd/LICENSE)
	$(call install-data,$(src_base)/NEWS,$(install_doc_dir)/xsd/NEWS)
	$(call install-data,$(src_base)/dist/README-UNIX,$(install_doc_dir)/xsd/README)


# Dist.
#
$(dist): $(out_base)/xsd/.dist           \
         $(out_base)/libxsd/.dist        \
         $(out_base)/examples/.dist      \
         $(out_base)/documentation/.dist
	$(call install-dir,$(src_base)/dist/examples/build,$(dist_prefix)/examples/build)
	$(call install-dir,$(src_base)/dist/examples/cxx,$(dist_prefix)/examples/cxx,-name makefile)
	$(call install-data,$(src_base)/dist/examples/makefile,$(dist_prefix)/examples/makefile)
	$(call install-data,$(src_base)/dist/README-UNIX,$(dist_prefix)/README)
	$(call install-data,$(src_base)/GPLv2,$(dist_prefix)/GPLv2)
	$(call install-data,$(src_base)/FLOSSE,$(dist_prefix)/FLOSSE)
	$(call install-data,$(src_base)/LICENSE,$(dist_prefix)/LICENSE)
	$(call install-data,$(src_base)/NEWS,$(dist_prefix)/NEWS)
	$(call install-data,$(src_base)/version,$(dist_prefix)/version)

$(dist-win): $(out_base)/xsd/.dist-win           \
             $(out_base)/libxsd/.dist-win        \
             $(out_base)/examples/.dist-win      \
             $(out_base)/documentation/.dist-win
	$(call install-dir,$(src_base)/dist/etc,$(dist_prefix)/etc)
	$(call install-dir,$(src_base)/dist/examples/build,$(dist_prefix)/examples/build)
	$(call install-dir,$(src_base)/dist/examples/cxx,$(dist_prefix)/examples/cxx)
	$(call install-data,$(src_base)/dist/examples/makefile,$(dist_prefix)/examples/makefile)
	$(call install-data,$(src_base)/dist/README-WINDOWS,$(dist_prefix)/README.txt)
	$(call message,,unix2dos $(dist_prefix)/README.txt)
	$(call install-data,$(src_base)/dist/README-UNIX,$(dist_prefix)/README-CYGWIN.txt)
	$(call message,,unix2dos $(dist_prefix)/README-CYGWIN.txt)
	$(call install-data,$(src_base)/GPLv2,$(dist_prefix)/GPLv2.txt)
	$(call message,,unix2dos $(dist_prefix)/GPLv2.txt)
	$(call install-data,$(src_base)/FLOSSE,$(dist_prefix)/FLOSSE.txt)
	$(call message,,unix2dos $(dist_prefix)/FLOSSE.txt)
	$(call install-data,$(src_base)/LICENSE,$(dist_prefix)/LICENSE.txt)
	$(call message,,unix2dos $(dist_prefix)/LICENSE.txt)
	$(call install-data,$(src_base)/NEWS,$(dist_prefix)/NEWS.txt)
	$(call message,,unix2dos $(dist_prefix)/NEWS.txt)
	$(call install-data,$(src_base)/version,$(dist_prefix)/version.txt)
	$(call message,,unix2dos $(dist_prefix)/version.txt)


# Clean.
#
$(clean): $(out_base)/xsd/.clean      \
          $(out_base)/tests/.clean    \
          $(out_base)/examples/.clean

$(cleandoc): $(out_base)/documentation/.cleandoc

$(call include,$(bld_root)/install.make)

$(call import,$(src_base)/xsd/makefile)
$(call import,$(src_base)/libxsd/makefile)
$(call import,$(src_base)/tests/makefile)
$(call import,$(src_base)/examples/makefile)
$(call import,$(src_base)/documentation/makefile)
