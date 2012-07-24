# file      : examples/cxx/parser/performance/makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
# license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))../../../../build/bootstrap.make

xsd := test.xsd
cxx := driver.cxx time.cxx

obj := $(addprefix $(out_base)/,$(cxx:.cxx=.o) $(xsd:.xsd=-pskel.o))
dep := $(obj:.o=.o.d)

driver   := $(out_base)/driver
install  := $(out_base)/.install
dist     := $(out_base)/.dist
dist-win := $(out_base)/.dist-win
clean    := $(out_base)/.clean


# Import.
#
$(call import,\
  $(scf_root)/import/libxerces-c/stub.make,\
  l: xerces_c.l,cpp-options: xerces_c.l.cpp-options)


# Build.
#
$(driver): $(obj) $(xerces_c.l)

$(obj) $(dep): cpp_options := -I$(src_root)/libxsd
$(obj) $(dep): $(xerces_c.l.cpp-options)

genf := $(xsd:.xsd=-pskel.hxx) $(xsd:.xsd=-pskel.ixx) $(xsd:.xsd=-pskel.cxx)
gen  := $(addprefix $(out_base)/,$(genf))

$(gen): xsd := $(out_root)/xsd/xsd
$(gen): xsd_options := --generate-inline
$(gen): $(out_root)/xsd/xsd

$(call include-dep,$(dep))

# Convenience alias for default target.
#
$(out_base)/: $(driver)


# Install & Dist.
#
dist-common := $(out_base)/.dist-common

$(install) $(dist) $(dist-win) $(dist-common): path := $(subst $(src_root)/,,$(src_base))

$(install):
	$(call install-data,$(src_base)/README,$(install_doc_dir)/xsd/$(path)/README)
	$(call install-data,$(src_base)/driver.cxx,$(install_doc_dir)/xsd/$(path)/driver.cxx)
	$(call install-data,$(src_base)/test.xsd,$(install_doc_dir)/xsd/$(path)/test.xsd)
	$(call install-data,$(src_base)/test-5k.xml,$(install_doc_dir)/xsd/$(path)/test-5k.xml)
	$(call install-data,$(src_base)/test-50k.xml,$(install_doc_dir)/xsd/$(path)/test-50k.xml)
	$(call install-data,$(src_base)/test-500k.xml,$(install_doc_dir)/xsd/$(path)/test-500k.xml)
	$(call install-data,$(src_base)/time.hxx,$(install_doc_dir)/xsd/$(path)/time.hxx)
	$(call install-data,$(src_base)/time.cxx,$(install_doc_dir)/xsd/$(path)/time.cxx)
	$(call install-data,$(src_base)/gen.cxx,$(install_doc_dir)/xsd/$(path)/gen.cxx)

$(dist-common):
	$(call install-data,$(src_base)/driver.cxx,$(dist_prefix)/$(path)/driver.cxx)
	$(call install-data,$(src_base)/test.xsd,$(dist_prefix)/$(path)/test.xsd)
	$(call install-data,$(src_base)/test-5k.xml,$(dist_prefix)/$(path)/test-5k.xml)
	$(call install-data,$(src_base)/test-50k.xml,$(dist_prefix)/$(path)/test-50k.xml)
	$(call install-data,$(src_base)/test-500k.xml,$(dist_prefix)/$(path)/test-500k.xml)
	$(call install-data,$(src_base)/time.hxx,$(dist_prefix)/$(path)/time.hxx)
	$(call install-data,$(src_base)/time.cxx,$(dist_prefix)/$(path)/time.cxx)
	$(call install-data,$(src_base)/gen.cxx,$(dist_prefix)/$(path)/gen.cxx)

$(dist): $(dist-common)
	$(call install-data,$(src_base)/README,$(dist_prefix)/$(path)/README)

$(dist-win): $(dist-common)
	$(call install-data,$(src_base)/README,$(dist_prefix)/$(path)/README.txt)
	$(call message,,unix2dos $(dist_prefix)/$(path)/README.txt)


# Clean.
#
$(clean): $(driver).o.clean                                 \
  $(addsuffix .cxx.clean,$(obj))                            \
  $(addsuffix .cxx.clean,$(dep))                            \
  $(addprefix $(out_base)/,$(xsd:.xsd=-pskel.cxx.xsd.clean))

# Generated .gitignore.
#
ifeq ($(out_base),$(src_base))
$(gen): | $(out_base)/.gitignore
$(driver): | $(out_base)/.gitignore

$(out_base)/.gitignore: files := driver $(genf)
$(clean): $(out_base)/.gitignore.clean

$(call include,$(bld_root)/git/gitignore.make)
endif

# How to.
#
$(call include,$(bld_root)/cxx/o-e.make)
$(call include,$(bld_root)/cxx/cxx-o.make)
$(call include,$(bld_root)/cxx/cxx-d.make)
$(call include,$(bld_root)/install.make)
$(call include,$(scf_root)/xsd/parser/xsd-cxx.make)


# Dependencies.
#
$(call import,$(src_root)/xsd/makefile)
