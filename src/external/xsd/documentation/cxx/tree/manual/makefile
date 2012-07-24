# file      : documentation/cxx/tree/manual/makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
# license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))../../../../build/bootstrap.make

default  := $(out_base)/
install  := $(out_base)/.install
dist     := $(out_base)/.dist
dist-win := $(out_base)/.dist-win
cleandoc := $(out_base)/.cleandoc

# Build.
#
$(default): $(out_base)/cxx-tree-manual.ps $(out_base)/cxx-tree-manual.pdf


$(out_base)/cxx-tree-manual.ps: $(src_base)/index.xhtml   \
                                $(src_base)/manual.html2ps \
                                | $(out_base)/.
	$(call message,html2ps $<,html2ps -f $(src_base)/manual.html2ps -o $@ $<)

$(out_base)/cxx-tree-manual.pdf: $(out_base)/cxx-tree-manual.ps | $(out_base)/.
	$(call message,ps2pdf $<,ps2pdf14 $< $@)

# Install & Dist.
#
$(install): path := $(subst $(src_root)/documentation/,,$(src_base))
$(dist): path := $(subst $(src_root)/,,$(src_base))

$(install): $(out_base)/cxx-tree-manual.ps $(out_base)/cxx-tree-manual.pdf
	$(call install-data,$(src_base)/index.xhtml,$(install_doc_dir)/xsd/$(path)/index.xhtml)
	$(call install-data,$(out_base)/cxx-tree-manual.ps,$(install_doc_dir)/xsd/$(path)/cxx-tree-manual.ps)
	$(call install-data,$(out_base)/cxx-tree-manual.pdf,$(install_doc_dir)/xsd/$(path)/cxx-tree-manual.pdf)

$(dist): $(out_base)/cxx-tree-manual.ps $(out_base)/cxx-tree-manual.pdf
	$(call install-data,$(src_base)/index.xhtml,$(dist_prefix)/$(path)/index.xhtml)
	$(call install-data,$(out_base)/cxx-tree-manual.ps,$(dist_prefix)/$(path)/cxx-tree-manual.ps)
	$(call install-data,$(out_base)/cxx-tree-manual.pdf,$(dist_prefix)/$(path)/cxx-tree-manual.pdf)

$(dist-win): $(dist)


# Clean
#
$(cleandoc):
	$(call message,rm $$1,rm -f $$1,$(out_base)/cxx-tree-manual.ps)
	$(call message,rm $$1,rm -f $$1,$(out_base)/cxx-tree-manual.pdf)

# How to.
#
$(call include,$(bld_root)/install.make)
