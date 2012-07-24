# file      : documentation/cxx/tree/makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
# license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))../../../build/bootstrap.make

docs := guide manual

default  := $(out_base)/
install  := $(out_base)/.install
dist     := $(out_base)/.dist
dist-win := $(out_base)/.dist-win
cleandoc := $(out_base)/.cleandoc

# Build.
#
$(default): $(addprefix $(out_base)/,$(addsuffix /,$(docs)))

# Install & Dist.
#
dist-common := $(out_base)/.dist-common

$(install): path := $(subst $(src_root)/documentation/,,$(src_base))
$(dist-common): path := $(subst $(src_root)/,,$(src_base))

$(install): $(addprefix $(out_base)/,$(addsuffix /.install,$(docs)))
	$(call install-dir,$(src_base)/dbxml,$(install_doc_dir)/xsd/$(path)/dbxml)
	$(call install-dir,$(src_base)/reference,$(install_doc_dir)/xsd/$(path)/reference)

$(dist-common):
	$(call install-dir,$(src_base)/dbxml,$(dist_prefix)/$(path)/dbxml)
	$(call install-dir,$(src_base)/reference,$(dist_prefix)/$(path)/reference)

$(dist): $(dist-common) $(addprefix $(out_base)/,$(addsuffix /.dist,$(docs)))
$(dist-win): $(dist-common) $(addprefix $(out_base)/,$(addsuffix /.dist-win,$(docs)))

# Clean.
#
$(cleandoc): $(addprefix $(out_base)/,$(addsuffix /.cleandoc,$(docs)))

$(call include,$(bld_root)/install.make)

$(foreach m,$(docs),$(call import,$(src_base)/$m/makefile))
