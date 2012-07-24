# file      : build/c/generic/configuration-sl.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(bld_root)/c/generic/configuration-sl-rules.make,$(dcf_root))

# Dynamic configuration.
#
c_generic_pic_option :=

$(call -include,$(dcf_root)/c/generic/configuration-sl-dynamic.make)

ifdef c_generic_pic_option

$(out_root)/%: c_generic_pic_option := $(c_generic_pic_option)
$(out_root)/%: c_generic_shared_option := $(c_generic_shared_option)

else
.NOTPARALLEL:
endif
