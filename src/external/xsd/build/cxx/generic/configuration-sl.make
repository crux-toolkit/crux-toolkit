# file      : build/cxx/generic/configuration-sl.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(bld_root)/cxx/generic/configuration-sl-rules.make,$(dcf_root))

# Dynamic configuration.
#
cxx_generic_pic_option :=

$(call -include,$(dcf_root)/cxx/generic/configuration-sl-dynamic.make)

ifdef cxx_generic_pic_option

$(out_root)/%: cxx_generic_pic_option := $(cxx_generic_pic_option)
$(out_root)/%: cxx_generic_shared_option := $(cxx_generic_shared_option)

else
.NOTPARALLEL:
endif
