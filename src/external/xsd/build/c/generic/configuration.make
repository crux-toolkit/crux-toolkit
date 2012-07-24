# file      : build/c/generic/configuration.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(bld_root)/c/generic/configuration-rules.make,$(dcf_root))

# Static configuration.
#
ifneq ($(bld_root),$(scf_root))
$(call -include,$(scf_root)/c/generic/configuration-static.make)
endif

# Dynamic configuration.
#
c_generic :=

$(call -include,$(dcf_root)/c/generic/configuration-dynamic.make)

ifndef c_generic
.NOTPARALLEL:
endif
