# file      : build/ld/configuration-lib.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(bld_root)/ld/configuration-lib-rules.make,$(dcf_root))

# Static configuration.
#
$(call include,$(bld_root)/ld/configuration-lib-static.make)

ifneq ($(bld_root),$(scf_root))
$(call -include,$(scf_root)/ld/configuration-lib-static.make)
endif

# Dynamic configuration.
#
ld_lib_type :=

$(call -include,$(dcf_root)/ld/configuration-lib-dynamic.make)

ifdef ld_lib_type

ifeq ($(ld_lib_type),archive)

ld_lib_suffix := $(ld_lib_archive_suffix)

else

ld_lib_suffix := $(ld_lib_shared_suffix)

endif

$(out_root)/%: ld_lib_type := $(ld_lib_type)
$(out_root)/%: ld_lib_prefix := $(ld_lib_prefix)
$(out_root)/%: ld_lib_suffix := $(ld_lib_suffix)

else

.NOTPARALLEL:

endif
