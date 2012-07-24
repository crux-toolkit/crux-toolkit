# file      : build/c/configuration.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(bld_root)/c/configuration-rules.make,$(dcf_root))

# Static configuration.
#
$(call include,$(bld_root)/c/configuration-static.make)

# Dynamic configuration.
#
c_id       :=
c_optimize :=
c_debug    :=
c_rpath    :=

c_pp_extra_options :=
c_extra_options    :=
c_ld_extra_options :=
c_extra_libs       :=
c_extra_lib_paths  :=

$(call -include,$(dcf_root)/c/configuration-dynamic.make)

ifdef c_id

$(out_root)/%: c_id       := $(c_id)
$(out_root)/%: c_optimize := $(c_optimize)
$(out_root)/%: c_debug    := $(c_debug)
$(out_root)/%: c_rpath    := $(c_rpath)

$(out_root)/%: c_pp_extra_options := $(c_pp_extra_options)
$(out_root)/%: c_extra_options    := $(c_extra_options)
$(out_root)/%: c_ld_extra_options := $(c_ld_extra_options)
$(out_root)/%: c_extra_libs       := $(c_extra_libs)
$(out_root)/%: c_extra_lib_paths  := $(c_extra_lib_paths)

else

.NOTPARALLEL:

endif
