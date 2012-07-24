# file      : build/cxx/configuration.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(bld_root)/cxx/configuration-rules.make,$(dcf_root))

# Static configuration.
#
$(call include,$(bld_root)/cxx/configuration-static.make)


# Dynamic configuration.
#
cxx_id       :=
cxx_optimize :=
cxx_debug    :=
cxx_rpath    :=

cxx_pp_extra_options :=
cxx_extra_options    :=
cxx_ld_extra_options :=
cxx_extra_libs       :=
cxx_extra_lib_paths  :=

$(call -include,$(dcf_root)/cxx/configuration-dynamic.make)

ifdef cxx_id

$(out_root)/%: cxx_id       := $(cxx_id)
$(out_root)/%: cxx_optimize := $(cxx_optimize)
$(out_root)/%: cxx_debug    := $(cxx_debug)
$(out_root)/%: cxx_rpath    := $(cxx_rpath)

$(out_root)/%: cxx_pp_extra_options := $(cxx_pp_extra_options)
$(out_root)/%: cxx_extra_options    := $(cxx_extra_options)
$(out_root)/%: cxx_ld_extra_options := $(cxx_ld_extra_options)
$(out_root)/%: cxx_extra_libs       := $(cxx_extra_libs)
$(out_root)/%: cxx_extra_lib_paths  := $(cxx_extra_lib_paths)

else

.NOTPARALLEL:

endif
