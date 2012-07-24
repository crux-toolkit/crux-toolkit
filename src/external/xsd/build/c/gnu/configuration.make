# file      : build/c/gnu/configuration.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(bld_root)/c/gnu/configuration-rules.make,$(dcf_root))

# Static configuration.
#
$(call include,$(bld_root)/c/gnu/configuration-static.make)

ifneq ($(bld_root),$(scf_root))
$(call -include,$(scf_root)/c/gnu/configuration-static.make)
endif

# Dynamic configuration.
#
c_gnu :=
c_gnu_optimization_options :=

$(call -include,$(dcf_root)/c/gnu/configuration-dynamic.make)

ifdef c_gnu

c_gnu_debugging_options := $(if $(findstring y,$(c_debug)),-g)

$(out_root)/%: c_gnu_debugging_options := $(c_gnu_debugging_options)
$(out_root)/%: c_gnu_optimization_options := $(c_gnu_optimization_options)

else

.NOTPARALLEL:

endif
