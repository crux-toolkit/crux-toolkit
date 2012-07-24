# file      : build/c/intel/configuration.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(bld_root)/c/intel/configuration-rules.make,$(dcf_root))

# Static configuration.
#
$(call include,$(bld_root)/c/intel/configuration-static.make)

ifneq ($(bld_root),$(scf_root))
$(call -include,$(scf_root)/c/intel/configuration-static.make)
endif

# Dynamic configuration.
#
c_intel :=
c_intel_optimization_options :=

$(call -include,$(dcf_root)/c/intel/configuration-dynamic.make)

ifdef c_intel

c_intel_debugging_options := $(if $(findstring y,$(c_debug)),-g)

$(out_root)/%: c_intel_debugging_options := $(c_intel_debugging_options)
$(out_root)/%: c_intel_optimization_options := $(c_intel_optimization_options)

else

.NOTPARALLEL:

endif
