# file      : build/cxx/intel/configuration.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(bld_root)/cxx/intel/configuration-rules.make,$(dcf_root))

# Static configuration.
#
$(call include,$(bld_root)/cxx/intel/configuration-static.make)

ifneq ($(bld_root),$(scf_root))
$(call -include,$(scf_root)/cxx/intel/configuration-static.make)
endif

# Dynamic configuration.
#
cxx_intel :=
cxx_intel_optimization_options :=

$(call -include,$(dcf_root)/cxx/intel/configuration-dynamic.make)

ifdef cxx_intel

cxx_intel_debugging_options := $(if $(findstring y,$(cxx_debug)),-g)

$(out_root)/%: cxx_intel_debugging_options := $(cxx_intel_debugging_options)
$(out_root)/%: cxx_intel_optimization_options := $(cxx_intel_optimization_options)

else

.NOTPARALLEL:

endif
