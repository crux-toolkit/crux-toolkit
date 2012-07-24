# file      : build/c/intel/o-e.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include,$(bld_root)/c/intel/configuration.make)

#@@ should it be lib%.so?
#
ifneq ($(c_extra_lib_paths),)
vpath %.so $(c_extra_lib_paths)
vpath %.a $(c_extra_lib_paths)
endif

ifneq ($(c_intel_libraries),)
vpath %.so $(c_intel_libraries)
vpath %.a $(c_intel_libraries)
endif

$(out_base)/%: ld := $(c_intel)
$(out_base)/%: ld_options := $(c_intel_optimization_options) $(c_intel_debugging_options) $(c_ld_extra_options)

$(out_base)/%: expand-l = $(if $(subst n,,$(c_rpath)),\
$(shell sed -e 's%^rpath:\(.*\)%-Wl,-rpath,\1%' $1),\
$(shell sed -e 's%^rpath:\(.*\)%%' $1))

ifeq ($(out_base),$(src_base))
$(out_base)/%: $(out_base)/%.o
else
$(out_base)/%: $(out_base)/%.o | $$(dir $$@).
endif
	$(call message,ld $@,$(ld) \
$(c_extra_options) $(ld_options) $(c_ld_extra_options) -o $@ \
$(foreach f,$^,$(if $(patsubst %.l,,$f),$f,$(call expand-l,$f))) $(c_extra_libs))

.PHONY: $(out_base)/%.o.clean

$(out_base)/%.o.clean:
	$(call message,rm $(basename $(basename $@)),rm -f $(basename $@) $(basename $(basename $@)))
