# file      : build/c/generic/o-e.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include,$(bld_root)/c/generic/configuration.make)

ifneq ($(c_extra_lib_paths),)
vpath %.so $(c_extra_lib_paths)
vpath %.a $(c_extra_lib_paths)
endif

$(out_base)/%: ld := $(c_generic)
$(out_base)/%: expand-l = $(shell sed -e 's%^rpath:\(.*\)%%' $1)

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
