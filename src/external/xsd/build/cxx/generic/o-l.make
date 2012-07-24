# file      : build/cxx/generic/o-l.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include,$(bld_root)/cxx/generic/configuration.make)

ifeq ($(ld_lib_type),shared)
$(call include,$(bld_root)/cxx/generic/configuration-sl.make)
endif

ifneq ($(cxx_extra_lib_paths),)
vpath %.so $(cxx_extra_lib_paths)
vpath %.a $(cxx_extra_lib_paths)
endif

ifdef ld_lib_type

.PHONY: $(out_base)/%.l.o.clean

ifeq ($(ld_lib_type),archive)

$(out_base)/%.l: ar := $(ld_lib_ar)
$(out_base)/%.l: ar_options ?= -rc

$(out_base)/%.l: ranlib := $(ld_lib_ranlib)
$(out_base)/%.l: ranlib_options ?=

ifeq ($(out_base),$(src_base))
$(out_base)/%.l:
else
$(out_base)/%.l: | $$(dir $$@).
endif
	$(call message,ar  $@,$(ar) $(ar_options) $(@D)/lib$(basename $(@F)).a $(filter %.o,$^))
	$(call message,,$(ranlib) $(ranlib_options) $(@D)/lib$(basename $(@F)).a)
	$(call message,,echo "$(@D)/lib$(basename $(@F)).a" >$@)
	$(call message,,echo "$(patsubst %.l,`cat %.l`,$(filter %.a %.so %.l,$^))" | xargs -n 1 echo >>$@)

$(out_base)/%.l.o.clean:
	$(call message,rm $$1,rm -f $$1 $(@D)/$(patsubst %.l.o.clean,lib%.a,$(@F)),$(basename $(basename $@)))

else

$(out_base)/%.l: ld := $(cxx_generic)
$(out_base)/%.l: c_pic_options := $(cxx_generic_pic_option)
$(out_base)/%.l: cxx_pic_options := $(cxx_generic_pic_option)
$(out_base)/%.l: comma_ := ,

$(out_base)/%.l: expand-l = $(shell sed -e 's%^rpath:\(.*\)%%' $1)

ifeq ($(out_base),$(src_base))
$(out_base)/%.l:
else
$(out_base)/%.l: | $$(dir $$@).
endif
	$(call message,ld  $@,$(ld) $(cxx_generic_shared_option) \
$(cxx_extra_options) $(ld_options) $(cxx_ld_extra_options) -o $(@D)/lib$(basename $(@F)).so \
$(foreach f,$^,$(if $(patsubst %.l,,$f),$f,$(call expand-l,$f))) $(cxx_extra_libs))
	$(call message,,echo "$(@D)/lib$(basename $(@F)).so" >$@)
	$(call message,,echo "rpath:$(@D)" >>$@)
	$(call message,,echo "$(patsubst %.l,`cat %.l`,$(filter %.a %.so %.l,$^))" | xargs -n 1 echo >>$@)

$(out_base)/%.l.o.clean:
	$(call message,rm $$1,rm -f $$1 $(@D)/$(patsubst %.l.o.clean,lib%.so,$(@F)),$(basename $(basename $@)))

endif
endif
