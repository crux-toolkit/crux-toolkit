# file      : build/c/gnu/o-l.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include,$(bld_root)/c/gnu/configuration.make)

ifneq ($(c_extra_lib_paths),)
vpath %.so $(c_extra_lib_paths)
vpath %.a $(c_extra_lib_paths)
endif

ifneq ($(c_gnu_libraries),)
vpath %.so $(c_gnu_libraries)
vpath %.a $(c_gnu_libraries)
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

$(out_base)/%.l: ld := $(c_gnu)
$(out_base)/%.l: ld_options := $(c_gnu_optimization_options) $(c_gnu_debugging_options)
$(out_base)/%.l: c_pic_options := -fPIC
$(out_base)/%.l: comma_ := ,

$(out_base)/%.l: expand-l = $(if $(subst n,,$(c_rpath)),\
$(shell sed -e 's%^rpath:\(.*\)%-Wl,-rpath,\1%' $1),\
$(shell sed -e 's%^rpath:\(.*\)%%' $1))

ifeq ($(out_base),$(src_base))
$(out_base)/%.l:
else
$(out_base)/%.l: | $$(dir $$@).
endif
	$(call message,ld  $@,$(ld) -shared \
$(c_extra_options) $(ld_options) $(c_ld_extra_options) \
-o $(@D)/lib$(basename $(@F)).so -Wl$(comma_)-soname=lib$(basename $(@F)).so \
$(foreach f,$^,$(if $(patsubst %.l,,$f),$f,$(call expand-l,$f))) $(c_extra_libs))
	$(call message,,echo "$(@D)/lib$(basename $(@F)).so" >$@)
	$(call message,,echo "rpath:$(@D)" >>$@)
	$(call message,,echo "$(patsubst %.l,`cat %.l`,$(filter %.a %.so %.l,$^))" | xargs -n 1 echo >>$@)

$(out_base)/%.l.o.clean:
	$(call message,rm $$1,rm -f $$1 $(@D)/$(patsubst %.l.o.clean,lib%.so,$(@F)),$(basename $(basename $@)))

endif
endif
