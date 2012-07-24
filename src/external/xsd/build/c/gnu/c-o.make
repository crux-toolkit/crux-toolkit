# file      : build/c/gnu/c-o.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(bld_root)/c/cpp-options.make,$(out_base))
$(call include,$(bld_root)/c/gnu/configuration.make)

$(out_base)/%.o: c := $(c_gnu)
$(out_base)/%.o: c_options := $(c_gnu_optimization_options) $(c_gnu_debugging_options)

#@@ wrong prefix
#
$(out_base)/%.o: expand-cpp-options-impl = \
$(if $1,$(shell sed -e 's%include: \(.*\)%\1%' -e t -e d $1))

$(out_base)/%.o: expand-cpp-options = \
$(call expand-cpp-options-impl,$(filter %.cpp-options,$1))

ifeq ($(out_base),$(src_base))
$(out_base)/%.o: $(src_base)/%.$(c_s_suffix)
else
$(out_base)/%.o: $(src_base)/%.$(c_s_suffix) | $$(dir $$@).
endif
	$(call message,c $<,$(c) \
$(cpp_options) $(call expand-cpp-options,$^) $(c_pp_extra_options) $(c_options) \
$(c_pic_options) $(subst y,-fexceptions,$(filter y,$(c_exceptions))) $(c_extra_options) -o $@ -c $<)


ifneq ($(out_base),$(src_base))

$(out_base)/%.o: $(out_base)/%.$(c_s_suffix) | $$(dir $$@).
	$(call message,c $<,$(c) \
$(cpp_options) $(call expand-cpp-options,$^) $(c_pp_extra_options) $(c_options) \
$(c_pic_options) $(subst y,-fexceptions,$(filter y,$(c_exceptions))) $(c_extra_options) -o $@ -c $<)

endif

.PHONY: $(out_base)/%.o.$(c_s_suffix).clean

$(out_base)/%.o.$(c_s_suffix).clean:
	$(call message,rm $$1,rm -f $$1,$(basename $(basename $@)))
