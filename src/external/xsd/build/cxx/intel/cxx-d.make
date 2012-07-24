# file      : build/cxx/intel/cxx-o.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(bld_root)/c/cpp-options.make,$(out_base))
$(call include,$(bld_root)/cxx/intel/configuration.make)

# Make will try to build dependecies (since they are ultimately included
# files) during configuartion phase without cxx_intel being discovered yet.
# This is also why dependecies should be included with -include.
#
ifdef cxx_intel

$(out_base)/%.o.d: cxx := $(cxx_intel)
$(out_base)/%.o.d: cpp_options ?=

#@@ This needs to be shared with cxx-o.
#
#@@ wrong prefix
#
$(out_base)/%.o.d: expand-cpp-options-impl = \
$(if $1,$(shell sed -e 's%include: \(.*\)%\1%' -e t -e d $1))

$(out_base)/%.o.d: expand-cpp-options = \
$(call expand-cpp-options-impl,$(filter %.cpp-options,$1))

.PRECIOUS: $(out_base)/%.o.d

ifeq ($(out_base),$(src_base))
$(out_base)/%.o.d: $(src_base)/%.$(cxx_s_suffix)
else
$(out_base)/%.o.d: $(src_base)/%.$(cxx_s_suffix) | $$(dir $$@).
endif
	$(call message,,if test -f $@; then mv $@ $@.old; fi && \
$(cxx) $(cpp_options) $(call expand-cpp-options,$^) $(cxx_pp_extra_options) \
$(cxx_extra_options) -M -MG -MP -MQ $@ -MQ $(basename $@) $< | \
$(bld_root)/c/intel/dep $(out_base) $(filter %.cpp-options,$^) >$@ && \
if test -f $@.old; then if cmp -s $@ $@.old; then rm -f $@.old && false; \
else rm -f $@.old && true; fi fi)

ifneq ($(out_base),$(src_base))

$(out_base)/%.o.d: $(out_base)/%.$(cxx_s_suffix) | $$(dir $$@).
	$(call message,,if test -f $@; then mv $@ $@.old; fi && \
$(cxx) $(cpp_options) $(call expand-cpp-options,$^) $(cxx_pp_extra_options) \
$(cxx_extra_options) -M -MG -MP -MQ $@ -MQ $(basename $@) $< | \
$(bld_root)/c/intel/dep $(out_base) $(filter %.cpp-options,$^) >$@ && \
if test -f $@.old; then if cmp -s $@ $@.old; then rm -f $@.old && false; \
else rm -f $@.old && true; fi fi)

endif

.PHONY: $(out_base)/%.o.d.$(cxx_s_suffix).clean

$(out_base)/%.o.d.$(cxx_s_suffix).clean:
	$(call message,rm $$1,rm -f $$1,$(basename $(basename $@)))

endif
