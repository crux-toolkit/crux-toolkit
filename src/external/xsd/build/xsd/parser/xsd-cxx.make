# file      : build/xsd/parser/xsd-cxx.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
# license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#@@ Need to use extensions from cxx config.
#

# C++/Parser mapping.
#
ifeq ($(xsd_parser_skel_suffix),)
xsd_parser_skel_suffix := -pskel
endif

xsd_parser_pattern :=                      \
$(out_base)/%$(xsd_parser_skel_suffix).cxx \
$(out_base)/%$(xsd_parser_skel_suffix).hxx \
$(out_base)/%$(xsd_parser_skel_suffix).ixx

ifneq ($(xsd_parser_impl_suffix),)
xsd_parser_pattern +=                      \
$(out_base)/%$(xsd_parser_impl_suffix).cxx \
$(out_base)/%$(xsd_parser_impl_suffix).hxx \
$(out_base)/%-driver.cxx
endif


$(xsd_parser_pattern): xsd := xsd
$(xsd_parser_pattern): xsd_command := cxx-parser

ops := --skel-file-suffix $(xsd_parser_skel_suffix)

ifneq ($(xsd_pimpl_suffix),)
ops += --impl-file-suffix $(xsd_parser_impl_suffix)
endif

$(xsd_parser_pattern): xsd_options := $(ops)


.PRECIOUS: $(xsd_parser_pattern)

ifeq ($(out_base),$(src_base))

$(xsd_parser_pattern): $(src_base)/%.xsd
	$(call message,xsd $<,$(xsd) $(xsd_command) $(xsd_options) --output-dir $(dir $@) $<)

else

$(xsd_parser_pattern): $(src_base)/%.xsd | $$(dir $$@).
	$(call message,xsd $<,$(xsd) $(xsd_command) $(xsd_options) --output-dir $(dir $@) $<)

$(xsd_parser_pattern): $(out_base)/%.xsd | $$(dir $$@).
	$(call message,xsd $<,$(xsd) $(xsd_command) $(xsd_options) --output-dir $(dir $@) $<)

endif


.PHONY: $(out_base)/%$(xsd_parser_skel_suffix).cxx.xsd.clean

$(out_base)/%$(xsd_parser_skel_suffix).cxx.xsd.clean:
	$(call message,rm $$1,rm -f $$1,$(@:.cxx.xsd.clean=.cxx))
	$(call message,rm $$1,rm -f $$1,$(@:.cxx.xsd.clean=.hxx))
	$(call message,rm $$1,rm -f $$1,$(@:.cxx.xsd.clean=.ixx))

ifneq ($(xsd_parser_impl_suffix),)
.PHONY: $(out_base)/%$(xsd_parser_impl_suffix).cxx.xsd.clean

$(out_base)/%$(xsd_parser_impl_suffix).cxx.xsd.clean:
	$(call message,rm $$1,rm -f $$1,$(@:.cxx.xsd.clean=.cxx))
	$(call message,rm $$1,rm -f $$1,$(@:.cxx.xsd.clean=.hxx))
	$(call message,rm $$1,rm -f $$1,$(out_base)/$*-driver.cxx)
endif

# Reset the config variables so they won't take effect in other places.
#
xsd_parser_skel_suffix :=
xsd_parser_impl_suffix :=
