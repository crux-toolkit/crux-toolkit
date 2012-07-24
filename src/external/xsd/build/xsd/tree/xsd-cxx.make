# file      : build/xsd/tree/xsd-cxx.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
# license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#@@ Need to use extensions from cxx config.
#

# C++/Tree mapping.
#
xsd_tree_pattern :=                     \
$(out_base)/%$(xsd_tree_suffix).cxx     \
$(out_base)/%$(xsd_tree_suffix).hxx     \
$(out_base)/%$(xsd_tree_suffix).ixx     \
$(out_base)/%$(xsd_tree_suffix)-fwd.hxx

$(xsd_tree_pattern): xsd := xsd
$(xsd_tree_pattern): xsd_command := cxx-tree
$(xsd_tree_pattern): xsd_options :=

.PRECIOUS: $(xsd_tree_pattern)

ifeq ($(out_base),$(src_base))

$(xsd_tree_pattern): $(src_base)/%.xsd
	$(call message,xsd $<,$(xsd) $(xsd_command) $(xsd_options) --output-dir $(dir $@) $<)

else

$(xsd_tree_pattern): $(src_base)/%.xsd | $$(dir $$@).
	$(call message,xsd $<,$(xsd) $(xsd_command) $(xsd_options) --output-dir $(dir $@) $<)

$(xsd_tree_pattern): $(out_base)/%.xsd | $$(dir $$@).
	$(call message,xsd $<,$(xsd) $(xsd_command) $(xsd_options) --output-dir $(dir $@) $<)

endif

.PHONY: $(out_base)/%$(xsd_tree_suffix).cxx.xsd.clean

$(out_base)/%$(xsd_tree_suffix).cxx.xsd.clean:
	$(call message,rm $(@:.cxx.xsd.clean=.cxx),rm -f $(@:.cxx.xsd.clean=.cxx))
	$(call message,rm $(@:.cxx.xsd.clean=.hxx),rm -f $(@:.cxx.xsd.clean=.hxx))
	$(call message,rm $(@:.cxx.xsd.clean=.ixx),rm -f $(@:.cxx.xsd.clean=.ixx))
	$(call message,rm $(@:.cxx.xsd.clean=-fwd.hxx),rm -f $(@:.cxx.xsd.clean=-fwd.hxx))

# Reset the config variables so they won't take effect in other places.
#
xsd_tree_suffix :=
