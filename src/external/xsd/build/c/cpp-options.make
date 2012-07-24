# file      : build/c/cpp-options.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(out_base)/%.cpp-options: c-cpp-options-expand-impl = $(if $1,cat $1 >>$2)
$(out_base)/%.cpp-options: c-cpp-options-expand = \
$(call c-cpp-options-expand-impl,$(filter %.cpp-options,$1),$2)


#@@ Asymmetry: I hide creation commands but show removals.
#
$(out_base)/%.cpp-options: c-cpp-options-force-impl = \
$(if $(subst $1,,$2)$(subst $2,,$1),FORCE)

$(out_base)/%.cpp-options: c-cpp-options-force = \
$(call c-cpp-options-force-impl,$(shell if test -f $1; then exec head -n 1 $1; fi),include: $2)


# @@ Need to incroporate prefix: into force-check. Then it will probably
#    have to be non-optional.
#

ifeq ($(src_base),$(out_base))
$(out_base)/%.cpp-options: $$(call c-cpp-options-force,$$@,$$(value))
else
$(out_base)/%.cpp-options: $$(call c-cpp-options-force,$$@,$$(value)) | $$(dir $$@).
endif
	@echo "include: $(value)" >$@
	$(if $(strip $(prefix)),@echo "prefix: $(prefix)" >>$@)
	@$(call c-cpp-options-expand,$^,$@)

# Clean.
#
$(out_base)/%.cpp-options.clean:
	$(call message,rm $(basename $@),rm -f $(basename $@))
