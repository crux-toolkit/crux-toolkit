# file      : build/system/configuration-rules.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file


$(dcf_root)/system/configuration-dynamic.make: | $(dcf_root)/system/.
	$(call message,,$(bld_root)/system/configure $@)

ifeq ($(.DEFAULT_GOAL),$(dcf_root)/system/configuration-dynamic.make)
.DEFAULT_GOAL :=
endif

ifndef %foreign%

$(dcf_root)/.disfigure::
	$(call message,rm $(dcf_root)/system/configuration-dynamic.make,\
rm -f $(dcf_root)/system/configuration-dynamic.make)

ifeq ($(.DEFAULT_GOAL),$(dcf_root)/.disfigure)
.DEFAULT_GOAL :=
endif

endif
