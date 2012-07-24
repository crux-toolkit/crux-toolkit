# file      : build/configuration-rules.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
# license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

$(dcf_root)/configuration-dynamic.make: | $(dcf_root)/.
	$(call message,,$(scf_root)/configure $@)

ifndef %foreign%

disfigure::
	$(call message,rm $$1,rm -f $$1,$(dcf_root)/configuration-dynamic.make)

endif

ifeq ($(.DEFAULT_GOAL),$(dcf_root)/configuration-dynamic.make)
.DEFAULT_GOAL :=
endif
