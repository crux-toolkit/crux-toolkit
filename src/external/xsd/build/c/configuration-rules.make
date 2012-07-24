# file      : build/c/configuration-rules.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(dcf_root)/c/configuration-dynamic.make: | $(dcf_root)/c/.
	$(call message,,$(bld_root)/c/configure $@ \
"$(origin c_pp_options)"  \
"$(origin c_options)"     \
"$(origin c_ld_options)"  \
"$(origin c_libs)")

ifndef %foreign%

$(dcf_root)/.disfigure::
	$(call message,rm $(dcf_root)/c/configuration-dynamic.make,\
rm -f $(dcf_root)/c/configuration-dynamic.make)

endif
