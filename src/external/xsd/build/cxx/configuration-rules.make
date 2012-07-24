# file      : build/cxx/configuration-rules.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(dcf_root)/cxx/configuration-dynamic.make: | $(dcf_root)/cxx/.
	$(call message,,$(bld_root)/cxx/configure $@ \
"$(origin cxx_pp_options)"  \
"$(origin cxx_options)"     \
"$(origin cxx_ld_options)"  \
"$(origin cxx_libs)")

ifndef %foreign%

$(dcf_root)/.disfigure::
	$(call message,rm $(dcf_root)/cxx/configuration-dynamic.make,\
rm -f $(dcf_root)/cxx/configuration-dynamic.make)

endif
