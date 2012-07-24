# file      : build/cxx/generic/configuration-sl-rules.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(dcf_root)/cxx/generic/configuration-sl-dynamic.make: | $(dcf_root)/cxx/generic/.
	$(call message,,$(bld_root)/cxx/generic/configure-sl $@)

ifndef %foreign%

$(dcf_root)/.disfigure::
	$(call message,rm $(dcf_root)/cxx/generic/configuration-sl-dynamic.make,\
rm -f $(dcf_root)/cxx/generic/configuration-sl-dynamic.make)

endif
