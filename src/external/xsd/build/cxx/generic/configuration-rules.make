# file      : build/cxx/generic/configuration-rules.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(dcf_root)/cxx/generic/configuration-dynamic.make: | $(dcf_root)/cxx/generic/.
	$(call message,,$(bld_root)/cxx/generic/configure $@)

ifndef %foreign%

$(dcf_root)/.disfigure::
	$(call message,rm $(dcf_root)/cxx/generic/configuration-dynamic.make,\
rm -f $(dcf_root)/cxx/generic/configuration-dynamic.make)

endif
