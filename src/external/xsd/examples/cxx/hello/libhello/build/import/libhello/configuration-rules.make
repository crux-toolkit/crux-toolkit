# file      : examples/cxx/hello/libhello/build/import/libhello/configuration-rules.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(dcf_root)/import/libhello/configuration-dynamic.make: | $(dcf_root)/import/libhello/.
	$(call message,,$(scf_root)/import/libhello/configure $@)

ifndef %foreign%

$(dcf_root)/.disfigure::
	$(call message,rm $(dcf_root)/import/libhello/configuration-dynamic.make,\
rm -f $(dcf_root)/import/libhello/configuration-dynamic.make)

endif
