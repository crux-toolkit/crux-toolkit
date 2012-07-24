# file      : build/import/libz/configuration-rules.make
# author    : Boris Kolpackov <boris@kolpackov.net>
# copyright : Copyright (c) 2009 Boris Kolpackov
# license   : GNU GPL v2; see accompanying LICENSE file

$(dcf_root)/import/libz/configuration-dynamic.make: | $(dcf_root)/import/libz/.
	$(call message,,$(scf_root)/import/libz/configure $@)

ifndef %foreign%

disfigure::
	$(call message,rm $(dcf_root)/import/libz/configuration-dynamic.make,\
rm -f $(dcf_root)/import/libz/configuration-dynamic.make)

endif
