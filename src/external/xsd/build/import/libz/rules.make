# file      : build/import/libz/rules.make
# author    : Boris Kolpackov <boris@kolpackov.net>
# copyright : Copyright (c) 2009 Boris Kolpackov
# license   : GNU GPL v2; see accompanying LICENSE file

$(dcf_root)/import/libz/%: root := $(libz_root)

ifeq ($(libz_type),archive)

$(dcf_root)/import/libz/z.l: $(libz_root)/libz.a
	@echo $< >$@
else

$(dcf_root)/import/libz/z.l: $(libz_root)/libz.so
	@echo $< >$@
	@echo rpath:$(root) >>$@
endif

$(dcf_root)/import/libz/z.l.cpp-options:
	@echo include: -I$(root) >$@

ifndef %foreign%

disfigure::
	$(call message,rm $(dcf_root)/import/libz/z.l,\
rm -f $(dcf_root)/import/libz/z.l)
	$(call message,,rm -f $(dcf_root)/import/libz/z.l.cpp-options)

endif
