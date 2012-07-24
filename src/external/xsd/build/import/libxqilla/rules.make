# file      : build/import/libxqilla/rules.make
# author    : Boris Kolpackov <boris@kolpackov.net>
# copyright : Copyright (c) 2005-2009 Boris Kolpackov
# license   : GNU GPL v2; see accompanying LICENSE file

$(dcf_root)/import/libxqilla/%: root := $(libxqilla_root)

ifeq ($(libxqilla_type),archive)

$(dcf_root)/import/libxqilla/xqilla.l: $(libxqilla_root)/.libs/libxqilla.a
	@echo $< >$@
else

$(dcf_root)/import/libxqilla/xqilla.l: $(libxqilla_root)/.libs/libxqilla.so
	@echo $< >$@
	@echo rpath:$(root)/.libs >>$@
endif

$(dcf_root)/import/libxqilla/xqilla.l.cpp-options:
	@echo include: -I$(root)/include >$@


ifndef %foreign%

disfigure::
	$(call message,rm $(dcf_root)/import/libxqilla/xqilla.l,\
rm -f $(dcf_root)/import/libxqilla/xqilla.l)
	$(call message,,rm -f $(dcf_root)/import/libxqilla/xqilla.l.cpp-options)

endif
