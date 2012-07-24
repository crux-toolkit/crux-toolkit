# file      : build/import/libdbxml/rules.make
# author    : Boris Kolpackov <boris@kolpackov.net>
# copyright : Copyright (c) 2005-2009 Boris Kolpackov
# license   : GNU GPL v2; see accompanying LICENSE file

#@@ .l construction is compiler-specific
#

$(dcf_root)/import/libdbxml/%: root := $(libdbxml_root)

ifeq ($(libdbxml_type),archive)

$(dcf_root)/import/libdbxml/dbxml.l:       \
$(libdbxml_root)/install/lib/libdbxml.a    \
$(libdbxml_root)/install/lib/libxqilla.a   \
$(libdbxml_root)/install/lib/libxerces-c.a \
$(libdbxml_root)/install/lib/libdb_cxx.a   \
$(libdbxml_root)/install/lib/libdb.a
	
else

$(dcf_root)/import/libdbxml/dbxml.l:        \
$(libdbxml_root)/install/lib/libdbxml.so    \
$(libdbxml_root)/install/lib/libxqilla.so   \
$(libdbxml_root)/install/lib/libxerces-c.so \
$(libdbxml_root)/install/lib/libdb_cxx.so   \
$(libdbxml_root)/install/lib/libdb.so
	@echo $^ | xargs -n 1 echo >$@	
	@echo rpath:$(root)/install/lib >>$@
endif

$(dcf_root)/import/libdbxml/dbxml.l.cpp-options:
	@echo include: -I$(root)/install/include >$@

ifndef %foreign%

disfigure::
	$(call message,rm $(dcf_root)/import/libdbxml/dbxml.l,\
rm -f $(dcf_root)/import/libdbxml/dbxml.l)
	$(call message,,rm -f $(dcf_root)/import/libdbxml/dbxml.l.cpp-options)

endif
