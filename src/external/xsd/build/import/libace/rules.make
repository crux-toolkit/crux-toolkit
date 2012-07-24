# file      : build/import/libace/rules.make
# author    : Boris Kolpackov <boris@kolpackov.net>
# copyright : Copyright (c) 2005-2009 Boris Kolpackov
# license   : GNU GPL v2; see accompanying LICENSE file

$(dcf_root)/import/libace/%: root := $(libace_root)

ifeq ($(libace_type),archive)

$(dcf_root)/import/libace/ace.l: $(libace_root)/lib/libACE.a
	@echo $< >$@
else

$(dcf_root)/import/libace/ace.l: $(libace_root)/lib/libACE.so
	@echo $< >$@
	@echo rpath:$(root)/lib >>$@
endif

$(dcf_root)/import/libace/ace.l.cpp-options:
	@echo include: -I$(root) >$@

ifndef %foreign%

disfigure::
	$(call message,rm $(dcf_root)/import/libace/ace.l,\
rm -f $(dcf_root)/import/libace/ace.l)
	$(call message,,rm -f $(dcf_root)/import/libace/ace.l.cpp-options)

endif
