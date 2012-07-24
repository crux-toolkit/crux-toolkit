# file      : build/import/libxerces-c/rules.make
# author    : Boris Kolpackov <boris@kolpackov.net>
# copyright : Copyright (c) 2005-2009 Boris Kolpackov
# license   : GNU GPL v2; see accompanying LICENSE file

$(dcf_root)/import/libxerces-c/%: root := $(libxerces_c_root)

ifneq ($(filter 3.%,$(libxerces_c_version)),)

# 3.x.y
#
ifeq ($(libxerces_c_type),archive)

$(dcf_root)/import/libxerces-c/xerces-c.l: $(libxerces_c_root)/src/.libs/libxerces-c.a
	@echo $< >$@
else

$(dcf_root)/import/libxerces-c/xerces-c.l: $(libxerces_c_root)/src/.libs/libxerces-c.so
	@echo $< >$@
	@echo rpath:$(root)/src/.libs >>$@
endif

$(dcf_root)/import/libxerces-c/xerces-c.l.cpp-options:
	@echo include: -I$(root)/src >$@
else

#  2.x.y
#
ifeq ($(libxerces_c_type),archive)

$(dcf_root)/import/libxerces-c/xerces-c.l: $(libxerces_c_root)/lib/libxerces-c.a
	@echo $< >$@
else

$(dcf_root)/import/libxerces-c/xerces-c.l: $(libxerces_c_root)/lib/libxerces-c.so
	@echo $< >$@
	@echo rpath:$(root)/lib >>$@
endif

$(dcf_root)/import/libxerces-c/xerces-c.l.cpp-options:
	@echo include: -I$(root)/include >$@
endif


ifndef %foreign%

disfigure::
	$(call message,rm $(dcf_root)/import/libxerces-c/xerces-c.l,\
rm -f $(dcf_root)/import/libxerces-c/xerces-c.l)
	$(call message,,rm -f $(dcf_root)/import/libxerces-c/xerces-c.l.cpp-options)

endif
