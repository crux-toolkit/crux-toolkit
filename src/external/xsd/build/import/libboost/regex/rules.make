# file      : build/import/libboost/regex/rules.make
# author    : Boris Kolpackov <boris@kolpackov.net>
# copyright : Copyright (c) 2005-2010 Boris Kolpackov
# license   : GNU GPL v2; see accompanying LICENSE file

$(dcf_root)/import/libboost/%: root := $(libboost_root)

$(dcf_root)/import/libboost/regex/regex.l: \
  | $(dcf_root)/import/libboost/regex/.

ifeq ($(libboost_type),archive)

ifeq ($(libboost_system),y)
$(dcf_root)/import/libboost/regex/regex.l: \
  $(libboost_root)/stage/lib/libboost_regex$(libboost_suffix).a \
  $(libboost_root)/stage/lib/libboost_system$(libboost_suffix).a
else
$(dcf_root)/import/libboost/regex/regex.l: \
  $(libboost_root)/stage/lib/libboost_regex$(libboost_suffix).a
endif
	@echo $^ >$@

else

ifeq ($(libboost_system),y)
$(dcf_root)/import/libboost/regex/regex.l: \
  $(libboost_root)/stage/lib/libboost_regex$(libboost_suffix).so \
  $(libboost_root)/stage/lib/libboost_system$(libboost_suffix).so
else
$(dcf_root)/import/libboost/regex/regex.l: \
  $(libboost_root)/stage/lib/libboost_regex$(libboost_suffix).so
endif
	@echo $^ >$@
	@echo rpath:$(root)/stage/lib >>$@

endif


$(dcf_root)/import/libboost/regex/regex.l.cpp-options: \
  | $(dcf_root)/import/libboost/regex/.
	@echo include: -I$(root) >$@

ifndef %foreign%

disfigure::
	$(call message,rm $(dcf_root)/import/libboost/regex/regex.l,\
rm -f $(dcf_root)/import/libboost/regex/regex.l)
	$(call message,,rm -f $(dcf_root)/import/libboost/regex/regex.l.cpp-options)

endif
