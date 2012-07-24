# file      : build/import/libxqilla/stub.make
# author    : Boris Kolpackov <boris@kolpackov.net>
# copyright : Copyright (c) 2005-2009 Boris Kolpackov
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(scf_root)/import/libxqilla/configuration-rules.make,$(dcf_root))

libxqilla_installed :=

$(call -include,$(dcf_root)/import/libxqilla/configuration-dynamic.make)

ifdef libxqilla_installed

ifeq ($(libxqilla_installed),y)

$(call export,l: -lxqilla,cpp-options: )

else

$(call include-once,$(scf_root)/import/libxqilla/rules.make,$(dcf_root))

$(call export,\
  l: $(dcf_root)/import/libxqilla/xqilla.l,\
  cpp-options: $(dcf_root)/import/libxqilla/xqilla.l.cpp-options)

endif

else

.NOTPARALLEL:

endif
