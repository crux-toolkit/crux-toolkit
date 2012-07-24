# file      : build/import/libz/stub.make
# author    : Boris Kolpackov <boris@kolpackov.net>
# copyright : Copyright (c) 2009 Boris Kolpackov
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(scf_root)/import/libz/configuration-rules.make,$(dcf_root))

libz_installed :=

$(call -include,$(dcf_root)/import/libz/configuration-dynamic.make)

ifdef libz_installed

ifeq ($(libz_installed),y)

$(call export,l: -lz,cpp-options: )

else

$(call include-once,$(scf_root)/import/libz/rules.make,$(dcf_root))

$(call export,\
  l: $(dcf_root)/import/libz/z.l,\
  cpp-options: $(dcf_root)/import/libz/z.l.cpp-options)

endif

else

.NOTPARALLEL:

endif
