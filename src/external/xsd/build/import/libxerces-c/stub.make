# file      : build/import/libxerces-c/stub.make
# author    : Boris Kolpackov <boris@kolpackov.net>
# copyright : Copyright (c) 2005-2009 Boris Kolpackov
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(scf_root)/import/libxerces-c/configuration-rules.make,$(dcf_root))

libxerces_c_installed :=

$(call -include,$(dcf_root)/import/libxerces-c/configuration-dynamic.make)

ifdef libxerces_c_installed

ifeq ($(libxerces_c_installed),y)

$(call export,l: -lxerces-c,cpp-options: )

else

$(call include-once,$(scf_root)/import/libxerces-c/rules.make,$(dcf_root))

$(call export,\
  l: $(dcf_root)/import/libxerces-c/xerces-c.l,\
  cpp-options: $(dcf_root)/import/libxerces-c/xerces-c.l.cpp-options)

endif

else

.NOTPARALLEL:

endif
