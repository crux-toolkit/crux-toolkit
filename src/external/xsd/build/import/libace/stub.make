# file      : build/import/libace/stub.make
# author    : Boris Kolpackov <boris@kolpackov.net>
# copyright : Copyright (c) 2005-2009 Boris Kolpackov
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(scf_root)/import/libace/configuration-rules.make,$(dcf_root))

libace_installed :=

$(call -include,$(dcf_root)/import/libace/configuration-dynamic.make)

ifdef libace_installed

ifeq ($(libace_installed),y)

$(call export,l: -lACE,cpp-options: )

else

$(call include-once,$(scf_root)/import/libace/rules.make,$(dcf_root))

$(call export,\
  l: $(dcf_root)/import/libace/ace.l,\
  cpp-options: $(dcf_root)/import/libace/ace.l.cpp-options)

endif

else

.NOTPARALLEL:

endif
