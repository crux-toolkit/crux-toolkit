# file      : build/import/libdbxml/stub.make
# author    : Boris Kolpackov <boris@kolpackov.net>
# copyright : Copyright (c) 2005-2009 Boris Kolpackov
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(scf_root)/import/libdbxml/configuration-rules.make,$(dcf_root))

libdbxml_installed :=

$(call -include,$(dcf_root)/import/libdbxml/configuration-dynamic.make)

ifdef libdbxml_installed

ifeq ($(libdbxml_installed),y)

$(call export,l: -ldbxml -lxqilla -lxerces-c -ldb_cxx -ldb,cpp-options: )

else

$(call include-once,$(scf_root)/import/libdbxml/rules.make,$(dcf_root))

$(call export,\
  l: $(dcf_root)/import/libdbxml/dbxml.l,\
  cpp-options: $(dcf_root)/import/libdbxml/dbxml.l.cpp-options)

endif

else

.NOTPARALLEL:

endif
