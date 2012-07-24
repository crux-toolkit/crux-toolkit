# file      : build/import/libboost/regex/stub.make
# author    : Boris Kolpackov <boris@kolpackov.net>
# copyright : Copyright (c) 2005-2010 Boris Kolpackov
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(scf_root)/import/libboost/configuration-rules.make,$(dcf_root))

libboost_installed :=

$(call -include,$(dcf_root)/import/libboost/configuration-dynamic.make)

ifdef libboost_installed

ifeq ($(libboost_installed),y)

ifeq ($(libboost_system),y)
$(call export,l: -lboost_regex$(libboost_suffix) -lboost_system$(libboost_suffix),cpp_options: )
else
$(call export,l: -lboost_regex$(libboost_suffix),cpp_options: )
endif

else

$(call include-once,$(scf_root)/import/libboost/regex/rules.make,$(dcf_root))

$(call export,\
  l: $(dcf_root)/import/libboost/regex/regex.l,\
  cpp-options: $(dcf_root)/import/libboost/regex/regex.l.cpp-options)

endif

else

.NOTPARALLEL:

endif
