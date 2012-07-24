# file      : examples/cxx/hello/libhello/build/import/libhello/stub.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(scf_root)/import/libhello/configuration-rules.make,$(dcf_root))

libhello_installed :=

$(call -include,$(dcf_root)/import/libhello/configuration-dynamic.make)

ifdef libhello_installed

ifeq ($(libhello_installed),y)

$(call export,l: -lhello,cpp-options:)

else

# Include export stub.
#
$(call include,$(scf_root)/export/libhello/stub.make)

endif

else

.NOTPARALLEL:

endif
