# file      : build/import/libcult/stub.make
# author    : Boris Kolpackov <boris@kolpackov.net>
# copyright : Copyright (c) 2005-2009 Boris Kolpackov
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(scf_root)/import/libcult/configuration-rules.make,$(dcf_root))

libcult_installed :=

$(call -include,$(dcf_root)/import/libcult/configuration-dynamic.make)

ifdef libcult_installed

ifeq ($(libcult_installed),y)

$(call export,l: -lcult,cpp-options: )

else

# Include export stub.
#
$(call include,$(scf_root)/export/libcult/stub.make)

endif

else

.NOTPARALLEL:

endif
