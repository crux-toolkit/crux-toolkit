# file      : build/import/libbackend-elements/stub.make
# author    : Boris Kolpackov <boris@kolpackov.net>
# copyright : Copyright (c) 2005-2009 Boris Kolpackov
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(scf_root)/import/libbackend-elements/configuration-rules.make,$(dcf_root))

libbackend_elements_installed :=

$(call -include,$(dcf_root)/import/libbackend-elements/configuration-dynamic.make)

ifdef libbackend_elements_installed

ifeq ($(libbackend_elements_installed),y)

#-lbackend-elements

$(call export,l: -lcult -lboost_regex,cpp_options: )

else

# Include export stub.
#
$(call include,$(scf_root)/export/libbackend-elements/stub.make)

endif

else

.NOTPARALLEL:

endif
