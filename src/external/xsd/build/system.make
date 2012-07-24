# file      : build/system.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2006-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include-once,$(bld_root)/system/configuration-rules.make,$(dcf_root))

$(call -include,$(dcf_root)/system/configuration-dynamic.make)

ifndef system_configuration

.NOTPARALLEL:

endif
