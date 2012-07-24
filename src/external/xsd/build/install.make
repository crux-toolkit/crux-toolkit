# file      : build/install.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include,$(bld_root)/install/configuration.make)
$(call include-once,$(bld_root)/install/install-functions.make,$(out_base))
