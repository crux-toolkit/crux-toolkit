# file      : build/c/generic/c-o.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(call include,$(bld_root)/c/generic/configuration.make)

# Make will try to build dependecies (since they are ultimately included
# files) during configuartion phase without c_generic being discovered
# yet. This is also why dependecies should be included with -include.
#
ifdef c_generic

.PHONY: $(out_base)/%.o.d.$(c_s_suffix).clean

$(out_base)/%.o.d.$(c_s_suffix).clean:
	@:

endif
