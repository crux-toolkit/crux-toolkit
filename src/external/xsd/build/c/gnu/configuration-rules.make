# file      : build/c/gnu/configuration-rules.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(dcf_root)/c/gnu/configuration-dynamic.make: | $(dcf_root)/c/gnu/.
	$(call message,,$(bld_root)/c/gnu/configure $@ $(c_optimize))

ifndef %foreign%

$(dcf_root)/.disfigure::
	$(call message,rm $(dcf_root)/c/gnu/configuration-dynamic.make,\
rm -f $(dcf_root)/c/gnu/configuration-dynamic.make)

endif
