# file      : build/ld/configuration-lib-rules.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(dcf_root)/ld/configuration-lib-dynamic.make: | $(dcf_root)/ld/.
	$(call message,,$(bld_root)/ld/configure-lib $@)


ifndef %foreign%

$(dcf_root)/.disfigure::
	$(call message,rm $(dcf_root)/ld/configuration-lib-dynamic.make,\
rm -f $(dcf_root)/ld/configuration-lib-dynamic.make)

endif
