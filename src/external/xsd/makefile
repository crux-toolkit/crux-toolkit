# file      : makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))build/bootstrap.make

# Help.
#
.PHONY: help

help:
	@echo
	@echo "Available targets:"
	@echo
	@echo "  help                                - print this message"
	@echo "  install [install_prefix=<prefix>]   - install build system"
	@echo
	@echo "You can also build various examples by executing make in "
	@echo "corresponding directories."
	@echo


# Install.
#
.PHONY: install $(out_base)/.install

install: $(out_base)/.install

$(out_base)/.install:
	$(call install-dir,$(bld_root),$(install_inc_dir)/build-$(shell \
sed -e 's%\(.*\)\.\(.*\)\.\(.*\)%\1.\2%' <$(src_base)/version))


# How to.
#
$(call include,$(bld_root)/install.make)
