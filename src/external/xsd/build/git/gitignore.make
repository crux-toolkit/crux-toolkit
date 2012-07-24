# file      : build/git/gitignore.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(out_base)/.gitignore:
	@$(bld_root)/git/gitignore $(files) >$@

.PHONY: $(out_base)/.gitignore.clean
$(out_base)/.gitignore.clean:
	@rm -f $(basename $@)
