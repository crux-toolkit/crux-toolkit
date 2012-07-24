# file      : build/clean.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

#@@ Maybe create file aliases.make for standard aliases
#   plus check for interactivity? What if some targets
#   are not defined?
#

ifdef %interactive%

.PHONY: clean

clean: $(out_base)/.clean

endif
