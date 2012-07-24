# file      : build/literals.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

literal_empty :=

define literal_newline


endef

literal_tab := $(literal_empty)	$(literal_empty)
literal_percent := %
