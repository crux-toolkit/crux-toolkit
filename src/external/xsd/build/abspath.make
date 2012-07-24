# file      : build/abspath.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

.PHONY: FORCE

%:: FORCE
	@echo $(abspath $@)

$(MAKEFILE_LIST):: ;@:

# Remove all default rules.
#
.SUFFIXES:

ifeq ($(filter -r,$(MAKEFLAGS)),)
MAKEFLAGS += -r
endif
