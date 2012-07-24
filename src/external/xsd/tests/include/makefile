# file      : build/include/makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))../../build/bootstrap.make

foo := one

$(call include-once,included,$(foo))
$(call include-once,included,$(foo))

foo := two

$(call include-once,included,$(foo))
$(call include-once,included,$(foo))

foo := one

$(call include-once,included,$(foo))

$(warning $(value %include_once_$(realpath included)%))

.PHONY: h
h:;@:
