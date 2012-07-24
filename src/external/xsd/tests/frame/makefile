# file      : build/frame/makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))../../build/bootstrap.make

out_dir := foo

%foo% := baz

$(warning $(out_dir))

# $(warning $(call frame-vars))

$(call frame-enter)

out_dir := foz
src_dir := baz

new := new var

$(call frame-enter)

out_dir := faz

$(warning $(out_dir))
$(warning $(new))

$(call frame-leave)

$(warning $(out_dir))
$(warning $(new))

$(call frame-leave)

$(warning $(out_dir))
$(warning $(new))

.PHONY: all
all:;@:
