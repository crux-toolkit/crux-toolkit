# file      : build/import/makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))../../build/bootstrap.make

foo := foz

$(call frame-push,foo bar)

bar := bar

foo := foo

$(call frame-pop)

# a := $(call frame-vars)

#$(warning $(filter-out $b,$a))

$(warning foo = $(foo))
$(warning bar = $(bar))

exclude := $(call exported-names,default: libhello,clean: libhello_clean)

$(call import-names,default: libhello,clean: libhello_clean)

$(warning $(default:))
$(warning $(clean:))
$(warning $(exclude))

$(call frame-push)

$(default:) := fofofof

$(warning $(libhello))

$(call frame-pop)


.PHONY: h
h:;@:
