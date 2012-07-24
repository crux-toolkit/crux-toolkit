# file      : documentation/cxx/makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
# license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))../../build/bootstrap.make

docs := parser tree

default  := $(out_base)/
install  := $(out_base)/.install
dist     := $(out_base)/.dist
dist-win := $(out_base)/.dist-win
cleandoc := $(out_base)/.cleandoc

$(default): $(addprefix $(out_base)/,$(addsuffix /,$(docs)))
$(install): $(addprefix $(out_base)/,$(addsuffix /.install,$(docs)))
$(dist): $(addprefix $(out_base)/,$(addsuffix /.dist,$(docs)))
$(dist-win): $(addprefix $(out_base)/,$(addsuffix /.dist-win,$(docs)))
$(cleandoc): $(addprefix $(out_base)/,$(addsuffix /.cleandoc,$(docs)))

$(foreach m,$(docs),$(call import,$(src_base)/$m/makefile))
