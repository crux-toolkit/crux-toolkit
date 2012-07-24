# file      : build/m4/m4.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2005-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

$(out_base)/%: m4 := m4
$(out_base)/%: m4_options +=

ifeq ($(out_base),$(src_base))
$(out_base)/%: $(src_base)/%.m4
else
$(out_base)/%: $(src_base)/%.m4 | $$(dir $$@).
endif
	$(call message,m4 $<,$(m4) $(m4_options) $< >$@)

ifneq ($(out_base),$(src_base))

$(out_base)/%: $(out_base)/%.m4 | $$(dir $$@).
	$(call message,m4 $<,$(m4) $(m4_options) $< >$@)

endif


# @@
# This is where things start breaking. Following standard logic I should
# make a $(out_base)/%.clean rule, i.e., "will clean anything" rule. If
# this rule happened to be before some other, more specialized rule, and
# that rule happened to rm some additional stuff (like %.o tries to rm
# .d file, which is also not quite correct...). In other word there
# doesn't seem to be a way to properly match "build" and "clean" rules.
# One idea is to make the "clean" rule depend on what "build" rule
# depends (%.m4 in our case) hoping that this way the rule won't match.
#
# There are two problems with this approach:
#
# 1. It is if not iff. However, since the rules come in pairs and make
#    pick the first implicit rule that matches, it is certain that if
#    make picked this "clean" rule it also picked corresponding "build"
#    rule.
#
# 2. The prerequisite (%.m4) can be an intermidiate file which itself
#    may not exist. We don't want make to build it just to clean it
#    or, even worse, to leave it laying around. I guess the only way
#    to work around this is to provide special do-nothing rules during
#    cleanup.
#
#
.PHONY: $(out_base)/%.m4.clean

$(out_base)/%.m4.clean:
	$(call message,rm $(@:.m4.clean=),rm -f $(@:.m4.clean=))
