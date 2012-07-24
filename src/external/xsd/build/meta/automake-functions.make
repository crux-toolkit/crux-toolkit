# file      : build/meta/automake-functions.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2009-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

# Process automake template and write output to $(dist_prefix)/<path>.
# Where path is computed as difference between src_base and src_root.
#
# Arguments:
#
# $1 - optional template path, if doesn't start with /, assume relative
#      to src_base default is $(src_base)/Makefile.am
#
$(out_base)/%: meta-automake = \
$(call meta-automake-body,$(if $1,$(if $(filter \
/%,$1),$1,$(src_base)/$1),$(src_base)/Makefile.am),$(subst \
$(src_root),,$(src_base)))

# $1 - template
# $2 - difference between src_base and src_root with leading '\'
#
$(out_base)/%: meta-automake-body = \
$(call message,meta $(dist_prefix)$2/$(notdir $1),$(bld_root)/meta/automake \
-o $(dist_prefix)$2/$(notdir $1) $1)
