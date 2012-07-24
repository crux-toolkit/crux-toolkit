# file      : build/meta/vctest-functions.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2009-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

# Process VC++ solution and test template and write output to
# $(dist_prefix)/<path>. Where path is computed as difference
# between src_base and src_root.
#
# Arguments:
#
# $1 - solution path, if doens't start with /, assume relative to
#      dist_prefix/<path>
# $2 - template path, if doesn't start with /, assume relative to src_base
# $3 - output name (optional)
#
$(out_base)/%: meta-vctest = \
$(call meta-vctest-body,$1,$(if $(filter /%,$2),$2,$(src_base)/$2),$(if \
$3,$3,$(notdir $2)),$(subst $(src_root),,$(src_base)))

# $1 - solution
# $2 - template
# $3 - output name
# $4 - difference between src_base and src_root with leading '\'
#
$(out_base)/%: meta-vctest-body = \
$(call message,meta $(dist_prefix)$4/$3,$(bld_root)/meta/vctest \
-r $(dist_prefix) -b $(dist_prefix)$4  -t $2 -o $(dist_prefix)$4/$3 \
$(if $(filter /%,$1),$1,$(dist_prefix)$4/$1))
