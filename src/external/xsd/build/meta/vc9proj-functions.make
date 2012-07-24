# file      : build/meta/vc9proj-functions.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2009-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

# Process VC++ project file template and write output to $(dist_prefix)/<path>.
# Where path is computed as difference between src_base and src_root.
#
# Arguments:
#
# $1   - template path, if doesn't start with /, assume relative to src_base
# $2   - output name (optional)
# $3-8 - optional pairs of additional varibales and values $3=$4, $5=$6, etc
#
#
$(out_base)/%: meta-vc9proj = \
$(call meta-vc9proj-body,$(if $(filter /%,$1),$1,$(src_base)/$1),$(if \
$2,$2,$(notdir $1)),$(subst $(src_root),,$(src_base)),$3,$4,$5,$6,$7,$8)

# $1   - template
# $2   - output name
# $3   - difference between src_base and src_root with leading '\'
# $4-9 - additional varibales
#
$(out_base)/%: meta-vc9proj-body = \
$(call message,meta $(dist_prefix)$3/$2,$(if $4,$4='$5'; export $4; )$(if \
$6,$6='$7'; export $6; )$(if $8,$8='$9'; export $8; )$(bld_root)/meta/vc9proj \
-o $(dist_prefix)$3/$2 $1)
