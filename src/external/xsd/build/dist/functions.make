# file      : build/dist/functions.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2009-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

dist_cmd  := $(bld_root)/install/install

dist_dir  := $(dist_cmd) -d -m 755
dist_data := $(dist_cmd) -p -m 644
dist_exec := $(dist_cmd) -p -m 755

$(out_root)/%: dist_dir  := $(dist_dir)
$(out_root)/%: dist_data := $(dist_data)
$(out_root)/%: dist_exec := $(dist_exec)

# Arguments:
#
# $1 - files; if not an absolute path, assume relative to src_base
# $2 - optional destination directory. If not specified,
#      $(dist_prefix)/<path> is used where path is computed as
#      difference between src_base and src_root.
#
$(out_base)/%: dist-data = \
$(call dist-data-body,$1,$(if $2,$2,$(dist_prefix)$(subst \
$(src_root),,$(src_base))))

$(out_base)/%: dist-exec = \
$(call dist-exec-body,$1,$(if $2,$2,$(dist_prefix)$(subst \
$(src_root),,$(src_base))))

# Arguments:
#
# $1 - directory
# $2 - optional find predicates
# $3 - optional destination directory. If not specified,
#      $(dist_prefix)/<path> is used where path is computed as
#      difference between src_base and src_root.
#
$(out_base)/%: dist-dir = \
$(call dist-dir-body,$1,$2,$(if $3,$3,$(dist_prefix)$(subst \
$(src_root),,$(src_base))))

$(out_base)/%: dist-data-body = \
$(call message,,$(dist_dir) $2)\
$(foreach d,$(sort $(dir $(addprefix $2/,$(foreach f,$1,$(if \
$(filter /%,$f), $(notdir $f), $f))))),$(literal_newline)\
$(literal_tab)$(call message,,$(dist_dir) $d))\
$(foreach f,$1,$(literal_newline)\
$(literal_tab)$(call message,dist $2/$(if $(filter /%,$f),$(notdir \
$f),$f),$(dist_data) $(if $(filter /%,$f),$f,$(src_base)/$f) $2/$(if \
$(filter /%,$f),$(notdir $f),$f)))

$(out_base)/%: dist-exec-body = \
$(call message,,$(dist_dir) $2)\
$(foreach d,$(sort $(dir $(addprefix $2/,$(foreach f,$1,$(if \
$(filter /%,$f), $(notdir $f), $f))))),$(literal_newline)\
$(literal_tab)$(call message,,$(dist_dir) $d))\
$(foreach f,$1,$(literal_newline)\
$(literal_tab)$(call message,dist $2/$(if $(filter /%,$f),$(notdir \
$f),$f),$(dist_exec) $(if $(filter /%,$f),$f,$(src_base)/$f) $2/$(if \
$(filter /%,$f),$(notdir $f),$f)))

$(out_base)/%: dist-dir-body = \
$(call message,dist $3/$1,find -L $(src_base)/$1 $2 -type f -print \
| xargs -n 1 $(bld_root)/run-if-arg dirname \
| sort -u \
| sed -e 's%$(src_base)/$1\(.*\)%$3/$1\1%' \
| xargs -n 1 $(bld_root)/run-if-arg "$(dist_dir)")$(literal_newline)\
$(literal_tab)$(call message,,\
find -L $(src_base)/$1 $2 -type f ! -perm -100 -print \
| sed -e 's%$(src_base)/$1\(.*\)%$(src_base)/$1\1 $3/$1\1%' \
| xargs -n 2 $(bld_root)/run-if-arg "$(dist_data)")$(literal_newline)\
$(literal_tab)$(call message,,\
find -L $(src_base)/$1 $2 -type f -perm -100 -print \
| sed -e 's%$(src_base)/$1\(.*\)%$(src_base)/$1\1 $3/$1\1%' \
| xargs -n 2 $(bld_root)/run-if-arg "$(dist_exec)")
