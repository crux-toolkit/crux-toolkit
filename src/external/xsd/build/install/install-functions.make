# file      : build/install/install-functions.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

#@@ could just call it functions.make
#

$(out_base)/%: install-exec = \
$(call message,,$(install_dir) $(dir $2))$(literal_newline)\
$(literal_tab)$(call message,install $2,$(install_exec) $1 $2)


$(out_base)/%: install-lib = \
$(call message,install $2,$(install_dir) $(dir $2))$(literal_newline)\
$(literal_tab)$(call message,,$(value $(if $(patsubst %.so,,$(shell\
 sed -e '1q' <$1)),install_data,install_exec)) $(shell sed -e '1q' <$1) $2)


$(out_base)/%: install-data = \
$(call message,install $2,$(install_dir) $(dir $2))$(literal_newline)\
$(literal_tab)$(call message,,$(install_data) $1 $2)


$(out_base)/%: install-dir = \
$(call message,install $2/,find -L $1 $3 -type f -print \
| xargs -n 1 $(bld_root)/run-if-arg dirname \
| sort -u \
| sed -e 's%$1\(.*\)%$2/\1%' \
| xargs -n 1 $(bld_root)/run-if-arg "$(install_dir)")$(literal_newline)\
$(literal_tab)$(call message,,\
find -L $1 $3 -type f ! -perm -100 -print \
| sed -e 's%$1\(.*\)%$1/\1 $2/\1%' \
| xargs -n 2 $(bld_root)/run-if-arg "$(install_data)")$(literal_newline)\
$(literal_tab)$(call message,,\
find -L $1 $3 -type f -perm -100 -print  \
| sed -e 's%$1\(.*\)%$1/\1 $2/\1%' \
| xargs -n 2 $(bld_root)/run-if-arg "$(install_exec)")
