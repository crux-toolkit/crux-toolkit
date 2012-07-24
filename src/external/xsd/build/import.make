# file      : build/import.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file


# This cannot be eliminated.
#
%frame_exclude% += import-set-value
define import-set-value
$(eval $$1 := $$2)
endef

%frame_exclude% += exported-names
define exported-names
$(foreach i,1 2 3 4 5 6 7,$(if $(value $i),$(lastword $(value $i))))
endef

%frame_exclude% += import-names
define import-names
$(foreach i,1 2 3 4 5 6 7,$(if $(lastword $(value $i)),\
$(call import-set-value,$(lastword $(value $i)),$(value $(firstword $(value $i))))))
endef

%frame_exclude% += import-body
define import-body
$(call frame-enter,$(call exported-names,$2,$3,$4,$5,$6,$7,$8))

$(call include,$1)

$(call import-names,$2,$3,$4,$5,$6,$7,$8)

$(call frame-leave)
endef


# Simple case where there are no names imported. Also note that we use
# include-once here since there is no reason to do it more than once.
# We key it onto out_root to make sure we don't process the same makefile
# more than once since that's what export stubs do.
#
%frame_exclude% += import-body-simple
define import-body-simple
$(call frame-enter)

$(call include-once,$1,$(out_root))

$(call frame-leave)
endef


#
#
%frame_exclude% += import
define import
$(eval $(if $(strip $2),\
$(call import-body,$1,$2,$3,$4,$5,$6,$7,$8),\
$(call import-body-simple,$1)))
endef


#
# export
#
%frame_exclude% += export-names
define export-names
$(foreach i,1 2 3 4 5 6 7,$(if $(value $i),\
$(call import-set-value,$(firstword $(value $i)),$(wordlist 2,$(words $(value $i)),$(value $i)))))
endef


%frame_exclude% += export
define export
$(call export-names,$1,$2,$3,$4,$5,$6,$7)
endef
