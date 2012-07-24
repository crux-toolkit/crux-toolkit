
include $(dir $(lastword $(MAKEFILE_LIST)))../../build/frame.make

define set-recursive
$1 = $2
endef

define set-simple
$1 := $2
endef

define save
$(eval $(if $(filter simple,$(flavor $1))\
,$(call set-simple,_$1,$(value $1))\
,$(call set-recursive,_$1,$(value $1))))
endef

define restore
$(eval $(if $(filter simple,$(flavor _$1))\
,$(call set-simple,$1,$(value _$1))\
,$(call set-recursive,$1,$(value _$1))))
endef

save = $(call frame-save,$1)
restore = $(call frame-restore,$1)

r = $(r2)

s := simple

$(call save,r)
$(call save,s)

#$(warning _r: $(value _r))
#$(warning _s: $(value _s))

r := foo
s = bar

$(call restore,r)
$(call restore,s)

r2 := recursive

$(warning r: $r)
$(warning s: $s)

.PHONY: all
all:;@:
