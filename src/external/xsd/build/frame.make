# file      : build/frame.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

%frame_exclude% += CURDIR SHELL MAKEFILE_LIST MAKEFLAGS
%frame_include% := #.DEFAULT_GOAL - for some reason it is already in the list


# return only vars with 'file' origin
#
%frame_exclude% += frame-vars-stage
define frame-vars-stage
$(foreach v,$1,$(if $(findstring file,$(origin $v)),$v))
endef


# $1 holds exclusion list
#
%frame_exclude% += frame-vars
define frame-vars
$(call frame-vars-stage,$(filter-out \%% $(%frame_exclude%) $1,$(.VARIABLES))) \
$(filter-out $1,$(%frame_include%))
endef


%frame% := frame%


%frame_exclude% += frame-set-recursive
define frame-set-recursive
define $1
$2
endef
endef

%frame_exclude% += frame-undefine
ifneq ($(filter undefine,$(.FEATURES)),)
define frame-undefine
$(eval undefine $$1)
endef
define frame-undefine-imm
$(eval undefine $1)
endef
else
define frame-undefine
$(eval $$1 :=)
endef
define frame-undefine-imm
$(eval $1 :=)
endef
endif

%frame_exclude% += frame-save
define frame-save
$(eval $(if $(filter simple,$(flavor $1)),\
%$1/$(%frame%) := $(value $1),\
$(call frame-set-recursive,%$1/$(%frame%),$(value $1))))
endef

%frame_exclude% += frame-restore
define frame-restore
$(eval $(if $(filter simple,$(flavor %$1/$(%frame%))),\
$1 := $(value %$1/$(%frame%)),\
$(call frame-set-recursive,$1,$(value %$1/$(%frame%)))))\
$(call frame-undefine-imm,%$1/$(%frame%))
endef

# Use debug messages to check for garbage being framed.
#

#$$(warning framing $$(value %vars_$(%frame%)))

%frame_exclude% += frame-enter-body
define frame-enter-body
%vars_$(%frame%) := $(call frame-vars,$1)
%excl_$(%frame%) := $1
$$(foreach v,$$(value %vars_$(%frame%)),$$(call frame-save,$$v))
%frame% := frame/$(%frame%)
endef

%frame_exclude% += frame-enter
define frame-enter
$(eval $(call frame-enter-body,$1))
endef


#$$(warning restoring $$(value %vars_$$(%frame%)))

%frame_exclude% += frame-leave-body
define frame-leave-body
%frame% := $(patsubst frame/%,%,$(%frame%))
$$(foreach v,$$(value %vars_$$(%frame%)),$$(call frame-restore,$$v))
$$(foreach v,\
$$(filter-out $$(value %vars_$$(%frame%)),$$(call frame-vars,$$(value %excl_$$(%frame%)))),\
$$(call frame-undefine,$$v))
$$(call frame-undefine,%vars_$$(%frame%))
$$(call frame-undefine,%excl_$$(%frame%))
endef

%frame_exclude% += frame-leave
define frame-leave
$(eval $(call frame-leave-body))
endef
