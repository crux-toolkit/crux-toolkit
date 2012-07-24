# file      : build/message.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

%frame_exclude% += build-message-expand-body
define build-message-expand-body
$(eval _1 = $1)$(call _1,$2,$3,$4,$5,$6,$7,$8,$9)
endef

%frame_exclude% += build-message-expand
define build-message-expand
$(call build-message-expand-body,$(subst #,\#,$1),$2,$3,$4,$5,$6,$7,$8,$9)
endef

%frame_exclude% += message

ifdef verbose

define message
$(call build-message-expand,$2,$3,$4,$5,$6,$7,$8,$9)
endef

else

define message
$(if $1,@echo $(call build-message-expand,$1,$3,$4,$5,$6,$7,$8,$9) && \
              $(call build-message-expand,$2,$3,$4,$5,$6,$7,$8,$9),\
        @$(call build-message-expand,$2,$3,$4,$5,$6,$7,$8,$9))
endef

endif
