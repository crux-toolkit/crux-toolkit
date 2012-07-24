# file      : build/configuration.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
# license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

$(call include-once,$(scf_root)/configuration-rules.make,$(dcf_root))

# Dynamic configuration.
#
xsd_with_zlib :=
xsd_with_ace :=
xsd_with_xdr :=
xsd_with_dbxml :=
xsd_with_xqilla :=
xsd_with_boost_date_time :=
xsd_with_boost_serialization :=

$(call -include,$(dcf_root)/configuration-dynamic.make)

ifdef xsd_with_zlib

$(out_root)/%: xsd_with_zlib := $(xsd_with_zlib)
$(out_root)/%: xsd_with_ace := $(xsd_with_ace)
$(out_root)/%: xsd_with_xdr := $(xsd_with_xdr)
$(out_root)/%: xsd_with_dbxml := $(xsd_with_dbxml)
$(out_root)/%: xsd_with_xqilla := $(xsd_with_xqilla)
$(out_root)/%: xsd_with_boost_date_time := $(xsd_with_boost_date_time)
$(out_root)/%: xsd_with_boost_serialization := $(xsd_with_boost_serialization)

else

.NOTPARALLEL:

endif
