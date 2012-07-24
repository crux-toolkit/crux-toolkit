# file      : build/cxx/configuration-static.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

cxx_h_suffix := hxx
cxx_t_suffix := txx
cxx_i_suffix := ixx
cxx_s_suffix := cxx

# Get user-supplied static configuration if any.
#
ifneq ($(bld_root),$(scf_root))
$(call -include,$(scf_root)/cxx/configuration-static.make)
endif