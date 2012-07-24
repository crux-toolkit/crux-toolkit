# file      : xsd/makefile
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
# license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

include $(dir $(lastword $(MAKEFILE_LIST)))../build/bootstrap.make

cxx_tun := xsd.cxx

cxx_tun += cxx/elements.cxx                            \
           cxx/literal-map.cxx

cxx_tun += cxx/parser/elements.cxx                     \
           cxx/parser/validator.cxx                    \
           cxx/parser/name-processor.cxx               \
           cxx/parser/type-processor.cxx               \
           cxx/parser/state-processor.cxx              \
           cxx/parser/generator.cxx                    \
           cxx/parser/parser-header.cxx                \
           cxx/parser/parser-inline.cxx                \
           cxx/parser/parser-source.cxx                \
           cxx/parser/parser-forward.cxx               \
	   cxx/parser/impl-header.cxx                  \
           cxx/parser/impl-source.cxx                  \
	   cxx/parser/driver-source.cxx                \
           cxx/parser/element-validation-source.cxx    \
           cxx/parser/attribute-validation-source.cxx  \
           cxx/parser/characters-validation-source.cxx

cxx_tun += cxx/tree/elements.cxx                      \
           cxx/tree/validator.cxx                     \
           cxx/tree/counter.cxx                       \
           cxx/tree/name-processor.cxx                \
           cxx/tree/polymorphism-processor.cxx        \
	   cxx/tree/default-value.cxx                 \
           cxx/tree/generator.cxx                     \
           cxx/tree/tree-forward.cxx                  \
           cxx/tree/tree-header.cxx                   \
           cxx/tree/tree-inline.cxx                   \
           cxx/tree/tree-source.cxx                   \
           cxx/tree/parser-header.cxx                 \
           cxx/tree/parser-source.cxx                 \
           cxx/tree/stream-header.cxx                 \
           cxx/tree/stream-source.cxx                 \
           cxx/tree/serialization-header.cxx          \
           cxx/tree/serialization-source.cxx          \
           cxx/tree/stream-insertion-header.cxx       \
           cxx/tree/stream-insertion-source.cxx       \
           cxx/tree/stream-extraction-source.cxx


# Type map
#
cxx_tun += type-map/lexer.cxx  \
           type-map/parser.cxx

# Processing
#
cxx_tun += processing/cardinality/processor.cxx \
           processing/inheritance/processor.cxx


cxx_obj := $(addprefix $(out_base)/,$(cxx_tun:.cxx=.o))
cxx_od  := $(cxx_obj:.o=.o.d)

xsd      := $(out_base)/xsd
install  := $(out_base)/.install
dist     := $(out_base)/.dist
dist-win := $(out_base)/.dist-win
clean    := $(out_base)/.clean

# Import.
#
$(call import,\
  $(scf_root)/import/libcult/stub.make,\
  l: cult.l,cpp-options: cult.l.cpp-options)

$(call import,\
  $(scf_root)/import/libbackend-elements/stub.make,\
  l: be.l,cpp-options: be.l.cpp-options)

$(call import,\
  $(scf_root)/import/libboost/regex/stub.make,\
  l: re.l,cpp-options: re.l.cpp-options)

$(call import,\
  $(scf_root)/import/libboost/filesystem/stub.make,\
  l: fs.l,cpp-options: fs.l.cpp-options)

$(call import,\
  $(scf_root)/import/libxsd-frontend/stub.make,\
  l: xsd_fe.l,cpp-options: xsd_fe.l.cpp-options)

$(call import,\
  $(scf_root)/import/libxerces-c/stub.make,\
  l: xerces_c.l,cpp-options: xerces_c.l.cpp-options)

# Build.
#
$(xsd): $(cxx_obj) $(xsd_fe.l) $(be.l) $(cult.l) $(fs.l) $(re.l) $(xerces_c.l)

$(cxx_obj) $(cxx_od): cpp_options := -I$(src_base)
$(cxx_obj) $(cxx_od):      \
  $(xsd_fe.l.cpp-options)  \
  $(be.l.cpp-options)      \
  $(cult.l.cpp-options)    \
  $(fs.l.cpp-options)      \
  $(re.l.cpp-options)      \
  $(xerces_c.l.cpp-options)

$(call include-dep,$(cxx_od))

# Alias for default target.
#
$(out_base)/: $(xsd)

# Install
#
$(install): $(xsd)
	$(call install-exec,$<,$(install_bin_dir)/xsd)

# Dist.
#
$(dist): $(xsd)
	$(call install-exec,$<,$(dist_prefix)/bin/xsd)

$(dist-win): $(xsd)
	$(call install-exec,$<,$(dist_prefix)/bin/xsd.exe)

# Clean.
#
$(clean):                            \
  $(xsd).o.clean                     \
  $(addsuffix .cxx.clean,$(cxx_obj)) \
  $(addsuffix .cxx.clean,$(cxx_od))

# Generated .gitignore.
#
ifeq ($(out_base),$(src_base))
$(xsd): | $(out_base)/.gitignore

$(out_base)/.gitignore: files := xsd
$(clean): $(out_base)/.gitignore.clean

$(call include,$(bld_root)/git/gitignore.make)
endif

# How to.
#
$(call include,$(bld_root)/cxx/o-e.make)
$(call include,$(bld_root)/cxx/cxx-o.make)
$(call include,$(bld_root)/cxx/cxx-d.make)
$(call include,$(bld_root)/install.make)
