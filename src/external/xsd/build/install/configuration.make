# file      : build/install/configuration.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

install_prefix       := /usr/local
install_data_prefix  := $(install_prefix)
install_exec_prefix  := $(install_prefix)

install_bin_dir      := $(install_exec_prefix)/bin
install_sbin_dir     := $(install_exec_prefix)/sbin
install_lib_dir      := $(install_exec_prefix)/lib

install_data_dir     := $(install_data_prefix)/share
install_inc_dir      := $(install_data_prefix)/include

install_doc_dir      := $(install_data_dir)/doc
install_man_dir      := $(install_data_dir)/man
install_info_dir     := $(install_data_dir)/info

$(out_root)/%: install_bin_dir  := $(install_bin_dir)
$(out_root)/%: install_sbin_dir := $(install_sbin_dir)
$(out_root)/%: install_lib_dir  := $(install_lib_dir)

$(out_root)/%: install_inc_dir  := $(install_inc_dir)
$(out_root)/%: install_data_dir := $(install_data_dir)

$(out_root)/%: install_doc_dir  := $(install_doc_dir)
$(out_root)/%: install_man_dir  := $(install_man_dir)
$(out_root)/%: install_info_dir := $(install_info_dir)

install_cmd  := $(bld_root)/install/install

install_dir  := $(install_cmd) -d -m 755
install_data := $(install_cmd) -p -m 644
install_exec := $(install_cmd) -p -m 755

$(out_root)/%: install_dir  := $(install_dir)
$(out_root)/%: install_data := $(install_data)
$(out_root)/%: install_exec := $(install_exec)


#@@ Installation process is very os-dependant. For example
#   on BSD there is libexec but there is no such thing on
#   GNU/Linux.
#
