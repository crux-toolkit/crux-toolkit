m4_divert(-1)
# file      : build/meta/automake.m4
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2009-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

m4_changequote([, ])

m4_include(__meta_base__[/common.m4])

m4_define([__path_impl__], [__env_impl__([$1])])
m4_define([__path__], [m4_equote()__path_impl__([$1])[]m4_dquote()])

m4_define([__file_impl__], [__env_impl__([$1])])
m4_define([__file__], [m4_equote()__file_impl__([$1])[]m4_dquote()])

# Disable comments and quoting.
#
m4_changecom([])
m4_changequote([])

m4_divert(0)m4_dnl
