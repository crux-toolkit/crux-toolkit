m4_divert(-1)
# file      : build/meta/vc9proj.m4
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2009-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

m4_changequote([, ])

m4_include(__meta_base__[/common.m4])
m4_include(__meta_base__[/windows-common.m4])

# header
#
m4_define([__header_entry_impl__],
[    <ClInclude Include="[$1]" />])

m4_define([__header_entry__],
[m4_equote()__header_entry_impl__([$1])[]m4_dquote()])


m4_define([__header_entries_impl__],
[m4_foreach_w([__f], __path_impl__([$1]), [
__header_entry_impl__(__f)])])

m4_define([__header_entries__],
[m4_equote()__header_entries_impl__([$1])[]m4_dquote()])

# header_filter
#
m4_define([__header_filter_entry_impl__],
[    <ClInclude Include="[$1]">
      <Filter>Header Files</Filter>
    </ClInclude>])

m4_define([__header_filter_entry__],
[m4_equote()__header_filter_entry_impl__([$1])[]m4_dquote()])


m4_define([__header_filter_entries_impl__],
[m4_foreach_w([__f], __path_impl__([$1]), [
__header_filter_entry_impl__(__f)])])

m4_define([__header_filter_entries__],
[m4_equote()__header_filter_entries_impl__([$1])[]m4_dquote()])

# source
#

# $1 - configuration
# $2 - directory
#
m4_define([__source_config_entry_impl__],
[      <ObjectFileName m4_dnl
Condition="'$(Configuration)|$(Platform)'==m4_dnl
'[]m4_patsubst([$1], ["], [])'">m4_dnl
$(IntDir)\[$2]\</ObjectFileName>])

m4_define([__source_entry_body__],
[m4_ifelse([$1], [$2],
[    <ClCompile Include="[$1]" />],
[    <ClCompile Include="[$1]">m4_dnl
m4_foreach([__c],
[__value_impl__([configurations])],
[
__source_config_entry_impl__(__c, [$2])])
    </ClCompile>])])

m4_define([__source_entry_impl__],
[__source_entry_body__([$1],
m4_patsubst([$1], [^\(.*\)\\\(.*\)$], [\1]))])

m4_define([__source_entry__],
[m4_equote()__source_entry_impl__([$1])[]m4_dquote()])

m4_define([__source_entries_impl__],
[m4_foreach_w([__f], __path_impl__([$1]), [
__source_entry_impl__(__f)])])

m4_define([__source_entries__],
[m4_equote()__source_entries_impl__([$1])[]m4_dquote()])

# source_filter
#
m4_define([__source_filter_entry_impl__],
[    <ClCompile Include="[$1]">
      <Filter>Source Files</Filter>
    </ClCompile>])

m4_define([__source_filter_entry__],
[m4_equote()__source_filter_entry_impl__([$1])[]m4_dquote()])


m4_define([__source_filter_entries_impl__],
[m4_foreach_w([__f], __path_impl__([$1]), [
__source_filter_entry_impl__(__f)])])

m4_define([__source_filter_entries__],
[m4_equote()__source_filter_entries_impl__([$1])[]m4_dquote()])

# custom_build
#
m4_define([__custom_build_entry_tag_impl__],
[      <[$1] Condition="'$(Configuration)|$(Platform)'==m4_dnl
'[]m4_patsubst([$2], ["], [])'">[$3]</[$1]>])

m4_define([__custom_build_entry_impl__],
[    <CustomBuild Include="[$1]">m4_dnl
m4_foreach([__c],
[__value_impl__([configurations])],
[
__custom_build_entry_tag_impl__([Message], __c, [$2])
__custom_build_entry_tag_impl__([Command], __c, [$3])
__custom_build_entry_tag_impl__([Outputs], __c, [$4];%(Outputs))])
    </CustomBuild>])

# (file, cmd-description, cmd, output;output)
#
m4_define([__custom_build_entry__],
[m4_equote()__custom_build_entry_impl__([$1],
[$2], [$3], [$4])[]m4_dquote()])

# Disable comments and quoting.
#
m4_changecom([])
m4_changequote([])

m4_divert(0)m4_dnl
