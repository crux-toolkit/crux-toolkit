m4_divert(-1)
# file      : build/meta/vcsln.m4
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2009-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

m4_changequote([, ])

m4_include(__meta_base__[/common.m4])
m4_include(__meta_base__[/windows-common.m4])

# solution_configuration
#
m4_define([__solution_configuration_entry_impl__], [		[$1] = [$1]])

m4_define([__solution_configurations_impl__],
  [m4_foreach([__c],
    [__value_impl__([configurations])],
    [__solution_configuration_entry_impl__(m4_patsubst(__c, ["], []))
])])

m4_define([__solution_configurations__],
  [m4_equote()__solution_configurations_impl__()[]m4_dquote()])

# project_configuration
#
m4_define([__project_configuration_entry_item_impl__],
[		{[$1]}.[$2].ActiveCfg = [$2]
		{[$1]}.[$2].Build.0 = [$2]])

m4_define([__project_configuration_entry_impl__],
  [m4_foreach(
    [__c],
    [__value_impl__([configurations])],
    [__project_configuration_entry_item_impl__([$1], m4_patsubst(__c, ["], []))
])])

m4_define([__project_configurations_impl__],
  [m4_foreach_w(
    [__u],
    [__value_impl__([project_uuids])],
    [__project_configuration_entry_impl__(__u)])])

m4_define([__project_configurations__],
  [m4_equote()__project_configurations_impl__()[]m4_dquote()])

# project_entry(name, file, uuid)
#
m4_define([__project_entry_impl__],
[Project("{8BC9CEB8-8B4A-11D0-8D11-00A0C91BC942}") = "[$1]", "[$2]", "{[$3]}"
EndProject])

m4_define([__project_entry__],
  [m4_equote()__project_entry_impl__([$1], [$2], [$3])[]m4_dquote()])

# projects
#
m4_define([__project_step_impl__],
[__project_entry_impl__(
  m4_substr([$1], 0, m4_ifelse(m4_index([$1], [ ]), [-1], [m4_len([$1])], [m4_index([$1], [ ])])),
  m4_substr([$2], 0, m4_ifelse(m4_index([$2], [ ]), [-1], [m4_len([$2])], [m4_index([$2], [ ])])),
  m4_substr([$3], 0, m4_ifelse(m4_index([$3], [ ]), [-1], [m4_len([$3])], [m4_index([$3], [ ])])))[]m4_dnl
m4_ifelse(m4_index([$1], [ ]), [-1],, [
__project_step_impl__(
  m4_substr([$1], m4_incr(m4_index([$1], [ ]))),
  m4_substr([$2], m4_incr(m4_index([$2], [ ]))),
  m4_substr([$3], m4_incr(m4_index([$3], [ ]))))])])

m4_define([__projects_impl__],
[__project_step_impl__(
  __value_impl__([project_names]),
  __path_impl__([project_files]),
  __value_impl__([project_uuids]))])

m4_define([__projects__],
  [m4_equote()__projects_impl__()[]m4_dquote()])

# Disable comments and quoting.
#
m4_changecom([])
m4_changequote([])

m4_divert(0)m4_dnl
