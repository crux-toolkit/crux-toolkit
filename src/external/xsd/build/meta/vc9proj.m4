m4_divert(-1)
# file      : build/meta/vc9proj.m4
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2009-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

m4_changequote([, ])

m4_include(__meta_base__[/common.m4])
m4_include(__meta_base__[/windows-common.m4])

# file
#
m4_define([__file_entry_impl__],[
			<File
				RelativePath=".\[$1]"
				>
			</File>])

m4_define([__file_entry__],
[m4_equote()__file_entry_impl__([$1])[]m4_dquote()])


m4_define([__file_entries_impl__],
[m4_foreach_w([__f], __path_impl__([$1]), [__file_entry_impl__(__f)])])

m4_define([__file_entries__],
[m4_equote()__file_entries_impl__([$1])[]m4_dquote()])

# file_custom_build
#
m4_define([__custom_build_entry_impl__],[
				<FileConfiguration
					Name=[$1]
					>
					<Tool
						Name="VCCustomBuildTool"
						Description="[$2]"
						CommandLine="[$3]"
						Outputs="[$4]"
					/>
				</FileConfiguration>])

m4_define([__file_entry_custom_build_impl__],[
			<File
				RelativePath=".\[$1]"
				>m4_dnl
m4_foreach([__c],
[__value_impl__([configurations])],
[__custom_build_entry_impl__(__c, [$2], [$3], [$4])])
			</File>])

# (file, cmd-description, cmd, output;output)
#
m4_define([__file_entry_custom_build__],
[m4_equote()__file_entry_custom_build_impl__([$1],
[$2], [$3], [$4])[]m4_dquote()])


# source
#

# $1 - configuration
# $2 - directory
#
m4_define([__source_config_entry_impl__],[
				<FileConfiguration
					Name=[$1]
					>
					<Tool
						Name="VCCLCompilerTool"
						ObjectFile="$(IntDir)\[$2]\"
					/>
				</FileConfiguration>])

m4_define([__source_config_impl__],
[m4_ifelse([$1], [$2], [],
[m4_foreach([__c],
[__value_impl__([configurations])],
[__source_config_entry_impl__(__c, [$1])])])])

m4_define([__source_entry_impl__],[
			<File
				RelativePath=".\[$1]"
				>m4_dnl
__source_config_impl__(m4_patsubst([$1], [^\(.*\)\\\(.*\)$], [\1]), [$1])
			</File>])

m4_define([__source_entry__],
[m4_equote()__source_entry_impl__([$1])[]m4_dquote()])

m4_define([__source_entries_impl__],
[m4_foreach_w([__f], __path_impl__([$1]), [__source_entry_impl__(__f)])])

m4_define([__source_entries__],
[m4_equote()__source_entries_impl__([$1])[]m4_dquote()])

# Disable comments and quoting.
#
m4_changecom([])
m4_changequote([])

m4_divert(0)m4_dnl
