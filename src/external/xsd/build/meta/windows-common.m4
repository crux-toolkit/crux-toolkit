# file      : build/meta/windows-common.m4
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2009-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

m4_define([__path_impl__], [m4_translit(__env_impl__([$1]),[/],[\])])
m4_define([__path__], [m4_equote()__path_impl__([$1])[]m4_dquote()])

m4_define([__file_impl__], [__env_impl__([$1])])
m4_define([__file__], [m4_equote()__file_impl__([$1])[]m4_dquote()])

# Convert POSIX shell quoting to Windows. Arguments that are in
# single quotes (') are converted to double quotes ("). Double
# quotes that are nested in single quotes are converted to \".
# Quote escaping is not recognized. Single-quoted arguments
# should be quoted as a whole, for example 'foo-bar' and not
# foo-'bar'.
#

m4_define([__shell_quotes_convert_impl__], [m4_patsubst([$1], ["], [\\"])])

# $1 string without the opening quote
# $2 index of the closing quote
#
m4_define([__shell_quotes_in_quote_impl__],
["__shell_quotes_convert_impl__(m4_substr([$1], [0], [$2]))"[]m4_dnl
 __shell_quotes_impl__(m4_substr([$1], m4_incr([$2])))])

# $1 string without the leading space
# $2 index of the next space or -1
#
m4_define([__shell_quotes_in_space_impl__],
[m4_substr([$1], [0], m4_ifelse([$2], [-1], [m4_len([$1])], [$2]))[]m4_dnl
m4_ifelse(
  [$2], [-1],
  [],
  [ __shell_quotes_impl__(m4_substr([$1], m4_incr([$2])))])])

#
#
m4_define([__shell_quotes_impl__],
[m4_ifelse(
  m4_index([$1], [']), [0],
  [__shell_quotes_in_quote_impl__(
    m4_substr([$1], [1]),
    m4_index(m4_substr([$1], [1]), [']))],
  [__shell_quotes_in_space_impl__([$1], m4_index([$1], [ ]))])])

m4_define([__shell_quotes__],
[m4_equote()__shell_quotes_impl__([$1])[]m4_dquote()])
