# file      : build/meta/common.m4
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2009-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

# m4_strip(STRING)
# ----------------
# Expands into STRING with tabs and spaces singled out into a single
# space, and removing leading and trailing spaces.
#
m4_define([m4_strip],
[m4_patsubst(m4_patsubst([ $1 ], [[	 ]+], [ ]), [^ \(.*\) $], [[[\1]]])])


# m4_normalize(STRING)
# --------------------
# Apply m4_strip to STRING.
#
m4_define([m4_normalize], [m4_strip([$1])])

# m4_split(STRING, [REGEXP])
# --------------------------
# Split STRING into an m4 list of quoted elements.  The elements are
# quoted with ] and ].  Beginning spaces and end spaces *are kept*.
# Use m4_strip to remove them.
#
# REGEXP specifies where to split.  Default is [\t ]+.
#
# If STRING is empty, the result is an empty list.
#

m4_define([m4_split],
[m4_ifelse([$1], [], [],
       [$2], [ ], [m4_ifelse(m4_index([$1], [ ]), [-1], [[[$1]]],
			 [_$0([$1], [$2], [, ])])],
       [$2], [], [_$0([$1], [[	 ]+], [, ])],
       [_$0([$1], [$2], [, ])])])

m4_define([_m4_split], [m4_patsubst([$1], [$2], [$3])])


# Simple foreach implementation. The list elements are expected to
# be fully expanded.
#
m4_define([m4_foreach], [m4_ifelse([$2], [], [],
  [m4_pushdef([$1])_$0([$1], [$3], [], $2)m4_popdef([$1])])])
m4_define([_m4_foreach], [m4_ifelse([$#], [3], [],
  [m4_define([$1], [$4])$2[]$0([$1], [$2],
    m4_shift(m4_shift(m4_shift($@))))])])

# m4_foreach_w(VARIABLE, LIST, EXPRESSION)
# ----------------------------------------
# Like m4_foreach, but the list is whitespace separated.
#
m4_define([m4_foreach_w],
[m4_foreach([$1], m4_split(m4_normalize([$2]), [ ]), [$3])])

# m4_strip_nl(STRING)
# ----------------------------------------
# Remove newlines
#
m4_define([m4_strip_nl], [m4_patsubst([$1], [
], [])])

# Enable/disable quoting.
#
m4_define([m4_equote], [m4_changequote([,])])
m4_define([m4_dquote], [m4_changequote([])])

#
#
m4_define([__env_impl__], [m4_esyscmd([echo -n $$1])])
m4_define([__env__], [m4_equote()__env_impl__([$1])[]m4_dquote()])

m4_define([__value_impl__], [__env_impl__([$1])])
m4_define([__value__], [m4_equote()__value_impl__([$1])[]m4_dquote()])

m4_define([__uuid_impl__],
  [m4_translit(m4_strip_nl(m4_esyscmd([uuidgen])), [a-z], [A-Z])])
m4_define([__uuid__], [m4_equote()__uuid_impl__([$1])[]m4_dquote()])

m4_define([__upcase_impl__], [m4_translit([$1], [a-z], [A-Z])])
m4_define([__upcase__], [m4_equote()__upcase_impl__([$1])[]m4_dquote()])

m4_define([__xml_impl__],
[m4_patsubst(
  m4_patsubst(
    m4_patsubst(
      [$1],
      [&],
      [&amp;]),
    [<],
    [&lt;]),
  ["],
  [&quot;])])

m4_define([__xml__], [m4_equote()__xml_impl__([$1])[]m4_dquote()])
