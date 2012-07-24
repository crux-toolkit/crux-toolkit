# file      : examples/cxx/multidir/build/bootstrap.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

project_name := multidir

include $(dir $(lastword $(MAKEFILE_LIST)))../../../../build/bootstrap.make
