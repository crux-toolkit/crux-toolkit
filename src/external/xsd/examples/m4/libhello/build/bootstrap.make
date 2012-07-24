# file      : examples/m4/hello/libhello/build/bootstrap.make
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

project_name := libhello

include $(dir $(lastword $(MAKEFILE_LIST)))../../../../build/bootstrap.make
