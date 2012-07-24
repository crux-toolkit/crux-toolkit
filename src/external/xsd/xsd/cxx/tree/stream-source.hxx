// file      : xsd/cxx/tree/stream-source.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#ifndef CXX_TREE_STREAM_SOURCE_HXX
#define CXX_TREE_STREAM_SOURCE_HXX

#include <cxx/tree/elements.hxx>

namespace CXX
{
  namespace Tree
  {
    Void
    generate_stream_source (Context&,
                            UnsignedLong first,
                            UnsignedLong last);
  }
}

#endif  // CXX_TREE_STREAM_SOURCE_HXX
