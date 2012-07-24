// file      : xsd/cxx/tree/serialization-source.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#ifndef CXX_TREE_SERIALIZATION_SOURCE_HXX
#define CXX_TREE_SERIALIZATION_SOURCE_HXX

#include <cxx/tree/elements.hxx>

namespace CXX
{
  namespace Tree
  {
    Void
    generate_serialization_source (Context&,
                                   UnsignedLong first,
                                   UnsignedLong last);
  }
}

#endif  // CXX_TREE_SERIALIZATION_SOURCE_HXX
