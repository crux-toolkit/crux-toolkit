// file      : xsde/cxx/tree/polymorphism-processor.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#ifndef CXX_TREE_POLYMORPHISM_PROCESSOR_HXX
#define CXX_TREE_POLYMORPHISM_PROCESSOR_HXX

#include <cult/types.hxx>

#include <xsd-frontend/semantic-graph.hxx>

#include <cxx/tree/cli.hxx>

#include <xsd.hxx>

namespace CXX
{
  namespace Tree
  {
    using namespace Cult::Types;

    class PolymorphismProcessor
    {
    public:
      Boolean
      process (CLI::Options const& options,
               XSDFrontend::SemanticGraph::Schema&,
               XSDFrontend::SemanticGraph::Path const& file,
               const WarningSet& disabled_warnings);
    };
  }
}

#endif // CXX_TREE_POLYMORPHISM_PROCESSOR_HXX
