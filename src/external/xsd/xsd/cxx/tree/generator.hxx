// file      : xsd/cxx/tree/generator.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#ifndef CXX_TREE_GENERATOR_HXX
#define CXX_TREE_GENERATOR_HXX

#include <cult/types.hxx>

#include <xsd-frontend/semantic-graph/elements.hxx> // Path
#include <xsd-frontend/semantic-graph/schema.hxx>

#include <xsd.hxx>

#include <cxx/literal-map.hxx>
#include <cxx/tree/cli.hxx>

namespace CXX
{
  namespace Tree
  {
    using namespace Cult::Types;

    class Generator
    {
    public:
      static Void
      usage ();

      static CLI::OptionsSpec
      options_spec ();

      struct Failed {};

      static UnsignedLong
      generate (CLI::Options const& options,
                XSDFrontend::SemanticGraph::Schema&,
                XSDFrontend::SemanticGraph::Path const& file,
                Boolean file_per_type,
                StringLiteralMap const&,
                const WarningSet& disabled_warnings,
                FileList& file_list,
                AutoUnlinks& unlinks);

    private:
      Generator ();
    };
  }
}

#endif // CXX_TREE_GENERATOR_HXX
