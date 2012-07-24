// file      : xsd/cxx/parser/validator.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#ifndef CXX_PARSER_VALIDATOR_HXX
#define CXX_PARSER_VALIDATOR_HXX

#include <cxx/parser/elements.hxx>
#include <cxx/parser/cli.hxx>

#include <xsd.hxx>

namespace CXX
{
  namespace Parser
  {
    using namespace Cult::Types;

    class Validator
    {
    public:
      Validator (); // Dummy ctor, helps with long symbols on HP-UX.

      Boolean
      validate (CLI::Options const& options,
                SemanticGraph::Schema&,
                SemanticGraph::Path const& tu,
                Boolean gen_driver,
                const WarningSet& disabled_warnings);
    };
  }
}

#endif  // CXX_PARSER_VALIDATOR_HXX
