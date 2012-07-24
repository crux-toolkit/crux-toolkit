// file      : xsd/type-map/parser.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2007-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#ifndef XSD_TYPE_MAP_PARSER_HXX
#define XSD_TYPE_MAP_PARSER_HXX

#include <cult/types.hxx>

#include <type-map/type-map.hxx>
#include <type-map/lexer.hxx>

namespace TypeMap
{
  using namespace Cult::Types;
  typedef WideString String;

  class Parser
  {
  public:
    Parser (Lexer&, String const& path);

    // Merge parsed namespaces.
    //
    Boolean
    parse (Namespaces&);

  private:
    Boolean
    namespace_ (Namespaces&);

    Boolean
    include (Namespace&);

    Boolean
    type (Lexer::Token, Namespace&);

  private:
    Lexer& lex_;
    String path_;
    std::wostream& e;
  };
}

#endif // XSD_TYPE_MAP_PARSER_HXX
