// file      : xsd/type-map/lexer.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2007-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#ifndef XSD_TYPE_MAP_LEXER_HXX
#define XSD_TYPE_MAP_LEXER_HXX

#include <locale>
#include <iosfwd>

#include <cult/types.hxx>

namespace TypeMap
{
  using namespace Cult::Types;
  typedef WideString String;

  class Lexer
  {
  public:
    class Token
    {
    public:
      enum Type
      {
        token,
        punct,
        eos
      };

      Token (Type type, String const& lexeme, UnsignedLong line)
          : type_ (type), lexeme_ (lexeme), line_ (line)
      {
      }

      Type
      type () const
      {
        return type_;
      }

      String const&
      lexeme () const
      {
        return lexeme_;
      }

      UnsignedLong
      line () const
      {
        return line_;
      }

    private:
      Type type_;
      String lexeme_;
      UnsignedLong line_;
    };

    Lexer (std::istream&, String const& path);

    struct Failed {};

    Token
    next ();

  private:
    std::locale locale_;
    std::istream& is_;
    String path_;
    UnsignedLong line_;
    String held_lexeme_;
    Boolean comment_;
  };

}

#endif // XSD_TYPE_MAP_LEXER_HXX

