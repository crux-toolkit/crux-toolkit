// file      : xsd/type-map/lexer.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2007-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <iostream>

#include <type-map/lexer.hxx>

using std::wcerr;
using std::endl;

namespace TypeMap
{
  Lexer::Lexer (std::istream& is, String const& path)
      : locale_ ("C"), is_ (is), path_ (path), line_ (1), comment_ (false)
  {
    is_.exceptions (std::ios_base::badbit);
  }

  Lexer::Token Lexer::
  next ()
  {
    if (held_lexeme_)
    {
      Token t (Token::punct, held_lexeme_, line_);
      held_lexeme_.clear ();
      return t;
    }

    typedef std::char_traits<char> Traits;
    typedef Traits::char_type CharType;
    typedef Traits::int_type IntType;

    IntType i;
    CharType c ('\0');
    NarrowString lexeme;

    // Skip all whitespaces including comments.
    //
    while (!is_.eof ())
    {
      i = is_.get ();

      if (i == Traits::eof ())
        break;

      c = Traits::to_char_type (i);

      if (comment_)
      {
        if (c == '\n')
          comment_ = false;
      }
      else
      {
        if (!(std::isspace (c, locale_) || c == '#'))
          break;

        if (c == '#')
          comment_ = true;
      }

      if (c == '\n')
        ++line_;
    }

    if (is_.eof ())
      return Token (Token::eos, L"<end-of-stream>", line_);

    Boolean quote (c == '"');

    if (!quote)
      lexeme += c;

    if (c != ';' && c != '{' && c != '}')
    {
      // Accumulate non-whitespace character sequence.
      //

      while (!is_.eof ())
      {
        i = is_.get ();

        if (i == Traits::eof ())
          break;

        c = Traits::to_char_type (i);

        if (!quote && c == '#')
        {
          comment_ = true;
          break;
        }

        if (std::isspace (c, locale_))
        {
          if (c == '\n')
            ++line_;

          if (!quote)
            break;
        }

        if (!quote && (c == ';' || c == '{' || c == '}'))
        {
          held_lexeme_ += c;
          break;
        }

        if (quote && c == '"')
          break;

        lexeme += c;
      }

      if (quote && c != '"')
      {
        wcerr << path_ << ":" << line_ << ": error: closing '\"' expected"
              << endl;

        throw Failed ();
      }
    }

    if (!quote && (lexeme == ";" || lexeme == "{" || lexeme == "}"))
    {
      return Token (Token::punct, lexeme, line_);
    }
    else
      return Token (Token::token, lexeme, line_);
  }
}
