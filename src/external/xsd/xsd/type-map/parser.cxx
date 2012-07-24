// file      : xsd/type-map/parser.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2007-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <iostream>

#include <backend-elements/regex.hxx>

#include <type-map/parser.hxx>

using std::endl;

namespace TypeMap
{
  typedef Lexer::Token Token;
  typedef BackendElements::Regex::Format<WideChar> Format;

  Parser::Parser (Lexer& lex, String const& path)
      : lex_ (lex), path_ (path), e (std::wcerr)
  {
  }

  Boolean Parser::
  parse (Namespaces& ns)
  {
    try
    {
      Namespace* global = 0;

      for (Token t (lex_.next ()); t.type () != Token::eos; t = lex_.next ())
      {
        String l (t.lexeme ());

        if (l == L"namespace")
        {
          global = 0;

          if (!namespace_ (ns))
            return false;
        }
        else if (l == L"include")
        {
          if (global == 0)
          {
            ns.push_back (Namespace (Pattern ()));
            global = &(*ns.rbegin ());
          }

          if (!include (*global))
            return false;
        }
        else if (l == L"type" || t.type () == Token::token)
        {
          // Type mapping can have 'type' specifier omitted.
          //
          if (l == L"type")
            t = lex_.next ();

          if (global == 0)
          {
            ns.push_back (Namespace (Pattern ()));
            global = &(*ns.rbegin ());
          }

          if (!type (t, *global))
            return false;
        }
        else
        {
          e << path_ << ":" << t.line () << ": unexpected '" << l << "'"
            << endl;

          return false;
        }
      }
    }
    catch (Lexer::Failed const&)
    {
      return false;
    }

    return true;
  }

  Boolean Parser::
  namespace_ (Namespaces& ns)
  {
    // First get XML namespace.
    //
    Token t (lex_.next ());

    Pattern xsd_name;

    try
    {
      xsd_name = t.lexeme ();
    }
    catch (Format const& ex)
    {
      e << path_ << ":" << t.line () << ": invalid namespace pattern: "
        << ex.description () << endl;
      return false;
    }

    if (t.type () != Token::token)
    {
      e << path_ << ":" << t.line () << ": expected XML namespace "
        << "instead of '" << xsd_name << "'" << endl;
      return false;
    }


    // See if we've got optional C++ mapping.
    //
    t = lex_.next ();

    Boolean has_cxx_name (false);
    String cxx_name;

    if (t.type () != Token::token)
    {
      if (t.lexeme () != L"{")
      {
        e << path_ << ":" << t.line () << ": expected C++ namespace or '{' "
          << "instead of '" << t.lexeme () << "'" << endl;
        return false;
      }
    }
    else
    {
      has_cxx_name = true;
      cxx_name = t.lexeme ();
    }

    // Swallow '{' if needed.
    //
    if (has_cxx_name)
    {
      t = lex_.next ();

      if (t.type () != Token::punct || t.lexeme () != L"{")
      {
        e << path_ << ":" << t.line () << ": expected '{' instead of '"
          << t.lexeme () << "'" << endl;
        return false;
      }
    }

    Namespace n (xsd_name, has_cxx_name, cxx_name);

    // Parse namespace body.
    //
    for (t = lex_.next ();; t = lex_.next ())
    {
      String l (t.lexeme ());

      if (l == L"include")
      {
        if (!include (n))
          return false;
      }
      else if (l == L"type" || t.type () == Token::token)
      {
        // Type mapping can have 'type' specifier omitted.
        //
        if (l == L"type")
          t = lex_.next ();

        if (!type (t, n))
          return false;
      }
      else if (t.type () == Token::punct && l == L"}")
      {
        break;
      }
      else
      {
        e << path_ << ":" << t.line () << ": unexpected '" << l << "'"
          << endl;
        return false;
      }
    }

    if (cxx_name || n.types_begin () != n.types_end () ||
        n.includes_begin () != n.includes_end ())
    {
      ns.push_back (n);
    }

    return true;
  }

  Boolean Parser::
  include (Namespace& n)
  {
    Token t (lex_.next ());

    String path (t.lexeme ());

    if (t.type () != Token::token)
    {
      e << path_ << ":" << t.line () << ": expected include path "
        << "instead of '" << path << "'" << endl;
      return false;
    }

    if (path && path[0] == L'<')
      n.includes_push_back (path);
    else
      n.includes_push_back (L'"' + path + L'"');

    t = lex_.next ();

    if (t.type () != Token::punct || t.lexeme () != L";")
    {
      e << path_ << ":" << t.line () << ": expected ';' after '"
        << path << "'" << endl;
      return false;
    }

    return true;
  }

  Boolean Parser::
  type (Token t, Namespace& n)
  {
    Pattern xsd_name;

    try
    {
      xsd_name = t.lexeme ();
    }
    catch (Format const& ex)
    {
      e << path_ << ":" << t.line () << ": invalid namespace pattern: "
        << ex.description () << endl;
      return false;
    }

    if (t.type () != Token::token)
    {
      e << path_ << ":" << t.line () << ": expected XML Schema type name "
        << "instead of '" << xsd_name << "'" << endl;
      return false;
    }

    t = lex_.next ();
    String cxx_ret_name (t.lexeme ());

    if (t.type () != Token::token)
    {
      e << path_ << ":" << t.line () << ": expected C++ type name "
        << "instead of '" << cxx_ret_name << "'" << endl;
      return false;
    }

    t = lex_.next ();

    String cxx_arg_name;

    // See if we've got optional argument type.
    //
    if (t.type () == Token::token)
    {
      cxx_arg_name = t.lexeme ();
      t = lex_.next ();
    }

    if (t.type () != Token::punct || t.lexeme () != L";")
    {
      e << path_ << ":" << t.line () << ": expected ';' after '"
        << cxx_arg_name << "'" << endl;
      return false;
    }

    n.types_push_back (xsd_name, cxx_ret_name, cxx_arg_name);

    return true;
  }
}
