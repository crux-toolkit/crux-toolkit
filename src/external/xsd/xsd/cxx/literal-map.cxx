// file      : xsd/cxx/literal-map.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <memory>  // std::auto_ptr
#include <cstddef> // std::size_t
#include <fstream>
#include <iostream>

#include <xercesc/util/XMLUni.hpp>
#include <xercesc/util/XercesVersion.hpp>

#include <xercesc/framework/LocalFileInputSource.hpp>

#include <xercesc/sax/Locator.hpp>
#include <xercesc/sax/SAXParseException.hpp>
#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>

#include <xsd-frontend/xml.hxx>

#include <cxx/literal-map.hxx>

using namespace std;
using namespace xercesc;
namespace XML = XSDFrontend::XML;

namespace CXX
{
  class Handler: public DefaultHandler
  {
  public:
    struct Failed {};

    Handler (String const& file, StringLiteralMap& map)
        : state_ (s_init), file_ (file), map_ (map)
    {
    }

    virtual void
    setDocumentLocator (const Locator* const l)
    {
      locator_ = l;
    }

    virtual Void
    startElement (const XMLCh* const,
                  const XMLCh* const lname,
                  const XMLCh* const,
                  const xercesc::Attributes&)
    {
      String n (XML::transcode (lname));

      if (n == L"string-literal-map" && state_ == s_init)
        state_ = s_map;
      else if (n == L"entry" && state_ == s_map)
      {
        str_seen_ = false;
        lit_seen_ = false;
        state_ = s_entry;
      }
      else if (n == L"string" && state_ == s_entry)
      {
        str_seen_ = true;
        str_.clear ();
        state_ = s_string;
      }
      else if (n == L"literal" && state_ == s_entry)
      {
        lit_seen_ = true;
        lit_.clear ();
        state_ = s_literal;
      }
      else
      {
        wcerr << file_ << ":" << line () << ":" << col () << ": error: "
              << "unexpected element '" << n << "'" << endl;
        throw Failed ();
      }
    }

    virtual Void
    endElement (const XMLCh* const,
                const XMLCh* const lname,
                const XMLCh* const)
    {
      String n (XML::transcode (lname));

      if (n == L"string-literal-map")
        state_ = s_init;
      else if (n == L"entry")
      {
        if (!str_seen_)
        {
          wcerr << file_ << ":" << line () << ":" << col () << ": error: "
                << "expected 'string' element" << endl;
          throw Failed ();
        }

        if (!lit_seen_)
        {
          wcerr << file_ << ":" << line () << ":" << col () << ": error: "
                << "expected 'literal' element" << endl;
          throw Failed ();
        }

        map_[str_] = lit_;
        state_ = s_map;
      }
      else if (n == L"string")
        state_ = s_entry;
      else if (n == L"literal")
        state_ = s_entry;
    }

#if _XERCES_VERSION >= 30000
    virtual Void
    characters (const XMLCh* const s, const XMLSize_t length)
#else
    virtual Void
    characters (const XMLCh* const s, const unsigned int length)
#endif
    {
      String str (XML::transcode (s, length));

      if (state_ == s_string)
        str_ += str;
      else if (state_ == s_literal)
        lit_ += str;
      else
      {
        for (Size i (0); i < str.size (); ++i)
        {
          WideChar c (str[i]);

          if (c != 0x20 && c != 0x0A && c != 0x0D && c != 0x09)
          {
            wcerr << file_ << ":" << line () << ":" << col () << ": error: "
                  << "unexpected character data" << endl;
            throw Failed ();
          }
        }
      }
    }

    // Error hanlding.
    //
    enum Severity {s_warning, s_error, s_fatal};

    virtual Void
    warning (const SAXParseException& e)
    {
      handle (e, s_warning);
    }

    virtual Void
    error (const SAXParseException& e)
    {
      handle (e, s_error);
    }

    virtual Void
    fatalError (const SAXParseException& e)
    {
      handle (e, s_fatal);
    }

    virtual Void
    resetErrors ()
    {
    }

    Void
    handle (const SAXParseException& e, Severity s)
    {
      wcerr << file_ << ":";

#if _XERCES_VERSION >= 30000
      wcerr << e.getLineNumber () << ":" << e.getColumnNumber () << ": ";
#else
      XMLSSize_t l (e.getLineNumber ());
      XMLSSize_t c (e.getColumnNumber ());
      wcerr << (l == -1 ? 0 : l) << ":" << (c == -1 ? 0 : c) << ": ";
#endif

      String msg (XML::transcode (e.getMessage ()));
      wcerr << (s == s_warning ? "warning: " : "error: ") << msg << endl;

      if (s != s_warning)
        throw Failed ();
    }

    size_t
    line () const
    {
      size_t r (0);

      if (locator_ != 0)
      {
#if _XERCES_VERSION >= 30000
        r = static_cast<size_t> (locator_->getLineNumber ());
#else
        XMLSSize_t l (locator_->getLineNumber ());
        r = l == -1 ? 0 : static_cast<size_t> (l);
#endif
      }

      return r;
    }

    size_t
    col () const
    {
      size_t r (0);

      if (locator_ != 0)
      {
#if _XERCES_VERSION >= 30000
        r = static_cast<size_t> (locator_->getColumnNumber ());
#else
        XMLSSize_t c (locator_->getColumnNumber ());
        r = c == -1 ? 0 : static_cast<size_t> (c);
#endif
      }

      return r;
    }

  private:
    const Locator* locator_;

    enum
    {
      s_init,
      s_map,
      s_entry,
      s_string,
      s_literal
    } state_;

    String file_;
    StringLiteralMap& map_;

    Boolean str_seen_;
    Boolean lit_seen_;

    String str_;
    String lit_;
  };

  bool
  read_literal_map (NarrowString const& file, StringLiteralMap& map)
  {
    try
    {
      // Try to open the file with fstream. This way we get to
      // report the error in a consistent manner.
      //
      {
        ifstream ifs (file.c_str ());
        if (!ifs.is_open ())
        {
          wcerr << file.c_str () << ": unable to open in read mode" << endl;
          return false;
        }
      }

      String wfile (file);

      LocalFileInputSource is (XML::XMLChString (wfile).c_str ());
      Handler h (wfile, map);

      auto_ptr<SAX2XMLReader> parser (
        XMLReaderFactory::createXMLReader ());

      parser->setFeature (XMLUni::fgSAX2CoreNameSpaces, true);
      parser->setFeature (XMLUni::fgSAX2CoreNameSpacePrefixes, true);
      parser->setFeature (XMLUni::fgSAX2CoreValidation, false);
      parser->setFeature (XMLUni::fgXercesSchema, false);
      parser->setFeature (XMLUni::fgXercesSchemaFullChecking, false);

      parser->setErrorHandler (&h);
      parser->setContentHandler (&h);

      parser->parse (is);
    }
    catch (Handler::Failed const&)
    {
      return false;
    }

    return true;
  }
}
