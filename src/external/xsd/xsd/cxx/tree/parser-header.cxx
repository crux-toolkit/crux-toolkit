// file      : xsd/cxx/tree/parser-header.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/tree/parser-header.hxx>

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

namespace CXX
{
  namespace Tree
  {
    namespace
    {
      struct ElementFunction: Traversal::Element,
                              GlobalElementBase,
                              Context
      {
        ElementFunction (Context& c)
            : GlobalElementBase (c), Context (c)
        {
        }

        virtual Void
        traverse (Type& e)
        {
          if (!doc_root_p (e))
            return;

          String const& name (eparser (e));
          String const& error_handler (error_handler_type);

          if (doxygen)
          {
            os << "/**" << endl
               << " * @name Parsing functions for the %" <<
	       comment (e.name ()) << " document root." << endl;

            if (e.annotated_p ())
            {
              os << " *" << endl;
              write_annotation (e.annotation ());
            }

            os << " */" << endl
               << "//@{" << endl
               << endl;
          }

          if (!doxygen)
          {
            os << "// Parse a URI or a local file." << endl
               << "//" << endl
               << endl;
          }

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Parse a URI or a local file." << endl
               << " *" << endl
               << " * @param uri A URI or a local file name." << endl
               << " * @param f Parsing flags." << endl
               << " * @param p Parsing properties. " << endl
               << " * @return A pointer to the root of the object model." << endl
               << " *" << endl
               << " * This function uses exceptions to report parsing errors." << endl
               << " */" << endl;
          }

          os << inst_exp
             << "::std::auto_ptr< " << type_name (e) << " >" << endl
             << name << " (const " << string_type << "& uri," << endl
             << flags_type << " f = 0," << endl
             << "const " << properties_type << "& p = " << properties_type << " ());"
             << endl;

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Parse a URI or a local file with an error handler." << endl
               << " *" << endl
               << " * @param uri A URI or a local file name." << endl
               << " * @param eh An error handler." << endl
               << " * @param f Parsing flags." << endl
               << " * @param p Parsing properties. " << endl
               << " * @return A pointer to the root of the object model." << endl
               << " *" << endl
               << " * This function reports parsing errors by calling the " <<
              "error handler." << endl
               << " */" << endl;
          }

          os << inst_exp
             << "::std::auto_ptr< " << type_name (e) << " >" << endl
             << name << " (const " << string_type << "& uri," << endl
             << error_handler << "& eh," << endl
             << flags_type << " f = 0," << endl
             << "const " << properties_type << "& p = " << properties_type << " ());"
             << endl;

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Parse a URI or a local file with a Xerces-C++ " <<
              "DOM error" << endl
               << " * handler." << endl
               << " *" << endl
               << " * @param uri A URI or a local file name." << endl
               << " * @param eh A Xerces-C++ DOM error handler." << endl
               << " * @param f Parsing flags." << endl
               << " * @param p Parsing properties. " << endl
               << " * @return A pointer to the root of the object model." << endl
               << " *" << endl
               << " * This function reports parsing errors by calling the " <<
              "error handler." << endl
               << " */" << endl;
          }

          os << inst_exp
             << "::std::auto_ptr< " << type_name (e) << " >" << endl
             << name << " (const " << string_type << "& uri," << endl
             << xerces_ns << "::DOMErrorHandler& eh," << endl
             << flags_type << " f = 0," << endl
             << "const " << properties_type << "& p = " << properties_type << " ());"
             << endl;

          if (!doxygen)
          {
            os << "// Parse std::istream." << endl
               << "//" << endl
               << endl;
          }

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Parse a standard input stream." << endl
               << " *" << endl
               << " * @param is A standrad input stream." << endl
               << " * @param f Parsing flags." << endl
               << " * @param p Parsing properties. " << endl
               << " * @return A pointer to the root of the object model." << endl
               << " *" << endl
               << " * This function uses exceptions to report parsing errors." << endl
               << " */" << endl;
          }

          os << inst_exp
             << "::std::auto_ptr< " << type_name (e) << " >" << endl
             << name << " (::std::istream& is," << endl
             << flags_type << " f = 0," << endl
             << "const " << properties_type << "& p = " << properties_type << " ());"
             << endl;

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Parse a standard input stream with an error handler." << endl
               << " *" << endl
               << " * @param is A standrad input stream." << endl
               << " * @param eh An error handler." << endl
               << " * @param f Parsing flags." << endl
               << " * @param p Parsing properties. " << endl
               << " * @return A pointer to the root of the object model." << endl
               << " *" << endl
               << " * This function reports parsing errors by calling the " <<
              "error handler." << endl
               << " */" << endl;
          }

          os << inst_exp
             << "::std::auto_ptr< " << type_name (e) << " >" << endl
             << name << " (::std::istream& is," << endl
             << error_handler << "& eh," << endl
             << flags_type << " f = 0," << endl
             << "const " << properties_type << "& p = " << properties_type << " ());"
             << endl;

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Parse a standard input stream with a " <<
              "Xerces-C++ DOM error" << endl
               << " * handler." << endl
               << " *" << endl
               << " * @param is A standrad input stream." << endl
               << " * @param eh A Xerces-C++ DOM error handler." << endl
               << " * @param f Parsing flags." << endl
               << " * @param p Parsing properties. " << endl
               << " * @return A pointer to the root of the object model." << endl
               << " *" << endl
               << " * This function reports parsing errors by calling the " <<
              "error handler." << endl
               << " */" << endl;
          }

          os << inst_exp
             << "::std::auto_ptr< " << type_name (e) << " >" << endl
             << name << " (::std::istream& is," << endl
             << xerces_ns << "::DOMErrorHandler& eh," << endl
             << flags_type << " f = 0," << endl
             << "const " << properties_type << "& p = " << properties_type << " ());"
             << endl;

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Parse a standard input stream with a resource id." << endl
               << " *" << endl
               << " * @param is A standrad input stream." << endl
               << " * @param id A resource id." << endl
               << " * @param f Parsing flags." << endl
               << " * @param p Parsing properties. " << endl
               << " * @return A pointer to the root of the object model." << endl
               << " *" << endl
               << " * The resource id is used to identify the document " <<
              "being parsed in" << endl
               << " * diagnostics as well as to resolve relative paths." << endl
               << " *" << endl
               << " * This function uses exceptions to report parsing errors." << endl
               << " */" << endl;
          }

          os << inst_exp
             << "::std::auto_ptr< " << type_name (e) << " >" << endl
             << name << " (::std::istream& is," << endl
             << "const " << string_type << "& id," << endl
             << flags_type << " f = 0," << endl
             << "const " << properties_type << "& p = " << properties_type << " ());"
             << endl;

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Parse a standard input stream with a resource " <<
              "id and an" << endl
               << " * error handler." << endl
               << " *" << endl
               << " * @param is A standrad input stream." << endl
               << " * @param id A resource id." << endl
               << " * @param eh An error handler." << endl
               << " * @param f Parsing flags." << endl
               << " * @param p Parsing properties. " << endl
               << " * @return A pointer to the root of the object model." << endl
               << " *" << endl
               << " * The resource id is used to identify the document " <<
              "being parsed in" << endl
               << " * diagnostics as well as to resolve relative paths." << endl
               << " *" << endl
               << " * This function reports parsing errors by calling the " <<
              "error handler." << endl
               << " */" << endl;
          }

          os << inst_exp
             << "::std::auto_ptr< " << type_name (e) << " >" << endl
             << name << " (::std::istream& is," << endl
             << "const " << string_type << "& id," << endl
             << error_handler << "& eh," << endl
             << flags_type << " f = 0," << endl
             << "const " << properties_type << "& p = " << properties_type << " ());"
             << endl;

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Parse a standard input stream with a resource " <<
              "id and a" << endl
               << " * Xerces-C++ DOM error handler." << endl
               << " *" << endl
               << " * @param is A standrad input stream." << endl
               << " * @param id A resource id." << endl
               << " * @param eh A Xerces-C++ DOM error handler." << endl
               << " * @param f Parsing flags." << endl
               << " * @param p Parsing properties. " << endl
               << " * @return A pointer to the root of the object model." << endl
               << " *" << endl
               << " * The resource id is used to identify the document " <<
              "being parsed in" << endl
               << " * diagnostics as well as to resolve relative paths." << endl
               << " *" << endl
               << " * This function reports parsing errors by calling the " <<
              "error handler." << endl
               << " */" << endl;
          }

          os << inst_exp
             << "::std::auto_ptr< " << type_name (e) << " >" << endl
             << name << " (::std::istream& is," << endl
             << "const " << string_type << "& id," << endl
             << xerces_ns << "::DOMErrorHandler& eh," << endl
             << flags_type << " f = 0," << endl
             << "const " << properties_type << "& p = " << properties_type << " ());"
             << endl;

          if (!doxygen)
          {
            os << "// Parse xercesc::InputSource." << endl
               << "//" << endl
               << endl;
          }

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Parse a Xerces-C++ input source." << endl
               << " *" << endl
               << " * @param is A Xerces-C++ input source." << endl
               << " * @param f Parsing flags." << endl
               << " * @param p Parsing properties. " << endl
               << " * @return A pointer to the root of the object model." << endl
               << " *" << endl
               << " * This function uses exceptions to report parsing errors." << endl
               << " */" << endl;
          }

          os << inst_exp
             << "::std::auto_ptr< " << type_name (e) << " >" << endl
             << name << " (" << xerces_ns << "::InputSource& is," << endl
             << flags_type << " f = 0," << endl
             << "const " << properties_type << "& p = " << properties_type << " ());"
             << endl;

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Parse a Xerces-C++ input source with an " <<
              "error handler." << endl
               << " *" << endl
               << " * @param is A Xerces-C++ input source." << endl
               << " * @param eh An error handler." << endl
               << " * @param f Parsing flags." << endl
               << " * @param p Parsing properties. " << endl
               << " * @return A pointer to the root of the object model." << endl
               << " *" << endl
               << " * This function reports parsing errors by calling the " <<
              "error handler." << endl
               << " */" << endl;
          }

          os << inst_exp
             << "::std::auto_ptr< " << type_name (e) << " >" << endl
             << name << " (" << xerces_ns << "::InputSource& is," << endl
             << error_handler << "& eh," << endl
             << flags_type << " f = 0," << endl
             << "const " << properties_type << "& p = " << properties_type << " ());"
             << endl;

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Parse a Xerces-C++ input source with a " <<
              "Xerces-C++ DOM" << endl
               << " * error handler." << endl
               << " *" << endl
               << " * @param is A Xerces-C++ input source." << endl
               << " * @param eh A Xerces-C++ DOM error handler." << endl
               << " * @param f Parsing flags." << endl
               << " * @param p Parsing properties. " << endl
               << " * @return A pointer to the root of the object model." << endl
               << " *" << endl
               << " * This function reports parsing errors by calling the " <<
              "error handler." << endl
               << " */" << endl;
          }

          os << inst_exp
             << "::std::auto_ptr< " << type_name (e) << " >" << endl
             << name << " (" << xerces_ns << "::InputSource& is," << endl
             << xerces_ns << "::DOMErrorHandler& eh," << endl
             << flags_type << " f = 0," << endl
             << "const " << properties_type << "& p = " << properties_type << " ());"
             << endl;

          if (!doxygen)
          {
            os << "// Parse xercesc::DOMDocument." << endl
               << "//" << endl
               << endl;
          }

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Parse a Xerces-C++ DOM document." << endl
               << " *" << endl
               << " * @param d A Xerces-C++ DOM document." << endl
               << " * @param f Parsing flags." << endl
               << " * @param p Parsing properties. " << endl
               << " * @return A pointer to the root of the object model." << endl
               << " */" << endl;
          }

          os << inst_exp
             << "::std::auto_ptr< " << type_name (e) << " >" << endl
             << name << " (const " << xerces_ns << "::DOMDocument& d," << endl
             << flags_type << " f = 0," << endl
             << "const " << properties_type << "& p = " << properties_type << " ());"
             << endl;

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Parse a Xerces-C++ DOM document." << endl
               << " *" << endl
               << " * @param d A pointer to the Xerces-C++ DOM document." << endl
               << " * @param f Parsing flags." << endl
               << " * @param p Parsing properties. " << endl
               << " * @return A pointer to the root of the object model." << endl
               << " *" << endl
               << " * This function is normally used together with the " <<
              "keep_dom and" << endl
               << " * own_dom parsing flags to assign ownership of the DOM " <<
              "document" << endl
               << " * to the object model." << endl
               << " */" << endl;
          }

          os << inst_exp
             << "::std::auto_ptr< " << type_name (e) << " >" << endl
             << name << " (" << dom_auto_ptr << "< " << xerces_ns <<
            "::DOMDocument >& d," << endl
             << flags_type << " f = 0," << endl
             << "const " << properties_type << "& p = " << properties_type << " ());"
             << endl;

          if (doxygen)
          {
            os << "//@}" << endl
               << endl;
          }
        }

      private:
        String
        type_name (Type& e)
        {
          std::wostringstream o;

          MemberTypeName type (*this, o);
          type.dispatch (e.type ());

          return o.str ();
        }
      };
    }


    Void
    generate_parser_header (Context& ctx)
    {
      ctx.os << "#include <iosfwd>" << endl
             << endl
             << "#include <xercesc/sax/InputSource.hpp>" << endl
               << "#include <xercesc/dom/DOMDocument.hpp>" << endl
             << "#include <xercesc/dom/DOMErrorHandler.hpp>" << endl
             << endl;

      Traversal::Schema schema;
      Traversal::Sources sources;
      Traversal::Names names_ns, names;
      Namespace ns (ctx);
      ElementFunction element (ctx);

      schema >> sources >> schema;
      schema >> names_ns >> ns >> names >> element;

      schema.dispatch (ctx.schema_root);
    }
  }
}
