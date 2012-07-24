// file      : xsd/cxx/tree/serialization-header.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/tree/serialization-header.hxx>

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

namespace CXX
{
  namespace Tree
  {
    namespace
    {
      struct List: Traversal::List, Context
      {
        List (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& l)
        {
          String name (ename (l));

          // If renamed name is empty then we do not need to generate
          // anything for this type.
          //
          if (renamed_type (l, name) && !name)
            return;

          // operator<< (xercesc::DOMElement)
          //
          os << inst_exp
             << "void" << endl
             << "operator<< (" << xerces_ns << "::DOMElement&, " <<
            "const " << name << "&);"
             << endl;

          // operator<< (xercesc::DOMAttr)
          //
          os << inst_exp
             << "void" << endl
             << "operator<< (" << xerces_ns << "::DOMAttr&, " <<
            "const " << name << "&);"
             << endl;

          // operator<< (list_stream)
          //
          os << inst_exp
             << "void" << endl
             << "operator<< (" << list_stream_type << "&," << endl
             << "const " << name << "&);"
             << endl;

        }
      };


      struct Union: Traversal::Union, Context
      {
        Union (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& u)
        {
          String name (ename (u));

          // If renamed name is empty then we do not need to generate
          // anything for this type.
          //
          if (renamed_type (u, name) && !name)
            return;

          // operator<< (xercesc::DOMElement)
          //
          os << inst_exp
             << "void" << endl
             << "operator<< (" << xerces_ns << "::DOMElement&, " <<
            "const " << name << "&);"
             << endl;

          // operator<< (xercesc::DOMAttr)
          //
          os << inst_exp
             << "void" << endl
             << "operator<< (" << xerces_ns << "::DOMAttr&, " <<
            "const " << name << "&);"
             << endl;

          // operator<< (list_stream)
          //
          os << inst_exp
             << "void" << endl
             << "operator<< (" << list_stream_type << "&," << endl
             << "const " << name << "&);"
             << endl;
        }
      };


      struct Enumeration: Traversal::Enumeration, Context
      {
        Enumeration (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& e)
        {
          String name (ename (e));

          // If renamed name is empty then we do not need to generate
          // anything for this type.
          //
          if (renamed_type (e, name) && !name)
            return;

          // operator<< (xercesc::DOMElement)
          //
          os << inst_exp
             << "void" << endl
             << "operator<< (" << xerces_ns << "::DOMElement&, " <<
            "const " << name << "&);"
             << endl;

          // operator<< (xercesc::DOMAttr)
          //
          os << inst_exp
             << "void" << endl
             << "operator<< (" << xerces_ns << "::DOMAttr&, " <<
            "const " << name << "&);"
             << endl;

          // operator<< (list_stream)
          //
          os << inst_exp
             << "void" << endl
             << "operator<< (" << list_stream_type << "&," << endl
             << "const " << name << "&);"
             << endl;
        }
      };

      struct Complex: Traversal::Complex, Context
      {
        Complex (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& c)
        {
          String name (ename (c));

          // If renamed name is empty then we do not need to generate
          // anything for this type.
          //
          if (renamed_type (c, name) && !name)
            return;

          // operator<< (xercesc::DOMElement)
          //
          os << inst_exp
             << "void" << endl
             << "operator<< (" << xerces_ns << "::DOMElement&, " <<
            "const " << name << "&);"
             << endl;

          Boolean simple (true);
          {
            IsSimpleType t (simple);
            t.dispatch (c);
          }

          if (simple)
          {
            // operator<< (xercesc::DOMAttr)
            //
            os << inst_exp
               << "void" << endl
               << "operator<< (" << xerces_ns << "::DOMAttr&, " <<
              "const " << name << "&);"
               << endl;

            // operator<< (list_stream)
            //
            os << inst_exp
               << "void" << endl
               << "operator<< (" << list_stream_type << "&," << endl
               << "const " << name << "&);"
               << endl;
          }
        }
      };

      struct ElementType: Traversal::Element,
                          GlobalElementBase,
                          Context
      {
        ElementType (Context& c)
            : GlobalElementBase (c), Context (c)
        {
        }

        virtual Void
        traverse (Type& e)
        {
          if (doc_root_p (e))
          {
            // operator<< (xercesc::DOMElement)
            //
            os << inst_exp
               << "void" << endl
               << "operator<< (" << xerces_ns << "::DOMElement&, " <<
              "const " << ename (e) << "&);"
               << endl;
          }
        }
      };

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

          String const& name (eserializer (e));
          String const& error_handler (error_handler_type);
          String const& namespace_infomap (namespace_infomap_type);

          if (doxygen)
          {
            os << "/**" << endl
               << " * @name Serialization functions for the %" <<
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
            os << "// Serialize to std::ostream." << endl
               << "//" << endl
               << endl;
          }

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Serialize to a standard output stream." << endl
               << " *" << endl
               << " * @param os A standrad output stream." << endl
               << " * @param x An object model to serialize." << endl
               << " * @param m A namespace information map." << endl
               << " * @param e A character encoding to produce XML in." << endl
               << " * @param f Serialization flags." << endl
               << " *" << endl
               << " * This function uses exceptions to report " <<
              "serialization errors." << endl
               << " */" << endl;
          }

          os << inst_exp
             << "void" << endl
             << name << " (::std::ostream& os," << endl
             << "const " << type_name (e) << "& x, " << endl
             << "const " << namespace_infomap << "& m = " <<
            namespace_infomap << " ()," << endl
             << "const " << string_type << "& e = " << L << "\"UTF-8\"," << endl
             << flags_type << " f = 0);"
             << endl;

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Serialize to a standard output stream with an " <<
              "error handler." << endl
               << " *" << endl
               << " * @param os A standrad output stream." << endl
               << " * @param x An object model to serialize." << endl
               << " * @param eh An error handler." << endl
               << " * @param m A namespace information map." << endl
               << " * @param e A character encoding to produce XML in." << endl
               << " * @param f Serialization flags." << endl
               << " *" << endl
               << " * This function reports serialization errors by " <<
              "calling the error" << endl
               << " * handler." << endl
               << " */" << endl;
          }

          os << inst_exp
             << "void" << endl
             << name << " (::std::ostream& os," << endl
             << "const " << type_name (e) << "& x, " << endl
             << error_handler << "& eh," << endl
             << "const " << namespace_infomap << "& m = "  <<
            namespace_infomap << " ()," << endl
             << "const " << string_type << "& e = " << L << "\"UTF-8\"," << endl
             << flags_type << " f = 0);"
             << endl;

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Serialize to a standard output stream with a " <<
              "Xerces-C++ DOM" << endl
               << " * error handler." << endl
               << " *" << endl
               << " * @param os A standrad output stream." << endl
               << " * @param x An object model to serialize." << endl
               << " * @param eh A Xerces-C++ DOM error handler." << endl
               << " * @param m A namespace information map." << endl
               << " * @param e A character encoding to produce XML in." << endl
               << " * @param f Serialization flags." << endl
               << " *" << endl
               << " * This function reports serialization errors by " <<
              "calling the error" << endl
               << " * handler." << endl
               << " */" << endl;
          }

          os << inst_exp
             << "void" << endl
             << name << " (::std::ostream& os," << endl
             << "const " << type_name (e) << "& x, " << endl
             << xerces_ns << "::DOMErrorHandler& eh," << endl
             << "const " << namespace_infomap << "& m = " <<
            namespace_infomap << " ()," << endl
             << "const " << string_type << "& e = " << L << "\"UTF-8\"," << endl
             << flags_type << " f = 0);"
             << endl;

          if (!doxygen)
          {
            os << "// Serialize to xercesc::XMLFormatTarget." << endl
               << "//" << endl
               << endl;
          }

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Serialize to a Xerces-C++ XML format target." << endl
               << " *" << endl
               << " * @param ft A Xerces-C++ XML format target." << endl
               << " * @param x An object model to serialize." << endl
               << " * @param m A namespace information map." << endl
               << " * @param e A character encoding to produce XML in." << endl
               << " * @param f Serialization flags." << endl
               << " *" << endl
               << " * This function uses exceptions to report " <<
              "serialization errors." << endl
               << " */" << endl;
          }

          os << inst_exp
             << "void" << endl
             << name << " (" << xerces_ns << "::XMLFormatTarget& ft," << endl
             << "const " << type_name (e) << "& x, " << endl
             << "const " << namespace_infomap << "& m = " <<
            namespace_infomap << " ()," << endl
             << "const " << string_type << "& e = " << L << "\"UTF-8\"," << endl
             << flags_type << " f = 0);"
             << endl;

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Serialize to a Xerces-C++ XML format target " <<
              "with an error" << endl
               << " * handler." << endl
               << " *" << endl
               << " * @param ft A Xerces-C++ XML format target." << endl
               << " * @param x An object model to serialize." << endl
               << " * @param eh An error handler." << endl
               << " * @param m A namespace information map." << endl
               << " * @param e A character encoding to produce XML in." << endl
               << " * @param f Serialization flags." << endl
               << " *" << endl
               << " * This function reports serialization errors by " <<
              "calling the error" << endl
               << " * handler." << endl
               << " */" << endl;
          }

          os << inst_exp
             << "void" << endl
             << name << " (" << xerces_ns << "::XMLFormatTarget& ft," << endl
             << "const " << type_name (e) << "& x, " << endl
             << error_handler << "& eh," << endl
             << "const " << namespace_infomap << "& m = " <<
            namespace_infomap << " ()," << endl
             << "const " << string_type << "& e = " << L << "\"UTF-8\"," << endl
             << flags_type << " f = 0);"
             << endl;

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Serialize to a Xerces-C++ XML format target " <<
              "with a" << endl
               << " * Xerces-C++ DOM error handler." << endl
               << " *" << endl
               << " * @param ft A Xerces-C++ XML format target." << endl
               << " * @param x An object model to serialize." << endl
               << " * @param eh A Xerces-C++ DOM error handler." << endl
               << " * @param m A namespace information map." << endl
               << " * @param e A character encoding to produce XML in." << endl
               << " * @param f Serialization flags." << endl
               << " *" << endl
               << " * This function reports serialization errors by " <<
              "calling the error" << endl
               << " * handler." << endl
               << " */" << endl;
          }

          os << inst_exp
             << "void" << endl
             << name << " (" << xerces_ns << "::XMLFormatTarget& ft," << endl
             << "const " << type_name (e) << "& x, " << endl
             << xerces_ns << "::DOMErrorHandler& eh," << endl
             << "const " << namespace_infomap << "& m = " <<
            namespace_infomap << " ()," << endl
             << "const " << string_type << "& e = " << L << "\"UTF-8\"," << endl
             << flags_type << " f = 0);"
             << endl;

          if (!doxygen)
          {
            os << "// Serialize to an existing xercesc::DOMDocument." << endl
               << "//" << endl
               << endl;
          }

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Serialize to an existing Xerces-C++ DOM document." << endl
               << " *" << endl
               << " * @param d A Xerces-C++ DOM document." << endl
               << " * @param x An object model to serialize." << endl
               << " * @param f Serialization flags." << endl
               << " *" << endl
               << " * Note that it is your responsibility to create the " <<
              "DOM document" << endl
               << " * with the correct root element as well as set the " <<
              "necessary" << endl
               << " * namespace mapping attributes." << endl
               << " */" << endl;
          }
          os << inst_exp
             << "void" << endl
             << name << " (" << xerces_ns << "::DOMDocument& d," << endl
             << "const " << type_name (e) << "& x," << endl
             << flags_type << " f = 0);"
             << endl;

          if (!doxygen)
          {
            os << "// Serialize to a new xercesc::DOMDocument." << endl
               << "//" << endl
               << endl;
          }

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Serialize to a new Xerces-C++ DOM document." << endl
               << " *" << endl
               << " * @param x An object model to serialize." << endl
               << " * @param m A namespace information map." << endl
               << " * @param f Serialization flags." << endl
               << " * @return A pointer to the new Xerces-C++ DOM document." << endl
               << " */" << endl;
          }

          os << inst_exp
             << dom_auto_ptr << "< " << xerces_ns <<
            "::DOMDocument >" << endl
             << name << " (const " << type_name (e) << "& x, " << endl
             << "const " << namespace_infomap << "& m = " <<
            namespace_infomap << " ()," << endl
             << flags_type << " f = 0);"
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
    generate_serialization_header (Context& ctx)
    {
      Boolean elemen_type (ctx.options.value<CLI::generate_element_type> ());

      if (!elemen_type)
        ctx.os << "#include <iosfwd>" << endl
               << endl
               << "#include <xercesc/dom/DOMDocument.hpp>" << endl
               << "#include <xercesc/dom/DOMErrorHandler.hpp>" << endl
               << "#include <xercesc/framework/XMLFormatter.hpp>" << endl
               << endl
               << "#include <xsd/cxx/xml/dom/auto-ptr.hxx>" << endl
               << endl;

      Traversal::Schema schema;

      Traversal::Sources sources;
      Traversal::Names names_ns, names;

      Namespace ns (ctx);

      List list (ctx);
      Union union_ (ctx);
      Complex complex (ctx);
      Enumeration enumeration (ctx);
      ElementType element_type (ctx);
      ElementFunction element_function (ctx);

      schema >> sources >> schema;
      schema >> names_ns >> ns >> names;

      names >> list;
      names >> union_;
      names >> complex;
      names >> enumeration;

      if (elemen_type)
        names >> element_type;
      else
        names >> element_function;

      schema.dispatch (ctx.schema_root);
    }
  }
}
