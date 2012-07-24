// file      : xsd/cxx/tree/parser-source.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/tree/parser-source.hxx>

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
          SemanticGraph::Type& t (e.type ());
          String type (type_name (e));
          String const& error_handler (error_handler_type);

          // Note that I am using fq-name in function calls because g++ gets
          // confused if the name is 'type'. (see tests/schema/anonymous)
          //

          // URI.
          //
          os << "::std::auto_ptr< " << type << " >" << endl
             << name << " (const " << string_type << "& u," << endl
             << flags_type << " f," << endl
             << "const " << properties_type << "& p)"
             << "{"
             << "::xsd::cxx::xml::auto_initializer i (" << endl
             << "(f & " << flags_type << "::dont_initialize) == 0," << endl
             << "(f & " << flags_type << "::keep_dom) == 0);"
             << endl
             << "::xsd::cxx::tree::error_handler< " << char_type << " > h;"
             << endl
             << dom_auto_ptr << "< " << xerces_ns <<
            "::DOMDocument > d (" << endl
             << "::xsd::cxx::xml::dom::parse< " << char_type << " > (" << endl
             << "u, h, p, f";

          if (options.value<CLI::disable_multi_import> ())
            os << " | ::xsd::cxx::xml::dom::no_muliple_imports";

            os << "));"
             << endl
             << "h.throw_if_failed< ::xsd::cxx::tree::parsing< " <<
            char_type << " > > ();"
             << endl
             << "::std::auto_ptr< " << type << " > r (" << endl
             << fq_name (e, "parser") << " (" << endl
             << "d, f | " << flags_type << "::own_dom, p));"
             << endl
             << "return r;"
             << "}";

          os << "::std::auto_ptr< " << type << " >" << endl
             << name << " (const " << string_type << "& u," << endl
             << error_handler << "& h," << endl
             << flags_type << " f," << endl
             << "const " << properties_type << "& p)"
             << "{"
             << "::xsd::cxx::xml::auto_initializer i (" << endl
             << "(f & " << flags_type << "::dont_initialize) == 0," << endl
             << "(f & " << flags_type << "::keep_dom) == 0);"
             << endl
             << dom_auto_ptr << "< " << xerces_ns <<
            "::DOMDocument > d (" << endl
             << "::xsd::cxx::xml::dom::parse< " << char_type << " > (" << endl
             << "u, h, p, f";

          if (options.value<CLI::disable_multi_import> ())
            os << " | ::xsd::cxx::xml::dom::no_muliple_imports";

          os << "));"
             << endl
             << "if (!d.get ())" << endl
             << "throw ::xsd::cxx::tree::parsing< " << char_type << " > ();"
             << endl
             << "::std::auto_ptr< " << type << " > r (" << endl
             << fq_name (e, "parser") << " (" << endl
             << "d, f | " << flags_type << "::own_dom, p));"
             << endl
             << "return r;"
             << "}";

          os << "::std::auto_ptr< " << type << " >" << endl
             << name << " (const " << string_type << "& u," << endl
             << xerces_ns << "::DOMErrorHandler& h," << endl
             << flags_type << " f," << endl
             << "const " << properties_type << "& p)"
             << "{"
             << dom_auto_ptr << "< " << xerces_ns <<
            "::DOMDocument > d (" << endl
             << "::xsd::cxx::xml::dom::parse< " << char_type << " > (" << endl
             << "u, h, p, f";

          if (options.value<CLI::disable_multi_import> ())
            os << " | ::xsd::cxx::xml::dom::no_muliple_imports";

          os << "));"
             << endl
             << "if (!d.get ())" << endl
             << "throw ::xsd::cxx::tree::parsing< " << char_type << " > ();"
             << endl
             << "::std::auto_ptr< " << type << " > r (" << endl
             << fq_name (e, "parser") << " (" << endl
             << "d, f | " << flags_type << "::own_dom, p));"
             << endl
             << "return r;"
             << "}";


          // istream
          //
          os << "::std::auto_ptr< " << type << " >" << endl
             << name << " (::std::istream& is," << endl
             << flags_type << " f," << endl
             << "const " << properties_type << "& p)"
             << "{"
             << "::xsd::cxx::xml::auto_initializer i (" << endl
             << "(f & " << flags_type << "::dont_initialize) == 0," << endl
             << "(f & " << flags_type << "::keep_dom) == 0);"
             << endl
             << "::xsd::cxx::xml::sax::std_input_source isrc (is);"
             << "return " << fq_name (e, "parser") << " (isrc, f, p);"
             << "}";

          os << "::std::auto_ptr< " << type << " >" << endl
             << name << " (::std::istream& is," << endl
             << error_handler << "& h," << endl
             << flags_type << " f," << endl
             << "const " << properties_type << "& p)"
             << "{"
             << "::xsd::cxx::xml::auto_initializer i (" << endl
             << "(f & " << flags_type << "::dont_initialize) == 0," << endl
             << "(f & " << flags_type << "::keep_dom) == 0);"
             << endl
             << "::xsd::cxx::xml::sax::std_input_source isrc (is);"
             << "return " << fq_name (e, "parser") << " (isrc, h, f, p);"
             << "}";

          os << "::std::auto_ptr< " << type << " >" << endl
             << name << " (::std::istream& is," << endl
             << xerces_ns << "::DOMErrorHandler& h," << endl
             << flags_type << " f," << endl
             << "const " << properties_type << "& p)"
             << "{"
             << "::xsd::cxx::xml::sax::std_input_source isrc (is);"
             << "return " << fq_name (e, "parser") << " (isrc, h, f, p);"
             << "}";

          os << "::std::auto_ptr< " << type << " >" << endl
             << name << " (::std::istream& is," << endl
             << "const " << string_type << "& sid," << endl
             << flags_type << " f," << endl
             << "const " << properties_type << "& p)"
             << "{"
             << "::xsd::cxx::xml::auto_initializer i (" << endl
             << "(f & " << flags_type << "::dont_initialize) == 0," << endl
             << "(f & " << flags_type << "::keep_dom) == 0);"
             << endl
             << "::xsd::cxx::xml::sax::std_input_source isrc (is, sid);"
             << "return " << fq_name (e, "parser") << " (isrc, f, p);"
             << "}";

          os << "::std::auto_ptr< " << type << " >" << endl
             << name << " (::std::istream& is," << endl
             << "const " << string_type << "& sid," << endl
             << error_handler << "& h," << endl
             << flags_type << " f," << endl
             << "const " << properties_type << "& p)"
             << "{"
             << "::xsd::cxx::xml::auto_initializer i (" << endl
             << "(f & " << flags_type << "::dont_initialize) == 0," << endl
             << "(f & " << flags_type << "::keep_dom) == 0);"
             << endl
             << "::xsd::cxx::xml::sax::std_input_source isrc (is, sid);"
             << "return " << fq_name (e, "parser") << " (isrc, h, f, p);"
             << "}";

          os << "::std::auto_ptr< " << type << " >" << endl
             << name << " (::std::istream& is," << endl
             << "const " << string_type << "& sid," << endl
             << xerces_ns << "::DOMErrorHandler& h," << endl
             << flags_type << " f," << endl
             << "const " << properties_type << "& p)"
             << "{"
             << "::xsd::cxx::xml::sax::std_input_source isrc (is, sid);"
             << "return " << fq_name (e, "parser") << " (isrc, h, f, p);"
             << "}";


          // InputSource.
          //
          os << "::std::auto_ptr< " << type << " >" << endl
             << name << " (" << xerces_ns << "::InputSource& i," << endl
             << flags_type << " f," << endl
             << "const " << properties_type << "& p)"
             << "{"
             << "::xsd::cxx::tree::error_handler< " << char_type << " > h;"
             << endl
             << dom_auto_ptr << "< " << xerces_ns <<
            "::DOMDocument > d (" << endl
             << "::xsd::cxx::xml::dom::parse< " << char_type << " > (" << endl
             << "i, h, p, f";

          if (options.value<CLI::disable_multi_import> ())
            os << " | ::xsd::cxx::xml::dom::no_muliple_imports";

          os << "));"
             << endl
             << "h.throw_if_failed< ::xsd::cxx::tree::parsing< " <<
            char_type << " > > ();"
             << endl
             << "::std::auto_ptr< " << type << " > r (" << endl
             << fq_name (e, "parser") << " (" << endl
             << "d, f | " << flags_type << "::own_dom, p));"
             << endl
             << "return r;"
             << "}";

          os << "::std::auto_ptr< " << type << " >" << endl
             << name << " (" << xerces_ns << "::InputSource& i," << endl
             << error_handler << "& h," << endl
             << flags_type << " f," << endl
             << "const " << properties_type << "& p)"
             << "{"
             << dom_auto_ptr << "< " << xerces_ns <<
            "::DOMDocument > d (" << endl
             << "::xsd::cxx::xml::dom::parse< " << char_type << " > (" << endl
             << "i, h, p, f";

          if (options.value<CLI::disable_multi_import> ())
            os << " | ::xsd::cxx::xml::dom::no_muliple_imports";

          os << "));"
             << endl
             << "if (!d.get ())" << endl
             << "throw ::xsd::cxx::tree::parsing< " << char_type << " > ();"
             << endl
             << "::std::auto_ptr< " << type << " > r (" << endl
             << fq_name (e, "parser") << " (" << endl
             << "d, f | " << flags_type << "::own_dom, p));"
             << endl
             << "return r;"
             << "}";


          os << "::std::auto_ptr< " << type << " >" << endl
             << name << " (" << xerces_ns << "::InputSource& i," << endl
             << xerces_ns << "::DOMErrorHandler& h," << endl
             << flags_type << " f," << endl
             << "const " << properties_type << "& p)"
             << "{"
             << dom_auto_ptr << "< " << xerces_ns <<
            "::DOMDocument > d (" << endl
             << "::xsd::cxx::xml::dom::parse< " << char_type << " > (" << endl
             << "i, h, p, f";

          if (options.value<CLI::disable_multi_import> ())
            os << " | ::xsd::cxx::xml::dom::no_muliple_imports";

          os << "));"
             << endl
             << "if (!d.get ())" << endl
             << "throw ::xsd::cxx::tree::parsing< " << char_type << " > ();"
             << endl
             << "::std::auto_ptr< " << type << " > r (" << endl
             << fq_name (e, "parser") << " (" << endl
             << "d, f | " << flags_type << "::own_dom, p));"
             << endl
             << "return r;"
             << "}";


          // DOM.
          //

          Boolean fund (false);
          {
            IsFundamentalType test (fund);
            test.dispatch (t);
          }

          // Check if we need to handle xsi:type and substitution groups.
          // If this element's type is anonymous then we don't need to do
          // anything.
          //
          Boolean poly (polymorphic && polymorphic_p (t) && !anonymous_p (t));

          // const DOMDocument&
          //
          os << "::std::auto_ptr< " << type << " >" << endl
             << name << " (const " << xerces_ns << "::DOMDocument& d," << endl
             << flags_type << " f," << endl
             << "const " << properties_type << "& p)"
             << "{"
             << "if (f & " << flags_type << "::keep_dom)"
             << "{"
             << dom_auto_ptr << "< " << xerces_ns <<
            "::DOMDocument > c (" << endl
             << "static_cast< " << xerces_ns <<
            "::DOMDocument* > (d.cloneNode (true)));"
             << endl
             << "::std::auto_ptr< " << type << " > r (" << endl
             << fq_name (e, "parser") << " (" << endl
             << "c, f | " << flags_type << "::own_dom, p));"
             << endl
             << "return r;"
             << "}"
             << "const " << xerces_ns << "::DOMElement& e (*d.getDocumentElement ());"
             << "const " << qname_type << " n (" << endl
             << "::xsd::cxx::xml::dom::name< " << char_type << " > (e));"
             << endl;

          if (poly)
          {
            // aCC cannot handle an inline call to type_factory_map_instance.
            //
            os << "::xsd::cxx::tree::type_factory_map< " << char_type <<
              " >& tfm (" << endl
               << "::xsd::cxx::tree::type_factory_map_instance< 0, " <<
              char_type << " > ());"
               << endl
               << "::std::auto_ptr< ::xsd::cxx::tree::type > tmp (" << endl
               << "tfm.create (" << endl
               << strlit (e.name ()) << "," << endl
               << strlit (e.namespace_().name ()) << "," << endl
               << "&::xsd::cxx::tree::factory_impl< " << type << " >," << endl
               << "true, true, e, n, f, 0));"
               << endl
               << "if (tmp.get () != 0)"
               << "{"
               << "::std::auto_ptr< " << type << " > r (" << endl
               << "dynamic_cast< " << type << "* > (tmp.get ()));"
               << endl
               << "if (r.get ())" << endl
               << "tmp.release ();"
               << "else" << endl
               << "throw ::xsd::cxx::tree::not_derived< " << char_type <<
              " > ();"
               << endl;
          }
          else
          {
            os << "if (n.name () == " << strlit (e.name ()) << " &&" << endl
               << "n.namespace_ () == " << strlit (e.namespace_().name ()) << ")"
               << "{";

            if (fund)
            {
              os << "::std::auto_ptr< " << type << " > r (" << endl
                 << "new " << type << " (" << endl
                 << "::xsd::cxx::tree::traits< " << type << ", " <<
                char_type;

              if (t.is_a<SemanticGraph::Fundamental::Double> ())
                os << ", ::xsd::cxx::tree::schema_type::double_";
              else if (t.is_a<SemanticGraph::Fundamental::Decimal> ())
                os << ", ::xsd::cxx::tree::schema_type::decimal";

              os << " >::create (" << endl
                 << "e, f, 0)));";
            }
            else
            {
              os << "::std::auto_ptr< " << type << " > r (" << endl
                 << "::xsd::cxx::tree::traits< " << type << ", " <<
                char_type << " >::create (" << endl
                 << "e, f, 0));";
            }
          }

          os << "return r;"
             << "}"
             << "throw ::xsd::cxx::tree::unexpected_element < " <<
            char_type << " > (" << endl
             << "n.name ()," << endl
             << "n.namespace_ ()," << endl
             << strlit (e.name ()) << "," << endl
             << strlit (e.namespace_().name ()) << ");"
             << "}";


          // dom::auto_ptr<DOMDocument>
          //
          os << "::std::auto_ptr< " << type << " >" << endl
             << name << " (" << dom_auto_ptr << "< " << xerces_ns <<
            "::DOMDocument >& d," << endl
             << flags_type << " f," << endl
             << "const " << properties_type << "&)"
             << "{"
             << dom_auto_ptr << "< " << xerces_ns <<
            "::DOMDocument > c (" << endl
             << "((f & " << flags_type << "::keep_dom) &&" << endl
             << "!(f & " << flags_type << "::own_dom))" << endl
             << "? static_cast< " << xerces_ns << "::DOMDocument* > (" <<
            "d->cloneNode (true))" << endl
             << ": 0);"
             << endl
             << xerces_ns << "::DOMDocument& doc (c.get () ? *c : *d);"
             << "const " << xerces_ns << "::DOMElement& e (" <<
            "*doc.getDocumentElement ());"
             << endl
             << "const " << qname_type << " n (" << endl
             << "::xsd::cxx::xml::dom::name< " << char_type << " > (e));"
             << endl
             << "if (f & " << flags_type << "::keep_dom)" << endl
             << "doc.setUserData (" << dom_node_key << "," << endl
             << "(c.get () ? &c : &d)," << endl
             << "0);"
             << endl;

          if (poly)
          {
            // aCC cannot handle an inline call to type_factory_map_instance.
            //
            os << "::xsd::cxx::tree::type_factory_map< " << char_type <<
              " >& tfm (" << endl
               << "::xsd::cxx::tree::type_factory_map_instance< 0, " <<
              char_type << " > ());"
               << endl
               << "::std::auto_ptr< ::xsd::cxx::tree::type > tmp (" << endl
               << "tfm.create (" << endl
               << strlit (e.name ()) << "," << endl
               << strlit (e.namespace_().name ()) << "," << endl
               << "&::xsd::cxx::tree::factory_impl< " << type << " >," << endl
               << "true, true, e, n, f, 0));"
               << endl
               << "if (tmp.get () != 0)"
               << "{";
          }
          else
          {
            os << "if (n.name () == " << strlit (e.name ()) << " &&" << endl
               << "n.namespace_ () == " << strlit (e.namespace_().name ()) << ")"
               << "{";

            if (fund)
            {
              os << "::std::auto_ptr< " << type << " > r (" << endl
                 << "new " << type << " (" << endl
                 << "::xsd::cxx::tree::traits< " << type << ", " <<
                char_type;

              if (t.is_a<SemanticGraph::Fundamental::Double> ())
                os << ", ::xsd::cxx::tree::schema_type::double_";
              else if (t.is_a<SemanticGraph::Fundamental::Decimal> ())
                os << ", ::xsd::cxx::tree::schema_type::decimal";

              os << " >::create (" << endl
                 << "e, f, 0)));";
            }
            else
            {
              os << "::std::auto_ptr< " << type << " > r (" << endl
                 << "::xsd::cxx::tree::traits< " << type << ", " <<
                char_type << " >::create (" << endl
                 << "e, f, 0));";
            }
          }

          if (poly)
          {
            os << endl
               << "::std::auto_ptr< " << type << " > r (" << endl
               << "dynamic_cast< " << type << "* > (tmp.get ()));"
               << endl
               << "if (r.get ())" << endl
               << "tmp.release ();"
               << "else" << endl
               << "throw ::xsd::cxx::tree::not_derived< " << char_type <<
              " > ();"
               << endl;
          }

          os << "return r;"
             << "}"
             << "throw ::xsd::cxx::tree::unexpected_element < " <<
            char_type << " > (" << endl
             << "n.name ()," << endl
             << "n.namespace_ ()," << endl
             << strlit (e.name ()) << "," << endl
             << strlit (e.namespace_().name ()) << ");"
             << "}";
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
    generate_parser_source (Context& ctx,
                            UnsignedLong first,
                            UnsignedLong last)
    {
      ctx.os << "#include <istream>" << endl
             << "#include <xsd/cxx/xml/sax/std-input-source.hxx>" << endl
             << "#include <xsd/cxx/tree/error-handler.hxx>" << endl
             << endl;

      Traversal::Schema schema;
      Traversal::Sources sources;
      Traversal::Names names_ns, names;
      Namespace ns (ctx, first, last);
      ElementFunction element (ctx);

      schema >> sources >> schema;
      schema >> names_ns >> ns >> names >> element;

      schema.dispatch (ctx.schema_root);
    }
  }
}
