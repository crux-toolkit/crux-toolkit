// file      : xsd/cxx/tree/serialization-source.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/tree/serialization-source.hxx>

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

namespace CXX
{
  namespace Tree
  {
    namespace
    {
      enum schema_type
      {
        st_other,
        st_double,
        st_decimal
      };

      enum schema_type
      schema_type (SemanticGraph::Type& t)
      {
        if (t.is_a<SemanticGraph::Fundamental::Double> ())
          return st_double;
        else if (t.is_a<SemanticGraph::Fundamental::Decimal> ())
          return st_decimal;
        else
          return st_other;
      }


      //
      //
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

          SemanticGraph::Type& item_type (l.argumented ().type ());
          String base (L"::xsd::cxx::tree::list< " +
                       item_type_name (item_type) + L", " + char_type);

          if (item_type.is_a<SemanticGraph::Fundamental::Double> ())
            base += L", ::xsd::cxx::tree::schema_type::double_";
          else if (item_type.is_a<SemanticGraph::Fundamental::Decimal> ())
            base += L", ::xsd::cxx::tree::schema_type::decimal";

          base += L" >";

          // operator<< (xercesc::DOMElement)
          //
          os << "void" << endl
             << "operator<< (" << xerces_ns << "::DOMElement& e, " <<
            "const " << name << "& i)"
             << "{"
             << "e << static_cast< const " << base << "& > (i);"
             << "}";

          // operator<< (xercesc::DOMAttr)
          //
          os << "void" << endl
             << "operator<< (" << xerces_ns << "::DOMAttr& a, " <<
            "const " << name << "& i)"
             << "{"
             << "a << static_cast< const " << base << "& > (i);"
             << "}";

          // operator<< (list_stream)
          //
          os << "void" << endl
             << "operator<< (" << list_stream_type << "& l," << endl
             << "const " << name << "& i)"
             << "{"
             << "l << static_cast< const " << base << "& > (i);"
             << "}";

          // Register with type factory map.
          //
          if (polymorphic && polymorphic_p (l) && !anonymous_p (l))
          {
            // Note that we are using the original type name.
            //
            String const& name (ename (l));

            os << "static" << endl
               << "const ::xsd::cxx::tree::type_serializer_initializer< 0, " <<
              char_type << ", " << name << " >" << endl
               << "_xsd_" << name << "_type_serializer_init (" << endl
               << strlit (l.name ()) << "," << endl
               << strlit (xml_ns_name (l)) << ");"
               << endl
               << endl;
          }
        }

      private:
        String
        item_type_name (SemanticGraph::Type& t)
        {
          std::wostringstream o;

          MemberTypeName type (*this, o);
          type.dispatch (t);

          return o.str ();
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

          String const& base (xs_string_type);

          // operator<< (xercesc::DOMElement)
          //
          os << "void" << endl
             << "operator<< (" << xerces_ns << "::DOMElement& e, " <<
            "const " << name << "& i)"
             << "{"
             << "e << static_cast< const " << base << "& > (i);"
             << "}";

          // operator<< (xercesc::DOMAttr)
          //
          os << "void" << endl
             << "operator<< (" << xerces_ns << "::DOMAttr& a, " <<
            "const " << name << "& i)"
             << "{"
             << "a << static_cast< const " << base << "& > (i);"
             << "}";

          // operator<< (list_stream)
          //
          os << "void" << endl
             << "operator<< (" << list_stream_type << "& l," << endl
             << "const " << name << "& i)"
             << "{"
             << "l << static_cast< const " << base << "& > (i);"
             << "}";

          // Register with type factory map.
          //
          if (polymorphic && polymorphic_p (u) && !anonymous_p (u))
          {
            // Note that we are using the original type name.
            //
            String const& name (ename (u));

            os << "static" << endl
               << "const ::xsd::cxx::tree::type_serializer_initializer< 0, " <<
              char_type << ", " << name << " >" << endl
               << "_xsd_" << name << "_type_serializer_init (" << endl
               << strlit (u.name ()) << "," << endl
               << strlit (xml_ns_name (u)) << ");"
               << endl
               << endl;
          }
        }
      };


      struct Enumeration: Traversal::Enumeration, Context
      {
        Enumeration (Context& c)
            : Context (c), base_ (c)
        {
          inherits_base_ >> base_;
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
          os << "void" << endl
             << "operator<< (" << xerces_ns << "::DOMElement& e, " <<
            "const " << name << "& i)"
             << "{"
             << "e << static_cast< const ";

          inherits (e, inherits_base_);

          os << "& > (i);"
             << "}";


          // operator<< (xercesc::DOMAttr)
          //
          os << "void" << endl
             << "operator<< (" << xerces_ns << "::DOMAttr& a, " <<
            "const " << name << "& i)"
             << "{"
             << "a << static_cast< const ";

          inherits (e, inherits_base_);

          os << "& > (i);"
             << "}";


          // operator<< (list_stream)
          //
          os << "void" << endl
             << "operator<< (" << list_stream_type << "& l," << endl
             << "const " << name << "& i)"
             << "{"
             << "l << static_cast< const ";

          inherits (e, inherits_base_);

          os << "& > (i);"
             << "}";


          // Register with type factory map.
          //
          if (polymorphic && polymorphic_p (e) && !anonymous_p (e))
          {
            // Note that we are using the original type name.
            //
            String const& name (ename (e));

            os << "static" << endl
               << "const ::xsd::cxx::tree::type_serializer_initializer< 0, " <<
              char_type << ", " << name << " >" << endl
               << "_xsd_" << name << "_type_serializer_init (" << endl
               << strlit (e.name ()) << "," << endl
               << strlit (xml_ns_name (e)) << ");"
               << endl
               << endl;
          }
        }

      private:
        Traversal::Inherits inherits_base_;
        BaseTypeName base_;
      };

      struct Element: Traversal::Element, Context
      {
        Element (Context& c, String const& scope_)
            : Context (c), scope (scope_)
        {
        }

        virtual Void
        traverse (Type& e)
        {
          if (skip (e))
            return;

          String const& aname (eaname (e));
          String ns (e.qualified_p () ? e.namespace_ ().name () : "");
          String type (scope + L"::" + etype (e));

          // Check if we need to handle xsi:type and substitution groups.
          // If this element's type is anonymous then we don't need to do
          // anything. Note that if the type is anonymous then it can't be
          // derived from which makes it impossible to substitute or
          // dynamically-type with xsi:type.
          //
          SemanticGraph::Type& t (e.type ());
          Boolean poly (polymorphic && polymorphic_p (t) && !anonymous_p (t));

          os << "// " << comment (e.name ()) << endl
             << "//" << endl;

          // aCC cannot handle an inline call to type_serializer_map_instance.
          //
          if (poly)
          {
            os << "{"
               << "::xsd::cxx::tree::type_serializer_map< " << char_type
               << " >& tsm (" << endl
               << "::xsd::cxx::tree::type_serializer_map_instance< 0, " <<
              char_type << " > ());"
               << endl;
          }

          if (max (e) != 1)
          {
            // sequence
            //
            os << "for (" << scope << "::" << econst_iterator (e) << endl
               << "b (i." << aname << " ().begin ()), " <<
              "n (i." << aname << " ().end ());" << endl
               << "b != n; ++b)"
               << "{";

            if (poly)
            {
              os << "if (typeid (" << type << ") == typeid (*b))"
                 << "{"
                 << xerces_ns << "::DOMElement& s (" << endl
                 << "::xsd::cxx::xml::dom::create_element (" << endl
                 << strlit (e.name ()) << "," << endl
                 << (ns ? strlit (ns) + L",\n" : L"")
                 << "e));"
                 << endl
                 << "s << *b;"
                 << "}"
                 << "else" << endl
                 << "tsm.serialize (" << endl
                 << strlit (e.name ()) << "," << endl
                 << strlit (ns) << "," << endl
                 << (e.global_p () ? "true" : "false") << ", " <<
                (e.qualified_p () ? "true" : "false") << ", e, *b);";
            }
            else
            {
              os << xerces_ns << "::DOMElement& s (" << endl
                 << "::xsd::cxx::xml::dom::create_element (" << endl
                 << strlit (e.name ()) << "," << endl
                 << (ns ? strlit (ns) + L",\n" : L"")
                 << "e));"
                 << endl;

              switch (schema_type (t))
              {
              case st_other:
                {
                  os << "s << *b;";
                  break;
                }
              case st_double:
                {
                  os << "s << " << as_double_type << " (*b);";
                  break;
                }
              case st_decimal:
                {
                  os << "s << " << as_decimal_type << " (*b);";
                  break;
                }
              }
            }

            os << "}";
          }
          else if (min (e) == 0)
          {
            // optional
            //
            os << "if (i." << aname << " ())"
               << "{";

            if (poly)
            {
              os << "const " << type << "& x (*i." << aname << " ());"
                 << "if (typeid (" << type << ") == typeid (x))"
                 << "{"
                 << xerces_ns << "::DOMElement& s (" << endl
                 << "::xsd::cxx::xml::dom::create_element (" << endl
                 << strlit (e.name ()) << "," << endl
                 << (ns ? strlit (ns) + L",\n" : L"")
                 << "e));"
                 << endl
                 << "s << x;"
                 << "}"
                 << "else" << endl
                 << "tsm.serialize (" << endl
                 << strlit (e.name ()) << "," << endl
                 << strlit (ns) << "," << endl
                 << (e.global_p () ? "true" : "false") << ", " <<
                (e.qualified_p () ? "true" : "false") << ", e, x);";
            }
            else
            {
              os << xerces_ns << "::DOMElement& s (" << endl
                 << "::xsd::cxx::xml::dom::create_element (" << endl
                 << strlit (e.name ()) << "," << endl
                 << (ns ? strlit (ns) + L",\n" : L"")
                 << "e));"
                 << endl;

              switch (schema_type (t))
              {
              case st_other:
                {
                  os << "s << *i." << aname << " ();";
                  break;
                }
              case st_double:
                {
                  os << "s << " << as_double_type << "(*i." << aname << " ());";
                  break;
                }
              case st_decimal:
                {
                  os << "s << " << as_decimal_type << "(*i." << aname << " ());";
                  break;
                }
              }
            }

            os << "}";
          }
          else
          {
            // one
            //
            if (poly)
            {
              os << "const " << type << "& x (i." << aname << " ());"
                 << "if (typeid (" << type << ") == typeid (x))"
                 << "{"
                 << xerces_ns << "::DOMElement& s (" << endl
                 << "::xsd::cxx::xml::dom::create_element (" << endl
                 << strlit (e.name ()) << "," << endl
                 << (ns ? strlit (ns) + L",\n" : L"")
                 << "e));"
                 << endl
                 << "s << x;"
                 << "}"
                 << "else" << endl
                 << "tsm.serialize (" << endl
                 << strlit (e.name ()) << "," << endl
                 << strlit (ns) << "," << endl
                 << (e.global_p () ? "true" : "false") << ", " <<
                (e.qualified_p () ? "true" : "false") << ", e, x);";
            }
            else
            {
              os << "{"
                 << xerces_ns << "::DOMElement& s (" << endl
                 << "::xsd::cxx::xml::dom::create_element (" << endl
                 << strlit (e.name ()) << "," << endl
                 << (ns ? strlit (ns) + L",\n" : L"")
                 << "e));"
                 << endl;

              switch (schema_type (t))
              {
              case st_other:
                {
                  os << "s << i." << aname << " ();";
                  break;
                }
              case st_double:
                {
                  os << "s << " << as_double_type << "(i." << aname << " ());";
                  break;
                }
              case st_decimal:
                {
                  os << "s << " << as_decimal_type << "(i." << aname << " ());";
                  break;
                }
              }

              os << "}";
            }
          }

          if (poly)
            os << "}";
        }

      private:
        String scope;
      };

      struct Any: Traversal::Any, Context
      {
        Any (Context& c, String const& scope_)
            : Context (c), scope (scope_)
        {
        }

        virtual Void
        traverse (Type& a)
        {
          String const& aname (eaname (a));

          os << "// " << ename (a) << endl
             << "//" << endl;

          if (max (a) != 1)
          {
            // sequence
            //
            os << "for (" << scope << "::" << econst_iterator (a) << endl
               << "b (i." << aname << " ().begin ()), " <<
              "n (i." << aname << " ().end ());" << endl
               << "b != n; ++b)"
               << "{"
               << "e.appendChild (" << endl
               << "e.getOwnerDocument ()->importNode (" << endl
               << "const_cast< " << xerces_ns <<
              "::DOMElement* > (&(*b)), true));"
               << "}";
          }
          else if (min (a) == 0)
          {
            // optional
            //
            os << "if (i." << aname << " ())"
               << "{"
               << "e.appendChild (" << endl
               << "e.getOwnerDocument ()->importNode (" << endl
               << "const_cast< " << xerces_ns << "::DOMElement* > (&(*i." <<
              aname << " ())), true));"
               << "}";
          }
          else
          {
            // one
            //
            os << "e.appendChild (" << endl
               << "e.getOwnerDocument ()->importNode (" << endl
               << "const_cast< " << xerces_ns << "::DOMElement* > (&(i." <<
              aname << " ())), true));"
               << endl;
          }
        }

      private:
        String scope;
      };

      struct Attribute: Traversal::Attribute, Context
      {
        Attribute (Context& c, String const& scope_)
            : Context (c), scope (scope_)
        {
        }

        virtual Void
        traverse (Type& a)
        {
          String const& aname (eaname (a));
          String ns (a.qualified_p () ? a.namespace_ ().name () : "");

          os << "// " << comment (a.name ()) << endl
             << "//" << endl;

          if (a.optional_p () && !a.default_p ())
          {
            os << "if (i." << aname << " ())"
               << "{"
               << xerces_ns << "::DOMAttr& a (" << endl
               << "::xsd::cxx::xml::dom::create_attribute (" << endl
               << strlit (a.name ()) << "," << endl
               << (ns ? strlit (ns) + L",\n" : L"")
               << "e));"
               << endl;

            switch (schema_type (a.type ()))
            {
            case st_other:
              {
                os << "a << *i." << aname << " ();";
                break;
              }
            case st_double:
              {
                os << "a << " << as_double_type << "(*i." << aname << " ());";
                break;
              }
            case st_decimal:
              {
                os << "a << " << as_decimal_type << "(*i." << aname << " ());";
                break;
              }
            }

            os << "}";
          }
          else
          {
            // Make sure we serialize required fixed attributes.
            //
            if (a.optional_p () &&
                options.value<CLI::omit_default_attributes> ())
            {
              os << "if (i." << aname << " () != " << scope <<
                "::" << edefault_value (a) << " ())";
            }

            os << "{"
               << xerces_ns << "::DOMAttr& a (" << endl
               << "::xsd::cxx::xml::dom::create_attribute (" << endl
               << strlit (a.name ()) << "," << endl
               << (ns ? strlit (ns) + L",\n" : L"")
               << "e));"
               << endl;

            switch (schema_type (a.type ()))
            {
            case st_other:
              {
                os << "a << i." << aname << " ();";
                break;
              }
            case st_double:
              {
                os << "a << " << as_double_type << "(i." << aname << " ());";
                break;
              }
            case st_decimal:
              {
                os << "a << " << as_decimal_type << "(i." << aname << " ());";
                break;
              }
            }

            os << "}";
          }
        }

      private:
        String scope;
      };

      struct AnyAttribute: Traversal::AnyAttribute, Context
      {
        AnyAttribute (Context& c, String const& scope_)
            : Context (c), scope (scope_)
        {
        }

        virtual Void
        traverse (Type& a)
        {
          String const& aname (eaname (a));

          os << "// " << ename (a) << endl
             << "//" << endl;

          os << "for (" << scope << "::" << econst_iterator (a) << endl
             << "b (i." << aname << " ().begin ()), " <<
            "n (i." << aname << " ().end ());" << endl
             << "b != n; ++b)"
             << "{"
             << xerces_ns << "::DOMAttr* a (" << endl
             << "static_cast< " << xerces_ns << "::DOMAttr* > (" << endl
             << "e.getOwnerDocument ()->importNode (" << endl
             << "const_cast< " << xerces_ns << "::DOMAttr* > (&(*b)), true)));"
             << endl
             << "if (a->getLocalName () == 0)" << endl
             << "e.setAttributeNode (a);"
             << "else" << endl
             << "e.setAttributeNodeNS (a);"
             << "}";
        }

      private:
        String scope;
      };

      struct Complex: Traversal::Complex, Context
      {
        Complex (Context& c)
            : Context (c), base_ (c)
        {
          inherits_ >> base_;
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
          os << "void" << endl
             << "operator<< (" << xerces_ns << "::DOMElement& e, " <<
            "const " << name << "& i)"
             << "{";

          if (c.inherits_p ())
          {
            os << "e << static_cast< const ";
            inherits (c, inherits_);
            os << "& > (i);"
               << endl;
          }
          else
            os << "e << static_cast< const " << any_type << "& > (i);"
               << endl;

          // Serialize anyAttribute content first so that is gets
          // overriden by schema-defined attributes.
          //
          if (options.value<CLI::generate_wildcard> ())
          {
            AnyAttribute any_attribute (*this, name);
            Traversal::Names names (any_attribute);

            Complex::names (c, names);
          }

          {
            Traversal::Names names;
            Any any (*this, name);
            Element element (*this, name);
            Attribute attribute (*this, name);

            names >> element;
            names >> attribute;

            if (options.value<CLI::generate_wildcard> ())
              names >> any;

            Complex::names (c, names);
          }

          os << "}";


          Boolean simple (true);
          {
            IsSimpleType t (simple);
            t.dispatch (c);
          }

          if (simple)
          {
            Boolean hb (c.inherits_p ());

            // operator<< (xercesc::DOMAttr)
            //
            os << "void" << endl
               << "operator<< (" << xerces_ns << "::DOMAttr&" <<
              (hb ? " a" : "") << ", " <<
              "const " << name << "&" << (hb ? " i" : "") << ")"
               << "{";

            if (hb)
            {
              os << "a << static_cast< const ";

              inherits (c, inherits_);

              os << "& > (i);";
            }

            os << "}";


            // operator<< (list_stream)
            //
            os << "void" << endl
               << "operator<< (" << list_stream_type << "&" <<
              (hb ? " l" : "") << "," << endl
               << "const " << name << "&" << (hb ? " i" : "") << ")"
               << "{";

            if (hb)
            {
              os << "l << static_cast< const ";

              inherits (c, inherits_);

              os << "& > (i);";
            }

            os << "}";
          }

          // Register with type factory map.
          //
          if (polymorphic && polymorphic_p (c) && !anonymous_p (c))
          {
            // Note that we are using the original type name.
            //
            String const& name (ename (c));

            os << "static" << endl
               << "const ::xsd::cxx::tree::type_serializer_initializer< 0, " <<
              char_type << ", " << name << " >" << endl
               << "_xsd_" << name << "_type_serializer_init (" << endl
               << strlit (c.name ()) << "," << endl
               << strlit (xml_ns_name (c)) << ");"
               << endl
               << endl;
          }
        }

      private:
        Traversal::Inherits inherits_;
        BaseTypeName base_;
      };


      // Generate substitution group map entry.
      //
      struct GlobalElement: Traversal::Element,
                            GlobalElementBase,
                            Context
      {
        GlobalElement (Context& c)
            : GlobalElementBase (c), Context (c), type_name_ (c)
        {
          belongs_ >> type_name_;
        }

        virtual Void
        traverse (Type& e)
        {
          if (polymorphic && e.substitutes_p ())
          {
            Type& r (e.substitutes ().root ());

            String const& name (ename (e));

            os << "static" << endl
               << "const ::xsd::cxx::tree::element_serializer_initializer< 0, " <<
              char_type << ", ";

            belongs (e, belongs_);

            os << " >" << endl
               << "_xsd_" << name << "_element_serializer_init (" << endl
               << strlit (r.name ()) << "," << endl
               << strlit (r.namespace_ ().name ()) << "," << endl
               << strlit (e.name ()) << "," << endl
               << strlit (e.namespace_ ().name ()) << ");"
               << endl
               << endl;
          }
        }

      private:
        Traversal::Belongs belongs_;
        MemberTypeName type_name_;

      };

      struct ElementType: Traversal::Element,
                          GlobalElementBase,
                          Context
      {
        ElementType (Context& c)
            : GlobalElementBase (c),
              Context (c),
              element_map_ (c.options.value<CLI::generate_element_map> ())
        {
        }

        virtual Void
        traverse (Type& e)
        {
          if (doc_root_p (e))
          {
            String const& name (ename (e));

            // operator<< (xercesc::DOMElement)
            //
            os << "void" << endl
               << "operator<< (" << xerces_ns << "::DOMElement& e, " <<
              "const " << name << "& i)"
               << "{"
               << "e << i." << eaname (e) << " ();"
               << "}";

            // Element map registration.
            //
            if (element_map_)
            {
              SemanticGraph::Context& ec (e.context ());
              String const& aname (ec.get<String> ("element-name"));
              String const& ans (ec.get<String> ("element-ns"));

              os << "static " << endl
                 << "const ::xsd::cxx::tree::serializer_init< " <<
                name << ", " << char_type << ", " << any_type << " >" << endl
                 << "_xsd_" << name << "_serializer_init (" <<
                name << "::" << aname << " (), " <<
                name << "::" << ans << " ());"
                 << endl;
            }
          }
        }

      private:
        Boolean element_map_;
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
          String ns (e.namespace_ ().name ());
          String const& error_handler (error_handler_type);
          String const& namespace_infomap (namespace_infomap_type);

          SemanticGraph::Type& type (e.type ());

          // Note that I am using fq-name in function calls because g++ gets
          // confused if the name is 'type'. (see tests/schema/anonymous)
          //

          // Check if we need to handle xsi:type and substitution groups.
          // If this element's type is anonymous then we don't need to do
          // anything.
          //
          Boolean poly (polymorphic &&
                        polymorphic_p (type) &&
                        !anonymous_p (type));

          // To std::ostream.
          //
          os << "void" << endl
             << name << " (::std::ostream& o," << endl
             << "const " << type_name (e) << "& s," << endl
             << "const " << namespace_infomap << "& m," << endl
             << "const " << string_type << "& e," << endl
             << flags_type << " f)"
             << "{"
             << "::xsd::cxx::xml::auto_initializer i (" << endl
             << "(f & " << flags_type << "::dont_initialize) == 0);"
             << endl
             << dom_auto_ptr << "< " << xerces_ns <<
            "::DOMDocument > d (" << endl
             << fq_name (e, "serializer") << " (s, m, f));"
             << endl
             << "::xsd::cxx::tree::error_handler< " << char_type << " > h;"
             << endl
             << "::xsd::cxx::xml::dom::ostream_format_target t (o);"
             << "if (!::xsd::cxx::xml::dom::serialize (t, *d, e, h, f))"
             << "{"
             << "h.throw_if_failed< ::xsd::cxx::tree::serialization< " <<
            char_type << " > > ();"
             << "}"
             << "}";

          os << "void" << endl
             << name << " (::std::ostream& o," << endl
             << "const " << type_name (e) << "& s," << endl
             << error_handler << "& h," << endl
             << "const " << namespace_infomap << "& m," << endl
             << "const " << string_type << "& e," << endl
             << flags_type << " f)"
             << "{"
             << "::xsd::cxx::xml::auto_initializer i (" << endl
             << "(f & " << flags_type << "::dont_initialize) == 0);"
             << endl
             << dom_auto_ptr << "< " << xerces_ns <<
            "::DOMDocument > d (" << endl
             << fq_name (e, "serializer") << " (s, m, f));"
             << "::xsd::cxx::xml::dom::ostream_format_target t (o);"
             << "if (!::xsd::cxx::xml::dom::serialize (t, *d, e, h, f))"
             << "{"
             << "throw ::xsd::cxx::tree::serialization< " <<
            char_type << " > ();"
             << "}"
             << "}";

          os << "void" << endl
             << name << " (::std::ostream& o," << endl
             << "const " << type_name (e) << "& s," << endl
             << xerces_ns << "::DOMErrorHandler& h," << endl
             << "const " << namespace_infomap << "& m," << endl
             << "const " << string_type << "& e," << endl
             << flags_type << " f)"
             << "{"
             << dom_auto_ptr << "< " << xerces_ns <<
            "::DOMDocument > d (" << endl
             << fq_name (e, "serializer") << " (s, m, f));"
             << "::xsd::cxx::xml::dom::ostream_format_target t (o);"
             << "if (!::xsd::cxx::xml::dom::serialize (t, *d, e, h, f))"
             << "{"
             << "throw ::xsd::cxx::tree::serialization< " <<
            char_type << " > ();"
             << "}"
             << "}";

          // To XMLFormatTarget.
          //
          os << "void" << endl
             << name << " (" << xerces_ns << "::XMLFormatTarget& t," << endl
             << "const " << type_name (e) << "& s," << endl
             << "const " << namespace_infomap << "& m," << endl
             << "const " << string_type << "& e," << endl
             << flags_type << " f)"
             << "{"
             << dom_auto_ptr << "< " << xerces_ns <<
            "::DOMDocument > d (" << endl
             << fq_name (e, "serializer") << " (s, m, f));"
             << endl
             << "::xsd::cxx::tree::error_handler< " << char_type << " > h;"
             << endl
             << "if (!::xsd::cxx::xml::dom::serialize (t, *d, e, h, f))"
             << "{"
             << "h.throw_if_failed< ::xsd::cxx::tree::serialization< " <<
            char_type << " > > ();"
             << "}"
             << "}";

          os << "void" << endl
             << name << " (" << xerces_ns << "::XMLFormatTarget& t," << endl
             << "const " << type_name (e) << "& s," << endl
             << error_handler << "& h," << endl
             << "const " << namespace_infomap << "& m," << endl
             << "const " << string_type << "& e," << endl
             << flags_type << " f)"
             << "{"
             << dom_auto_ptr << "< " << xerces_ns <<
            "::DOMDocument > d (" << endl
             << fq_name (e, "serializer") << " (s, m, f));"
             << "if (!::xsd::cxx::xml::dom::serialize (t, *d, e, h, f))"
             << "{"
             << "throw ::xsd::cxx::tree::serialization< " <<
            char_type << " > ();"
             << "}"
             << "}";

          os << "void" << endl
             << name << " (" << xerces_ns << "::XMLFormatTarget& t," << endl
             << "const " << type_name (e) << "& s," << endl
             << xerces_ns << "::DOMErrorHandler& h," << endl
             << "const " << namespace_infomap << "& m," << endl
             << "const " << string_type << "& e," << endl
             << flags_type << " f)"
             << "{"
             << dom_auto_ptr << "< " << xerces_ns <<
            "::DOMDocument > d (" << endl
             << fq_name (e, "serializer") << " (s, m, f));"
             << "if (!::xsd::cxx::xml::dom::serialize (t, *d, e, h, f))"
             << "{"
             << "throw ::xsd::cxx::tree::serialization< " <<
            char_type << " > ();"
             << "}"
             << "}";

          // To an existing DOM instance.
          //
          os << "void" << endl
             << name << " (" << xerces_ns << "::DOMDocument& d," << endl
             << "const " << type_name (e) << "& s," << endl
             << flags_type << ")"
             << "{"
             << xerces_ns << "::DOMElement& e (*d.getDocumentElement ());"
             << "const " << qname_type << " n (" << endl
             << "::xsd::cxx::xml::dom::name< " << char_type << " > (e));"
             << endl;

          if (poly)
          {
            os << "if (typeid (" << type_name (e) << ") == typeid (s))"
               << "{";
          }

          os << "if (n.name () == " << strlit (e.name ()) << " &&" << endl
             << "n.namespace_ () == " << strlit (ns) << ")"
             << "{";

          switch (schema_type (type))
          {
          case st_other:
            {
              os << "e << s;";
              break;
            }
          case st_double:
            {
              os << "e << " << as_double_type << "(s);";
              break;
            }
          case st_decimal:
            {
              os << "e << " << as_decimal_type << "(s);";
              break;
            }
          }

          os << "}"
             << "else"
             << "{"
             << "throw ::xsd::cxx::tree::unexpected_element < " <<
            char_type << " > (" << endl
             << "n.name ()," << endl
             << "n.namespace_ ()," << endl
             << strlit (e.name ()) << "," << endl
             << strlit (ns) << ");"
             << "}";

          if (poly)
          {
            // aCC cannot handle an inline call to
            // type_serializer_map_instance.
            //
            os << "}"
               << "else"
               << "{"
               << "::xsd::cxx::tree::type_serializer_map< " << char_type
               << " >& tsm (" << endl
               << "::xsd::cxx::tree::type_serializer_map_instance< 0, " <<
              char_type << " > ());"
               << endl
               << "tsm.serialize (" << endl
               << strlit (e.name ()) << "," << endl
               << strlit (e.namespace_().name ()) << "," << endl
               << "e, n, s);"
               << "}";
          }

          os << "}";


          // To a new DOM instance.
          //
          os << dom_auto_ptr << "< " << xerces_ns << "::DOMDocument >" << endl
             << name << " (const " << type_name (e) << "& s," << endl
             << "const " << namespace_infomap << "& m," << endl
             << flags_type << " f)"
             << "{";

          if (poly)
          {
            // aCC cannot handle an inline call to
            // type_serializer_map_instance as well as the direct
            // auto_ptr assignment.
            //
            os << dom_auto_ptr << "< " << xerces_ns << "::DOMDocument > d;"
               << endl
               << "if (typeid (" << type_name (e) << ") == typeid (s))"
               << "{"
               << dom_auto_ptr << "< " << xerces_ns <<
              "::DOMDocument > r (" << endl
               << "::xsd::cxx::xml::dom::serialize< " <<
              char_type << " > (" << endl
               << strlit (e.name ()) << "," << endl
               << strlit (ns) << "," << endl
               << "m, f));"
               << "d = r;"
               << "}"
               << "else"
               << "{"
               << "::xsd::cxx::tree::type_serializer_map< " << char_type
               << " >& tsm (" << endl
               << "::xsd::cxx::tree::type_serializer_map_instance< 0, " <<
              char_type << " > ());"
               << endl
               << dom_auto_ptr << "< " << xerces_ns <<
              "::DOMDocument > r (" << endl
               << "tsm.serialize (" << endl
               << strlit (e.name ()) << "," << endl
               << strlit (e.namespace_().name ()) << "," << endl
               << "m, s, f));"
               << "d = r;"
               << "}";
          }
          else
          {
            os << dom_auto_ptr << "< " << xerces_ns <<
              "::DOMDocument > d (" << endl
               << "::xsd::cxx::xml::dom::serialize< " <<
              char_type << " > (" << endl
               << strlit (e.name ()) << "," << endl
               << strlit (ns) << "," << endl
               << "m, f));"
               << endl;
          }

          os << fq_name (e, "serializer") << " (*d, s, f);"
             << "return d;"
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
    generate_serialization_source (Context& ctx,
                                   UnsignedLong first,
                                   UnsignedLong last)
    {
      Boolean elemen_type (ctx.options.value<CLI::generate_element_type> ());

      if (!elemen_type)
        ctx.os << "#include <ostream>" << endl
               << "#include <xsd/cxx/tree/error-handler.hxx>" << endl;

      ctx.os << "#include <xsd/cxx/xml/dom/serialization-source.hxx>" << endl
             << endl;

      if (ctx.polymorphic)
      {
        ctx.os << "#include <xsd/cxx/tree/type-serializer-map.hxx>" << endl
               << endl;

        Boolean import_maps (ctx.options.value<CLI::import_maps> ());
        Boolean export_maps (ctx.options.value<CLI::export_maps> ());

        if (import_maps || export_maps)
        {
          ctx.os << "#ifndef XSD_NO_EXPORT" << endl
                 << endl
                 << "namespace xsd"
                 << "{"
                 << "namespace cxx"
                 << "{"
                 << "namespace tree"
                 << "{"
                 << "#ifdef _MSC_VER" << endl;

          if (export_maps)
            ctx.os << "template struct __declspec (dllexport) " <<
              "type_serializer_plate< 0, " << ctx.char_type << " >;";

          if (import_maps)
            ctx.os << "template struct __declspec (dllimport) " <<
              "type_serializer_plate< 0, " << ctx.char_type << " >;";

          ctx.os << "#elif defined(__GNUC__) && __GNUC__ >= 4" << endl
                 << "template struct __attribute__ ((visibility(\"default\"))) " <<
            "type_serializer_plate< 0, " << ctx.char_type << " >;";

          ctx.os << "#elif defined(XSD_MAP_VISIBILITY)" << endl
                 << "template struct XSD_MAP_VISIBILITY " <<
            "type_serializer_plate< 0, " << ctx.char_type << " >;";

          ctx.os << "#endif" << endl
                 << "}"  // tree
                 << "}"  // cxx
                 << "}"  // xsd
                 << "#endif // XSD_NO_EXPORT" << endl
                 << endl;
        }

        ctx.os << "namespace _xsd"
               << "{"
               << "static" << endl
               << "const ::xsd::cxx::tree::type_serializer_plate< 0, " <<
          ctx.char_type << " >" << endl
               << "type_serializer_plate_init;"
               << "}";
      }

      Traversal::Schema schema;

      Traversal::Sources sources;
      Traversal::Names names_ns, names;

      Namespace ns (ctx, first, last);

      List list (ctx);
      Union union_ (ctx);
      Complex complex (ctx);
      Enumeration enumeration (ctx);
      GlobalElement element (ctx);
      ElementType element_type (ctx);
      ElementFunction element_function (ctx);


      schema >> sources >> schema;
      schema >> names_ns >> ns >> names;

      names >> list;
      names >> union_;
      names >> complex;
      names >> enumeration;
      names >> element;

      if (elemen_type)
        names >> element_type;
      else
        names >> element_function;

      schema.dispatch (ctx.schema_root);
    }
  }
}
