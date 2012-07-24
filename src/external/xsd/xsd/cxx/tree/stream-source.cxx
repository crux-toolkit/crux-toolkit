// file      : xsd/cxx/tree/stream-source.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/tree/stream-source.hxx>

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

          SemanticGraph::Type& item_type (l.argumented ().type ());
          String base (L"::xsd::cxx::tree::list< " +
                       item_type_name (item_type) + L", " + char_type);

          if (item_type.is_a<SemanticGraph::Fundamental::Double> ())
            base += L", ::xsd::cxx::tree::schema_type::double_";
          else if (item_type.is_a<SemanticGraph::Fundamental::Decimal> ())
            base += L", ::xsd::cxx::tree::schema_type::decimal";

          base += L" >";

          os << std_ostream_type << "&" << endl
             << "operator<< (" << std_ostream_type << "& o, " <<
            "const " << name << "& i)"
             << "{"
             << "return o << static_cast< const " << base << "& > (i);"
             << "}";

          // Register with ostream map.
          //
          if (polymorphic && polymorphic_p (l) && !anonymous_p (l))
          {
            // Note that we are using the original type name.
            //
            String const& name (ename (l));

            os << "static" << endl
               << "const ::xsd::cxx::tree::std_ostream_initializer< 0, " <<
              char_type << ", " << name << " >" << endl
               << "_xsd_" << name << "_std_ostream_init;"
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

          os << std_ostream_type << "&" << endl
             << "operator<< (" << std_ostream_type << "& o, " <<
            "const " << name << "& i)"
             << "{"
             << "return o << static_cast< const " << xs_string_type << "& > (i);"
             << "}";

          // Register with ostream map.
          //
          if (polymorphic && polymorphic_p (u) && !anonymous_p (u))
          {
            // Note that we are using the original type name.
            //
            String const& name (ename (u));

            os << "static" << endl
               << "const ::xsd::cxx::tree::std_ostream_initializer< 0, " <<
              char_type << ", " << name << " >" << endl
               << "_xsd_" << name << "_std_ostream_init;"
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

          Boolean string_based (false);
          {
            IsStringBasedType t (string_based);
            t.dispatch (e);
          }

          Boolean enum_based (false);
          if (string_based)
          {
            SemanticGraph::Enumeration* be (0);
            IsEnumBasedType t (be);
            t.dispatch (e);

            enum_based = (be != 0);
          }

          // If we are based on an enum then the value type is just an
          // alias and we don't need to generate this operator again.
          //
          if (string_based && !enum_based)
          {
            os << std_ostream_type << "&" << endl
               << "operator<< (" << std_ostream_type << "& o, " <<
              name << "::" << evalue (e) << " i)"
               << "{"
               << "return o << " << name << "::_xsd_" << name <<
              "_literals_[i];"
               << "}";
          }

          os << std_ostream_type << "&" << endl
             << "operator<< (" << std_ostream_type << "& o, " <<
            "const " << name << "& i)"
             << "{"
             << "return o << static_cast< const ";

          inherits (e, inherits_base_);

          os << "& > (i);"
             << "}";

          // Register with ostream map.
          //
          if (polymorphic && polymorphic_p (e) && !anonymous_p (e))
          {
            // Note that we are using the original type name.
            //
            String const& name (ename (e));

            os << "static" << endl
               << "const ::xsd::cxx::tree::std_ostream_initializer< 0, " <<
              char_type << ", " << name << " >" << endl
               << "_xsd_" << name << "_std_ostream_init;"
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

          // Check if we need to handle xsi:type and substitution groups.
          // If this element's type is anonymous then we don't need to do
          // anything. Note that if the type is anonymous then it can't be
          // derived from which makes it impossible to substitute or
          // dynamically-type with xsi:type.
          //
          SemanticGraph::Type& t (e.type ());
          Boolean poly (polymorphic && polymorphic_p (t) && !anonymous_p (t));

          // aCC cannot handle an inline call to std_ostream_map_instance.
          //
          if (poly)
          {
            os << "{"
               << "::xsd::cxx::tree::std_ostream_map< " << char_type
               << " >& om (" << endl
               << "::xsd::cxx::tree::std_ostream_map_instance< 0, " <<
              char_type << " > ());"
               << endl;
          }

          if (max (e) != 1)
          {
            // sequence
            //
            os << "for (" << scope << "::" << econst_iterator (e) << endl
               << "b (i." << aname << " ().begin ()), " <<
              "e (i." << aname << " ().end ());" << endl
               << "b != e; ++b)"
               << "{"
               << "o << ::std::endl << " << strlit (e.name () + L": ");

            if (!poly)
              os << " << *b;";
            else
              os << ";"
                 << "om.insert (o, *b);";

            os << "}";
          }
          else if (min (e) == 0)
          {
            // optional
            //

            os << "if (i." << aname << " ())"
               << "{"
               << "o << ::std::endl << " << strlit (e.name () + L": ");

            if (!poly)
              os << " << *i." << aname << " ();";
            else
              os << ";"
                 << "om.insert (o, *i." << aname << " ());";

            os << "}";
          }
          else
          {
            // one
            //
            os << "o << ::std::endl << " << strlit (e.name () + L": ");

            if (!poly)
              os << " << i." << aname << " ();";
            else
              os << ";"
                 << "om.insert (o, i." << aname << " ());";
          }

          if (poly)
            os << "}";
        }

      private:
        String scope;
      };

      struct Attribute: Traversal::Attribute, Context
      {
        Attribute (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& a)
        {
          String const& aname (eaname (a));

          if (a.optional_p () && !a.default_p ())
          {
            os << "if (i." << aname << " ())"
               << "{"
               << "o << ::std::endl << " << strlit (a.name () + L": ") <<
              " << *i." << aname << " ();"
               << "}";
          }
          else
          {
            os << "o << ::std::endl << " << strlit (a.name () + L": ") <<
              " << i." << aname << " ();";
          }
        }
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

          //
          //
          Boolean has_body (has<Traversal::Member> (c) || c.inherits_p ());

          os << std_ostream_type << "&" << endl
             << "operator<< (" << std_ostream_type << "& o, " <<
            "const " << name << "&" << (has_body ? " i" : "") << ")"
             << "{";

          if (c.inherits_p ())
          {
            os << "o << static_cast< const ";

            inherits (c, inherits_);

            os << "& > (i);"
               << endl;
          }

          {
            Traversal::Names names_member;
            Element element (*this, name);
            Attribute attribute (*this);

            names_member >> element;
            names_member >> attribute;

            names (c, names_member);
          }

          os << "return o;"
             << "}";

          // Register with ostream map.
          //
          if (polymorphic && polymorphic_p (c) && !anonymous_p (c))
          {
            // Note that we are using the original type name.
            //
            String const& name (ename (c));

            os << "static" << endl
               << "const ::xsd::cxx::tree::std_ostream_initializer< 0, " <<
              char_type << ", " << name << " >" << endl
               << "_xsd_" << name << "_std_ostream_init;"
               << endl;
          }
        }

      private:
        Traversal::Inherits inherits_;
        BaseTypeName base_;
      };
    }

    Void
    generate_stream_source (Context& ctx,
                            UnsignedLong first,
                            UnsignedLong last)
    {
      String c (ctx.char_type);

      ctx.os << "#include <ostream>" << endl
             << endl;

      if (ctx.polymorphic)
      {
        ctx.os << "#include <xsd/cxx/tree/std-ostream-map.hxx>" << endl
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
              "std_ostream_plate< 0, " << ctx.char_type << " >;";

          if (import_maps)
            ctx.os << "template struct __declspec (dllimport) " <<
              "std_ostream_plate< 0, " << ctx.char_type << " >;";

          ctx.os << "#elif defined(__GNUC__) && __GNUC__ >= 4" << endl
                 << "template struct __attribute__ ((visibility(\"default\"))) " <<
            "std_ostream_plate< 0, " << ctx.char_type << " >;";

          ctx.os << "#elif defined(XSD_MAP_VISIBILITY)" << endl
                 << "template struct XSD_MAP_VISIBILITY " <<
            "std_ostream_plate< 0, " << ctx.char_type << " >;";

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
               << "const ::xsd::cxx::tree::std_ostream_plate< 0, " <<
          ctx.char_type << " >" << endl
               << "std_ostream_plate_init;"
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

      schema >> sources >> schema;
      schema >> names_ns >> ns >> names;

      names >> list;
      names >> union_;
      names >> complex;
      names >> enumeration;

      schema.dispatch (ctx.schema_root);
    }
  }
}
