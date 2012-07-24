// file      : xsd/cxx/tree/stream-insertion-source.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/tree/stream-insertion-source.hxx>

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

namespace CXX
{
  namespace Tree
  {
    namespace
    {
      typedef Containers::Vector<NarrowString> Streams;

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

          UnsignedLong n (0);
          Streams const& st (options.value<CLI::generate_insertion> ());
          for (Streams::ConstIterator i (st.begin ()); i != st.end (); ++i)
          {
            String stream_type ("::xsd::cxx::tree::ostream< " + *i + " >");

            os << stream_type << "&" << endl
               << "operator<< (" << stream_type << "& s," << endl
               << "const " << name << "& x)"
               << "{"
               << "return s << static_cast< const " << base << "& > (x);"
               << "}";

            // Register with type map.
            //
            if (polymorphic && polymorphic_p (l) && !anonymous_p (l))
            {
              // Note that we are using the original type name.
              //
              String const& name (ename (l));

              os << "static" << endl
                 << "const ::xsd::cxx::tree::stream_insertion_initializer< " <<
                "0, " << i->c_str () << ", " <<
                char_type << ", " << name << " >" << endl
                 << "_xsd_" << name << "_stream_insertion_init_" <<
                n++ << " (" << endl
                 << strlit (l.name ()) << "," << endl
                 << strlit (xml_ns_name (l)) << ");"
                 << endl;
            }
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

          UnsignedLong n (0);
          Streams const& st (options.value<CLI::generate_insertion> ());
          for (Streams::ConstIterator i (st.begin ()); i != st.end (); ++i)
          {
            String stream_type ("::xsd::cxx::tree::ostream< " + *i + " >");

            os << stream_type << "&" << endl
               << "operator<< (" << stream_type << "& s," << endl
               << "const " << name << "& x)"
               << "{"
               << "return s << static_cast< const " << base << "& > (x);"
               << "}";

            // Register with type map.
            //
            if (polymorphic && polymorphic_p (u) && !anonymous_p (u))
            {
              // Note that we are using the original type name.
              //
              String const& name (ename (u));

              os << "static" << endl
                 << "const ::xsd::cxx::tree::stream_insertion_initializer< " <<
                "0, " << i->c_str () << ", " <<
                char_type << ", " << name << " >" << endl
                 << "_xsd_" << name << "_stream_insertion_init_" <<
                n++ << " (" << endl
                 << strlit (u.name ()) << "," << endl
                 << strlit (xml_ns_name (u)) << ");"
                 << endl;
            }
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

          UnsignedLong n (0);
          Streams const& st (options.value<CLI::generate_insertion> ());
          for (Streams::ConstIterator i (st.begin ()); i != st.end (); ++i)
          {
            String stream_type ("::xsd::cxx::tree::ostream< " + *i + " >");

            os << stream_type << "&" << endl
               << "operator<< (" << stream_type << "& s," << endl
               << "const " << name << "& x)"
               << "{"
               << "return s << static_cast< const ";

            inherits (e, inherits_base_);

            os << "& > (x);"
               << "}";

            // Register with type map.
            //
            if (polymorphic && polymorphic_p (e) && !anonymous_p (e))
            {
              // Note that we are using the original type name.
              //
              String const& name (ename (e));

              os << "static" << endl
                 << "const ::xsd::cxx::tree::stream_insertion_initializer< " <<
                "0, " << i->c_str () << ", " <<
                char_type << ", " << name << " >" << endl
                 << "_xsd_" << name << "_stream_insertion_init_" <<
                n++ << " (" << endl
                 << strlit (e.name ()) << "," << endl
                 << strlit (xml_ns_name (e)) << ");"
                 << endl;
            }
          }
        }

      private:
        Traversal::Inherits inherits_base_;
        BaseTypeName base_;
      };

      struct Element: Traversal::Element, Context
      {
        Element (Context& c, String const& scope_, String const& stream_)
            : Context (c), scope (scope_), stream (stream_)
        {
        }

        virtual Void
        traverse (Type& e)
        {
          if (skip (e)) return;

          String const& aname (eaname (e));
          SemanticGraph::Type& t (e.type ());
          String type (scope + L"::" + etype (e));

          // Figure out if we need to generate polymorphic code. If this
          // elemen's type is anonymous then we don't need to do anything.
          // Note that if the type is anonymous then it can't be derived
          // from which makes it impossible to substitute or dynamically-
          // type with xsi:type.
          //
          Boolean poly (polymorphic && polymorphic_p (t) && !anonymous_p (t));

          if (max (e) != 1)
          {
            // sequence
            //
            os << "{"
               << "const " << scope << "::" << econtainer (e) << "& c (" <<
              "x." << aname << " ());"
               << "s << ::xsd::cxx::tree::ostream_common::as_size< " <<
              "::std::size_t > (c.size ());";

            os << "for (" << scope << "::" << econst_iterator (e) << endl
               << "i (c.begin ()), e (c.end ());" << endl
               << "i != e; ++i)"
               << "{";

            if (poly)
            {
              os << "bool d (typeid (" << type << ") != typeid (*i));"
                 << "s << d;"
                 << "if (!d)" << endl
                 << "s << *i;"
                 << "else" << endl
                 << "::xsd::cxx::tree::stream_insertion_map_instance< 0, " <<
                stream << ", " << char_type << " > ().insert (s, *i);";
            }
            else
              os << "s << *i;";

            os << "}" // for
               << "}";
          }
          else if (min (e) == 0)
          {
            // optional
            //
            os << "{"
               << "bool p (x." << aname << " ());"
               << "s << p;"
               << "if (p)";

            if (poly)
            {
              os << "{"
                 << "const " << type << "& i (*x." << aname << " ());"
                 << "bool d (typeid (" << type << ") != typeid (i));"
                 << "s << d;"
                 << "if (!d)" << endl
                 << "s << i;"
                 << "else" << endl
                 << "::xsd::cxx::tree::stream_insertion_map_instance< 0, " <<
                stream << ", " << char_type << " > ().insert (s, i);"
                 << "}";
            }
            else
              os << endl
                 << "s << *x." << aname << " ();";

            os << "}";
          }
          else
          {
            // one
            //
            if (poly)
            {
              os << "{"
                 << "const " << type << "& i (x." << aname << " ());"
                 << "bool d (typeid (" << type << ") != typeid (i));"
                 << "s << d;"
                 << "if (!d)" << endl
                 << "s << i;"
                 << "else" << endl
                 << "::xsd::cxx::tree::stream_insertion_map_instance< 0, " <<
                stream << ", " << char_type << " > ().insert (s, i);"
                 << "}";
            }
            else
              os << "s << x." << aname << " ();";
          }
        }

      private:
        String scope;
        String stream;
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
            os << "{"
               << "bool p (x." << aname << " ());"
               << "s << p;"
               << "if (p)" << endl
               << "s << *x." << aname << " ();"
               << "}";
          }
          else
          {
            os << "s << x." << aname << " ();";
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

          Boolean has_body (has<Traversal::Member> (c) || c.inherits_p ());

          UnsignedLong n (0);
          Streams const& st (options.value<CLI::generate_insertion> ());
          for (Streams::ConstIterator i (st.begin ()); i != st.end (); ++i)
          {
            String stream_type ("::xsd::cxx::tree::ostream< " + *i + " >");

            os << stream_type << "&" << endl
               << "operator<< (" << stream_type << "& s," << endl
               << "const " << name << "&" << (has_body ? " x" : "") << ")"
               << "{";

            if (c.inherits_p ())
            {
              os << "s << static_cast< const ";

              inherits (c, inherits_);

              os << "& > (x);";
            }

            {
              Traversal::Names names_member;
              Element element (*this, name, *i);
              Attribute attribute (*this);

              names_member >> element;
              names_member >> attribute;

              names (c, names_member);
            }

            os << "return s;"
               << "}";


            // Register with type map.
            //
            if (polymorphic && polymorphic_p (c) && !anonymous_p (c))
            {
              // Note that we are using the original type name.
              //
              String const& name (ename (c));

              os << "static" << endl
                 << "const ::xsd::cxx::tree::stream_insertion_initializer< " <<
                "0, " << i->c_str () << ", " <<
                char_type << ", " << name << " >" << endl
                 << "_xsd_" << name << "_stream_insertion_init_" <<
                n++ << " (" << endl
                 << strlit (c.name ()) << "," << endl
                 << strlit (xml_ns_name (c)) << ");"
                 << endl;
            }
          }
        }

      private:
        Traversal::Inherits inherits_;
        BaseTypeName base_;
      };
    }

    Void
    generate_stream_insertion_source (Context& ctx)
    {
      if (ctx.polymorphic)
      {
        Streams const& st (ctx.options.value<CLI::generate_insertion> ());

        ctx.os << "#include <xsd/cxx/tree/stream-insertion-map.hxx>" << endl
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
                 << "{";

          for (Streams::ConstIterator i (st.begin ()); i != st.end (); ++i)
          {
            String stream (*i);

            ctx.os << "#ifdef _MSC_VER" << endl;

            if (export_maps)
              ctx.os << "template struct __declspec (dllexport) " <<
                "stream_insertion_plate< 0, " << stream << ", " <<
                ctx.char_type << " >;";

            if (import_maps)
              ctx.os << "template struct __declspec (dllimport) " <<
                "stream_insertion_plate< 0, " << stream << ", " <<
                ctx.char_type << " >;";

            ctx.os << "#elif defined(__GNUC__) && __GNUC__ >= 4" << endl
                   << "template struct __attribute__ ((visibility(\"default\"))) " <<
              "stream_insertion_plate< 0, " << stream << ", " <<
              ctx.char_type << " >;"
                   << "#elif defined(XSD_MAP_VISIBILITY)" << endl
                   << "template struct XSD_MAP_VISIBILITY " <<
              "stream_insertion_plate< 0, " << stream << ", " <<
              ctx.char_type << " >;"
                   << "#endif" << endl;
          }

          ctx.os << "}"  // tree
                 << "}"  // cxx
                 << "}"  // xsd
                 << "#endif // XSD_NO_EXPORT" << endl
                 << endl;
        }

        ctx.os << "namespace _xsd"
               << "{";

        UnsignedLong n (0);
        for (Streams::ConstIterator i (st.begin ()); i != st.end (); ++i)
        {
          String stream (*i);

          ctx.os << "static" << endl
                 << "const ::xsd::cxx::tree::stream_insertion_plate< 0, " <<
            stream << ", " << ctx.char_type << " >" << endl
                 << "stream_insertion_plate_init_" << n++ << ";";
        }

        ctx.os << "}";
      }

      Traversal::Schema schema;

      Traversal::Sources sources;
      Traversal::Names names_ns, names;

      Namespace ns (ctx);

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
