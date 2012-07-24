// file      : xsd/cxx/tree/tree-forward.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/tree/tree-forward.hxx>
#include <cxx/tree/fundamental-header.hxx>

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
          String const& name (ename (l));

          if (String custom = custom_type (l))
          {
            String new_name;
            renamed_type (l, new_name);

            if (new_name)
              os << "class " << new_name << ";";

            if (custom == name)
              os << "class " << name << ";";
            else
              os << "typedef " << custom << " " << name << ";";
          }
          else
            os << "class " << name << ";";
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
          String const& name (ename (u));

          if (String custom = custom_type (u))
          {
            String new_name;
            renamed_type (u, new_name);

            if (new_name)
              os << "class " << new_name << ";";

            if (custom == name)
              os << "class " << name << ";";
            else
              os << "typedef " << custom << " " << name << ";";
          }
          else
            os << "class " << name << ";";
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
          String const& name (ename (e));

          if (String custom = custom_type (e))
          {
            String new_name;
            renamed_type (e, new_name);

            if (new_name)
              os << "class " << new_name << ";";

            if (custom == name)
              os << "class " << name << ";";
            else
              os << "typedef " << custom << " " << name << ";";
          }
          else
            os << "class " << name << ";";
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
          String const& name (ename (c));

          if (String custom = custom_type (c))
          {
            String new_name;
            renamed_type (c, new_name);

            if (new_name)
              os << "class " << new_name << ";";

            if (custom == name)
              os << "class " << name << ";";
            else
              os << "typedef " << custom << " " << name << ";";
          }
          else
            os << "class " << name << ";";
        }
      };
    }

    Void
    generate_forward (Context& ctx)
    {
      NarrowString xml_schema (ctx.options.value<CLI::extern_xml_schema> ());

      // Inlcude or Emit fundamental types.
      //
      if (xml_schema)
      {
        String name (ctx.hxx_expr->merge (xml_schema));

        ctx.os << "#include " << ctx.process_include_path (name) << endl
               << endl;
      }
      else
      {
        if (ctx.char_type == L"char" && ctx.char_encoding != L"custom")
        {
          ctx.os << "#include <xsd/cxx/xml/char-" << ctx.char_encoding <<
            ".hxx>" << endl
                 << endl;
        }

        ctx.os << "#include <xsd/cxx/tree/exceptions.hxx>" << endl
               << "#include <xsd/cxx/tree/elements.hxx>" << endl
               << "#include <xsd/cxx/tree/types.hxx>" << endl
               << endl;

        if (!ctx.options.value<CLI::suppress_parsing> () ||
            ctx.options.value<CLI::generate_serialization> ())
        {
          ctx.os << "#include <xsd/cxx/xml/error-handler.hxx>" << endl
                 << endl;
        }

        if (!ctx.options.value<CLI::suppress_parsing> () ||
            ctx.options.value<CLI::generate_serialization> ())
        {
          ctx.os << "#include <xsd/cxx/xml/dom/auto-ptr.hxx>" << endl
                 << endl;
        }

        Boolean element_map (ctx.options.value<CLI::generate_element_map> ());

        if (element_map)
          ctx.os << "#include <xsd/cxx/tree/element-map.hxx>" << endl
                 << endl;

        // I need to include all the "optional" headers here (instead of
        // later in the individual generators for each feature because
        // those headers provide implementation for the fundamental types.
        //
        if (!ctx.options.value<CLI::suppress_parsing> ())
        {
          ctx.os << "#include <xsd/cxx/tree/parsing.hxx>" << endl;

          Traversal::Schema schema, xsd;
          Traversal::Implies implies;
          Traversal::Names names;
          Traversal::Namespace ns;
          Traversal::Names ns_names;
          FundIncludes type (ctx, "parsing");

          schema >> implies >> xsd >> names >> ns >> ns_names >> type;

          schema.dispatch (ctx.schema_root);

          if (element_map)
            ctx.os << "#include <xsd/cxx/tree/parsing/element-map.txx>" <<
              endl;

          ctx.os << endl;
        }

        if (ctx.options.value<CLI::generate_serialization> ())
        {
          ctx.os << "#include <xsd/cxx/xml/dom/serialization-header.hxx>" << endl
                 << "#include <xsd/cxx/tree/serialization.hxx>" << endl;

          Traversal::Schema schema, xsd;
          Traversal::Implies implies;
          Traversal::Names names;
          Traversal::Namespace ns;
          Traversal::Names ns_names;
          FundIncludes type (ctx, "serialization");

          schema >> implies >> xsd >> names >> ns >> ns_names >> type;

          schema.dispatch (ctx.schema_root);

          if (element_map)
            ctx.os << "#include <xsd/cxx/tree/serialization/element-map.txx>" <<
              endl;

          ctx.os << endl;
        }

        if (ctx.options.value<CLI::generate_ostream> ())
        {
          ctx.os << "#include <xsd/cxx/tree/std-ostream-operators.hxx>" << endl
                 << endl;
        }

        typedef Containers::Vector<NarrowString> Streams;

        Streams const& ist (ctx.options.value<CLI::generate_insertion> ());
        if (!ist.empty ())
        {
          for (Streams::ConstIterator i (ist.begin ()); i != ist.end (); ++i)
          {
            if (*i == "ACE_OutputCDR")
              ctx.os << "#include <xsd/cxx/tree/ace-cdr-stream-insertion.hxx>"
                     << endl;
            else if (*i == "XDR")
              ctx.os << "#include <xsd/cxx/tree/xdr-stream-insertion.hxx>"
                     << endl;
          }

          ctx.os << "#include <xsd/cxx/tree/stream-insertion.hxx>" << endl
                 << endl;
        }

        Streams const& est (ctx.options.value<CLI::generate_extraction> ());
        if (!est.empty ())
        {
          for (Streams::ConstIterator i (est.begin ()); i != est.end (); ++i)
          {
            if (*i == "ACE_InputCDR")
              ctx.os << "#include <xsd/cxx/tree/ace-cdr-stream-extraction.hxx>"
                     << endl;
            else if (*i == "XDR")
              ctx.os << "#include <xsd/cxx/tree/xdr-stream-extraction.hxx>"
                     << endl;
          }

          ctx.os << "#include <xsd/cxx/tree/stream-extraction.hxx>" << endl
                 << endl;
        }


        Traversal::Schema schema, xsd;
        Traversal::Implies implies;
        Traversal::Names names;
        FundamentalNamespace ns (ctx);

        schema >> implies >> xsd >> names >> ns;

        schema.dispatch (ctx.schema_root);
      }

      // First emit header includes.
      //
      if (ctx.options.value<CLI::generate_forward> ())
      {
        Traversal::Schema schema;
        Includes includes (ctx, Includes::forward);

        schema >> includes;

        schema.dispatch (ctx.schema_root);
      }

      ctx.os << "// Forward declarations." << endl
             << "//" << endl;

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

      ctx.os << endl;
    }
  }
}
