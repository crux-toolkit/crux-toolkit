// file      : xsd/cxx/parser/type-processor.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cult/containers/set.hxx>

#include <cxx/parser/elements.hxx>
#include <cxx/parser/type-processor.hxx>

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

namespace CXX
{
  namespace Parser
  {
    TypeProcessor::
    TypeProcessor ()
    {
      // Dummy ctor, helps with long symbols on HP-UX.
    }

    namespace
    {
      //
      //
      struct Type: Traversal::Type
      {
        Type (SemanticGraph::Schema& schema,
              TypeMap::Namespaces& type_map,
              Boolean add_includes)
            : schema_ (schema),
              type_map_ (type_map),
              add_includes_ (add_includes)
        {
        }

        virtual Void
        traverse (SemanticGraph::Type& type)
        {
          using TypeMap::Namespace;
          using TypeMap::Namespaces;

          SemanticGraph::Context& tc (type.context ());

          // There are two situations where we may try to process the
          // same type more than once. The first is when the type is
          // used in several element declarations in the same schema.
          // The second situation only occurs when we are in the file-
          // per-type mode. In this case the type was processed as part
          // of another root schema. In the second case, while the ret
          // and arg types are assumed to be the same, we need to re-
          // match the type in order to add include directives to the
          // new root schema.
          //
          Boolean set (true);

          if (tc.count ("ret-type"))
          {
            SemanticGraph::Schema* s (
              tc.get<SemanticGraph::Schema*> ("root-schema"));

            if (&schema_ == s)
              return;

            set = false;
          }

          SemanticGraph::Namespace& ns (
            dynamic_cast<SemanticGraph::Namespace&> (type.scope ()));

          String ns_name (ns.name ());
          String t_name (type.name ());

          // std::wcerr << "traversing: " << ns_name << "#" << t_name << endl;

          for (Namespaces::ConstIterator n (type_map_.begin ());
               n != type_map_.end (); ++n)
          {
            // Check if the namespace matches.
            //
            Boolean ns_match;

            if (!n->xsd_name ().empty ())
            {
              ns_match = n->xsd_name ().match (ns_name);
            }
            else
              ns_match = ns_name.empty ();

            // std::wcerr << "considering ns expr: " << n->xsd_name () << endl;

            if (ns_match)
            {
              // Namespace matched. See if there is a type that matches.
              //
              for (Namespace::TypesIterator t (n->types_begin ());
                   t != n->types_end (); ++t)
              {
                if (t->xsd_name ().match (t_name))
                {
                  if (set)
                  {
                    // Got a match. See if the namespace has the C++
                    // namespace mapping.
                    //
                    String cxx_ns;

                    if (n->has_cxx_name ())
                    {
                      if (!n->xsd_name ().empty ())
                      {
                        cxx_ns = n->xsd_name ().merge (
                          n->cxx_name (), ns_name, true);
                      }
                      else
                        cxx_ns = n->cxx_name ();

                      cxx_ns += L"::";
                    }

                    // Figure out ret and arg type names.
                    //
                    String ret_type (cxx_ns);

                    ret_type += t->xsd_name ().merge (
                      t->cxx_ret_name (), t_name, true);

                    String arg_type;

                    if (t->cxx_arg_name ())
                    {
                      arg_type = cxx_ns;
                      arg_type += t->xsd_name ().merge (
                        t->cxx_arg_name (), t_name, true);
                    }
                    else
                    {
                      if (ret_type == L"void")
                        arg_type = ret_type;
                      else
                      {
                        WideChar last (ret_type[ret_type.size () - 1]);

                        // If it is already a pointer or reference then use
                        // it as is.
                        //
                        if (last == L'*' || last == L'&')
                          arg_type = ret_type;
                        else
                          arg_type = L"const " + ret_type + L"&";
                      }
                    }

                    tc.set ("ret-type", ret_type);
                    tc.set ("arg-type", arg_type);
                  }

                  tc.set ("root-schema", &schema_);

                  //std::wcerr << t_name << " -> " << ret_type << endl;

                  // See of we need to add any includes to the translations
                  // unit.
                  //
                  if (add_includes_)
                  {
                    if (n->includes_begin () != n->includes_end ())
                    {
                      typedef Cult::Containers::Set<String> Includes;

                      if (!schema_.context ().count ("includes"))
                        schema_.context ().set ("includes", Includes ());

                      Includes& is (
                        schema_.context ().get<Includes> ("includes"));

                      for (Namespace::IncludesIterator i (n->includes_begin ());
                           i != n->includes_end (); ++i)
                      {
                        is.insert (*i);
                      }
                    }
                  }

                  return;
                }
              }
            }
          }
        }

      private:
        SemanticGraph::Schema& schema_;
        TypeMap::Namespaces& type_map_;
        Boolean add_includes_;
      };


      //
      //
      struct GlobalType: Traversal::Type,
                         Traversal::List,
                         Traversal::Complex,
                         Traversal::Enumeration
      {
        GlobalType (SemanticGraph::Schema& schema,
                    TypeMap::Namespaces& type_map,
                    Boolean add_includes)
            : type_ (schema, type_map, add_includes)
        {
          inherits_ >> type_;
          names_ >> instance_ >> belongs_ >> type_;
          argumented_ >> type_;
        }

        virtual Void
        traverse (SemanticGraph::Type& t)
        {
          type_.traverse (t);
        }

        virtual Void
        traverse (SemanticGraph::List& l)
        {
          type_.traverse (l);
          Traversal::List::argumented (l, argumented_);
        }

        virtual Void
        traverse (SemanticGraph::Complex& c)
        {
          type_.traverse (c);
          Complex::inherits (c, inherits_);
          Complex::names (c, names_);
        }

        virtual Void
        traverse (SemanticGraph::Enumeration& e)
        {
          type_.traverse (e);
          Complex::inherits (e, inherits_);
        }

      private:
        Parser::Type type_;
        Traversal::Names names_;
        Traversal::Instance instance_;
        Traversal::Inherits inherits_;
        Traversal::Belongs belongs_;
        Traversal::Argumented argumented_;
      };

      Void
      process_impl (CLI::Options const& options,
                    XSDFrontend::SemanticGraph::Schema& tu,
                    Boolean gen_driver,
                    TypeMap::Namespaces& type_map)
      {
        if (tu.names_begin ()->named ().name () ==
            L"http://www.w3.org/2001/XMLSchema")
        {
          // XML Schema namespace.
          //
          Traversal::Schema schema;

          Traversal::Names schema_names;
          Traversal::Namespace ns;
          Traversal::Names ns_names;
          GlobalType global_type (tu, type_map, true);

          schema >> schema_names >> ns >> ns_names >> global_type;

          schema.dispatch (tu);
        }
        else
        {
          // If --extern-xml-schema is specified, then we don't want
          // includes from the XML Schema type map.
          //
          Boolean extern_xml_schema (
            options.value<CLI::extern_xml_schema> ());

          //
          //
          Traversal::Schema schema;
          Traversal::Schema xs_schema;
          Traversal::Sources sources;
          Traversal::Implies implies;

          schema >> sources >> schema;
          schema >> implies >> xs_schema;

          Traversal::Names schema_names;
          Traversal::Namespace ns;
          Traversal::Names ns_names;
          GlobalType global_type (tu, type_map, true);

          schema >> schema_names >> ns >> ns_names >> global_type;

          Traversal::Names xs_schema_names;
          Traversal::Namespace xs_ns;
          Traversal::Names xs_ns_names;
          GlobalType xs_global_type (tu, type_map, !extern_xml_schema);

          xs_schema >> xs_schema_names >> xs_ns >> xs_ns_names >>
            xs_global_type;

          schema.dispatch (tu);

          // If we are generating the test driver, make sure the root
          // element type is processed.
          //
          if (gen_driver && options.value<CLI::generate_test_driver> ())
          {
            // Figure out the root element. Validator should have made sure
            // it is unique.
            //
            SemanticGraph::Element* root (0);
            {
              Traversal::Schema schema;
              Traversal::Sources sources;

              schema >> sources >> schema;

              Traversal::Names schema_names;
              Traversal::Namespace ns;
              Traversal::Names ns_names;
              RootElement root_element (options, root);

              schema >> schema_names >> ns >> ns_names >> root_element;

              schema.dispatch (tu);
            }

            global_type.dispatch (root->type ());
          }
        }
      }
    }

    Void TypeProcessor::
    process (CLI::Options const& options,
             XSDFrontend::SemanticGraph::Schema& s,
             Boolean gen_driver,
             TypeMap::Namespaces& tm)
    {
      process_impl (options, s, gen_driver, tm);
    }
  }
}
