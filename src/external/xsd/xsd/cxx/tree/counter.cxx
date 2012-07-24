// file      : xsd/cxx/tree/counter.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/tree/counter.hxx>

#include <iostream>

namespace CXX
{
  namespace Tree
  {
    namespace
    {
      struct Member: Traversal::Member
      {
        Member (UnsignedLong& complexity)
            : complexity_ (complexity)
        {
        }

        virtual Void
        traverse (Type&)
        {
          complexity_++;
        }

        UnsignedLong& complexity_;
      };

      struct Any: Traversal::Any, Traversal::AnyAttribute
      {
        Any (UnsignedLong& complexity)
            : complexity_ (complexity)
        {
        }

        virtual Void
        traverse (SemanticGraph::Any&)
        {
          complexity_++;
        }

        virtual Void
        traverse (SemanticGraph::AnyAttribute&)
        {
          complexity_++;
        }

        UnsignedLong& complexity_;
      };

      struct TypeBase: Traversal::List,
                       Traversal::Union,
                       Traversal::Enumeration,
                       Traversal::Complex,
                       Context
      {
        TypeBase (Context& c, UnsignedLong& complexity)
            : Context (c), complexity_ (complexity)
        {
        }

        virtual Void
        traverse (SemanticGraph::List&)
        {
          complexity_++;
        }

        virtual Void
        traverse (SemanticGraph::Union&)
        {
          complexity_++;
        }

        virtual Void
        traverse (SemanticGraph::Enumeration& e)
        {
          Boolean string_based (false);
          {
            IsStringBasedType t (string_based);
            t.dispatch (e);
          }

          complexity_ += (string_based ? 1 : 2);
        }

        virtual Void
        traverse (SemanticGraph::Complex& c)
        {
          complexity_++; // One for the type itself.

          // Plus some for each member.
          //
          Any any (complexity_);
          Member member (complexity_);
          Traversal::Names names;

          names >> member;

          if (options.value<CLI::generate_wildcard> ())
            names >> any;

          Complex::names (c, names);
        }

      private:
        UnsignedLong& complexity_;
      };


      //
      //
      struct GlobalType: Traversal::Type, Context

      {
        GlobalType (Context& c, Counts& counts)
            : Context (c), counts_ (counts)
        {
        }

        virtual Void
        traverse (SemanticGraph::Type& t)
        {
          counts_.global_types++;

          UnsignedLong complexity (0);
          TypeBase type (*this, complexity);
          type.dispatch (t);

          counts_.complexity_total += complexity;
          counts_.complexity.push_back (complexity);
        }

      private:
        Counts& counts_;
      };

      //
      //
      struct GlobalElement: Traversal::Element,
                            GlobalElementBase,
                            Context
      {
        GlobalElement (Context& c, Counts& counts)
            : GlobalElementBase (c),
              Context (c),
              counts_ (counts),
              last_ (0)
        {
        }

        ~GlobalElement ()
        {
          if (last_ != 0)
          {
            last_->context ().set ("last", true);
            count_last ();
          }
        }

        virtual Void
        traverse (Type& e)
        {
          // Check if the previous element we saw needs to be generated.
          //
          if (last_ != 0)
            count_last ();

          last_ = &e;

          if (counts_.global_elements == 0)
            e.context ().set ("first", true);

          counts_.global_elements++;
        }

      private:
        Void
        count_last ()
        {
          if (generate_p (*last_))
          {
            counts_.generated_global_elements++;

            UnsignedLong complexity (0);

            if (doc_root_p (*last_))
            {
              if (options.value<CLI::generate_element_type> ())
              {
                complexity += 1; // For c-tors and d-tor.

                if (!options.value<CLI::suppress_parsing> ())
                  complexity += 1;

                if (options.value<CLI::generate_serialization> ())
                  complexity += 1;
              }
              else
              {
                if (!options.value<CLI::suppress_parsing> ())
                  complexity += 6; // 13 parsing functions.

                if (options.value<CLI::generate_serialization> ())
                  complexity += 4; // 8 serialization functions.
              }
            }

            if (complexity == 0)
            {
              // This element must be a substitution group members. For
              // such elements we are only generating an entry in a map.
              // We will assign it a complexity of 1 so that we don't
              // end up with the total complexity that is less than the
              // number of elements and types.
              //
              complexity = 1;
            }

            counts_.complexity_total += complexity;
            counts_.complexity.push_back (complexity);
          }
        }

      private:
        Counts& counts_;
        SemanticGraph::Element* last_;
      };
    }

    Counter::
    Counter ()
    {
    }

    Counts Counter::
    count (CLI::Options const& options,
           SemanticGraph::Schema& tu,
           SemanticGraph::Path const& path)
    {
      Counts counts;
      Context ctx (std::wcerr, tu, path, options, counts, false, 0, 0, 0, 0);

      Traversal::Schema schema;
      Traversal::Sources sources;

      schema >> sources >> schema;

      Traversal::Names schema_names;
      Traversal::Namespace ns;
      Traversal::Names ns_names;
      GlobalType global_type (ctx, counts);
      GlobalElement global_element (ctx, counts);

      schema >> schema_names >> ns >> ns_names;

      ns_names >> global_element;
      ns_names >> global_type;

      schema.dispatch (tu);

      return counts;
    }
  }
}
