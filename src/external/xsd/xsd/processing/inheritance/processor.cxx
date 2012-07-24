// file      : processing/inheritance/processor.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <processing/inheritance/processor.hxx>

#include <elements.hxx>

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

#include <cult/containers/set.hxx>

#include <iostream>
using std::wcerr;
using std::endl;

namespace Processing
{
  using namespace Cult;

  namespace SemanticGraph = XSDFrontend::SemanticGraph;
  namespace Traversal = XSDFrontend::Traversal;

  typedef WideString String;

  namespace Inheritance
  {
    namespace
    {
      struct Dep
      {
        Dep (SemanticGraph::Type& t,
             SemanticGraph::Member* m = 0,
             String const& xpath = L"")
            : type (t), member (m), member_xpath (xpath)
        {
        }

        SemanticGraph::Type& type;
        SemanticGraph::Member* member; // Member if type is anonymous.
        String member_xpath;
      };

      inline Boolean
      operator< (Dep const& a, Dep const& b)
      {
        return &a.type < &b.type;
      }

      typedef Containers::Set<Dep> DepSet;
      typedef Containers::Set<SemanticGraph::Type*> TypeSet;


      String
      xpath (SemanticGraph::Nameable& n)
      {
        if (dynamic_cast<SemanticGraph::Namespace*> (&n) != 0)
          return L"<namespace-level>"; // There is a bug if you see this.

        if (n.named_p ())
        {
          SemanticGraph::Scope& scope (n.scope ());

          if (dynamic_cast<SemanticGraph::Namespace*> (&scope) != 0)
            return n.name ();

          return xpath (scope) + L"/" + n.name ();
        }
        else
        {
          return L"(anonymous type for " +
            n.context ().get<String> ("instance-name") + L")";
        }
      }


      // Calculate the list of dependencies for this complex
      // type.
      //
      struct ComplexType: Traversal::Complex,
                          Traversal::Member
      {
        ComplexType (DepSet& dep_set)
            : dep_set_ (dep_set), last_ (0)
        {
          *this >> names_ >> *this;
        }

        virtual Void
        traverse (SemanticGraph::Complex& c)
        {
          using SemanticGraph::Complex;

          if (c.inherits_p ())
            dep_set_.insert (Dep (c.inherits ().base (), last_, last_xpath_));

          types_seen_.insert (&c);

          // Go after anonymous types.
          //
          names (c);
        }

        virtual Void
        traverse (SemanticGraph::Member& m)
        {
          SemanticGraph::Type& t (m.type ());

          if (!t.named_p () && types_seen_.find (&t) == types_seen_.end ())
          {
            FrontendElements::Context& ctx (t.context ());

            last_xpath_ = xpath (m);

            String prev_xpath;

            if (ctx.count ("instance-name"))
              prev_xpath = ctx.get<String> ("instance-name");

            ctx.set ("instance-name", last_xpath_);

            last_ = &m;
            dispatch (t);

            if (prev_xpath)
              ctx.set ("instance-name", prev_xpath);
            else
              ctx.remove ("instance-name");
          }
        }

      private:
        DepSet& dep_set_;
        TypeSet types_seen_;

        SemanticGraph::Member* last_;
        String last_xpath_;

        Traversal::Names names_;
      };


      //
      //
      template <typename N, typename A>
      struct NodeArgs
      {
        NodeArgs (N& node, A arg)
            : node_ (node), arg_ (arg)
        {
        }

        operator N& () const
        {
          return node_;
        }

        template <typename E>
        Void
        add_edge_left (E& e)
        {
          node_.add_edge_left (e, arg_);
        }

        template <typename E>
        Void
        add_edge_right (E& e)
        {
          node_.add_edge_right (e, arg_);
        }

      private:
        N& node_;
        A arg_;
      };


      //
      //
      struct Global: Traversal::Type,
                     Traversal::Complex,
                     Traversal::Element
      {
        Global (SemanticGraph::Schema& root,
                SemanticGraph::Schema& schema,
                Boolean& failed)
            : root_ (root), schema_ (schema), failed_ (failed)
        {
        }

        virtual Void
        traverse (SemanticGraph::Type& t)
        {
          if (t.named_p ())
            types_seen_.insert (&t);
        }

        virtual Void
        traverse (SemanticGraph::Complex& c)
        {
          check_dep (c, c);
          types_seen_.insert (&c);
        };

        virtual Void
        traverse (SemanticGraph::Element& e)
        {
          SemanticGraph::Type& t (e.type ());

          if (!t.named_p ())
          {
            t.context ().set ("instance-name", xpath (e));
            check_dep (e, t);
            t.context ().remove ("instance-name");
          }
        };

      private:
        Void
        check_dep (SemanticGraph::Nameable& global,
                   SemanticGraph::Type& type)
        {
          using SemanticGraph::Type;
          using SemanticGraph::Scope;
          using SemanticGraph::Names;
          using SemanticGraph::Schema;

          DepSet prereqs;

          // Calculate our prerequisistes.
          //
          {
            ComplexType complex (prereqs);
            complex.dispatch (type);
          }

          for (DepSet::ConstIterator i (prereqs.begin ());
               i != prereqs.end (); ++i)
          {
            Dep const& dep (*i);
            Type& t (dep.type);

            // We won't be able to generate compilable code in case of a
            // dependency on ourselves (e.g., a member element with
            // anonymous type that inherits from us).
            //
            if (&t == &type)
            {
              assert (dep.member != 0);

              SemanticGraph::Member& m (*dep.member);

              wcerr << t.file () << ":" << t.line () << ":" << t.column ()
                    << " error: nested anonymous type for '"
                    << dep.member_xpath << "' cyclicly inherits from '"
                    << t.name () << "'" << endl;

              wcerr << t.file () << ":" << t.line () << ":" << t.column ()
                    << " error: unable to generate valid code for such "
                    << "cyclic inheritance" << endl;

              wcerr << m.file () << ":" << m.line () << ":" << m.column ()
                    << " info: '" << m.name () << "' element is declared here"
                    << endl;

              wcerr << t.file () << ":" << t.line () << ":" << t.column ()
                    << ": info: consider explicitly naming this type "
                    << "or remove the --preserve-anonymous option"
                    << endl;

              failed_ = true;
              continue;
            }

            if (types_seen_.find (&t) == types_seen_.end ())
            {
              Scope& scope (t.scope ());
              Schema& schema (dynamic_cast<Schema&> (scope.scope ()));

              // Don't worry about types that are in included/imported
              // schemas.
              //
              if (&schema != &schema_ && !sources_p (schema_, schema))
                continue;

              if (t.context ().count ("seen"))
              {
                wcerr << t.file () << ":" << t.line () << ":" << t.column ()
                      << " error: nested anonymous type in '" << t.name ()
                      << "' or '" << type.name () << "' inherits from one of "
                      << "these types and makes them mutually dependant"
                      << endl;

                wcerr << t.file () << ":" << t.line () << ":" << t.column ()
                      << " error: unable to generate valid code for such "
                      << "cyclic dependency" << endl;

                wcerr << type.file () << ":" << type.line () << ":"
                      << type.column () << " info: '" << type.name ()
                      << "' type is defined here"
                      << endl;

                wcerr << t.file () << ":" << t.line () << ":" << t.column ()
                      << ": info: consider explicitly naming the anonymous "
                      << "type or remove the --preserve-anonymous option"
                      << endl;

                failed_ = true;
                continue;
              }


              //wcerr << "type '" << t.name () << "' needs to be moved " <<
              //  "before " << (global.is_a<Type> () ? "type" : "element") <<
              //  " '" << global.name () << "'" << endl;


              // Delete current Names edge.
              //
              String name (t.name ());
              {
                Names& n (t.named ());
                root_.delete_edge (scope, t, n);
              }

              // Insert a new Names edge before global.
              //
              {
                // Convert to the insert-after call.
                //
                Scope::NamesIterator i (scope.find (global.named ()));

                if (i == scope.names_begin ())
                  i = scope.names_end ();
                else
                  --i;

                NodeArgs<Scope, Scope::NamesIterator> na (scope, i);
                root_.new_edge<Names> (na, t, name);
              }

              // Recursively process the moved type.
              //
              global.context ().set ("seen", true);
              dispatch (t);
              global.context ().remove ("seen");
            }
          }
        }

      private:
        // Return true if root sources s.
        //
        Boolean
        sources_p (SemanticGraph::Schema& root, SemanticGraph::Schema& s)
        {
          using SemanticGraph::Schema;
          using SemanticGraph::Sources;

          for (Schema::UsesIterator i (root.uses_begin ());
               i != root.uses_end (); ++i)
          {
            if (i->is_a<Sources> ())
            {
              if (&i->schema () == &s || sources_p (i->schema (), s))
                return true;
            }
          }

          return false;
        }

      private:
        SemanticGraph::Schema& root_;
        SemanticGraph::Schema& schema_;
        TypeSet types_seen_;
        Boolean& failed_;
      };


      // Go into included/imported schemas while making sure we don't
      // process the same stuff more than once.
      //
      struct Uses: Traversal::Includes, Traversal::Imports
      {
        Uses (SemanticGraph::Schema& root, Boolean& failed)
            : root_ (root), failed_ (failed)
        {
        }

        virtual Void
        traverse (SemanticGraph::Includes& i)
        {
          traverse (i.schema ());
        }

        virtual Void
        traverse (SemanticGraph::Imports& i)
        {
          traverse (i.schema ());
        }

      private:
        Void
        traverse (SemanticGraph::Schema& s)
        {
          if (!s.context ().count ("processing-inheritance-seen"))
          {
            Traversal::Schema schema;
            Traversal::Sources sources;

            schema >> sources >> schema;
            schema >> *this;

            Traversal::Names schema_names;
            Traversal::Namespace ns;
            Traversal::Names ns_names;

            schema >> schema_names >> ns >> ns_names;

            Global global (root_, s, failed_);

            ns_names >> global;

            s.context ().set ("processing-inheritance-seen", true);
            schema.dispatch (s);
          }
        }

      private:
        SemanticGraph::Schema& root_;
        Boolean& failed_;
      };
    }

    Void Processor::
    process (SemanticGraph::Schema& tu, SemanticGraph::Path const&)
    {
      Boolean failed (false);

      // We need to process include/imported schemas since other
      // parts of the process, for example, name processors can
      // rely on the order of types in the schema.
      //
      Traversal::Schema schema;
      Traversal::Sources sources;
      Uses uses (tu, failed);

      schema >> sources >> schema;
      schema >> uses;

      Traversal::Names schema_names;
      Traversal::Namespace ns;
      Traversal::Names ns_names;

      schema >> schema_names >> ns >> ns_names;

      Global global (tu, tu, failed);

      ns_names >> global;

      // Some twisted schemas do recusive self-inclusion.
      //
      tu.context ().set ("processing-inheritance-seen", true);

      schema.dispatch (tu);

      if (failed)
        throw Failed ();
    }
  }
}
