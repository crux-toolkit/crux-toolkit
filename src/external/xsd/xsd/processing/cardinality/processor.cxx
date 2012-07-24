// file      : processing/cardinality/processor.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <processing/cardinality/processor.hxx>

#include <elements.hxx>

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

#include <cult/containers/map.hxx>

namespace Processing
{
  using namespace Cult;

  namespace SemanticGraph = XSDFrontend::SemanticGraph;
  namespace Traversal = XSDFrontend::Traversal;

  typedef WideString String;

  namespace Cardinality
  {
    namespace
    {
      //
      //
      struct ElementInfo
      {
        ElementInfo ()
            : min (0), max (0), e_ (0)
        {
        }

        ElementInfo (SemanticGraph::Element& e)
            : min (1), max (1), e_ (&e)
        {
        }

        ElementInfo (SemanticGraph::Element& e,
                     UnsignedLong min_, UnsignedLong max_)
            : min (min_), max (max_), e_ (&e)
        {
        }

        SemanticGraph::Element&
        element ()
        {
          assert (e_ != 0);
          return *e_;
        }

      public:
        UnsignedLong min, max;

      private:
        SemanticGraph::Element* e_;
      };

      typedef Cult::Containers::Map<String, ElementInfo> ElementInfoMap;

      //
      //
      struct AnyInfo
      {
        AnyInfo ()
            : min (0), max (0), a_ (0)
        {
        }

        AnyInfo (SemanticGraph::Any& a)
            : min (1), max (1), a_ (&a)
        {
        }

        AnyInfo (SemanticGraph::Any& a,
                 UnsignedLong min_, UnsignedLong max_)
            : min (min_), max (max_), a_ (&a)
        {
        }

        SemanticGraph::Any&
        any ()
        {
          assert (a_ != 0);
          return *a_;
        }

      public:
        UnsignedLong min, max;

      private:
        SemanticGraph::Any* a_;
      };

      typedef Cult::Containers::Map<String, AnyInfo> AnyInfoMap;

      //
      //
      struct Particle: Traversal::All,
                       Traversal::Choice,
                       Traversal::Sequence,
                       Traversal::Element,
                       Traversal::Any
      {
        virtual Void
        traverse (SemanticGraph::All& a)
        {
          traverse_sequence (a);
        }

        virtual Void
        traverse (SemanticGraph::Choice& c)
        {
          using SemanticGraph::Compositor;

          // Go over all particles we contain and add them to the map.
          //
          for (Compositor::ContainsIterator ci (c.contains_begin ());
               ci != c.contains_end (); ++ci)
          {
            Particle t;
            t.dispatch (ci->particle ());

            // Handle elements.
            //
            if (ci == c.contains_begin ())
              el_map = t.el_map; // First arm.
            else
            {
              // For elements that are in the map but not in this
              // arm of choice, we need to set min to 0 while for
              // those that are we need to choose minimum between
              // the two for min and maximum for max.
              //
              for (ElementInfoMap::Iterator i (el_map.begin ());
                   i != el_map.end (); ++i)
              {
                String const& name (i->first);
                ElementInfo& ei (i->second);

                ElementInfoMap::Iterator j (t.el_map.find (name));

                if (j == t.el_map.end ())
                  ei.min = 0;
                else
                {
                  ei.min = j->second.min < ei.min ? j->second.min : ei.min;
                  ei.max = j->second.max > ei.max ? j->second.max : ei.max;
                }
              }

              // Now elements that are in this arm of choice but are
              // not in the map, we need to add to the map and set their
              // min to 0.
              //
              for (ElementInfoMap::Iterator i (t.el_map.begin ());
                   i != t.el_map.end (); ++i)
              {
                String const& name (i->first);
                ElementInfo& ei (i->second);

                ElementInfoMap::Iterator j (el_map.find (name));

                if (j == el_map.end ())
                  el_map[name] = ElementInfo (ei.element (), 0, ei.max);
              }
            }

            // Handle wildcards. Since each wildcard is treated as unique,
            // we need to copy them from each arm of choice and set min to
            // 0.
            //
            for (AnyInfoMap::Iterator i (t.any_map.begin ());
                 i != t.any_map.end (); ++i)
            {
              String const& name (i->first);
              AnyInfo& ai (i->second);

              assert (any_map.find (name) == any_map.end ());

              any_map[name] = AnyInfo (ai.any (), 0, ai.max);
            }
          }

          // Choice's min and max.
          //
          UnsignedLong cmin (c.min ()), cmax (c.max ());

          // Iterate over elements and wildcards in the maps and multiply
          // their cardinality by cmin and cmax.
          //
          for (ElementInfoMap::Iterator i (el_map.begin ());
               i != el_map.end (); ++i)
          {
            i->second.min *= cmin;
            i->second.max *= cmax;
          }

          for (AnyInfoMap::Iterator i (any_map.begin ());
               i != any_map.end (); ++i)
          {
            i->second.min *= cmin; // Not really necessary since min == 0.
            i->second.max *= cmax;
          }
        }

        virtual Void
        traverse (SemanticGraph::Sequence& s)
        {
          traverse_sequence (s);
        }

        Void
        traverse_sequence (SemanticGraph::Compositor& c)
        {
          using SemanticGraph::Compositor;

          // Sequence's min and max.
          //
          UnsignedLong smin (c.min ()), smax (c.max ());

          // Go over all particles we contain and add them to the map.
          //
          for (Compositor::ContainsIterator ci (c.contains_begin ());
               ci != c.contains_end (); ++ci)
          {
            Particle t;
            t.dispatch (ci->particle ());

            // Handle elements.
            //
            for (ElementInfoMap::Iterator i (t.el_map.begin ());
                 i != t.el_map.end (); ++i)
            {
              String const& name (i->first);
              ElementInfo& ei (i->second);
              UnsignedLong min (ei.min * smin);
              UnsignedLong max (ei.max * smax);
              ElementInfoMap::Iterator j (el_map.find (name));

              if (j != el_map.end ())
              {
                // Add i's cardinality to j
                //
                j->second.min += min;
                j->second.max = (j->second.max == 0 || max == 0) ?
                  0 : (j->second.max + max);
              }
              else
                el_map[name] = ElementInfo (ei.element (), min, max);
            }

            // Handle wildcards.
            //
            for (AnyInfoMap::Iterator i (t.any_map.begin ());
                 i != t.any_map.end (); ++i)
            {
              String const& name (i->first);
              AnyInfo& ai (i->second);
              UnsignedLong min (ai.min * smin);
              UnsignedLong max (ai.max * smax);

              assert (any_map.find (name) == any_map.end ());

              any_map[name] = AnyInfo (ai.any (), min, max);
            }
          }
        }

        virtual Void
        traverse (SemanticGraph::Element& e)
        {
          SemanticGraph::ContainsParticle& cp (e.contained_particle ());

          String name (e.qualified_p ()
                       ? e.namespace_ ().name () + L" " + e.name ()
                       : e.name ());

          el_map[name] = ElementInfo (e, cp.min (), cp.max ());
        }

        virtual Void
        traverse (SemanticGraph::Any& a)
        {
          SemanticGraph::ContainsParticle& cp (a.contained_particle ());

          any_map[a.name ()] = AnyInfo (a, cp.min (), cp.max ());
        }

      public:
        AnyInfoMap any_map;
        ElementInfoMap el_map;
      };


      //
      //
      struct Complex: Traversal::Complex
      {
        virtual Void
        traverse (Type& c)
        {
          if (c.contains_compositor_p ())
          {
            Particle t;
            t.dispatch (c.contains_compositor ().compositor ());

            for (ElementInfoMap::Iterator i (t.el_map.begin ());
                 i != t.el_map.end (); ++i)
            {
              ElementInfo& ei (i->second);
              FrontendElements::Context& ctx (ei.element ().context ());

              ctx.set ("min", ei.min);
              ctx.set ("max", ei.max);
            }

            for (AnyInfoMap::Iterator i (t.any_map.begin ());
                 i != t.any_map.end (); ++i)
            {
              AnyInfo& ai (i->second);
              FrontendElements::Context& ctx (ai.any ().context ());

              ctx.set ("min", ai.min);
              ctx.set ("max", ai.max);
            }
          }

          // Traverse attributes and anonymous types (via elements).
          //
          Complex::names (c);
        }
      };


      //
      //
      struct Attribute: Traversal::Attribute
      {
        virtual Void
        traverse (Type& a)
        {
          FrontendElements::Context& ctx (a.context ());

          ctx.set ("min", a.optional_p () ? 0UL : 1UL);
          ctx.set ("max", 1UL);
        }
      };

      // Go into implied/included/imported schemas while making sure
      // we don't process the same stuff more than once.
      //
      struct Uses: Traversal::Uses
      {
        virtual Void
        traverse (Type& u)
        {
          SemanticGraph::Schema& s (u.schema ());

          if (!s.context ().count ("processing-cardinality-seen"))
          {
            s.context ().set ("processing-cardinality-seen", true);
            Traversal::Uses::traverse (u);
          }
        }
      };
    }

    Void Processor::
    process (SemanticGraph::Schema& tu, SemanticGraph::Path const&)
    {
      Traversal::Schema schema;
      Uses uses;

      schema >> uses >> schema;

      Traversal::Names schema_names;
      Traversal::Namespace ns;
      Traversal::Names ns_names;

      schema >> schema_names >> ns >> ns_names;

      Complex complex_type;
      AnonymousBase anonymous (complex_type);

      ns_names >> complex_type;
      ns_names >> anonymous;

      Attribute attribute;
      Traversal::Names names;

      complex_type >> names;

      names >> attribute;
      names >> anonymous;

      // Some twisted schemas do recusive inclusions.
      //
      tu.context ().set ("processing-cardinality-seen", true);

      schema.dispatch (tu);
    }
  }
}
