// file      : xsd/cxx/tree/tree-inline.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

#include <cxx/tree/tree-inline.hxx>
#include <cxx/tree/default-value.hxx>

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
          String item_name (item_type_name (item_type));
          String base_type (L"::xsd::cxx::tree::list< " + item_name +
                            L", " + char_type);

          if (item_type.is_a<SemanticGraph::Fundamental::Double> ())
            base_type += L", ::xsd::cxx::tree::schema_type::double_";
          else if (item_type.is_a<SemanticGraph::Fundamental::Decimal> ())
            base_type += L", ::xsd::cxx::tree::schema_type::decimal";

          base_type += L" >";

          os << "// " << name << endl
             << "//" << endl
             << endl;

          // c-tor ()
          //
          os << inl
             << name << "::" << endl
             << name << " ()" << endl
             << ": " << base_type << " ( " << flags_type << " (0), this)"
             << "{"
             << "}";

          // c-tor (size_type, const X& x)
          //
          String size_type (name != L"size_type"
                            ? String (L"size_type")
                            : base_type + L"::size_type");

          os << inl
             << name << "::" << endl
             << name << " (" << size_type << " n, const " << item_name <<
            "& x)" << endl
             << ": " << base_type << " (n, x, this)"
             << "{"
             << "}";

          // copy c-tor ()
          //
          os << inl
             << name << "::" << endl
             << name << " (const " << name << "& o," << endl
             << flags_type << " f," << endl
             << container << "* c)" << endl
             << ": " << any_simple_type << " (o, f, c)," << endl
             << "  " << base_type << " (o, f, this)"
             << "{"
             << "}";
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

          os << "// " << name << endl
             << "//" << endl
             << endl;

          if (options.value<CLI::generate_default_ctor> ())
          {
            // c-tor ()
            //
            os << inl
               << name << "::" << endl
               << name << " ()" << endl
               << ": " << base << " ()"
               << "{"
               << "}";
          }

          // c-tor (const char*)
          //
          os << inl
             << name << "::" << endl
             << name << " (const " << char_type << "* s)" << endl
             << ": " << base << " (s)"
             << "{"
             << "}";

          // c-tor (string const&)
          //
          os << inl
             << name << "::" << endl
             << name << " (const " << string_type << "& s)" << endl
             << ": " << base << " (s)"
             << "{"
             << "}";

          // copy c-tor ()
          //
          os << inl
             << name << "::" << endl
             << name << " (const " << name << "& o," << endl
             << flags_type << " f," << endl
             << container << "* c)" << endl
             << ": " << base << " (o, f, c)"
             << "{"
             << "}";
        }
      };

      // Enumeration
      //

      // Generate a sequence of explicit c-tor calls until we reach
      // one of the fundamental string types that can be constructed
      // from char literals.
      //
      struct CtorCallSequence: Traversal::Complex,
                               Traversal::Fundamental::Type,
                               Context
      {
        CtorCallSequence (Context& c, String const& arg)
            : Context (c), arg_ (arg), base_type_name_ (c)
        {
        }

        virtual Void
        traverse (SemanticGraph::Complex& c)
        {
          // This type should be ultimately string based.
          //
          assert (c.inherits_p ());

          os << ename (c) << " (" << endl;

          dispatch (c.inherits ().base ());

          os << ")";
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Type& t)
        {
          base_type_name_.dispatch (t);

          os << " (" << arg_ << ")";
        }

      private:
        String arg_;
        BaseTypeName base_type_name_;
      };

      struct Enumeration: Traversal::Enumeration, Context
      {
        Enumeration (Context& c)
            : Context (c), member_ (c)
        {
          inherits_member_ >> member_;
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

          String value;
          if (string_based)
            value = evalue (e);

          // Get to the ultimate base and see if is a fundamental type.
          //
          Boolean fund_based (false);
          SemanticGraph::Type& ult_base (ultimate_base (e));
          {
            IsFundamentalType t (fund_based);
            t.dispatch (ult_base);
          }

          //
          //
          String base; // base type name
          {
            std::wostringstream o;

            BaseTypeName base_type (*this, o);
            Traversal::Inherits inherits_type (base_type);

            inherits (e, inherits_type);
            base = o.str ();
          }

          os << "// " << name << endl
             << "// " << endl
             << endl;

          // default c-tor
          //
          if (options.value<CLI::generate_default_ctor> ())
          {
            os << inl
               << name << "::" <<  endl
               << name << " ()" << endl
               << ": " << base << " ()"
               << "{"
               << "}";
          }

          // c-tor (value)
          //
          if (string_based)
          {
            os << inl
               << name << "::" <<  endl
               << name << " (" << value << " v)" << endl
               << ": ";

            // If we are enum-based then we can just use the corresponding
            // base c-tor directly. Otherwise we will use the from-string
            // c-tor.
            //
            if (enum_based)
              os << base << " (v)";
            else
            {
              CtorCallSequence t (*this, L"_xsd_" + name + L"_literals_[v]");
              t.dispatch (e.inherits ().base ());
            }

            os << "{"
               << "}";
          }

          // c-tor (const char*)
          //
          if (string_based)
          {
            os << inl
               << name << "::" <<  endl
               << name << " (const " << char_type << "* v)" << endl
               << ": " << base << " (v)"
               << "{"
               << "}";
          }

          // c-tor (const std::string&)
          //
          if (string_based)
          {
            os << inl
               << name << "::" <<  endl
               << name << " (const " << string_type << "& v)" << endl
               << ": " << base << " (v)"
               << "{"
               << "}";
          }

          // c-tor (fundamental)
          //
          if (fund_based)
          {
            os << inl
               << name << "::" <<  endl
               << name << " (";

            member_.dispatch (ult_base);

            os << " v)"
               << ": " << base << " (v)"
               << "{"
               << "}";
          }

          // c-tor (base)
          //
          // If the ultimate is also our immediate base and it is a
          // fundamental type then this c-tor clashes with c-tor
          // (fundamental) above.
          //
          if (!fund_based || &ult_base != &e.inherits ().base ())
          {
            os << inl
               << name << "::" <<  endl
               << name << " (const ";

          inherits (e, inherits_member_);

          os << "& v)" << endl
             << ": " << base << " (v)"
             << "{"
             << "}";
          }

          // copy c-tor
          //
          os << inl
             << name << "::" <<  endl
             << name << " (const " << name << "& v," << endl
             << flags_type << " f," << endl
             << container << "* c)" << endl
             << ": " << base << " (v, f, c)"
             << "{"
             << "}";

          // operato= (value)
          //
          if (string_based)
          {
            os << inl
               << name << "& " << name << "::" << endl
               << "operator= (" << value << " v)"
               << "{"
               << "static_cast< " << base << "& > (*this) = ";

            // If we are enum-based then we can just use the corresponding
            // base assignment directly. Otherwise we will use the from-
            // base assignment and a from-string base c-tor.
            //
            if (enum_based)
              os << "v";
            else
            {
              os << endl;

              CtorCallSequence t (*this, L"_xsd_" + name + L"_literals_[v]");
              t.dispatch (e.inherits ().base ());
            }

            os << ";"
               << endl
               << "return *this;"
               << "}";
          }

          os << endl;
        }

      private:
        Traversal::Inherits inherits_member_;
        MemberTypeName member_;
      };

      struct Member: Traversal::Member, Context
      {
        Member (Context& c, String const& scope)
            : Context (c), scope_ (scope), lit_value_ (c)
        {
        }

        virtual Void
        traverse (Type& m)
        {
          if (skip (m))
            return;

          String const& aname (eaname (m));
          String const& mname (emname (m));
          String const& member (emember (m));

          Boolean fund (false);
          {
            IsFundamentalType t (fund);
            t.dispatch (m.type ());
          }

          Boolean def_attr (m.default_p () &&
                            m.is_a<SemanticGraph::Attribute> ());

          if (max (m) != 1)
          {
            // sequence
            //
            String container (econtainer (m));
            String q_container (scope_ + L"::" + container);

            // container const&
            // name () const;
            //
            os << inl
               << "const " << q_container << "& " << scope_ << "::" << endl
               << aname << " () const"
               << "{"
               << "return this->" << member << ";"
               << "}";

            // container&
            // name ();
            //
            os << inl
               << q_container << "& " << scope_ << "::" << endl
               << aname << " ()"
               << "{"
               << "return this->" << member << ";"
               << "}";

            // void
            // name (container const&);
            //
            os << inl
               << "void " << scope_ << "::" << endl
               << mname << " (const " << container << "& s)"
               << "{"
               << "this->" << member << " = s;"
               << "}";
          }
          else if (min (m) == 0 && !def_attr)
          {
            // optional
            //
            String type (etype (m));
            String container (econtainer (m));
            String q_container (scope_ + L"::" + container);

            // container const&
            // name () const;
            //
            os << inl
               << "const " << q_container << "& " << scope_ << "::" << endl
               << aname << " () const"
               << "{"
               << "return this->" << member << ";"
               << "}";

            // container&
            // name ();
            //
            os << inl
               << q_container << "& " << scope_ << "::" << endl
               << aname << " ()"
               << "{"
               << "return this->" << member << ";"
               << "}";

            // void
            // name (type const&);
            //
            os << inl
               << "void " << scope_ << "::" << endl
               << mname << " (const " << type << "& x)"
               << "{"
               << "this->" << member << ".set (x);"
               << "}";

            // void
            // name (container const&);
            //
            os << inl
               << "void " << scope_ << "::" << endl
               << mname << " (const " << container << "& x)"
               << "{"
               << "this->" << member << " = x;"
               << "}";

            // void
            // name (auto_ptr<type>);
            //
            if (!fund)
              os << inl
                 << "void " << scope_ << "::" << endl
                 << mname << " (::std::auto_ptr< " << type << " > x)"
                 << "{"
                 << "this->" << member << ".set (x);"
                 << "}";
          }
          else
          {
            // one
            //
            String type (etype (m));
            String q_type (scope_ + L"::" + type);

            // type const&
            // name () const;
            //
            os << inl
               << "const " << q_type << "& " << scope_ << "::" << endl
               << aname << " () const"
               << "{"
               << "return this->" << member << ".get ();"
               << "}";

            // Do not generate modifiers for fixed attributes.
            //
            if (!(def_attr && m.fixed_p ()))
            {
              // type&
              // name ();
              //
              os << inl
                 << q_type << "& " << scope_ << "::" << endl
                 << aname << " ()"
                 << "{"
                 << "return this->" << member << ".get ();"
                 << "}";

              // void
              // name (type const&);
              //
              os << inl
                 << "void " << scope_ << "::" << endl
                 << mname << " (const " << type << "& x)"
                 << "{"
                 << "this->" << member << ".set (x);"
                 << "}";

              // void
              // name (auto_ptr<type>);
              //
              if (!fund)
                os << inl
                   << "void " << scope_ << "::" << endl
                   << mname << " (::std::auto_ptr< " << type << " > x)"
                   << "{"
                   << "this->" << member << ".set (x);"
                   << "}";

              // auto_ptr<type>
              // detach_name ();
              //
              if (detach && !fund)
                os << inl
                   << "::std::auto_ptr< " << q_type << " > " <<
                  scope_ << "::" << endl
                   << edname (m) << " ()"
                   << "{"
                   << "return this->" << member << ".detach ();"
                   << "}";
            }
          }


          // default_value
          //
          if (m.default_p ())
          {
            Boolean simple (true);

            if (m.is_a<SemanticGraph::Element> ())
            {
              IsSimpleType test (simple);
              test.dispatch (m.type ());
            }

            if (simple)
            {
              String lit (lit_value_.dispatch (m.type (), m.value ()));

              os << inl;

              if (lit)
                os << scope_ << "::" << etype (m) << " ";
              else
                os << "const " << scope_ << "::" << etype (m) << "& ";

              os << scope_ << "::" << endl
                 << edefault_value (m) << " ()"
                 << "{";

              if (lit)
                os << "return " << etype (m) << " (" << lit << ");";
              else
                os << "return " << edefault_value_member (m) << ";";

              os << "}";
            }
          }
        }

      private:
        String scope_;
        LiteralValue lit_value_;
      };

      struct Any: Traversal::Any,
                  Traversal::AnyAttribute,
                  Context
      {
        Any (Context& c, String const& scope)
            : Context (c), scope_ (scope)
        {
        }

        virtual Void
        traverse (SemanticGraph::Any& a)
        {
          String const& aname (eaname (a));
          String const& mname (emname (a));

          String const& member (emember (a));

          if (max (a) != 1)
          {
            // sequence
            //
            String container (econtainer (a));
            String q_container (scope_ + L"::" + container);

            // container const&
            // name () const;
            //
            os << inl
               << "const " << q_container << "& " << scope_ << "::" << endl
               << aname << " () const"
               << "{"
               << "return this->" << member << ";"
               << "}";

            // container&
            // name ();
            //
            os << inl
               << q_container << "& " << scope_ << "::" << endl
               << aname << " ()"
               << "{"
               << "return this->" << member << ";"
               << "}";

            // void
            // name (container const&);
            //
            os << inl
               << "void " << scope_ << "::" << endl
               << mname << " (const " << container << "& s)"
               << "{"
               << "this->" << member << " = s;"
               << "}";
          }
          else if (min (a) == 0)
          {
            // optional
            //
            String container (econtainer (a));
            String q_container (scope_ + L"::" + container);

            // container const&
            // name () const;
            //
            os << inl
               << "const " << q_container << "& " << scope_ << "::" << endl
               << aname << " () const"
               << "{"
               << "return this->" << member << ";"
               << "}";

            // container&
            // name ();
            //
            os << inl
               << q_container << "& " << scope_ << "::" << endl
               << aname << " ()"
               << "{"
               << "return this->" << member << ";"
               << "}";

            // void
            // name (type const&);
            //
            os << inl
               << "void " << scope_ << "::" << endl
               << mname << " (const " << xerces_ns << "::DOMElement& e)"
               << "{"
               << "this->" << member << ".set (e);"
               << "}";

            // void
            // name (type*);
            //
            os << inl
               << "void " << scope_ << "::" << endl
               << mname << " (" << xerces_ns << "::DOMElement* e)"
               << "{"
               << "this->" << member << ".set (e);"
               << "}";

            // void
            // name (container const&);
            //
            os << inl
               << "void " << scope_ << "::" << endl
               << mname << " (const " << container << "& x)"
               << "{"
               << "this->" << member << " = x;"
               << "}";
          }
          else
          {
            // one
            //

            // type const&
            // name () const;
            //
            os << inl
               << "const " << xerces_ns << "::DOMElement& " <<
              scope_ << "::" << endl
               << aname << " () const"
               << "{"
               << "return this->" << member << ".get ();"
               << "}";

            // type&
            // name ();
            //
            os << inl
               << xerces_ns << "::DOMElement& " << scope_ << "::" << endl
               << aname << " ()"
               << "{"
               << "return this->" << member << ".get ();"
               << "}";

            // void
            // name (type const&);
            //
            os << inl
               << "void " << scope_ << "::" << endl
               << mname << " (const " << xerces_ns << "::DOMElement& e)"
               << "{"
               << "this->" << member << ".set (e);"
               << "}";

            // void
            // name (type*);
            //
            os << inl
               << "void " << scope_ << "::" << endl
               << mname << " (" << xerces_ns << "::DOMElement* e)"
               << "{"
               << "this->" << member << ".set (e);"
               << "}";
          }
        }

        virtual Void
        traverse (SemanticGraph::AnyAttribute& a)
        {
          String const& aname (eaname (a));
          String const& mname (emname (a));

          String const& member (emember (a));
          String container (econtainer (a));
          String q_container (scope_ + L"::" + container);

          // container const&
          // name () const;
          //
          os << inl
             << "const " << q_container << "& " << scope_ << "::" << endl
             << aname << " () const"
             << "{"
             << "return this->" << member << ";"
             << "}";

          // container&
          // name ();
          //
          os << inl
             << q_container << "& " << scope_ << "::" << endl
             << aname << " ()"
             << "{"
             << "return this->" << member << ";"
             << "}";

          // void
          // name (container const&);
          //
          os << inl
             << "void " << scope_ << "::" << endl
             << mname << " (const " << container << "& s)"
             << "{"
             << "this->" << member << " = s;"
             << "}";
        }

      private:
        String scope_;
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

          os << "// " << name << endl
             << "// " << endl
             << endl;

          // Generate accessors and modifiers.
          //
          Any any (*this, name);
          Member member (*this, name);
          Traversal::Names names;

          if (options.value<CLI::generate_wildcard> ())
            names >> any;

          names >> member;

          Complex::names (c, names);

          // dom_document accessors.
          //
          if (edom_document_member_p (c))
          {
            os << inl
               << "const " << xerces_ns << "::DOMDocument& " <<
              name << "::" << endl
               << edom_document (c) << " () const"
               << "{"
               << "return *" << edom_document_member (c) << ";"
               << "}";

            os << inl
               << xerces_ns << "::DOMDocument& " << name << "::" << endl
               << edom_document (c) << " ()"
               << "{"
               << "return *" << edom_document_member (c) << ";"
               << "}";
          }

          os << endl;
        }
      };


      struct GlobalElement: Traversal::Element,
                            GlobalElementBase,
                            Context
      {
        GlobalElement (Context& c)
            : GlobalElementBase (c), Context (c)
        {
        }

        virtual Void
        traverse (Type& e)
        {
          if (!doc_root_p (e))
            return;

          Boolean fund (false);
          {
            IsFundamentalType test (fund);
            test.dispatch (e.type ());
          }

          String const& name (ename (e));

          os << "// " << name << endl
             << "// " << endl
             << endl;

          // Accessors/modifiers.
          //
          String const& type (etype (e));
          String const& aname (eaname (e));
          String const& mname (emname (e));
          String const& member (emember (e));

          // type const&
          // name () const;
          //
          os << inl
             << "const " << name << "::" << type << "& " << name << "::" << endl
             << aname << " () const"
             << "{"
             << "return this->" << member << ".get ();"
             << "}";

          // type&
          // name ();
          //
          os << inl
             << name << "::" << type << "& " << name << "::" << endl
             << aname << " ()"
             << "{"
             << "return this->" << member << ".get ();"
             << "}";

          // void
          // name (type const&);
          //
          os << inl
             << "void " << name << "::" << endl
             << mname << " (const " << type << "& x)"
             << "{"
             << "return this->" << member << ".set (x);"
             << "}";

          // void
          // name (auto_ptr<type>);
          //
          if (!fund)
          {
            os << inl
               << "void " << name << "::" << endl
               << mname << " (::std::auto_ptr< " << type << " > p)"
               << "{"
               << "return this->" << member << ".set (p);"
               << "}";
          }

          // auto_ptr<type>
          // detach_name ();
          //
          if (detach && !fund)
            os << inl
               << "::std::auto_ptr< " << name << "::" << type << " > " <<
              name << "::" << endl
               << edname (e) << " ()"
               << "{"
               << "return this->" << member << ".detach ();"
               << "}";
        }
      };
    }

    Void
    generate_tree_inline (Context& ctx, UnsignedLong first, UnsignedLong last)
    {
      // Generate includes.
      //
      if (ctx.options.value<CLI::generate_inline> ())
      {
        Traversal::Schema schema;
        Includes includes (ctx, Includes::inline_);

        schema >> includes;

        schema.dispatch (ctx.schema_root);
      }
      else
      {
        // Emit "weak" header includes that are used in the file-per-type
        // compilation model.
        //
        Traversal::Schema schema;
        Includes includes (ctx, Includes::source);

        schema >> includes;

        schema.dispatch (ctx.schema_root);
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

      schema >> sources >> schema;
      schema >> names_ns >> ns >> names;

      names >> list;
      names >> union_;
      names >> complex;
      names >> enumeration;

      if (ctx.options.value<CLI::generate_element_type> ())
        names >> element;

      schema.dispatch (ctx.schema_root);
    }
  }
}
