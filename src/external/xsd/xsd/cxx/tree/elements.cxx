// file      : xsd/cxx/tree/elements.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/tree/elements.hxx>

namespace CXX
{
  namespace Tree
  {
    // Context
    //
    Void Context::
    update_ns_scope () // Keeping this function first helps HP-UX
    {                  // (long symbols).
      ns_scope.clear ();

      Boolean first (true);

      for (NamespaceStack::Iterator i (ns_scope_stack.begin ());
           i != ns_scope_stack.end ();
           ++i)
      {
        // We only qualify names until the namespace level.
        //
        if (first)
          first = false;
        else
          ns_scope += L"::";

        ns_scope += *i;
      }
    }

    Context::
    Context (std::wostream& o,
             SemanticGraph::Schema& root,
             SemanticGraph::Path const& path,
             CLI::Options const& ops,
             Counts const& counts_,
             Boolean generate_xml_schema__,
             StringLiteralMap const* map,
             Regex const* fe,
             Regex const* he,
             Regex const* ie)
        : CXX::Context (o,
                        root,
                        path,
                        map,
                        ops.value<CLI::char_type> (),
                        ops.value<CLI::char_encoding> (),
                        ops.value<CLI::include_with_brackets> (),
                        ops.value<CLI::include_prefix> (),
                        ops.value<CLI::export_symbol> (),
                        ops.value<CLI::namespace_map> (),
                        ops.value<CLI::namespace_regex> (),
                        ops.value<CLI::namespace_regex_trace> (),
                        ops.value<CLI::include_regex> (),
                        ops.value<CLI::include_regex_trace> (),
                        ops.value<CLI::generate_inline> (),
                        ops.value<CLI::reserved_name> ()),
          options (ops),
          counts (counts_),
          any_type (any_type_),
          any_simple_type (any_simple_type_),
          element_type (element_type_),
          container (container_),
          flags_type (flags_type_),
          qname_type (qname_type_),
          xs_string_type (xs_string_type_),
          properties_type (properties_type_),
          error_handler_type (error_handler_type_),
          list_stream_type (list_stream_type_),
          namespace_infomap_type (namespace_infomap_type_),
          parser_type (parser_type_),
          std_ostream_type (std_ostream_type_),
          ostream_type (ostream_type_),
          istream_type (istream_type_),
          xerces_ns (xerces_ns_),
          dom_auto_ptr (dom_auto_ptr_),
          dom_node_key (dom_node_key_),
          as_double_type (as_double_type_),
          as_decimal_type (as_decimal_type_),
          generate_xml_schema (generate_xml_schema_),
          doxygen (doxygen_),
          polymorphic (ops.value<CLI::generate_polymorphic> ()),
          polymorphic_all (ops.value<CLI::polymorphic_type_all> ()),
          detach (ops.value<CLI::generate_detach> ()),
          fwd_expr (fe),
          hxx_expr (he),
          ixx_expr (ie),
          ns_scope (ns_scope_),
          regex_custom_type_map (regex_custom_type_map_),
          direct_custom_type_map (direct_custom_type_map_),
          qname_type_ (L"::xsd::cxx::xml::qualified_name< " + char_type + L" >"),
          parser_type_ (L"::xsd::cxx::xml::dom::parser< " + char_type + L" >"),
          generate_xml_schema_ (generate_xml_schema__),
          doxygen_ (ops.value<CLI::generate_doxygen> ()),
          ns_scope_stack (ns_scope_stack_),
          cxx_uq_id_expr_ (L"^[a-zA-Z_]\\w*$"),
          cxx_uq_id_expr (cxx_uq_id_expr_)
    {
      SemanticGraph::Namespace& xs (xs_ns ());
      SemanticGraph::Context& xsc (xs.context ());

      // Cache some often-used names from the XML Schema namespace
      // if names have already been processed.
      //
      if (xsc.count ("container"))
      {
        String xs_name (ns_name (xs));

        any_type = fq_name (xs.find ("anyType").first->named ());
        any_simple_type = fq_name (xs.find ("anySimpleType").first->named ());
        xs_string_type = fq_name (xs.find ("string").first->named ());

        container = xs_name + L"::" + xsc.get<String> ("container");
        flags_type = xs_name + L"::" + xsc.get<String> ("flags");

        if (ops.value<CLI::generate_element_type> ())
          element_type = xs_name + L"::" + xsc.get<String> ("element-type");

        properties_type = xs_name + L"::" + xsc.get<String> ("properties");

        if (!ops.value<CLI::suppress_parsing> () ||
            ops.value<CLI::generate_serialization> ())
        {
          error_handler_type = xs_name + L"::" +
            xsc.get<String> ("error-handler");
        }

        dom_auto_ptr_ = xs_name + L"::dom::auto_ptr";
        dom_node_key_ = xs_name + L"::dom::" +
          xsc.get<String> ("tree-node-key");

        if (ops.value<CLI::generate_serialization> ())
        {
          as_double_type_ = xs_name + L"::" +
            xsc.get<String> ("as-double");

          as_decimal_type_ = xs_name + L"::" +
            xsc.get<String> ("as-decimal");

          list_stream_type  = xs_name + L"::" +
            xsc.get<String> ("list-stream");

          namespace_infomap_type  = xs_name + L"::" +
            xsc.get<String> ("namespace-infomap");
        }

        // istream and ostream are templates and for now use the same
        // names regardless of the naming convention.
        //
        if (!ops.value<CLI::generate_extraction> ().empty ())
          istream_type = xs_name + L"::istream";

        if (!ops.value<CLI::generate_insertion> ().empty ())
          ostream_type = xs_name + L"::ostream";
      }

      // Xerces-C++ namespace. IntelliSense for some reason does not like
      // it fully-qualified (maybe because it's a namespace alias).
      //
      if (ops.value<CLI::generate_intellisense> ())
        xerces_ns = "xercesc";
      else
        xerces_ns = "::xercesc";

      //
      //
      if (char_type == L"char")
        std_ostream_type_ = L"::std::ostream";
      else if (char_type == L"wchar_t")
        std_ostream_type_ = L"::std::wostream";
      else
        std_ostream_type_ = L"::std::basic_ostream< " + char_type + L" >";

      // Custom type mapping.
      //
      typedef Containers::Vector<NarrowString> Vector;

      // Direct custom type mapping.
      //
      {
        Vector const& v (ops.value<CLI::custom_type> ());

        for (Vector::ConstIterator i (v.begin ()), e (v.end ()); i != e; ++i)
        {
          String s (*i);

          if (s.empty ())
            throw InvalidCustomTypeMapping (s, "mapping string is empty");

          // Split the string in two parts at the last '='.
          //
          Size pos (s.rfind ('='));

          // If no delimiter found then both type and base are empty.
          //
          if (pos == String::npos)
          {
            direct_custom_type_map[s].type.clear ();
            direct_custom_type_map[s].base.clear ();
            continue;
          }

          String name (s, 0, pos);
          String rest (s, pos + 1);

          // See if we've got the base part after '/'.
          //
          pos = rest.rfind ('/');

          String type, base;

          if (pos != String::npos)
          {
            type.assign (rest, 0, pos);
            base.assign (rest, pos + 1, String::npos);
          }
          else
            type = rest;

          // type can be a potentially-qualified template-id. base is
          // an unqualified C++ name.
          //

          if (!base.empty () && !cxx_uq_id_expr.match (base))
            throw InvalidCustomTypeMapping (s, "invalid C++ identifier");

          direct_custom_type_map[name].type = type;
          direct_custom_type_map[name].base = base;
        }
      }

      // Regex custom type mapping.
      //
      {
        Vector const& v (ops.value<CLI::custom_type_regex> ());

        for (Vector::ConstIterator i (v.begin ()), e (v.end ()); i != e; ++i)
        {
          String s (*i);

          if (s.empty ())
            throw InvalidCustomTypeMapping (s, "mapping string is empty");

          WideChar delimiter (s[0]);

          // First get pattern.
          //
          Size pos (s.find (delimiter, 1));

          if (pos == String::npos)
            throw InvalidCustomTypeMapping (
              s, "missing pattern-substitution separator");

          String pat (s, 1, pos - 1);
          String rest (s, pos + 1);

          String type, base;

          // See if we've got type and base.
          //
          if (!rest.empty ())
          {
            pos = rest.find (delimiter);

            if (pos == String::npos)
              throw InvalidCustomTypeMapping (
                s, "missing pattern-substitution separator");

            type.assign (rest, 0, pos);
            rest = String (rest, pos + 1);

            if (!rest.empty ())
            {
              pos = rest.find (delimiter);

              if (pos == String::npos)
                throw InvalidCustomTypeMapping (
                  s, "missing pattern-substitution separator");

              base.assign (rest, 0, pos);
              rest = String (rest, pos + 1);

              if (!rest.empty ())
                throw InvalidCustomTypeMapping (s, "invalid format");
            }
          }

          regex_custom_type_map.push_back (
            RegexCustomTypeMapInfo (pat, type, base));
        }
      }
    }

    Context::
    Context (Context& c)
        : CXX::Context (c),
          options (c.options),
          counts (c.counts),
          any_type (c.any_type),
          any_simple_type (c.any_simple_type),
          element_type (c.element_type),
          container (c.container),
          flags_type (c.flags_type),
          qname_type (c.qname_type),
          xs_string_type (c.xs_string_type),
          properties_type (c.properties_type),
          error_handler_type (c.error_handler_type),
          list_stream_type (c.list_stream_type),
          namespace_infomap_type (c.namespace_infomap_type),
          parser_type (c.parser_type),
          std_ostream_type (c.std_ostream_type),
          ostream_type (c.ostream_type),
          istream_type (c.istream_type),
          xerces_ns (c.xerces_ns),
          dom_auto_ptr (c.dom_auto_ptr),
          dom_node_key (c.dom_node_key),
          as_double_type (c.as_double_type),
          as_decimal_type (c.as_decimal_type),
          generate_xml_schema (c.generate_xml_schema),
          doxygen (c.doxygen),
          polymorphic (c.polymorphic),
          polymorphic_all (c.polymorphic_all),
          detach (c.detach),
          fwd_expr (c.fwd_expr),
          hxx_expr (c.hxx_expr),
          ixx_expr (c.ixx_expr),
          ns_scope (c.ns_scope),
          regex_custom_type_map (c.regex_custom_type_map),
          direct_custom_type_map (c.direct_custom_type_map),
          ns_scope_stack (c.ns_scope_stack),
          cxx_uq_id_expr (c.cxx_uq_id_expr)
    {
    }

    Context::
    Context (Context& c, std::wostream& o)
        : CXX::Context (c, o),
          options (c.options),
          counts (c.counts),
          any_type (c.any_type),
          any_simple_type (c.any_simple_type),
          element_type (c.element_type),
          container (c.container),
          flags_type (c.flags_type),
          qname_type (c.qname_type),
          xs_string_type (c.xs_string_type),
          properties_type (c.properties_type),
          error_handler_type (c.error_handler_type),
          list_stream_type (c.list_stream_type),
          namespace_infomap_type (c.namespace_infomap_type),
          parser_type (c.parser_type),
          std_ostream_type (c.std_ostream_type),
          ostream_type (c.ostream_type),
          istream_type (c.istream_type),
          xerces_ns (c.xerces_ns),
          dom_auto_ptr (c.dom_auto_ptr),
          dom_node_key (c.dom_node_key),
          as_double_type (c.as_double_type),
          as_decimal_type (c.as_decimal_type),
          generate_xml_schema (c.generate_xml_schema),
          doxygen (c.doxygen),
          polymorphic (c.polymorphic),
          polymorphic_all (c.polymorphic_all),
          detach (c.detach),
          fwd_expr (c.fwd_expr),
          hxx_expr (c.hxx_expr),
          ixx_expr (c.ixx_expr),
          ns_scope (c.ns_scope),
          regex_custom_type_map (c.regex_custom_type_map),
          direct_custom_type_map (c.direct_custom_type_map),
          ns_scope_stack (c.ns_scope_stack),
          cxx_uq_id_expr (c.cxx_uq_id_expr)
    {
    }

    Boolean Context::
    custom_type (SemanticGraph::Type const& t, String& r) const
    {
      String const& name (t.name ());

      // First search the direct mapping.
      //
      {
        DirectCustomTypeMap::ConstIterator i (
          direct_custom_type_map.find (name));

        if (i != direct_custom_type_map.end ())
        {
          r = i->second.type;
          return true;
        }
      }


      // Second search the regex mapping.
      //
      for (RegexCustomTypeMap::ConstIterator
             i (regex_custom_type_map.begin ()),
             e (regex_custom_type_map.end ());
           i != e; ++i)
      {
        if (i->pat.match (name))
        {
          // Empty type sub tells us to use the original name.
          //
          if (i->type_sub.empty ())
          {
            r.clear ();
            return true;
          }

          r = i->pat.merge (i->type_sub, name);
          return true;
        }
      }

      return false;
    }

    String Context::
    custom_type (SemanticGraph::Type const& t) const
    {
      String r;
      if (custom_type (t, r))
      {
        // Empty type name tells us to use the original name.
        //
        if (r.empty ())
          r = ename (t);
      }

      return r;
    }

    Boolean Context::
    renamed_type (SemanticGraph::Type const& t, String& r) const
    {
      String const& name (t.name ());

      // First search the direct mapping.
      //
      {
        DirectCustomTypeMap::ConstIterator i (
          direct_custom_type_map.find (name));

        if (i != direct_custom_type_map.end ())
        {
          r = i->second.base;
          return true;
        }
      }


      // Second search the regex mapping.
      //
      for (RegexCustomTypeMap::ConstIterator
             i (regex_custom_type_map.begin ()),
             e (regex_custom_type_map.end ());
           i != e; ++i)
      {
        if (i->pat.match (name))
        {
          if (!i->base_sub.empty ())
          {
            r = i->pat.merge (i->base_sub, name);
          }
          else
            r.clear ();

          return true;
        }
      }

      return false;
    }

    Void Context::
    write_annotation (SemanticGraph::Annotation& a)
    {
      String const& doc (a.documentation ());
      WideChar const* s (doc.c_str ());
      Size size (doc.size ());

      // Remove leading and trailing whitespaces.
      //
      while (*s == WideChar (0x20) || *s == WideChar (0x0A) ||
             *s == WideChar (0x0D) || *s == WideChar (0x09))
      {
        s++;
        size--;
      }

      if (size != 0)
      {
        WideChar const* e (s + size - 1);

        while (e > s &&
               (*e == WideChar (0x20) || *e == WideChar (0x0A) ||
                *e == WideChar (0x0D) || *e == WideChar (0x09)))
          --e;

        size = s <= e ? e - s + 1 : 0;
      }

      if (size != 0)
      {
        os << " * ";

        // Go over the data, forcing newline after 80 chars and adding
        // ' * ' after each new line.
        //
        WideChar const* last_space (0);
        WideChar const* b (s);
        WideChar const* e (s);
        Boolean after_newline (false);
        Boolean rogue (false);

        for (; e < s + size; ++e)
        {
          UnsignedLong u (unicode_char (e)); // May advance e.

          // We are going to treat \v and \f as rogue here even though
          // they can be present in C++ source code.
          //
          if (u > 127 || (u < 32 && u != '\t' && u != '\n'))
            rogue = true;

          if (u == ' ' || u == '\t')
          {
            if (after_newline)
            {
              if (e == b)
                b++; // Skip leading spaces after newline.

              continue;
            }
            else
              last_space = e;
          }
          else if (after_newline)
          {
            os << " * ";
            after_newline = false;
          }

          if (u == '\n')
          {
            write_rogue_text (b, e - b + 1, rogue);

            b = e + 1;
            last_space = 0;
            after_newline = true;
            rogue = false;
            continue;
          }

          if (e - b >= 70 && last_space != 0)
          {
            write_rogue_text (b, last_space - b, rogue);
            os << endl;

            b = last_space + 1;
            last_space = 0;
            after_newline = true;
            // Cannot reset rogue since we don't output the whole string.
          }
        }

        if (e != b)
          write_rogue_text (b, e - b, rogue);

        if (!after_newline)
          os << endl;
      }
    }

    Void Context::
    write_rogue_text (WideChar const* s, Size size, Boolean rogue)
    {
      if (!rogue)
        os.write (s, size);
      else
      {
        for (WideChar const* p (s); p < s + size; ++p)
        {
          UnsignedLong u (unicode_char (p)); // May advance p.

          // We are going to treat \v and \f as rogue here even though
          // they can be present in C++ source code.
          //
          if (u > 127 || (u < 32 && u != '\t' && u != '\n'))
            os.put ('?');
          else
            os.put (static_cast<WideChar> (u));
        }
      }
    }

    Boolean Context::
    polymorphic_p (SemanticGraph::Type& t)
    {
      if (polymorphic_all)
      {
        Boolean fund (false);
        IsFundamentalType test (fund);
        test.dispatch (t);
        return !fund;
      }
      else
        return t.context ().get<Boolean> ("polymorphic");
    }

    // GenerateDefautCtor
    //
    GenerateDefaultCtor::
    GenerateDefaultCtor (Context& c, Boolean& generate, Boolean no_base)
        : Context (c), generate_ (generate), no_base_ (no_base)
    {
      *this >> inherits_ >> *this;
      *this >> names_ >> *this;
    }

    Void GenerateDefaultCtor::
    traverse (SemanticGraph::Complex& c)
    {
      // Make sure we figure out if we have any required members before
      // we base our decision on the base type.
      //
      Complex::names (c, names_);

      if (!generate_)
        Complex::inherits (c, inherits_);
    }

    Void GenerateDefaultCtor::
    traverse (SemanticGraph::Type&)
    {
      if (!no_base_)
        generate_ = true;
    }

    Void GenerateDefaultCtor::
    traverse (SemanticGraph::Enumeration&)
    {
      if (!no_base_)
        generate_ = true;
    }

    Void GenerateDefaultCtor::
    traverse (SemanticGraph::Element& e)
    {
      if (!skip (e) && min (e) == 1 && max (e) == 1)
        generate_ = true;
    }

    Void GenerateDefaultCtor::
    traverse (SemanticGraph::Attribute& a)
    {
      if (min (a) == 1 && !a.fixed_p ())
        generate_ = true;
    }

    Void GenerateDefaultCtor::
    traverse (SemanticGraph::Any& a)
    {
      if (options.value<CLI::generate_wildcard> () &&
          min (a) == 1 && max (a) == 1)
        generate_ = true;
    }


    // GenerateFromBaseCtor
    //
    GenerateFromBaseCtor::
    GenerateFromBaseCtor (Context& c, Boolean& generate)
        : generate_ (generate),
          custom_ (false),
          traverser_ (c, generate, custom_)
    {
      inherits_ >> traverser_;
    }

    Void GenerateFromBaseCtor::
    traverse (SemanticGraph::Complex& c)
    {
      inherits (c, inherits_);

      if (!generate_ && custom_)
      {
        // We have a customized type in the hierarchy. In this case we
        // want to generate the c-tor unless base and ultimate-base are
        // the same (see CtorArgs).
        //
        SemanticGraph::Type& b (c.inherits ().base ());
        generate_ = b.is_a<SemanticGraph::Complex> () &&
          !b.is_a<SemanticGraph::Enumeration> ();
      }
    }

    GenerateFromBaseCtor::Traverser::
    Traverser (Context& c, Boolean& generate, Boolean& custom)
        : Context (c), generate_ (generate), custom_ (custom)
    {
      *this >> inherits_ >> *this;
      *this >> names_ >> *this;
    }

    Void GenerateFromBaseCtor::Traverser::
    traverse (SemanticGraph::Type& t)
    {
      if (!custom_)
      {
        String tmp;
        custom_ = custom_type (t, tmp);
      }
    }

    Void GenerateFromBaseCtor::Traverser::
    traverse (SemanticGraph::Complex& c)
    {
      names (c, names_);

      if (!generate_)
        inherits (c, inherits_);

      if (!generate_)
        traverse (static_cast<SemanticGraph::Type&> (c));
    }

    Void GenerateFromBaseCtor::Traverser::
    traverse (SemanticGraph::Element& e)
    {
      if (!skip (e) && min (e) == 1 && max (e) == 1)
        generate_ = true;
    }

    Void GenerateFromBaseCtor::Traverser::
    traverse (SemanticGraph::Attribute& a)
    {
      if (min (a) == 1 && !a.fixed_p ())
        generate_ = true;
    }

    Void GenerateFromBaseCtor::Traverser::
    traverse (SemanticGraph::Any& a)
    {
      if (options.value<CLI::generate_wildcard> () &&
          min (a) == 1 && max (a) == 1)
        generate_ = true;
    }

    // HasComplexNonOptArgs
    //
    HasComplexPolyNonOptArgs::
    HasComplexPolyNonOptArgs (Context& c,
                              Boolean base,
                              Boolean& complex,
                              Boolean& poly,
                              Boolean& clash)
        : Context (c),
          complex_ (complex),
          poly_ (poly),
          clash_ (clash)
    {
      if (base)
        *this >> inherits_ >> *this;

      *this >> names_ >> *this;
    }

    Void HasComplexPolyNonOptArgs::
    traverse (SemanticGraph::Complex& c)
    {
      // No optimizations: need to check every arg for clashes.
      //
      inherits (c, inherits_);
      names (c, names_);
    }

    Void HasComplexPolyNonOptArgs::
    traverse (SemanticGraph::Element& e)
    {
      if (!skip (e) && min (e) == 1 && max (e) == 1)
      {
        Boolean poly (polymorphic && polymorphic_p (e.type ()));

        Boolean simple (true);
        IsSimpleType t (simple);
        t.dispatch (e.type ());

        if (poly)
          poly_ = true;

        if (!simple)
          complex_ = true;

        if (poly && simple)
          clash_ = false;
      }
    }

    // FromBaseCtorArg
    //
    FromBaseCtorArg::
    FromBaseCtorArg (Context& c, ArgType at, Boolean arg)
        : Context (c), arg_type_ (at), arg_ (arg)
    {
    }

    Void FromBaseCtorArg::
    traverse (SemanticGraph::Any& a)
    {
      if (!options.value<CLI::generate_wildcard> ())
        return;

      if (min (a) == 1 && max (a) == 1)
      {
        String const& name (ename (a));

        os << "," << endl
           << "const " << xerces_ns << "::DOMElement&";

        if (arg_)
          os << " " << name;
      }
    }

    Void FromBaseCtorArg::
    traverse (SemanticGraph::Element& e)
    {
      if (skip (e))
        return;

      if (min (e) == 1 && max (e) == 1)
      {
        String const& name (ename (e));

        os << "," << endl;

        Boolean auto_ptr (false);

        switch (arg_type_)
        {
        case arg_complex_auto_ptr:
          {
            Boolean simple (true);
            IsSimpleType t (simple);
            t.dispatch (e.type ());
            auto_ptr = !simple;
            break;
          }
        case arg_poly_auto_ptr:
          {
            auto_ptr = polymorphic && polymorphic_p (e.type ());
            break;
          }
        case arg_type:
          break;
        }

        if (auto_ptr)
          os << "::std::auto_ptr< " << etype (e) << " >&";
        else
          os << "const " << etype (e) << "&";

        if (arg_)
          os << " " << name;
      }
    }

    Void FromBaseCtorArg::
    traverse (SemanticGraph::Attribute& a)
    {
      // Note that we are not going to include attributes with
      // default or required fixed values here. Instead we are
      // going to default-initialize them.
      //
      if (min (a) == 1 && !a.fixed_p ())
      {
        String const& name (ename (a));

        os << "," << endl
           << "const " << etype (a) << "&";

        if (arg_)
          os << " " << name;
      }
    }

    // CtorArgs
    //
    CtorArgs::
    CtorArgs (Context& c, ArgType at)
        : Context (c),
          arg_type_ (at),
          base_arg_ (0),
          first_ (true),
          member_name_ (c)
    {
      *this >> inherits_ >> *this;
      *this >> names_ >> *this;
    }

    CtorArgs::
    CtorArgs (Context& c, ArgType at, String& base_arg)
        : Context (c),
          arg_type_ (at),
          base_arg_ (&base_arg),
          first_ (true),
          member_name_ (c)
    {
      *this >> inherits_ >> *this;
      *this >> names_ >> *this;
    }

    Void CtorArgs::
    traverse (SemanticGraph::Type& t)
    {
      os << comma () << "const ";

      member_name_.dispatch (t);

      os << "&";

      if (base_arg_ != 0)
      {
        *base_arg_ = L"_xsd_" + ename (t) + L"_base";

        os << " " << *base_arg_;
      }
    }

    Void CtorArgs::
    traverse (SemanticGraph::Enumeration& e)
    {
      os << comma () << "const ";

      member_name_.traverse (e);

      os << "&";

      if (base_arg_ != 0)
      {
        *base_arg_ = L"_xsd_" + ename (e) + L"_base";

        os << " " << *base_arg_;
      }
    }

    Void CtorArgs::
    traverse (SemanticGraph::Any& a)
    {
      if (!options.value<CLI::generate_wildcard> ())
        return;

      if (min (a) == 1 && max (a) == 1)
      {
        os << comma () << "const " << xerces_ns << "::DOMElement&";

        if (base_arg_ != 0)
          os << " " << ename (a);
      }
    }

    Void CtorArgs::
    traverse (SemanticGraph::Element& e)
    {
      if (skip (e))
        return;

      if (min (e) == 1 && max (e) == 1)
      {
        Boolean auto_ptr (false);

        switch (arg_type_)
        {
        case arg_complex_auto_ptr:
          {
            Boolean simple (true);
            IsSimpleType t (simple);
            t.dispatch (e.type ());
            auto_ptr = !simple;
            break;
          }
        case arg_poly_auto_ptr:
          {
            auto_ptr = polymorphic && polymorphic_p (e.type ());
            break;
          }
        case arg_type:
          break;
        }

        if (auto_ptr)
          os << comma () << "::std::auto_ptr< " << etype (e) << " >&";
        else
          os << comma () << "const " << etype (e) << "&";

        if (base_arg_ != 0)
          os << " " << ename (e);
      }
    }

    Void CtorArgs::
    traverse (SemanticGraph::Attribute& a)
    {
      // Note that we are not going to include attributes with
      // default or required fixed values here. Instead we are
      // going to default-initialize them.
      //
      if (min (a) == 1 && !a.fixed_p ())
      {
        os << comma () << "const " << etype (a) << "&";

        if (base_arg_ != 0)
          os << " " << ename (a);
      }
    }

    String CtorArgs::
    comma ()
    {
      Boolean tmp (first_);
      first_ = false;
      return tmp ? "" : ",\n";
    }


    // CtorArgsWithoutBase
    //
    CtorArgsWithoutBase::
    CtorArgsWithoutBase (Context& c, ArgType at, Boolean arg, Boolean first)
        : Context (c), arg_type_ (at), arg_ (arg), first_ (first)
    {
      *this >> inherits_ >> *this;
      *this >> names_ >> *this;
    }

    Void CtorArgsWithoutBase::
    traverse (SemanticGraph::Any& a)
    {
      if (!options.value<CLI::generate_wildcard> ())
        return;

      if (min (a) == 1 && max (a) == 1)
      {
        os << comma () << "const " << xerces_ns << "::DOMElement&";

        if (arg_)
          os << " " << ename (a);
      }
    }

    Void CtorArgsWithoutBase::
    traverse (SemanticGraph::Element& e)
    {
      if (skip (e))
        return;

      if (min (e) == 1 && max (e) == 1)
      {
        Boolean auto_ptr (false);

        switch (arg_type_)
        {
        case arg_complex_auto_ptr:
          {
            Boolean simple (true);
            IsSimpleType t (simple);
            t.dispatch (e.type ());
            auto_ptr = !simple;
            break;
          }
        case arg_poly_auto_ptr:
          {
            auto_ptr = polymorphic && polymorphic_p (e.type ());
            break;
          }
        case arg_type:
          break;
        }

        if (auto_ptr)
          os << comma () << "::std::auto_ptr< " << etype (e) << " >&";
        else
          os << comma () << "const " << etype (e) << "&";

        if (arg_)
          os << " " << ename (e);
      }
    }

    Void CtorArgsWithoutBase::
    traverse (SemanticGraph::Attribute& a)
    {
      // Note that we are not going to include attributes with
      // default or required fixed values here. Instead we are
      // going to default-initialize them.
      //
      if (min (a) == 1 && !a.fixed_p ())
      {
        os << comma () << "const " << etype (a) << "&";

        if (arg_)
          os << " " << ename (a);
      }
    }

    String CtorArgsWithoutBase::
    comma ()
    {
      Boolean tmp (first_);
      first_ = false;
      return tmp ? "" : ",\n";
    }

    // GlobalElementBase
    //
    Boolean GlobalElementBase::
    generate_p (SemanticGraph::Element& e)
    {
      if (e.substitutes_p () && ctx_.polymorphic)
        return true;

      if (!doc_root_p (e))
        return false;

      // If we are not generating element types nor parsing/serialization
      // code then we won't generate anything from it.
      //
      if (!ctx_.options.value<CLI::generate_element_type> () &&
          ctx_.options.value<CLI::suppress_parsing> () &&
          !ctx_.options.value<CLI::generate_serialization> ())
        return false;

      return true;
    }

    Boolean GlobalElementBase::
    doc_root_p (SemanticGraph::Element& e)
    {
      if (!ctx_.options.value<CLI::root_element_first> () &&
          !ctx_.options.value<CLI::root_element_last> () &&
          !ctx_.options.value<CLI::root_element_all> () &&
          !ctx_.options.value<CLI::root_element_none> () &&
          ctx_.options.value<CLI::root_element> ().empty ())
        return true; // By default treat them all.

      if (ctx_.options.value<CLI::root_element_none> ())
        return false;

      if (ctx_.options.value<CLI::root_element_all> ())
        return true;

      if (ctx_.options.value<CLI::root_element_first> () &&
          e.context ().count ("first") != 0)
        return true;

      if (ctx_.options.value<CLI::root_element_last> () &&
          e.context ().count ("last") != 0)
        return true;

      typedef Cult::Containers::Vector<NarrowString> Names;
      Names const& names (ctx_.options.value<CLI::root_element> ());

      // Hopefully nobody will specify more than a handful of names ;-).
      //
      for (Names::ConstIterator i (names.begin ()); i != names.end (); ++i)
      {
        String name (*i);

        if (e.name () == name)
          return true;
      }

      return false;
    }

    // Namespace
    //
    Namespace::
    Namespace (Context& c,
               UnsignedLong first,
               UnsignedLong last)
        : CXX::Namespace (c, *this),
          GlobalElementBase (c),
          ctx_ (c),
          first_ (first),
          last_ (last),
          count_ (0)
    {
    }

    Void Namespace::
    traverse (Type& ns)
    {
      using SemanticGraph::Element;

      if (first_ > last_)
        CXX::Namespace::traverse (ns);
      else
      {
        Boolean opened (false);

        for (Type::NamesIterator i (ns.names_begin ());
             i != ns.names_end (); ++i)
        {
          SemanticGraph::Nameable& n (i->named ());

          if (n.is_a<SemanticGraph::Type> () ||
              (n.is_a<Element> () && generate_p (dynamic_cast<Element&> (n))))
          {
            if (count_ >= first_ && count_ <= last_)
            {
              if (!opened)
              {
                opened = true;
                pre (ns);
              }

              edge_traverser ().dispatch (*i);
            }

            ++count_;
          }
        }

        if (opened)
          post (ns);
      }
    }

    Void Namespace::
    enter (Type&, String const& name, Boolean)
    {
      ctx_.enter_ns_scope (name);
    }

    Void Namespace::
    leave ()
    {
      ctx_.leave_ns_scope ();
    }

    // Includes
    //
    Void TypeForward::
    traverse (SemanticGraph::Type& t)
    {
      String const& name (ename (t));

      if (String custom = custom_type (t))
      {
        String new_name;
        renamed_type (t, new_name);

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

    Void Includes::
    traverse_ (SemanticGraph::Uses& u)
    {
      // Support for weak (forward) inclusion used in the file-per-type
      // compilation model.
      //
      Type t (type_);
      Boolean weak (u.context ().count ("weak"));
      SemanticGraph::Schema& s (u.schema ());

      if (weak && t == header)
      {
        // Generate forward declarations.
        //
        if (forward_)
          t = forward;
        else
        {
          schema_.dispatch (s);
          return;
        }
      }

      if (t == source && !weak)
        return;

      SemanticGraph::Path path (
        s.context ().count ("renamed")
        ? s.context ().get<SemanticGraph::Path> ("renamed")
        : u.path ());

      // Try to use the portable representation of the path. If that
      // fails, fall back to the native representation.
      //
      NarrowString path_str;
      try
      {
        path_str = path.string ();
      }
      catch (SemanticGraph::InvalidPath const&)
      {
#if !defined(BOOST_FILESYSTEM_VERSION) || BOOST_FILESYSTEM_VERSION == 2
        path_str = path.native_file_string ();
#else
        path_str = path.string ();
#endif
      }

      String inc_path;

      switch (t)
      {
      case forward:
        {
          inc_path = ctx_.fwd_expr->merge (path_str);
          break;
        }
      case header:
      case source:
        {
          inc_path = ctx_.hxx_expr->merge (path_str);
          break;
        }
      case inline_:
        {
          if (weak)
          {
            inc_path = ctx_.hxx_expr->merge (path_str);
            ctx_.os << "#include " << ctx_.process_include_path (inc_path)
                    << endl;
          }

          inc_path = ctx_.ixx_expr->merge (path_str);
          break;
        }
      }

      ctx_.os << "#include " << ctx_.process_include_path (inc_path) << endl
              << endl;
    }
  }
}
