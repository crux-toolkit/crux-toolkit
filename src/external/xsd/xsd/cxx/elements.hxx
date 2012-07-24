// file      : xsd/cxx/elements.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#ifndef CXX_ELEMENTS_HXX
#define CXX_ELEMENTS_HXX

#include <ostream>

#include <cult/types.hxx>
#include <cult/containers/set.hxx>
#include <cult/containers/map.hxx>
#include <cult/containers/vector.hxx>

#include <backend-elements/regex.hxx>

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

#include <elements.hxx>
#include <cxx/literal-map.hxx>

namespace CXX
{
  using std::endl;
  typedef WideString String;


  // On some platforms std::toupper can be something other than a
  // function with C++ linkage.
  //
  wchar_t
  upcase (wchar_t c);


  // Exceptions.
  //

  struct UnrepresentableCharacter
  {
    UnrepresentableCharacter (String const& str, Size pos)
        : str_ (str), pos_ (pos)
    {
    }

    String const&
    string () const
    {
      return str_;
    }

    Size
    position () const
    {
      return pos_;
    }

  private:
    String str_;
    Size pos_;
  };

  struct NoNamespaceMapping
  {
    NoNamespaceMapping (SemanticGraph::Path const& file,
                        UnsignedLong line,
                        UnsignedLong column,
                        String const& ns)
        : file_ (file),
          line_ (line),
          column_ (column),
          ns_ (ns)
    {
    }


    SemanticGraph::Path const&
    file () const
    {
      return file_;
    }

    UnsignedLong
    line () const
    {
      return line_;
    }

    UnsignedLong
    column () const
    {
      return column_;
    }

    String const&
    ns () const
    {
      return ns_;
    }

  private:
    SemanticGraph::Path file_;
    UnsignedLong line_;
    UnsignedLong column_;
    String ns_;
  };

  struct InvalidNamespaceMapping
  {
    InvalidNamespaceMapping (String const& mapping,
                             String const& reason)
        : mapping_ (mapping), reason_ (reason)
    {
    }

    String const&
    mapping () const
    {
      return mapping_;
    }

    String const&
    reason () const
    {
      return reason_;
    }

  private:
    String mapping_;
    String reason_;
  };

  //
  //
  class Context
  {
  public:
    typedef BackendElements::Regex::Pattern<WideChar> RegexPat;
    typedef BackendElements::Regex::Expression<WideChar> Regex;
    typedef Cult::Containers::Vector<Regex> RegexMapping;
    typedef Cult::Containers::Map<String, String> MapMapping;
    typedef Cult::Containers::Map<String, String> MappingCache;

    typedef Cult::Containers::Map<String, String> ReservedNameMap;
    typedef Cult::Containers::Set<String> KeywordSet;

  public:
    Context (std::wostream& o,
             SemanticGraph::Schema& root,
             SemanticGraph::Path const& path,
             StringLiteralMap const* custom_literals_map,
             NarrowString const& char_type__,
             NarrowString const& char_encoding__,
             Boolean include_with_brackets__,
             NarrowString const& include_prefix__,
             NarrowString const& esymbol,
             Containers::Vector<NarrowString> const& nsm,
             Containers::Vector<NarrowString> const& nsr,
             Boolean trace_namespace_regex_,
             Containers::Vector<NarrowString> const& include_regex,
             Boolean trace_include_regex_,
             Boolean inline_,
             Containers::Vector<NarrowString> const& reserved_name);

  protected:
    Context (Context& c)
        : os (c.os),
          schema_root (c.schema_root),
          schema_path (c.schema_path),
          char_type (c.char_type),
          char_encoding (c.char_encoding),
          L (c.L),
          string_type (c.string_type),
          string_literal_map (c.string_literal_map),
          include_with_brackets (c.include_with_brackets),
          include_prefix (c.include_prefix),
          type_exp (c.type_exp),
          inst_exp (c.inst_exp),
          inl (c.inl),
          ns_mapping_cache (c.ns_mapping_cache),
          xs_ns_ (c.xs_ns_),
          cxx_id_expr (c.cxx_id_expr),
          trace_namespace_regex (c.trace_namespace_regex),
          urn_mapping (c.urn_mapping),
          nsr_mapping (c.nsr_mapping),
          nsm_mapping (c.nsm_mapping),
          include_mapping (c.include_mapping),
          trace_include_regex (c.trace_include_regex),
          reserved_name_map (c.reserved_name_map),
          keyword_set (c.keyword_set)
    {
    }

    Context (Context& c, std::wostream& o)
        : os (o),
          schema_root (c.schema_root),
          schema_path (c.schema_path),
          char_type (c.char_type),
          char_encoding (c.char_encoding),
          L (c.L),
          string_type (c.string_type),
          string_literal_map (c.string_literal_map),
          include_with_brackets (c.include_with_brackets),
          include_prefix (c.include_prefix),
          type_exp (c.type_exp),
          inst_exp (c.inst_exp),
          inl (c.inl),
          ns_mapping_cache (c.ns_mapping_cache),
          xs_ns_ (c.xs_ns_),
          cxx_id_expr (c.cxx_id_expr),
          trace_namespace_regex (c.trace_namespace_regex),
          urn_mapping (c.urn_mapping),
          nsr_mapping (c.nsr_mapping),
          nsm_mapping (c.nsm_mapping),
          include_mapping (c.include_mapping),
          trace_include_regex (c.trace_include_regex),
          reserved_name_map (c.reserved_name_map),
          keyword_set (c.keyword_set)
    {
    }

  public:
    static String
    unclash (String const& name, String const& new_name)
    {
      return name == new_name ? (new_name + L'_') : new_name;
    }

  public:
    // Return UTF-32 character starting at this position. Position is
    // advanced by 1 if this Unicode character takes more than one
    // underlying character.
    //
    static UnsignedLong
    unicode_char (String const& str, Size& pos);

    static UnsignedLong
    unicode_char (WideChar const*& p);

    // Escape C++ keywords and illegal characters.
    //
    String
    escape (String const&);

    // Create a string literal so that it can be used in C++ source
    // code. It includes "".
    //
    String
    strlit (String const&);

    // Escape the string so that it can be used in C++ comment.
    //
    String
    comment (String const&);

    // Translate XML namespace name to a C++ identifier.
    //
    String
    ns_name (SemanticGraph::Namespace&);

    // XML Schema namespace.
    //
    SemanticGraph::Namespace&
    xs_ns ();

    // C++ namespace for XML Schema.
    //
    String
    xs_ns_name ();

    //
    //
    SemanticGraph::Namespace&
    namespace_ (SemanticGraph::Nameable& n);

    // Original XML namespace name.
    //
    String
    xml_ns_name (SemanticGraph::Nameable& ns);


    // Fully-qualified C++ name.
    //
    String
    fq_name (SemanticGraph::Nameable& n, Char const* name_key = "name");

  public:
    static SemanticGraph::Type&
    ultimate_base (SemanticGraph::Complex&);

  public:
    String
    process_include_path (String const&) const;

  public:
    static Boolean
    skip (SemanticGraph::Member& m)
    {
      // "Subsequent" local element.
      //
      return !m.scope ().is_a<SemanticGraph::Namespace> () &&
        m.context ().count ("min") == 0;
    }

    static UnsignedLong
    min (SemanticGraph::Member const& m)
    {
      return m.context ().get<UnsignedLong> ("min");
    }

    static UnsignedLong
    min (SemanticGraph::Any const& a)
    {
      return a.context ().get<UnsignedLong> ("min");
    }

    static UnsignedLong
    max (SemanticGraph::Member const& m)
    {
      return m.context ().get<UnsignedLong> ("max");
    }

    static UnsignedLong
    max (SemanticGraph::Any const& a)
    {
      return a.context ().get<UnsignedLong> ("max");
    }

  public:
    // Get escaped name.
    //
    static String const&
    ename (SemanticGraph::Nameable const& n)
    {
      return n.context ().get<String> ("name");
    }

  public:
    std::wostream& os;

    SemanticGraph::Schema& schema_root;
    SemanticGraph::Path const& schema_path;

    String& char_type;
    String& char_encoding;
    String& L;                  // string literal prefix
    String& string_type;
    StringLiteralMap const* string_literal_map;

    Boolean& include_with_brackets;
    String& include_prefix;

    String& type_exp;
    String& inst_exp;
    String& inl;

  public:
    MappingCache& ns_mapping_cache;

  private:
    SemanticGraph::Path const schema_path_;

    SemanticGraph::Namespace* xs_ns_;

    String char_type_;
    String char_encoding_;
    String L_;
    String string_type_;

    Boolean include_with_brackets_;
    String include_prefix_;

    String type_exp_;
    String inst_exp_;
    String inl_;

  private:
    RegexPat const cxx_id_expr_;
    RegexPat const& cxx_id_expr;
    Boolean trace_namespace_regex;
    Regex urn_mapping_;
    RegexMapping nsr_mapping_;
    MapMapping nsm_mapping_;
    Regex const& urn_mapping;
    RegexMapping const& nsr_mapping;
    MapMapping const& nsm_mapping;
    MappingCache ns_mapping_cache_;

    RegexMapping include_mapping_;
    RegexMapping const& include_mapping;
    Boolean trace_include_regex;

    ReservedNameMap const& reserved_name_map;
    ReservedNameMap reserved_name_map_;

    KeywordSet const& keyword_set;
    KeywordSet keyword_set_;
  };

  inline UnsignedLong Context::
  unicode_char (String const& str, Size& pos)
  {
    if (sizeof (WideChar) == 4)
    {
      return str[pos];
    }
    else if (sizeof (WideChar) == 2)
    {
      WideChar x (str[pos]);

      if (x < 0xD800 || x > 0xDBFF)
        return x;
      else
        return ((x - 0xD800) << 10) + (str[++pos] - 0xDC00) + 0x10000;
    }
    else
      return 0;
  }

  inline UnsignedLong Context::
  unicode_char (WideChar const*& p)
  {
    if (sizeof (WideChar) == 4)
    {
      return *p;
    }
    else if (sizeof (WideChar) == 2)
    {
      WideChar x (*p);

      if (x < 0xD800 || x > 0xDBFF)
        return x;
      else
        return ((x - 0xD800) << 10) + (*(++p) - 0xDC00) + 0x10000;
    }
    else
      return 0;
  }

  // Usual namespace mapping.
  //
  struct Namespace: Traversal::Namespace
  {
    struct ScopeTracker
    {
      // First scope name if always empty (global scope). The last flag
      // signals the last scope.
      //
      virtual Void
      enter (Type&, String const& name, Boolean last) = 0;

      virtual Void
      leave () = 0;
    };


    Namespace (Context& c)
        : ctx_ (c), st_ (0)
    {
    }

    Namespace (Context& c, ScopeTracker& st)
        : ctx_ (c), st_ (&st)
    {
    }

    virtual Void
    pre (Type&);

    virtual Void
    post (Type&);

  private:
    Context& ctx_;
    ScopeTracker* st_;
  };


  //
  //
  template <typename X>
  struct Has : X
  {
    Has (Boolean& result)
        : result_ (result)
    {
    }

    virtual Void
    traverse (typename X::Type&)
    {
      result_ = true;
    }

  private:
    Boolean& result_;
  };

  // Checks if scope 'Y' names any of 'X'
  //
  template <typename X, typename Y>
  Boolean
  has (Y& y)
  {
    using SemanticGraph::Scope;

    Boolean result (false);
    Has<X> t (result);

    for (Scope::NamesIterator i (y.names_begin ()), e (y.names_end ());
         !result && i != e; ++i)
      t.dispatch (i->named ());

    return result;
  }

  // Checks if the compositor has any particle of 'X'
  //
  template <typename X>
  Boolean
  has_particle (SemanticGraph::Compositor& y)
  {
    using SemanticGraph::Compositor;

    Boolean result (false);
    Has<X> t (result);

    for (Compositor::ContainsIterator i (y.contains_begin ()),
           e (y.contains_end ()); !result && i != e; ++i)
    {
      SemanticGraph::Particle& p (i->particle ());

      t.dispatch (p);

      if (!result && p.is_a<Compositor> ())
        result = has_particle<X> (dynamic_cast<Compositor&> (p));
    }

    return result;
  }

  // Specialization for Complex
  //
  template <typename X>
  Boolean
  has_particle (SemanticGraph::Complex& c)
  {
    return c.contains_compositor_p () &&
      has_particle<X> (c.contains_compositor ().compositor ());
  }

  // Fundamental type mapping helper.
  //
  struct Fundamental: Traversal::Fundamental::Type,
                      Traversal::Fundamental::String,
                      Traversal::Fundamental::NormalizedString,
                      Traversal::Fundamental::Token,
                      Traversal::Fundamental::Name,
                      Traversal::Fundamental::NameToken,
                      Traversal::Fundamental::NCName,
                      Traversal::Fundamental::Id,
                      Traversal::Fundamental::IdRef
  {
    virtual Void
    fundamental_type (SemanticGraph::Fundamental::Type& t) = 0;

    virtual Void
    fundamental_template (SemanticGraph::Fundamental::Type& t) = 0;

    virtual Void
    traverse (SemanticGraph::Fundamental::Type& t)
    {
      fundamental_type (t);
    }

    virtual Void
    traverse (SemanticGraph::Fundamental::String& t)
    {
      fundamental_template (t);
    }

    virtual Void
    traverse (SemanticGraph::Fundamental::NormalizedString& t)
    {
      fundamental_template (t);
    }

    virtual Void
    traverse (SemanticGraph::Fundamental::Token& t)
    {
      fundamental_template (t);
    }

    virtual Void
    traverse (SemanticGraph::Fundamental::Name& t)
    {
      fundamental_template (t);
    }

    virtual Void
    traverse (SemanticGraph::Fundamental::NameToken& t)
    {
      fundamental_template (t);
    }

    virtual Void
    traverse (SemanticGraph::Fundamental::NCName& t)
    {
      fundamental_template (t);
    }

    virtual Void
    traverse (SemanticGraph::Fundamental::Id& t)
    {
      fundamental_template (t);
    }

    virtual Void
    traverse (SemanticGraph::Fundamental::IdRef& t)
    {
      fundamental_template (t);
    }
  };
}

#endif  // CXX_TREE_ELEMENTS_HXX
