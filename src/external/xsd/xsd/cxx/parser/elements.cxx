// file      : xsd/cxx/parser/elements.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/parser/elements.hxx>

namespace CXX
{
  namespace Parser
  {
    // Keep this symbol first to help HP-UX linker (long symbols?).
    //
    Content::Value Context::
    content (SemanticGraph::Complex& c)
    {
      using namespace SemanticGraph;

      if (c.mixed_p ())
        return Content::mixed;

      if (c.inherits_p ())
      {
        Type& base (c.inherits ().base ());

        if (Complex* cb = dynamic_cast<Complex*> (&base))
          return content (*cb);

        if (base.is_a<AnyType> ())
          return Content::complex;

        // Everyhting else (built-in type and AnySimpleType) is simple
        // content.
        //
        return Content::simple;
      }
      else
        return Content::complex;
    }

    Context::
    Context (std::wostream& o,
             SemanticGraph::Schema& root,
             SemanticGraph::Path const& path,
             CLI::Options const& ops,
             StringLiteralMap const* map,
             Regex const* he,
             Regex const* ie,
             Regex const* hie)
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
          xml_parser (xml_parser_),
          simple_base (simple_base_),
          complex_base (complex_base_),
          list_base (list_base_),
          cout_inst (cout_inst_),
          cerr_inst (cerr_inst_),
          parser_map (parser_map_),
          std_string_type (std_string_type_),
          validation (validation_),
          polymorphic (polymorphic_),
          hxx_expr (he),
          ixx_expr (ie),
          hxx_impl_expr (hie),
          xml_parser_ (ops.value<CLI::xml_parser> ()),
          validation_ ((ops.value<CLI::xml_parser> () == "expat" ||
                        ops.value<CLI::generate_validation> ()) &&
                       !ops.value<CLI::suppress_validation> ()),
          polymorphic_ (ops.value<CLI::generate_polymorphic> ())
    {
      if (char_type == L"char")
        std_string_type = L"::std::string";
      else if (char_type == L"wchar_t")
        std_string_type = L"::std::wstring";
      else
        std_string_type = L"::std::basic_string< " + char_type + L" >";

      String xs_ns (xs_ns_name ());

      string_type = xs_ns + L"::ro_string";
      simple_base = xs_ns + L"::simple_content";
      complex_base = xs_ns + L"::complex_content";
      list_base = xs_ns + L"::list_base";

      cout_inst = (char_type == L"char" ? L"std::cout" : L"std::wcout");
      cerr_inst = (char_type == L"char" ? L"std::cerr" : L"std::wcerr");

      if (polymorphic)
        parser_map_ = xs_ns + L"::parser_map";
    }

    Context::
    Context (Context& c)
        : CXX::Context (c),
          options (c.options),
          xml_parser (c.xml_parser),
          simple_base (c.simple_base),
          complex_base (c.complex_base),
          list_base (c.list_base),
          cout_inst (c.cout_inst),
          cerr_inst (c.cerr_inst),
          parser_map (c.parser_map),
          std_string_type (c.std_string_type),
          validation (c.validation),
          polymorphic (c.polymorphic),
          hxx_expr (c.hxx_expr),
          ixx_expr (c.ixx_expr),
          hxx_impl_expr (c.hxx_impl_expr)
    {
    }

    Context::
    Context (Context& c, std::wostream& o)
        : CXX::Context (c, o),
          options (c.options),
          xml_parser (c.xml_parser),
          simple_base (c.simple_base),
          complex_base (c.complex_base),
          list_base (c.list_base),
          cout_inst (c.cout_inst),
          cerr_inst (c.cerr_inst),
          parser_map (c.parser_map),
          std_string_type (c.std_string_type),
          validation (c.validation),
          polymorphic (c.polymorphic),
          hxx_expr (c.hxx_expr),
          ixx_expr (c.ixx_expr),
          hxx_impl_expr (c.hxx_impl_expr)
    {
    }

    Boolean Context::
    anonymous (SemanticGraph::Type& t)
    {
      return t.context ().count ("anonymous");
    }

    String const& Context::
    ret_type (SemanticGraph::Type& t)
    {
      return t.context ().get<String> ("ret-type");
    }

    String const& Context::
    arg_type (SemanticGraph::Type& t)
    {
      return t.context ().get<String> ("arg-type");
    }

    String const& Context::
    post_name (SemanticGraph::Type& t)
    {
      return t.context ().get<String> ("post");
    }

    String const& Context::
    eparser (SemanticGraph::Member& m)
    {
      return m.context ().get<String> ("parser");
    }

    String const& Context::
    emember (SemanticGraph::Member& m)
    {
      return m.context ().get<String> ("member");
    }

    String const& Context::
    emember_cache (SemanticGraph::Member& m)
    {
      return m.context ().get<String> ("member-cache");
    }

    String const& Context::
    emember_map (SemanticGraph::Member& m)
    {
      return m.context ().get<String> ("member-map");
    }

    String const& Context::
    eimpl (SemanticGraph::Type& t)
    {
      return t.context ().get<String> ("impl");
    }

    // Includes
    //
    Void TypeForward::
    traverse (SemanticGraph::Type& t)
    {
      os << "class " << t.context ().get<String> (name_key_) << ";";
    }

    Void Includes::
    traverse_ (SemanticGraph::Uses& u)
    {
      // Support for weak (forward) inclusion used in the file-per-type
      // compilation model.
      //
      SemanticGraph::Schema& s (u.schema ());
      Boolean weak (u.context ().count ("weak"));

      if (weak && (type_ == header || type_ == impl_header))
      {
        // Generate forward declarations. We don't really need them
        // in the impl files.
        //
        if (type_ == header)
          schema_.dispatch (s);

        return;
      }

      if (type_ == source && !weak)
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

      switch (type_)
      {
      case header:
      case source:
        {
          inc_path = ctx_.hxx_expr->merge (path_str);
          break;
        }
      case impl_header:
        {
          inc_path = ctx_.hxx_impl_expr->merge (path_str);
          break;
        }
      }

      ctx_.os << "#include " << ctx_.process_include_path (inc_path) << endl
              << endl;
    }
  }
}
