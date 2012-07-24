// file      : xsd/cxx/parser/cli.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#ifndef CXX_PARSER_CLI_HXX
#define CXX_PARSER_CLI_HXX

#include <cult/types.hxx>

#include <cult/containers/vector.hxx>

#include <cult/cli/options.hxx>
#include <cult/cli/options-spec.hxx>

namespace CXX
{
  namespace Parser
  {
    namespace CLI
    {
      using namespace Cult::Types;

      typedef Char const Key[];

      extern Key type_map;
      extern Key char_encoding;
      extern Key char_type;
      extern Key output_dir;
      extern Key xml_parser;
      extern Key generate_inline;
      extern Key generate_validation;
      extern Key suppress_validation;
      extern Key generate_polymorphic;
      extern Key generate_noop_impl;
      extern Key generate_print_impl;
      extern Key generate_test_driver;
      extern Key force_overwrite;
      extern Key root_element_first;
      extern Key root_element_last;
      extern Key root_element;
      extern Key generate_xml_schema;
      extern Key extern_xml_schema;
      extern Key skel_type_suffix;
      extern Key skel_file_suffix;
      extern Key impl_type_suffix;
      extern Key impl_file_suffix;
      extern Key namespace_map;
      extern Key namespace_regex;
      extern Key namespace_regex_trace;
      extern Key reserved_name;
      extern Key include_with_brackets;
      extern Key include_prefix;
      extern Key include_regex;
      extern Key include_regex_trace;
      extern Key guard_prefix;
      extern Key hxx_suffix;
      extern Key ixx_suffix;
      extern Key cxx_suffix;
      extern Key hxx_regex;
      extern Key ixx_regex;
      extern Key cxx_regex;
      extern Key hxx_prologue;
      extern Key ixx_prologue;
      extern Key cxx_prologue;
      extern Key prologue;
      extern Key hxx_epilogue;
      extern Key ixx_epilogue;
      extern Key cxx_epilogue;
      extern Key epilogue;
      extern Key hxx_prologue_file;
      extern Key ixx_prologue_file;
      extern Key cxx_prologue_file;
      extern Key prologue_file;
      extern Key hxx_epilogue_file;
      extern Key ixx_epilogue_file;
      extern Key cxx_epilogue_file;
      extern Key epilogue_file;
      extern Key export_symbol;
      extern Key export_maps;
      extern Key import_maps;
      extern Key show_anonymous;
      extern Key show_sloc;
      extern Key proprietary_license;

      typedef Cult::CLI::Options<
        type_map,                 Cult::Containers::Vector<NarrowString>,
        char_type,                NarrowString,
        char_encoding,            NarrowString,
        output_dir,               NarrowString,
        xml_parser,               NarrowString,
        generate_inline,          Boolean,
        generate_validation,      Boolean,
        suppress_validation,      Boolean,
        generate_polymorphic,     Boolean,
        generate_noop_impl,       Boolean,
        generate_print_impl,      Boolean,
        generate_test_driver,     Boolean,
        force_overwrite,          Boolean,
        root_element_first,       Boolean,
        root_element_last,        Boolean,
        root_element,             NarrowString,
        generate_xml_schema,      Boolean,
        extern_xml_schema,        NarrowString,
        skel_type_suffix,         NarrowString,
        skel_file_suffix,         NarrowString,
        impl_type_suffix,         NarrowString,
        impl_file_suffix,         NarrowString,
        namespace_map,            Cult::Containers::Vector<NarrowString>,
        namespace_regex,          Cult::Containers::Vector<NarrowString>,
        namespace_regex_trace,    Boolean,
        reserved_name,            Cult::Containers::Vector<NarrowString>,
        include_with_brackets,    Boolean,
        include_prefix,           NarrowString,
        include_regex,            Cult::Containers::Vector<NarrowString>,
        include_regex_trace,      Boolean,
        guard_prefix,             NarrowString,
        hxx_suffix,               NarrowString,
        ixx_suffix,               NarrowString,
        cxx_suffix,               NarrowString,
        hxx_regex,                NarrowString,
        ixx_regex,                NarrowString,
        cxx_regex,                NarrowString,
        hxx_prologue,             Cult::Containers::Vector<NarrowString>,
        ixx_prologue,             Cult::Containers::Vector<NarrowString>,
        cxx_prologue,             Cult::Containers::Vector<NarrowString>,
        prologue,                 Cult::Containers::Vector<NarrowString>,
        hxx_epilogue,             Cult::Containers::Vector<NarrowString>,
        ixx_epilogue,             Cult::Containers::Vector<NarrowString>,
        cxx_epilogue,             Cult::Containers::Vector<NarrowString>,
        epilogue,                 Cult::Containers::Vector<NarrowString>,
        hxx_prologue_file,        NarrowString,
        ixx_prologue_file,        NarrowString,
        cxx_prologue_file,        NarrowString,
        prologue_file,            NarrowString,
        hxx_epilogue_file,        NarrowString,
        ixx_epilogue_file,        NarrowString,
        cxx_epilogue_file,        NarrowString,
        epilogue_file,            NarrowString,
        export_symbol,            NarrowString,
        export_maps,              Boolean,
        import_maps,              Boolean,
        show_anonymous,           Boolean,
        show_sloc,                Boolean,
        proprietary_license,      Boolean

        > Options;

      struct OptionsSpec: Cult::CLI::OptionsSpec<Options> {};
    }
  }
}

#endif // CXX_PARSER_CLI_HXX
