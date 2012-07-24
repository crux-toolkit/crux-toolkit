// file      : xsd/cxx/tree/cli.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#ifndef CXX_TREE_CLI_HXX
#define CXX_TREE_CLI_HXX

#include <cult/types.hxx>

#include <cult/containers/vector.hxx>

#include <cult/cli/options.hxx>
#include <cult/cli/options-spec.hxx>

namespace CXX
{
  namespace Tree
  {
    namespace CLI
    {
      using namespace Cult::Types;

      typedef Char const Key[];

      extern Key char_type;
      extern Key char_encoding;
      extern Key output_dir;
      extern Key generate_polymorphic;
      extern Key polymorphic_type;
      extern Key polymorphic_type_all;
      extern Key generate_serialization;
      extern Key generate_inline;
      extern Key generate_ostream;
      extern Key generate_doxygen;
      extern Key generate_comparison;
      extern Key generate_default_ctor;
      extern Key generate_from_base_ctor;
      extern Key generate_detach;
      extern Key generate_wildcard;
      extern Key generate_insertion;
      extern Key generate_extraction;
      extern Key generate_forward;
      extern Key generate_xml_schema;
      extern Key extern_xml_schema;
      extern Key suppress_parsing;
      extern Key generate_element_type;
      extern Key generate_element_map;
      extern Key generate_intellisense;
      extern Key omit_default_attributes;
      extern Key namespace_map;
      extern Key namespace_regex;
      extern Key namespace_regex_trace;
      extern Key reserved_name;
      extern Key type_naming;
      extern Key function_naming;
      extern Key type_regex;
      extern Key accessor_regex;
      extern Key one_accessor_regex;
      extern Key opt_accessor_regex;
      extern Key seq_accessor_regex;
      extern Key modifier_regex;
      extern Key one_modifier_regex;
      extern Key opt_modifier_regex;
      extern Key seq_modifier_regex;
      extern Key parser_regex;
      extern Key serializer_regex;
      extern Key enumerator_regex;
      extern Key element_type_regex;
      extern Key name_regex_trace;
      extern Key include_with_brackets;
      extern Key include_prefix;
      extern Key include_regex;
      extern Key include_regex_trace;
      extern Key guard_prefix;
      extern Key root_element_first;
      extern Key root_element_last;
      extern Key root_element_all;
      extern Key root_element_none;
      extern Key root_element;
      extern Key custom_type;
      extern Key custom_type_regex;
      extern Key hxx_suffix;
      extern Key ixx_suffix;
      extern Key cxx_suffix;
      extern Key fwd_suffix;
      extern Key hxx_regex;
      extern Key ixx_regex;
      extern Key cxx_regex;
      extern Key fwd_regex;
      extern Key hxx_prologue;
      extern Key ixx_prologue;
      extern Key cxx_prologue;
      extern Key fwd_prologue;
      extern Key prologue;
      extern Key hxx_epilogue;
      extern Key ixx_epilogue;
      extern Key cxx_epilogue;
      extern Key fwd_epilogue;
      extern Key epilogue;
      extern Key hxx_prologue_file;
      extern Key ixx_prologue_file;
      extern Key cxx_prologue_file;
      extern Key fwd_prologue_file;
      extern Key prologue_file;
      extern Key hxx_epilogue_file;
      extern Key ixx_epilogue_file;
      extern Key cxx_epilogue_file;
      extern Key fwd_epilogue_file;
      extern Key epilogue_file;
      extern Key parts;
      extern Key parts_suffix;
      extern Key export_symbol;
      extern Key export_xml_schema;
      extern Key export_maps;
      extern Key import_maps;
      extern Key show_anonymous;
      extern Key show_sloc;
      extern Key proprietary_license;
      extern Key disable_multi_import; // Undocumented.


      typedef Cult::CLI::Options<

        char_type,                NarrowString,
        char_encoding,            NarrowString,
        output_dir,               NarrowString,
        generate_polymorphic,     Boolean,
        polymorphic_type,         Cult::Containers::Vector<NarrowString>,
        polymorphic_type_all,     Boolean,
        generate_serialization,   Boolean,
        generate_inline,          Boolean,
        generate_ostream,         Boolean,
        generate_doxygen,         Boolean,
        generate_comparison,      Boolean,
        generate_default_ctor,    Boolean,
        generate_from_base_ctor,  Boolean,
        generate_detach,          Boolean,
        generate_wildcard,        Boolean,
        generate_insertion,       Cult::Containers::Vector<NarrowString>,
        generate_extraction,      Cult::Containers::Vector<NarrowString>,
        generate_forward,         Boolean,
        generate_xml_schema,      Boolean,
        extern_xml_schema,        NarrowString,
        suppress_parsing,         Boolean,
        generate_element_type,    Boolean,
        generate_element_map,     Boolean,
        generate_intellisense,    Boolean,
        omit_default_attributes,  Boolean,
        namespace_map,            Cult::Containers::Vector<NarrowString>,
        namespace_regex,          Cult::Containers::Vector<NarrowString>,
        namespace_regex_trace,    Boolean,
        reserved_name,            Cult::Containers::Vector<NarrowString>,
        type_naming,              NarrowString,
        function_naming,          NarrowString,
        type_regex,               Cult::Containers::Vector<NarrowString>,
        accessor_regex,           Cult::Containers::Vector<NarrowString>,
        one_accessor_regex,       Cult::Containers::Vector<NarrowString>,
        opt_accessor_regex,       Cult::Containers::Vector<NarrowString>,
        seq_accessor_regex,       Cult::Containers::Vector<NarrowString>,
        modifier_regex,           Cult::Containers::Vector<NarrowString>,
        one_modifier_regex,       Cult::Containers::Vector<NarrowString>,
        opt_modifier_regex,       Cult::Containers::Vector<NarrowString>,
        seq_modifier_regex,       Cult::Containers::Vector<NarrowString>,
        parser_regex,             Cult::Containers::Vector<NarrowString>,
        serializer_regex,         Cult::Containers::Vector<NarrowString>,
        enumerator_regex,         Cult::Containers::Vector<NarrowString>,
        element_type_regex,       Cult::Containers::Vector<NarrowString>,
        name_regex_trace,         Boolean,
        include_with_brackets,    Boolean,
        include_prefix,           NarrowString,
        include_regex,            Cult::Containers::Vector<NarrowString>,
        include_regex_trace,      Boolean,
        guard_prefix,             NarrowString,
        root_element_first,       Boolean,
        root_element_last,        Boolean,
        root_element_all,         Boolean,
        root_element_none,        Boolean,
        root_element,             Cult::Containers::Vector<NarrowString>,
        custom_type,              Cult::Containers::Vector<NarrowString>,
        custom_type_regex,        Cult::Containers::Vector<NarrowString>,
        hxx_suffix,               NarrowString,
        ixx_suffix,               NarrowString,
        cxx_suffix,               NarrowString,
        fwd_suffix,               NarrowString,
        hxx_regex,                NarrowString,
        ixx_regex,                NarrowString,
        cxx_regex,                NarrowString,
        fwd_regex,                NarrowString,
        hxx_prologue,             Cult::Containers::Vector<NarrowString>,
        ixx_prologue,             Cult::Containers::Vector<NarrowString>,
        cxx_prologue,             Cult::Containers::Vector<NarrowString>,
        fwd_prologue,             Cult::Containers::Vector<NarrowString>,
        prologue,                 Cult::Containers::Vector<NarrowString>,
        hxx_epilogue,             Cult::Containers::Vector<NarrowString>,
        ixx_epilogue,             Cult::Containers::Vector<NarrowString>,
        cxx_epilogue,             Cult::Containers::Vector<NarrowString>,
        fwd_epilogue,             Cult::Containers::Vector<NarrowString>,
        epilogue,                 Cult::Containers::Vector<NarrowString>,
        hxx_prologue_file,        NarrowString,
        ixx_prologue_file,        NarrowString,
        cxx_prologue_file,        NarrowString,
        fwd_prologue_file,        NarrowString,
        prologue_file,            NarrowString,
        hxx_epilogue_file,        NarrowString,
        ixx_epilogue_file,        NarrowString,
        cxx_epilogue_file,        NarrowString,
        fwd_epilogue_file,        NarrowString,
        epilogue_file,            NarrowString,
        parts,                    UnsignedLong,
        parts_suffix,             NarrowString,
        export_symbol,            NarrowString,
        export_xml_schema,        Boolean,
        export_maps,              Boolean,
        import_maps,              Boolean,
        show_anonymous,           Boolean,
        show_sloc,                Boolean,
        proprietary_license,      Boolean,
        disable_multi_import,     Boolean

        > Options;

      struct OptionsSpec: Cult::CLI::OptionsSpec<Options> {};
    }
  }
}

#endif // CXX_TREE_CLI_HXX
