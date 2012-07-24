// file      : xsd/xsd.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cult/types.hxx>

#include <cult/trace/log.hxx>

#include <cult/containers/map.hxx>
#include <cult/containers/vector.hxx>

#include <cult/cli/exceptions.hxx>
#include <cult/cli/file-arguments.hxx>
#include <cult/cli/options.hxx>
#include <cult/cli/options-spec.hxx>
#include <cult/cli/options-parser.hxx>

#include <xsd-frontend/parser.hxx>
#include <xsd-frontend/transformations/anonymous.hxx>
#include <xsd-frontend/transformations/enum-synthesis.cxx>
#include <xsd-frontend/transformations/restriction.hxx>
#include <xsd-frontend/transformations/schema-per-type.hxx>
#include <xsd-frontend/transformations/simplifier.hxx>

#include <backend-elements/regex.hxx>
#include <backend-elements/indentation/clip.hxx>

#include <cxx/tree/generator.hxx>
#include <cxx/parser/generator.hxx>

#include <processing/cardinality/processor.hxx>
#include <processing/inheritance/processor.hxx>

#include <iostream>
#include <boost/filesystem/fstream.hpp>

#include <xercesc/util/PlatformUtils.hpp>

#include <xsd.hxx>
#include <usage.hxx>

#include "../libxsd/xsd/cxx/version.hxx"

using namespace Cult::Types;

typedef Cult::Containers::Vector<NarrowString> NarrowStrings;

namespace SemanticGraph = XSDFrontend::SemanticGraph;
namespace Transformations = XSDFrontend::Transformations;

using std::wcerr;
using std::endl;

namespace CLI
{
  using namespace Cult::CLI;

  typedef Char const Key[];

  extern Key help                = "help";
  extern Key version             = "version";
  extern Key proprietary_license = "proprietary-license";

  typedef Cult::CLI::Options
  <
    help,                Boolean,
    version,             Boolean,
    proprietary_license, Boolean
  >
  HelpOptions;

  struct HelpOptionsSpec: Cult::CLI::OptionsSpec<HelpOptions> {};


  extern Key disable_warning         = "disable-warning";
  extern Key sloc_limit              = "sloc-limit";
  extern Key morph_anonymous         = "morph-anonymous";
  extern Key preserve_anonymous      = "preserve-anonymous";
  extern Key anonymous_regex         = "anonymous-regex";
  extern Key anonymous_regex_trace   = "anonymous-regex-trace";
  extern Key location_map            = "location-map";
  extern Key location_regex          = "location-regex";
  extern Key location_regex_trace    = "location-regex-trace";
  extern Key custom_literals         = "custom-literals";
  extern Key file_per_type           = "file-per-type";
  extern Key type_file_regex         = "type-file-regex";
  extern Key type_file_regex_trace   = "type-file-regex-trace";
  extern Key schema_file_regex       = "schema-file-regex";
  extern Key schema_file_regex_trace = "schema-file-regex-trace";
  extern Key file_list               = "file-list";
  extern Key file_list_prologue      = "file-list-prologue";
  extern Key file_list_epilogue      = "file-list-epilogue";
  extern Key file_list_delim         = "file-list-delim";
  extern Key disable_multi_import    = "disable-multi-import"; // Undocumented.
  extern Key disable_full_check      = "disable-full-check";   // Undocumented.


  typedef Cult::CLI::Options
  <
    disable_warning,         Cult::Containers::Vector<NarrowString>,
    sloc_limit,              UnsignedLong,
    morph_anonymous,         Boolean,
    preserve_anonymous,      Boolean,
    anonymous_regex,         NarrowStrings,
    anonymous_regex_trace,   Boolean,
    location_map,            NarrowStrings,
    location_regex,          NarrowStrings,
    location_regex_trace,    Boolean,
    custom_literals,         NarrowString,
    file_per_type,           Boolean,
    type_file_regex,         NarrowStrings,
    type_file_regex_trace,   Boolean,
    schema_file_regex,       NarrowStrings,
    schema_file_regex_trace, Boolean,
    file_list,               NarrowString,
    file_list_prologue,      NarrowString,
    file_list_epilogue,      NarrowString,
    file_list_delim,         NarrowString,
    disable_multi_import,    Boolean,
    disable_full_check,      Boolean
  >
  CommonOptions;

  struct CommonOptionsSpec: Cult::CLI::OptionsSpec<CommonOptions> {};
}

//
//
struct LocationTranslator: XSDFrontend::LocationTranslator
{
  struct Failed {};

  LocationTranslator (NarrowStrings const& map,
                      NarrowStrings const& regex,
                      Boolean trace);

  virtual NarrowString
  translate (NarrowString const&);

private:
  typedef Cult::Containers::Map<NarrowString, NarrowString> Map;

  typedef BackendElements::Regex::Expression<Char> Regex;
  typedef BackendElements::Regex::Format<Char> RegexFormat;
  typedef Cult::Containers::Vector<Regex> RegexVector;

  typedef Cult::Containers::Map<NarrowString, NarrowString> Cache;

  Map map_;
  RegexVector regex_;
  Cache cache_;
  Boolean trace_;
};

//
//
struct AnonymousNameTranslator: Transformations::AnonymousNameTranslator
{
  struct Failed {};

  AnonymousNameTranslator (NarrowStrings const& regex, Boolean trace);

  virtual WideString
  translate (WideString const& file,
             WideString const& ns,
             WideString const& name,
             WideString const& xpath);

private:
  typedef BackendElements::Regex::Expression<WideChar> Regex;
  typedef BackendElements::Regex::Format<WideChar> RegexFormat;
  typedef Cult::Containers::Vector<Regex> RegexVector;

  RegexVector regex_;
  Boolean trace_;

};

//
//
struct SchemaPerTypeTranslator: Transformations::SchemaPerTypeTranslator
{
  struct Failed {};

  SchemaPerTypeTranslator (NarrowStrings const& type_regex,
                           Boolean type_trace,
                           NarrowStrings const& schema_regex,
                           Boolean schema_trace);

  virtual WideString
  translate_type (WideString const& ns, WideString const& name);

  virtual NarrowString
  translate_schema (NarrowString const& file);

private:
  typedef BackendElements::Regex::Expression<WideChar> TypeRegex;
  typedef BackendElements::Regex::Format<WideChar> TypeRegexFormat;
  typedef Cult::Containers::Vector<TypeRegex> TypeRegexVector;

  TypeRegexVector type_regex_;
  Boolean type_trace_;

  typedef BackendElements::Regex::Expression<Char> SchemaRegex;
  typedef BackendElements::Regex::Format<Char> SchemaRegexFormat;
  typedef Cult::Containers::Vector<SchemaRegex> SchemaRegexVector;

  SchemaRegexVector schema_regex_;
  Boolean schema_trace_;
};

//
//
struct XercesInitializer
{
  XercesInitializer ()
  {
    xercesc::XMLPlatformUtils::Initialize ();
  }

  ~XercesInitializer ()
  {
    xercesc::XMLPlatformUtils::Terminate ();
  }
};

// Expand the \n escape sequence.
//
Void
expand_nl (NarrowString& s);


Int
main (Int argc, Char* argv[])
{
  std::wostream& e (wcerr);

  Cult::Trace::Log::instance ().level (0);

  try
  {
    CLI::FileArguments args (argc, argv, "--options-file");

    CLI::HelpOptions help_options (
      CLI::parse (CLI::HelpOptionsSpec (), args, CLI::UnknownMode::stop));

    NarrowString cmd;

    if (args.size () > 1)
    {
      cmd = args[1];
      args.erase (1);
    }

    if (help_options.value<CLI::version> () || cmd == "version")
    {
      e << "CodeSynthesis XSD XML Schema to C++ compiler " <<
        XSD_STR_VERSION << endl
        << "Copyright (C) 2005-2010 Code Synthesis Tools CC" << endl;

      if (!help_options.value<CLI::proprietary_license> () &&
          cmd == "version")
      {
        // Parse the options after the command to detect trailing
        // --proprietary-license.
        //
        help_options = CLI::parse (
          CLI::HelpOptionsSpec (), args, CLI::UnknownMode::stop);
      }

      if (help_options.value<CLI::proprietary_license> ())
      {
        e << "The compiler was invoked in the Proprietary License mode. You "
          << "should have\nreceived a proprietary license from Code Synthesis "
          << "Tools CC that entitles\nyou to use it in this mode." << endl;
      }
      else
      {
        e << "This is free software; see the source for copying conditions. "
          << "There is NO\nwarranty; not even for MERCHANTABILITY or FITNESS "
          << "FOR A PARTICULAR PURPOSE." << endl;
      }

      return 0;
    }

    if (help_options.value<CLI::help> () || cmd == "help")
    {
      if (cmd == "help" && args.size () > 1)
      {
        NarrowString arg (args[1]);

        if (arg == "cxx-tree")
        {
          e << "Usage: " << args[0] << " cxx-tree [options] file [file ...]"
            << endl
            << "Options:" << endl;

          CXX::Tree::Generator::usage ();
        }
        else if (arg == "cxx-parser")
        {
          e << "Usage: " << args[0] << " cxx-parser [options] file [file ...]"
            << endl
            << "Options:" << endl;

          CXX::Parser::Generator::usage ();
        }
        else
        {
          e << "error: unknown command '" << arg.c_str () << "'" << endl
            << "info: try '" << args[0] << " help' for the list of commands"
            << endl;

          return 1;
        }

        ::CLI::Indent::Clip< ::CLI::OptionsUsage, WideChar> clip (e);

        // Disable warning option.
        //
        e << "--disable-warning <warn>" << endl
          << " Disable printing warning with id <warn>. If 'all'\n"
          << " is specified for the warning id then all warnings\n"
          << " are disabled."
          << endl;

        // Anonymous morphing options.
        //
        e << "--preserve-anonymous" << endl
          << " Preserve anonymous types. By default anonymous\n"
          << " types are automatically named with names derived\n"
          << " from the enclosing elements/attributes."
          << endl;

        e << "--anonymous-regex <regex>" << endl
          << " Add the provided regular expression to the list of\n"
          << " regular expressions used to derive names for\n"
          << " anonymous types from the names of the enclosing\n"
          << " attributes/elements."
          << endl;

        e << "--anonymous-regex-trace" << endl
          << " Trace the process of applying regular expressions\n"
          << " specified with the --anonymous-regex option."
          << endl;

        // Location mapping options.
        //
        e << "--location-map <ol>=<nl>" << endl
          << " Map the original schema location <ol> that is\n"
          << " specified in the XML Schema include or import\n"
          << " elements to new schema location <nl>. Repeat\n"
          << " this option to map more than one schema location."
          << endl;

        e << "--location-regex <regex>" << endl
          << " Add <regex> to the list of regular expressions\n"
          << " used to map schema locations that are specified\n"
          << " in the XML Schema include or import elements."
          << endl;

        e << "--location-regex-trace" << endl
          << " Trace the process of applying regular expressions\n"
          << " specified with the --location-regex option."
          << endl;

        // File-per-type compilation mode options.
        //
        e << "--file-per-type" << endl
          << " Generate a separate set of C++ files for each\n"
          << " type defined in XML Schema."
          << endl;

        e << "--type-file-regex <regex>" << endl
          << " Add the provided regular expression to the list of\n"
          << " regular expressions used to translate type names\n"
          << " to file names when the --file-per-type option is\n"
          << " specified."
          << endl;

        e << "--type-file-regex-trace" << endl
          << " Trace the process of applying regular expressions\n"
          << " specified with the --type-file-regex option."
          << endl;

        e << "--schema-file-regex <regex>" << endl
          << " Add the provided regular expression to the list\n"
          << " of regular expressions used to translate schema\n"
          << " file names when the --file-per-type option is\n"
          << " specified."
          << endl;

        e << "--schema-file-regex-trace" << endl
          << " Trace the process of applying regular expressions\n"
          << " specified with the --schema-file-regex option."
          << endl;

        // File list options.
        //
        e << "--file-list <file>" << endl
          << " Write a list of generated C++ files to <file>."
          << endl;

        e << "--file-list-prologue <p>" << endl
          << " Insert <p> at the beginning of the file list. All\n"
          << " occurrences of the \\n character sequence in <p>\n"
          << " are replaced with new lines."
          << endl;

        e << "--file-list-prologue <e>" << endl
          << " Insert <e> at the end of the file list. All\n"
          << " occurrences of the \\n character sequence in <e>\n"
          << " are replaced with new lines."
          << endl;

        e << "--file-list-delim <d>" << endl
          << " Delimit file names written to the file list with\n"
          << " <d> instead of new lines. All occurrences of the\n"
          << " \\n character sequence in <d> are replaced with\n"
          << " new lines."
          << endl;
      }
      else
      {
        e << "Usage: " << args[0] << " <cmd> ..." << endl
          << "Commands:" << endl;

        e << "  help            Print usage information and exit. Use\n"
          << "                  'help <cmd>' for command-specific options."
          << endl;

        e << "  version         Print version and exit."
          << endl;

        e << "  cxx-tree        Generate the C++/Tree mapping."
          << endl;

        e << "  cxx-parser      Generate the C++/Parser mapping."
          << endl;
      }

      return 0;
    }

    if (cmd.empty ())
    {
      e << "error: no command specified" << endl
        << "info: try '" << args[0] << " help' for usage information" << endl;

      return 1;
    }

    if (cmd != "cxx-tree" && cmd != "cxx-parser")
    {
      e << "error: unknown command '" << cmd.c_str () << "'" << endl
        << "info: try '" << args[0] << " help' for the list of commands"
        << endl;

      return 1;
    }

    // We need to parse command line options before we can get to
    // the arguments.
    //
    CLI::CommonOptionsSpec common_spec;
    common_spec.option<CLI::file_list_delim> ().default_value ("\n");

    CLI::CommonOptions common_ops (
      CLI::parse (
        common_spec,
        args,
        CLI::UnknownMode::skip,
        CLI::UnknownMode::skip));

    WarningSet disabled_w;
    {
      typedef Cult::Containers::Vector<NarrowString> Warnings;
      Warnings const& w (common_ops.value<CLI::disable_warning> ());

      for (Warnings::ConstIterator i (w.begin ()); i != w.end (); ++i)
        disabled_w.insert (*i);
    }

    Boolean disabled_w_all (disabled_w.find ("all") != disabled_w.end ());

    if (common_ops.value<CLI::morph_anonymous> () &&
        !disabled_w_all && disabled_w.find ("D001") == disabled_w.end ())
    {
      e << "warning D001: the --morph-anonymous option is on by default and "
        << "no longer required"
        << endl;
    }

    Evptr<CXX::Tree::CLI::Options> tree_ops;
    Evptr<CXX::Parser::CLI::Options> parser_ops;

    Boolean show_sloc (false);

    if (cmd == "cxx-tree")
    {
      tree_ops = new CXX::Tree::CLI::Options (
        CLI::parse (CXX::Tree::Generator::options_spec (), args));

      tree_ops->value<CXX::Tree::CLI::disable_multi_import> () =
        common_ops.value<CLI::disable_multi_import> ();

      show_sloc = tree_ops->value<CXX::Tree::CLI::show_sloc> ();
    }
    else if (cmd == "cxx-parser")
    {
      parser_ops = new CXX::Parser::CLI::Options (
        CLI::parse (CXX::Parser::Generator::options_spec (), args));

      show_sloc = parser_ops->value<CXX::Parser::CLI::show_sloc> ();
    }


    if (args.size () < 2)
    {
      e << "error: no input file specified" << endl;
      return 1;
    }

    Boolean fpt (common_ops.value<CLI::file_per_type> ());

    if (cmd == "cxx-tree" || cmd == "cxx-parser")
    {
      Boolean gen (false), use (false);

      if (cmd == "cxx-tree")
      {
        gen = tree_ops->value<CXX::Tree::CLI::generate_xml_schema> ();
        use = tree_ops->value<CXX::Tree::CLI::extern_xml_schema> ();
      }
      else if (cmd == "cxx-parser")
      {
        gen = parser_ops->value<CXX::Parser::CLI::generate_xml_schema> ();
        use = parser_ops->value<CXX::Parser::CLI::extern_xml_schema> ();
      }

      // Things get complicated when we are compiling several schemas at
      // once (non-file-per-type mode) and use the --generate-xml-schema/
      // --extern-xml-schema options. The only way we can figure out which
      // file corresponds to XML Schema is if the --extern-xml-schema option
      // is also present. So we are going to require it for this case,
      // especially since it generally makes sense.
      //
      if (!fpt)
      {
        if (args.size () > 2 && gen && !use)
        {
          e << "error: --extern-xml-schema is required when compiling more "
            << "than one schema and --generate-xml-schema is specified"
            << endl;

          return 1;
        }

        if (args.size () == 2 && gen && use)
        {
          e << "error: --generate-xml-schema and --extern-xml-schema are "
            << "mutually exclusive when compiling a single schema" << endl;

          return 1;
        }
      }
      else
      {
        // The --file-per-type and --generate-xml-schema options are
        // incompatible. It also makes sense to use --file-per-type
        // and --extern-xml-schema.
        //
        if (gen)
        {
          e << "error: --file-per-type and --generate-xml-schema are "
            << "incompatible" << endl
            << "info: use --generate-xml-schema in a separate invocation "
            << "of the compiler" << endl;

          return 1;
        }

        if (!use &&
            !disabled_w_all && disabled_w.find ("D002") == disabled_w.end ())
        {
          e << "warning D002: --extern-xml-schema is recommended when "
            << "--file-per-type is specified to reduce generated code size"
            << endl;
        }
      }
    }

    //
    //
    FileList file_list;
    AutoUnlinks unlinks;
    UnsignedLong sloc (0);

    LocationTranslator loc_translator (
      common_ops.value<CLI::location_map> (),
      common_ops.value<CLI::location_regex> (),
      common_ops.value<CLI::location_regex_trace> ());

    AnonymousNameTranslator anon_translator (
      common_ops.value<CLI::anonymous_regex> (),
      common_ops.value<CLI::anonymous_regex_trace> ());

    // Load custom string literals, if any.
    //
    CXX::StringLiteralMap string_literal_map;

    if (NarrowString file = common_ops.value<CLI::custom_literals> ())
    {
      XercesInitializer xerces_init;

      if (!CXX::read_literal_map (file, string_literal_map))
      {
        // Diagnostics has already been issued.
        //
        return 1;
      }
    }

    if (!fpt)
    {
      // File-per-schema compilation mode.
      //

      for (Size i (1); i < args.size (); ++i)
      {
        // Parse schema.
        //
        SemanticGraph::Path tu;

        try
        {
#if !defined(BOOST_FILESYSTEM_VERSION) || BOOST_FILESYSTEM_VERSION == 2
          tu = SemanticGraph::Path (args[i], boost::filesystem::native);
#else
          tu = SemanticGraph::Path (args[i]);
#endif
        }
        catch (SemanticGraph::InvalidPath const&)
        {
          e << "error: '" << args[i] << "' is not a valid "
            << "filesystem path" << endl;

          return 1;
        }

        XSDFrontend::Parser parser (
          cmd != "cxx-tree",
          !common_ops.value<CLI::disable_multi_import> (),
          !common_ops.value<CLI::disable_full_check> (),
          loc_translator,
          disabled_w);

        Evptr<SemanticGraph::Schema> schema;

        if (cmd == "cxx-tree" || cmd == "cxx-parser")
        {
          // See if we are generating code for the XML Schema namespace.
          // We could be compiling several schemas at once in which case
          // handling of the --generate-xml-schema option gets tricky: we
          // will need to rely on the presence of the --extern-xml-schema
          // to tell us which (fake) schema file corresponds to XML Schema.
          //
          Boolean gen_xml_schema (false);

          if (cmd == "cxx-tree")
          {
            gen_xml_schema =
              tree_ops->value<CXX::Tree::CLI::generate_xml_schema> ();

            if (gen_xml_schema)
            {
              if (NarrowString name =
                  tree_ops->value<CXX::Tree::CLI::extern_xml_schema> ())
              {
#if !defined(BOOST_FILESYSTEM_VERSION) || BOOST_FILESYSTEM_VERSION == 2
                if (tu.native_file_string () != name)
#else
                if (tu.string () != name)
#endif
                  gen_xml_schema = false;
              }
            }
          }
          else if (cmd == "cxx-parser")
          {
            gen_xml_schema =
              parser_ops->value<CXX::Parser::CLI::generate_xml_schema> ();

            if (gen_xml_schema)
            {
              if (NarrowString name =
                  parser_ops->value<CXX::Parser::CLI::extern_xml_schema> ())
              {
#if !defined(BOOST_FILESYSTEM_VERSION) || BOOST_FILESYSTEM_VERSION == 2
                if (tu.native_file_string () != name)
#else
                if (tu.string () != name)
#endif
                  gen_xml_schema = false;
              }
            }
          }

          if (gen_xml_schema)
            schema = parser.xml_schema (tu);
          else
            schema = parser.parse (tu);
        }
        else
          schema = parser.parse (tu);

        // Morph anonymous types.
        //
        if (!common_ops.value<CLI::preserve_anonymous> ())
        {
          try
          {
            Transformations::Anonymous trans (anon_translator);
            trans.transform (*schema, tu, true);
          }
          catch (Transformations::Anonymous::Failed const&)
          {
            return 1; // Diagnostic has already been issued.
          }
        }

        // Synthesize enumerations from unions.
        //
        if (cmd == "cxx-tree")
        {
          Transformations::EnumSynthesis trans;
          trans.transform (*schema, tu);
        }

        // Simplify the schema graph.
        //
        if (cmd == "cxx-parser")
        {
          Transformations::Simplifier trans;
          trans.transform (*schema, tu);
        }

        // Try to rearrange definitions so that there is no forward
        // inheritance.
        //
        try
        {
          Processing::Inheritance::Processor proc;
          proc.process (*schema, tu);
        }
        catch (Processing::Inheritance::Processor::Failed const&)
        {
          return 1; // Diagnostic has already been issued.
        }

        // Normalize and annotate complex content restrictions.
        //
        if (cmd == "cxx-parser")
        {
          try
          {
            Transformations::Restriction trans;
            trans.transform (*schema, tu);
          }
          catch (Transformations::Restriction::Failed const&)
          {
            return 1; // Diagnostic has already been issued.
          }
        }

        // Calculate cardinality.
        //
        {
          Processing::Cardinality::Processor proc;
          proc.process (*schema, tu);
        }

        // Generate mapping.
        //
        if (cmd == "cxx-tree")
        {
          try
          {
            sloc += CXX::Tree::Generator::generate (
              *tree_ops,
              *schema,
              tu,
              false,
              string_literal_map,
              disabled_w,
              file_list,
              unlinks);
          }
          catch (CXX::Tree::Generator::Failed const&)
          {
            // Diagnostic has already been issued.
            //
            return 1;
          }
        }
        else if (cmd == "cxx-parser")
        {
          try
          {
            sloc += CXX::Parser::Generator::generate (
              *parser_ops,
              *schema,
              tu,
              false,
              string_literal_map,
              true,
              disabled_w,
              file_list,
              unlinks);
          }
          catch (CXX::Parser::Generator::Failed const&)
          {
            // Diagnostic has already been issued.
            //
            return 1;
          }
        }
      }
    }
    else
    {
      // File-per-type compilation mode.
      //
      SemanticGraph::Paths paths;

      for (Size i (1); i < args.size (); ++i)
      {
        try
        {
#if !defined(BOOST_FILESYSTEM_VERSION) || BOOST_FILESYSTEM_VERSION == 2
          paths.push_back (
            SemanticGraph::Path (args[i], boost::filesystem::native));
#else
          paths.push_back (
            SemanticGraph::Path (args[i]));
#endif
        }
        catch (SemanticGraph::InvalidPath const&)
        {
          e << "error: '" << args[i] << "' is not a valid "
            << "filesystem path" << endl;

          return 1;
        }
      }

      if (cmd == "cxx-parser" &&
          paths.size () > 1 &&
          parser_ops->value<CXX::Parser::CLI::generate_test_driver> ())
      {
        e << "info: generating test driver for the first schema only: '" <<
          paths[0] << "'" << endl;
      }

      XSDFrontend::Parser parser (
        cmd != "cxx-tree",
        !common_ops.value<CLI::disable_multi_import> (),
        !common_ops.value<CLI::disable_full_check> (),
        loc_translator,
        disabled_w);

      Evptr<SemanticGraph::Schema> schema (parser.parse (paths));

      // Morph anonymous types.
      //
      if (!common_ops.value<CLI::preserve_anonymous> ())
      {
        try
        {
          Transformations::Anonymous trans (anon_translator);
          trans.transform (*schema, "", false);
        }
        catch (Transformations::Anonymous::Failed const&)
        {
          return 1; // Diagnostic has already been issued.
        }
      }

      // Synthesize enumerations from unions.
      //
      if (cmd == "cxx-tree")
      {
        Transformations::EnumSynthesis trans;
        trans.transform (*schema, "");
      }

      // Simplify the schema graph.
      //
      if (cmd == "cxx-parser")
      {
        Transformations::Simplifier trans;
        trans.transform (*schema, "");
      }

      // Normalize and annotate complex content restrictions.
      //
      if (cmd == "cxx-parser")
      {
        try
        {
          Transformations::Restriction trans;
          trans.transform (*schema, "");
        }
        catch (Transformations::Restriction::Failed const&)
        {
          return 1; // Diagnostic has already been issued.
        }
      }

      // Calculate cardinality.
      //
      {
        Processing::Cardinality::Processor proc;
        proc.process (*schema, "");
      }

      // Rearrange the graph so that each type is in a seperate
      // schema file.
      //
      typedef Cult::Containers::Vector<SemanticGraph::Schema*> Schemas;

      SchemaPerTypeTranslator type_translator (
        common_ops.value<CLI::type_file_regex> (),
        common_ops.value<CLI::type_file_regex_trace> (),
        common_ops.value<CLI::schema_file_regex> (),
        common_ops.value<CLI::schema_file_regex_trace> ());

      Transformations::SchemaPerType trans (type_translator);
      Schemas schemas (trans.transform (*schema));

      // Generate code.
      //
      for (Schemas::Iterator b (schemas.begin ()), i (b), e (schemas.end ());
           i != e; ++i)
      {
        SemanticGraph::Schema& s (**i);
        SemanticGraph::Path path (
          s.context ().count ("renamed")
          ? s.context ().get<SemanticGraph::Path> ("renamed")
          : s.used_begin ()->path ());

        if (cmd == "cxx-tree")
        {
          try
          {
            sloc += CXX::Tree::Generator::generate (
              *tree_ops,
              s,
              path,
              true,
              string_literal_map,
              disabled_w,
              file_list,
              unlinks);
          }
          catch (CXX::Tree::Generator::Failed const&)
          {
            // Diagnostic has already been issued.
            //
            return 1;
          }
        }
        else if (cmd == "cxx-parser")
        {
          try
          {
            // Only generate driver for the first schema.
            //
            sloc += CXX::Parser::Generator::generate (
              *parser_ops,
              s,
              path,
              true,
              string_literal_map,
              i == b,
              disabled_w,
              file_list,
              unlinks);
          }
          catch (CXX::Parser::Generator::Failed const&)
          {
            // Diagnostic has already been issued.
            //
            return 1;
          }
        }
      }
    }

    // See if we need to produce the file list.
    //
    if (NarrowString fl = common_ops.value<CLI::file_list> ())
    {
      typedef boost::filesystem::ofstream OutputFileStream;

      try
      {
        OutputFileStream ofs;
#if !defined(BOOST_FILESYSTEM_VERSION) || BOOST_FILESYSTEM_VERSION == 2
        SemanticGraph::Path path (fl);
#else
        SemanticGraph::Path path (fl.c_str());
#endif

#if !defined(BOOST_FILESYSTEM_VERSION) || BOOST_FILESYSTEM_VERSION == 2
        ofs.open (fl, std::ios_base::out);
#else
        ofs.open (fl.c_str(), std::ios_base::out);
#endif

        if (!ofs.is_open ())
        {
          wcerr << path << ": error: unable to open in write mode" << endl;
          return 1;
        }

        NarrowString d (common_ops.value<CLI::file_list_delim> ());
        expand_nl (d);

        if (NarrowString p = common_ops.value<CLI::file_list_prologue> ())
        {
          expand_nl (p);
          ofs << p;
        }

        for (FileList::Iterator i (file_list.begin ()), e (file_list.end ());
             i != e;)
        {
          ofs << *i;

          if (++i != e)
            ofs << d;
        }

        if (NarrowString e = common_ops.value<CLI::file_list_epilogue> ())
        {
          expand_nl (e);
          ofs << e;
        }
      }
      catch (SemanticGraph::InvalidPath const&)
      {
        wcerr << "error: '" << fl.c_str () << "' is not a valid "
              << "filesystem path" << endl;
        return 1;
      }
    }

    if (show_sloc)
      e << "total: " << sloc << endl;

    if (UnsignedLong sloc_limit = common_ops.value<CLI::sloc_limit> ())
    {
      if (sloc_limit < sloc)
      {
        e << "error: SLOC limit of " << sloc_limit
          << " lines has been exceeded" << endl;

        return 1;
      }
    }

    unlinks.cancel ();

    return 0;
  }
  catch (LocationTranslator::Failed const&)
  {
    // Diagnostic has already been issued.
  }
  catch (AnonymousNameTranslator::Failed const&)
  {
    // Diagnostic has already been issued.
  }
  catch (SchemaPerTypeTranslator::Failed const&)
  {
    // Diagnostic has already been issued.
  }
  catch (Transformations::SchemaPerType::Failed const&)
  {
    // Diagnostic has already been issued.
  }
  catch (XSDFrontend::InvalidSchema const&)
  {
    // Diagnostic has already been issued.
  }
  catch (CLI::UnexpectedOption const& e)
  {
    wcerr << "error: unknown option '" << e.option ().c_str () << "'" << endl
          << "info: try '" << argv[0] << " help' for usage information"
          << endl;
  }
  catch (CLI::OptionFormat const& e)
  {
    wcerr << "error: value for option '" << e.option ().c_str ()
          << "' is invalid or missing" << endl
          << "info: try '" << argv[0] << " help' for usage information"
          << endl;
  }
  catch (CLI::OptionFile const& e)
  {
    if (e.value ())
      wcerr << "error: " << e.value ().c_str () << ": "
            << e.description ().c_str () << endl;
    else
      wcerr << "error: missing --options-file argument" << endl;
  }

  return 1;
}

// LocationTranslator
//

LocationTranslator::
LocationTranslator (NarrowStrings const& map,
                    NarrowStrings const& regex,
                    Boolean trace)
    : trace_ (trace)
{
  // Map.
  //
  for (NarrowStrings::ConstIterator i (map.begin ()); i != map.end (); ++i)
  {
    // Split the string in two parts at the last '='.
    //
    Size pos (i->rfind ('='));

    if (pos == NarrowString::npos)
    {
      wcerr << "error: invalid location map: '" << i->c_str () <<
        "': delimiter ('=') not found" << endl;

      throw Failed ();
    }

    map_[NarrowString (*i, 0, pos)] = NarrowString (*i, pos + 1);
  }

  // Regex.
  //
  for (NarrowStrings::ConstIterator i (regex.begin ()); i != regex.end (); ++i)
  {
    try
    {
      regex_.push_back (Regex (*i));
    }
    catch (RegexFormat const& e)
    {
      wcerr << "error: invalid location regex: '" <<
        e.expression ().c_str () << "': " <<
        e.description ().c_str () << endl;

      throw Failed ();
    }
  }
}

NarrowString LocationTranslator::
translate (NarrowString const& l)
{
  // First check the cache.
  //
  Cache::ConstIterator ci (cache_.find (l));

  if (ci != cache_.end ())
    return ci->second;

  // Then check the direct map.
  //
  Map::ConstIterator mi (map_.find (l));

  if (mi != map_.end ())
  {
    cache_[l] = mi->second;
    return mi->second;
  }

  // Finally try regex.
  //
  if (trace_)
    wcerr << "location: '" << l.c_str () << "'" << endl;

  for (RegexVector::ReverseIterator i (regex_.rbegin ());
       i != regex_.rend (); ++i)
  {
    if (trace_)
      wcerr << "try: '" << i->pattern () << "' : ";

    if (i->match (l))
    {
      NarrowString r (i->merge (l));

      if (trace_)
        wcerr << "'" << r.c_str () << "' : +" << endl;

      cache_[l] = r;
      return r;
    }

    if (trace_)
      wcerr << '-' << endl;
  }

  // No match - return the original location.
  //
  cache_[l] = l;
  return l;
}

// AnonymousNameTranslator
//

AnonymousNameTranslator::
AnonymousNameTranslator (NarrowStrings const& regex, Boolean trace)
    : trace_ (trace)
{
  for (NarrowStrings::ConstIterator i (regex.begin ()); i != regex.end (); ++i)
  {
    try
    {
      regex_.push_back (Regex (*i));
    }
    catch (RegexFormat const& e)
    {
      wcerr << "error: invalid anonymous type regex: '" <<
        e.expression () << "': " << e.description () << endl;

      throw Failed ();
    }
  }
}

WideString AnonymousNameTranslator::
translate (WideString const& file,
           WideString const& ns,
           WideString const& name,
           WideString const& xpath)
{
  WideString s (file + L' ' + ns + L' ' + xpath);

  if (trace_)
    wcerr << "anonymous type: '" << s << "'" << endl;

  for (RegexVector::ReverseIterator i (regex_.rbegin ());
       i != regex_.rend (); ++i)
  {
    if (trace_)
      wcerr << "try: '" << i->pattern () << "' : ";

    if (i->match (s))
    {
      WideString r (i->merge (s));

      if (trace_)
        wcerr << "'" << r << "' : +" << endl;

      return r;
    }

    if (trace_)
      wcerr << '-' << endl;
  }

  // No match - return the name.
  //
  return name;
}

// SchemaPerTypeTranslator
//

SchemaPerTypeTranslator::
SchemaPerTypeTranslator (NarrowStrings const& type_regex,
                         Boolean type_trace,
                         NarrowStrings const& schema_regex,
                         Boolean schema_trace)
    : type_trace_ (type_trace), schema_trace_ (schema_trace)
{
  for (NarrowStrings::ConstIterator i (type_regex.begin ());
       i != type_regex.end (); ++i)
  {
    try
    {
      type_regex_.push_back (TypeRegex (*i));
    }
    catch (TypeRegexFormat const& e)
    {
      wcerr << "error: invalid type file regex: '" <<
        e.expression () << "': " << e.description () << endl;

      throw Failed ();
    }
  }

  for (NarrowStrings::ConstIterator i (schema_regex.begin ());
       i != schema_regex.end (); ++i)
  {
    try
    {
      schema_regex_.push_back (SchemaRegex (*i));
    }
    catch (SchemaRegexFormat const& e)
    {
      wcerr << "error: invalid type file regex: '" <<
        e.expression ().c_str () << "': " << e.description ().c_str () << endl;

      throw Failed ();
    }
  }
}

WideString SchemaPerTypeTranslator::
translate_type (WideString const& ns, WideString const& name)
{
  WideString s (ns + L' ' + name);

  if (type_trace_)
    wcerr << "type: '" << s << "'" << endl;

  for (TypeRegexVector::ReverseIterator i (type_regex_.rbegin ());
       i != type_regex_.rend (); ++i)
  {
    if (type_trace_)
      wcerr << "try: '" << i->pattern () << "' : ";

    if (i->match (s))
    {
      WideString r (i->merge (s));

      if (type_trace_)
        wcerr << "'" << r << "' : +" << endl;

      return r;
    }

    if (type_trace_)
      wcerr << '-' << endl;
  }

  // No match - return empty string.
  //
  return L"";
}

NarrowString SchemaPerTypeTranslator::
translate_schema (NarrowString const& file)
{
  if (schema_trace_)
    wcerr << "schema: '" << file.c_str () << "'" << endl;

  for (SchemaRegexVector::ReverseIterator i (schema_regex_.rbegin ());
       i != schema_regex_.rend (); ++i)
  {
    if (schema_trace_)
      wcerr << "try: '" << i->pattern () << "' : ";

    if (i->match (file))
    {
      NarrowString r (i->merge (file));

      if (schema_trace_)
        wcerr << "'" << r.c_str () << "' : +" << endl;

      return r;
    }

    if (schema_trace_)
      wcerr << '-' << endl;
  }

  // No match - return empty string.
  //
  return "";
}

//
//
Void
expand_nl (NarrowString& s)
{
  for (Size i (0); i < s.size ();)
  {
    if (s[i] == '\\' && (i + 1) < s.size () && s[i + 1] == 'n')
    {
      NarrowString tmp (s, 0, i);
      tmp += '\n';
      tmp.append (s.c_str () + i + 2);
      s = tmp;
    }
    else
      ++i;
  }
}
