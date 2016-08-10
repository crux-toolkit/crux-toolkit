#include "ArgParser.h"
#include "StringUtils.h"

using namespace std;

ArgParser::ArgParser() {
}

void ArgParser::Parse(int argc, char** argv, const vector<string>& args) {
  args_.clear();
  options_.clear();

  vector<ArgSpec> argSpecs = ArgStringsToArgSpecs(args);
  int multi = -1; // Which argument accepts multiple values?
  for (size_t i = 0; i < argSpecs.size(); i++) {
    if (!argSpecs[i].IsMulti()) {
      continue;
    } else if (multi != -1) {
      throw runtime_error("Cannot have multiple arguments that accept "
                          "multiple values.");
    }
    multi = i;
  }

  vector<string> parsedArgs;
  // Start at 1 to skip the application name
  for (int i = 1; i < argc; i++) {
    string arg(argv[i]);
    if (StringUtils::StartsWith(arg, "--")) {
      // This is a option name
      string option = arg.substr(2);
      if (++i == argc) {
        throw ArgParserException("No value found for option '" + option + "'");
      }
      options_[option] = argv[i];
    } else {
      // This is an argument
      parsedArgs.push_back(arg);
    }
  }

  int excess_args = parsedArgs.size() - argSpecs.size();
  if (multi == -1 && excess_args != 0) {
    throw ArgParserException("Expected " + StringUtils::ToString(argSpecs.size()) +
                             " arguments, but found " +
                             StringUtils::ToString(parsedArgs.size()),
                             argc == 1);
  } else if (multi != -1 && excess_args < 0) {
    throw ArgParserException("Expected at least " +
                             StringUtils::ToString(argSpecs.size()) +
                             " arguments, but found " +
                             StringUtils::ToString(parsedArgs.size()),
                             argc == 1);
  }

  int curParsedArg = -1;
  for (vector<ArgSpec>::const_iterator i = argSpecs.begin(); i != argSpecs.end(); i++) {
    string argName = i->GetName();
    args_[argName] = vector<string>(1, parsedArgs[++curParsedArg]);
    if (i->IsMulti()) {
      for (int j = 0; j < excess_args; j++) {
        args_[argName].push_back(parsedArgs[++curParsedArg]);
      }
    }
  }
}

ArgParser::~ArgParser() {
}

const map< string, vector<string> >& ArgParser::GetArgs() const {
  return args_;
}

const map<string, string>& ArgParser::GetOptions() const {
  return options_;
}

const string ArgParser::GetArg(const string& name) const {
  map< string, vector<string> >::const_iterator i = args_.find(name);
  return (i != args_.end()) ? i->second.front() : "";
}

const vector<string> ArgParser::GetArgMulti(const string& name) const {
  map< string, vector<string> >::const_iterator i = args_.find(name);
  return (i != args_.end()) ? i->second : vector<string>();
}

const string ArgParser::GetOption(const string& name) const {
  map<string, string>::const_iterator i = options_.find(name);
  return (i != options_.end()) ? i->second : "";
}

ArgParser::ArgSpec::ArgSpec(const string& name, bool multi)
  : name_(name), multi_(multi) {
}

string ArgParser::ArgSpec::GetName() const {
  return name_;
}

bool ArgParser::ArgSpec::IsMulti() const {
  return multi_;
}

vector<ArgParser::ArgSpec> ArgParser::ArgStringsToArgSpecs(const vector<string>& args) {
  vector<ArgSpec> specs;
  for (vector<string>::const_iterator i = args.begin(); i != args.end(); i++) {
    string arg = *i;
    size_t plus = arg.find('+');
    if (plus == i->length() - 1) {
      specs.push_back(ArgSpec(arg.substr(0, plus), true));
    } else {
      specs.push_back(ArgSpec(arg, false));
    }
  }
  return specs;
}

ArgParserException::ArgParserException(const string& what, bool fullUsage)
  : runtime_error(what), fullUsage_(fullUsage) {
}

bool ArgParserException::ShowFullUsage() const {
  return fullUsage_;
}

