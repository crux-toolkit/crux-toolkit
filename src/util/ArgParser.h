#ifndef ARGPARSER_H
#define ARGPARSER_H

#include <map>
#include <stdexcept>
#include <string>
#include <vector>

class ArgParser {
 public:
  ArgParser();
  ~ArgParser();

  void Parse(int argc, char** argv, const std::vector<std::string>& args);
  const std::map< std::string, std::vector<std::string> >& GetArgs() const;
  const std::map<std::string, std::string>& GetOptions() const;
  const std::string GetArg(const std::string& name) const;
  const std::vector<std::string> GetArgMulti(const std::string& name) const;
  const std::string GetOption(const std::string& name) const;

 protected:
  class ArgSpec {
   public:
    explicit ArgSpec(const std::string& name, bool multi = false);
    std::string GetName() const;
    bool IsMulti() const;
   protected:
    std::string name_;
    bool multi_;
  };

  static std::vector<ArgSpec> ArgStringsToArgSpecs(const std::vector<std::string>& args);

  std::map< std::string, std::vector<std::string> > args_;
  std::map<std::string, std::string> options_;
};

class ArgParserException : public std::runtime_error {
 public:
  explicit ArgParserException(const std::string& what, bool fullUsage = false);
  bool ShowFullUsage() const;
 private:
  bool fullUsage_;
};

#endif

