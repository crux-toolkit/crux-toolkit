#ifndef CREATEDOCS_H
#define CREATEDOCS_H

#include <map>
#include <ostream>
#include <string>
#include <vector>

#include "CruxApplicationList.h"

class CreateDocs : public CruxApplication {
 public:
  CreateDocs();
  virtual ~CreateDocs();

  virtual int main(int argc, char** argv);
  virtual std::string getName() const;
  virtual std::string getDescription() const;
  virtual std::vector<std::string> getArgs() const;
  virtual std::vector<std::string> getOptions() const;
  virtual std::vector< std::pair<std::string, std::string> > getOutputs() const;
  virtual bool needsOutputDirectory() const;
  virtual bool hidden() const;

 protected:
  void checkParams(const CruxApplicationList* apps = NULL);

  void makeParamTable(const CruxApplicationList* apps);
  void generateToolHtml(std::ostream* outStream, const CruxApplication* application);

  void makeReplacements(
    std::string* templateStr,
    const std::map<std::string, std::string>& replacements
  );

  static const std::string TOOL_TEMPLATE;
  static const std::string TOOL_INPUT_TEMPLATE;
  static const std::string TOOL_OUTPUT_TEMPLATE;
  static const std::string TOOL_OPTION_CATEGORY_TEMPLATE;
  static const std::string TOOL_NO_OPTIONS_TEMPLATE;
  static const std::string TOOL_OPTION_TEMPLATE;
};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

