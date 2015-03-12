#include <algorithm>
#include <cctype>
#include <functional>
#include <fstream>
#include <iostream>

#include "CreateDocs.h"
#include "CruxApplicationList.h"
#include "util/FileUtils.h"
#include "util/Params.h"

#include "qranker-barista/Barista.h"
#include "ComputeQValues.h"
#include "CruxBullseyeApplication.h"
#include "CruxHardklorApplication.h"
#include "ExtractColumns.h"
#include "ExtractRows.h"
#include "GeneratePeptides.h"
#include "GetMs2Spectrum.h"
#include "PercolatorApplication.h"
#include "PredictPeptideIons.h"
#include "PrintProcessedSpectra.h"
#include "qranker-barista/QRanker.h"
#include "ReadTideIndex.h"
#include "xlink/SearchForXLinks.h"
#include "SortColumn.h"
#include "SpectralCounts.h"
#include "StatColumn.h"
#include "TideIndexApplication.h"
#include "TideSearchApplication.h"
#include "CometApplication.h"

using namespace std;

CreateDocs::CreateDocs() {
  htmlEntities_['&'] = "&amp;";
  htmlEntities_['"'] = "&quot;";
  htmlEntities_['\''] = "&#039;";
  htmlEntities_['<'] = "&lt;";
  htmlEntities_['>'] = "&gt;";
}

CreateDocs::~CreateDocs() {
}

int CreateDocs::main(int argc, char** argv) {
  initialize(argc, argv);

  CruxApplicationList apps("crux");
  apps.add(new Barista());
  apps.add(new CometApplication());
  apps.add(new ComputeQValues());
  apps.add(new CreateDocs());
  apps.add(new CruxBullseyeApplication());
  apps.add(new CruxHardklorApplication());
  apps.add(new ExtractColumns());
  apps.add(new ExtractRows());
  apps.add(new GeneratePeptides());
  apps.add(new GetMs2Spectrum());
  apps.add(new PercolatorApplication());
  apps.add(new PredictPeptideIons());
  apps.add(new PrintProcessedSpectra());
  apps.add(new QRanker());
  apps.add(new ReadTideIndex());
  apps.add(new SearchForXLinks());
  apps.add(new SortColumn());
  apps.add(new SpectralCounts());
  apps.add(new StatColumn());
  apps.add(new TideIndexApplication());
  apps.add(new TideSearchApplication());

  string targetApp = Params::GetString("tool name");
  if (targetApp == "list") {
    for (vector<CruxApplication*>::const_iterator i = apps.begin(); i != apps.end(); i++) {
      cout << (*i)->getName() << endl;
    }
  } else if (targetApp == "default-params") {
    Params::Write(&cout, true);
  } else {
    CruxApplication* app = apps.find(targetApp);
    if (app == NULL) {
      carp(CARP_FATAL, "Invalid application '%s'", targetApp.c_str());
    }
    readParams(&apps); // This is unnecessary, but checks for problems
    generateToolHtml(&cout, app);
  }

  return 0;
}

void CreateDocs::readParams(const CruxApplicationList* apps) {
  carp(CARP_INFO, "Running parameter validity checks...");
  for (map<string, Param*>::const_iterator i = Params::BeginAll();
       i != Params::EndAll();
       i++) {
    Param* param = i->second;
    string name = param->GetName();
    bool isArg = param->IsArgument();
    // Check for unused if passed an application list
    if (apps) {
      vector<CruxApplication*> appsUsing;
      int appArgCount = 0;
      int appOptionCount = 0;
      for (vector<CruxApplication*>::const_iterator j = apps->begin();
           j != apps->end();
           j++) {
        vector<string> appArgs = (*j)->getArgs();
        for (vector<string>::iterator k = appArgs.begin(); k != appArgs.end(); k++) {
          if (StringUtils::EndsWith(*k, "+")) {
            *k = k->substr(k->length() - 1);
          }
        }
        vector<string> appOptions = (*j)->getOptions();
        bool appArg =
          std::find(appArgs.begin(), appArgs.end(), name) != appArgs.end();
        bool appOption =
          std::find(appOptions.begin(), appOptions.end(), name) != appOptions.end();
        if (!appArg && !appOption) {
          continue;
        }
        string appName = (*j)->getName();
        appsUsing.push_back(*j);
        if (appArg) {
          ++appArgCount;
          if (!isArg) {
            carp(CARP_WARNING, "'%s' is an option, but is listed as an argument for '%s'",
                 name.c_str(), appName.c_str());
          }
        }
        if (appOption) {
          ++appOptionCount;
          if (isArg) {
            carp(CARP_WARNING, "'%s' is an argument, but is listed as an option for '%s'",
                 name.c_str(), appName.c_str());
          }
          if (!Params::IsVisible(name)) {
            carp(CARP_WARNING, "'%s' is marked as hidden, but is listed as an option for '%s'",
                 name.c_str(), appName.c_str());
          }
        }
      }
      if (appArgCount > 0 && appOptionCount > 0) {
        carp(CARP_WARNING, "'%s' is both an option and an argument", name.c_str());
      }
      if (!appsUsing.empty()) {
        stringstream ss;
        for (vector<CruxApplication*>::const_iterator j = appsUsing.begin();
             j != appsUsing.end();
             j++) {
          if (j > appsUsing.begin()) {
            ss << ", ";
          }
          ss << (*j)->getName();
        }
        carp(CARP_DEBUG, "'%s' is used by: %s", name.c_str(), ss.str().c_str());
      } else {
        carp(CARP_WARNING, "No applications are using '%s'", name.c_str());
      }
    }
  }
}

void CreateDocs::generateToolHtml(
  ostream* outStream,
  const CruxApplication* application
) {
  carp(CARP_INFO, "Generating documentation for '%s'", application->getName().c_str());
  string doc = getToolTemplate();
  string inputTemplate = getToolInputTemplate();
  string outputTemplate = getToolOutputTemplate();
  string categoryTemplate = getToolOptionCategoryTemplate();
  string optionTemplate = getToolOptionTemplate();

  string appName = application->getName();
  string appDescription = application->getDescription();
  vector<string> args = application->getArgs();
  map<string, string> out = application->getOutputs();
  vector<string> opts = application->getOptions();

  // Build usage and input strings
  string usage = "crux " + appName + " [options]";
  string inputs;
  for (vector<string>::const_iterator i = args.begin(); i != args.end(); i++) {
    bool multiArg = StringUtils::EndsWith(*i, "+");
    string argName = multiArg ? i->substr(0, i->length() - 1) : *i;
    usage += " &lt;" + argName + "&gt;";
    if (multiArg) {
      usage += "+";
    }

    if (!Params::Exists(argName)) {
      carp(CARP_FATAL, "Invalid argument '%s' for application '%s'",
           argName.c_str(), appName.c_str());
    }
    string single = inputTemplate;
    map<string, string> replaceMap;
    replaceMap["#INPUTNAME#"] = argName;
    replaceMap["#INPUTDESCRIPTION#"] =
      Params::ProcessHtmlDocTags(htmlEscape(Params::GetUsage(argName)), true);
    makeReplacements(&single, replaceMap);
    inputs += single;
  }
  // Build outputs string
  string outputs;
  for (map<string, string>::const_iterator i = out.begin(); i != out.end(); i++) {
    string single = outputTemplate;
    map<string, string> replaceMap;
    replaceMap["#OUTPUTNAME#"] = i->first;
    replaceMap["#OUTPUTDESCRIPTION#"] = Params::ProcessHtmlDocTags(i->second, true);
    makeReplacements(&single, replaceMap);
    outputs += single;
  }
  // Build options string
  string options;
  for (int i = opts.size() - 1; i >= 0; i--) {
    if (!Params::IsVisible(opts[i])) {
      opts.erase(opts.begin() + i);
    }
  }
  vector< pair< string, vector<string> > > categories = Params::GroupByCategory(opts); 
  for (vector< pair< string, vector<string> > >::const_iterator i = categories.begin();
       i != categories.end();
       i++) {
    string categoryName = i->first;
    if (categoryName.empty()) {
      // Give category to all uncategorized options
      categoryName = appName + " options";
    }
    string optionSubset;
    const vector<string>& items = i->second;
    for (vector<string>::const_iterator j = items.begin(); j != items.end(); j++) {
      if (!Params::Exists(*j)) {
        carp(CARP_FATAL, "Invalid option '%s' for application '%s'",
             j->c_str(), appName.c_str());
      }
      string single = optionTemplate;
      map<string, string> replaceMap;
      replaceMap["#OPTIONNAME#"] = *j;
      replaceMap["#OPTIONDESCRIPTION#"] =
        Params::ProcessHtmlDocTags(htmlEscape(Params::GetUsage(*j)), true);
      replaceMap["#OPTIONTYPE#"] = Params::GetType(*j);
      replaceMap["#OPTIONDEFAULT#"] = Params::GetStringDefault(*j);
      makeReplacements(&single, replaceMap);
      optionSubset += single;
    }

    string single = categoryTemplate;
    map<string, string> replaceMap;
    replaceMap["#CATEGORYNAME#"] = categoryName;
    replaceMap["#CATEGORYOPTIONS#"] = optionSubset;
    makeReplacements(&single, replaceMap);
    options += single;
  }

  map<string, string> replacements;
  replacements["#TOOLNAME#"] = appName;
  replacements["#TOOLDESCRIPTION#"] = Params::ProcessHtmlDocTags(appDescription, true);
  replacements["#TOOLUSAGE#"] = usage;
  replacements["#TOOLINPUTS#"] = inputs;
  replacements["#TOOLOUTPUTS#"] = outputs;
  replacements["#TOOLOPTIONS#"] = options;
  makeReplacements(&doc, replacements);
  *outStream << doc;
}

void CreateDocs::makeReplacements(
  string* templateStr,
  const map<string, string>& replacements
) {
  const string OPEN_TAG = "<!--";
  const string CLOSE_TAG = "-->";

  unsigned idx, end_idx = 0;
  while ((idx = templateStr->find(OPEN_TAG, end_idx)) != string::npos) {
    if (idx > templateStr->length()) {
      break;
    }
    idx += OPEN_TAG.length();
    end_idx = templateStr->find(CLOSE_TAG, idx);
    if (end_idx == string::npos) {
      break;
    }
    string comment = templateStr->substr(idx, end_idx - idx);
    comment.erase(comment.begin(),
                  find_if(comment.begin(), comment.end(),
                  not1(ptr_fun<int, int>(isspace))));
    comment.erase(find_if(comment.rbegin(), comment.rend(),
                  not1(ptr_fun<int, int>(isspace))).base(), comment.end());
    map<string, string>::const_iterator iter = replacements.find(comment);
    if (iter == replacements.end()) {
      continue;
    }
    idx -= OPEN_TAG.length();
    end_idx += CLOSE_TAG.length();
    templateStr->replace(idx, end_idx - idx, iter->second);
    end_idx = idx + iter->second.length();
  }
}

string CreateDocs::htmlEscape(string s) {
  for (size_t i = 0; i < s.length(); i++) {
    for (map<char, string>::const_iterator j = htmlEntities_.begin();
         j != htmlEntities_.end();
         j++) {
      if (s[i] == j->first) {
        s.replace(i, 1, j->second);
        i += j->second.length() - 1;
        break;
      }
    }
  }
  return s;
}

string CreateDocs::getName() const {
  return "create-docs";
}

string CreateDocs::getDescription() const {
  return "[[html:<p>]]Creates documentation.[[html:</p>]]";
}

vector<string> CreateDocs::getArgs() const {
  string arr[] = {
    "tool name"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector<string> CreateDocs::getOptions() const {
  string arr[] = {
    "doc-template",
    "doc-input-template",
    "doc-output-template",
    "doc-option-category-template",
    "doc-option-template"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

map<string, string> CreateDocs::getOutputs() const {
  map<string, string> outputs;
  outputs["stdout"] =
    "HTML documentation for the specified tool.";
  return outputs;
}

bool CreateDocs::needsOutputDirectory() const {
  return false;
}

bool CreateDocs::hidden() const {
  return true;
}

string CreateDocs::getToolTemplate() {
  string path = Params::GetString("doc-template");
  if (!path.empty()) {
    return FileUtils::Read(path);
  }
  return
  "<!DOCTYPE HTML>\n"
  "<html>\n"
  "<head>\n"
  "<meta charset=\"UTF-8\">\n"
  "<title>crux <!-- #TOOLNAME# --></title>\n"
  "<script type=\"text/javascript\"\n"
  "  src=\"http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML\">\n"
  "</script>\n"
  "<script type=\"text/javascript\">\n"
  "  MathJax.Hub.Config({jax: [\"input/TeX\",\"output/HTML-CSS\"], displayAlign: \"left\"});\n"
  "</script>\n"
  "</head>\n"
  "<body>\n"
  "<h1><!-- #TOOLNAME# --></h1>\n"
  "<h2>Usage:</h2>\n"
  "<p><code><!-- #TOOLUSAGE# --></code></p>\n"
  "<h2>Description:</h2>\n"
  "<!-- #TOOLDESCRIPTION# -->\n"
  "<h2>Input:</h2>\n"
  "<ul>\n"
  "<!-- #TOOLINPUTS# --></ul>\n"
  "<h2>Output:</h2>\n"
  "<ul>\n"
  "<!-- #TOOLOUTPUTS# --></ul>\n"
  "<h2>Options:</h2>\n"
  "<ul>\n"
  "<!-- #TOOLOPTIONS# -->\n"
  "</ul>\n"
  "<hr>\n"
  "<a href=\"/\">Home</a>\n"
  "</body>\n"
  "</html>\n";
}

string CreateDocs::getToolInputTemplate() {
  string path = Params::GetString("doc-input-template");
  if (!path.empty()) {
    return FileUtils::Read(path);
  }
  return "  <li>&lt;<!-- #INPUTNAME# -->&gt; - <!-- #INPUTDESCRIPTION# --></li>\n";
}

string CreateDocs::getToolOutputTemplate() {
  string path = Params::GetString("doc-output-template");
  if (!path.empty()) {
    return FileUtils::Read(path);
  }
  return "  <li><!-- #OUTPUTNAME# --> - <!-- #OUTPUTDESCRIPTION# --></li>\n";
}

string CreateDocs::getToolOptionCategoryTemplate() {
  string path = Params::GetString("doc-option-category-template");
  if (!path.empty()) {
    return FileUtils::Read(path);
  }
  return
  "<li>\n"
  "<h3><!-- #CATEGORYNAME# --></h3>\n"
  "<ul>\n"
  "<!-- #CATEGORYOPTIONS# --></ul>\n"
  "</li>\n";
}

string CreateDocs::getToolOptionTemplate() {
  string path = Params::GetString("doc-option-template");
  if (!path.empty()) {
    return FileUtils::Read(path);
  }
  return
  "  <li><!-- #OPTIONNAME# --> &lt;<!-- #OPTIONTYPE# -->&gt; - "
  "<!-- #OPTIONDESCRIPTION# --></li>\n";
}

