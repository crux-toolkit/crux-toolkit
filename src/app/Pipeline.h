#ifndef PIPELINE_H
#define PIPELINE_H

#include "CruxApplication.h"

class PipelineApplication : public CruxApplication {
 public:
  PipelineApplication();
  virtual ~PipelineApplication();

  virtual int main(int argc, char** argv);
  virtual std::string getName() const;
  virtual std::string getDescription() const;
  virtual std::vector<std::string> getArgs() const;
  virtual std::vector<std::string> getOptions() const;
  virtual std::vector< std::pair<std::string, std::string> > getOutputs() const;
  virtual std::string getFileStem() const;
  virtual COMMAND_T getCommand() const;
  virtual bool needsOutputDirectory() const;
  virtual bool hidden() const;
  virtual void processParams();

 private:
  std::vector<CruxApplication*> apps_;

  static void checkParams();
  static std::vector<std::string> getExpectedResultsFiles(
    CruxApplication* app,
    const std::vector<std::string>& spectra
  );
  int runBullseye(CruxApplication* app,
                  std::vector<std::string>* spectra);
  int runSearch(CruxApplication* app,
                const std::vector<std::string>& spectra,
                const std::string& database,
                std::vector<std::string>* resultsFiles);
  int runPostProcessor(CruxApplication* app,
                       const std::vector<std::string>& resultsFiles);
};

#endif
