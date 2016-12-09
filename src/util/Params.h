#ifndef PARAMS_H
#define PARAMS_H

#include <ostream>
#include <limits>
#include <map>
#include <set>
#include <vector>

#ifdef _MSC_VER
// Turn off Microsoft min and max macros
#undef max
#undef min
#endif

class Param;

class Params {
 public:
  Params();
  ~Params();

  // Get the value of the parameter
  static bool GetBool(const std::string& name);
  static int GetInt(const std::string& name);
  static double GetDouble(const std::string& name);
  static std::string GetString(const std::string& name);

  // Get the default value of the parameter
  static bool GetBoolDefault(const std::string& name);
  static int GetIntDefault(const std::string& name);
  static double GetDoubleDefault(const std::string& name);
  static std::string GetStringDefault(const std::string& name);

  // Get the values of an argument that accepts multiple values
  static const std::vector<std::string>& GetStrings(const std::string& name);

  // Get the usage statement for the parameter
  static std::string GetUsage(const std::string& name);

  // Get whether the parameter has been modified from its original value or not
  static std::string GetFileNotes(const std::string& name);

  // Get whether the parameter is visible to the user or not
  static bool IsVisible(const std::string& name);

  // Get whether the parameter is an argument or an option
  static bool IsArgument(const std::string& name);

  // Get the accepted values for the parameter
  static std::string GetAcceptedValues(const std::string& name);

  // Get whether the parameter has been modified from its original value or not
  static bool IsDefault(const std::string& name);

  // Get whether the parameter exists or not
  static bool Exists(const std::string& name);

  // Set the value of a parameter after it has been initialized
  // Throws exception if parameter doesn't exist or parameters have been finalized
  static void Set(const std::string& name, bool value);
  static void Set(const std::string& name, int value);
  static void Set(const std::string& name, double value);
  static void Set(const std::string& name, const char* value);
  static void Set(const std::string& name, const std::string& value);

  // Add the value of an argument
  // Throws exception if parameter has an invalid value, already exists as a non-argument,
  // or parameters have been finalized
  static void AddArgValue(const std::string& name, const std::string& value);

  // Lock parameters and prevent them from being modified
  static void Finalize();

  // Write all contents of the ordered parameter list to file
  static void Write(std::ostream* out, bool defaults = false);

  // Get iterators for the beginning and end of the entire parameter list
  static std::map<std::string, Param*>::const_iterator BeginAll();
  static std::map<std::string, Param*>::const_iterator EndAll();

  // Get iterators for the beginning and end of the ordered parameter list
  static std::vector<const Param*>::const_iterator Begin();
  static std::vector<const Param*>::const_iterator End();

  // Process tags in a string
  static std::string ProcessHtmlDocTags(std::string s, bool html = false);

  // Group the given options by category
  static std::vector< std::pair< std::string, std::vector<std::string> > > GroupByCategory(
    const std::vector<std::string>& options);

 private:
  struct ParamCategory {
    explicit ParamCategory(std::string name): Name(name), Items(std::set<const Param*>()) {}
    std::string Name;
    std::set<const Param*> Items;
  };

  // Group parameters by category
  void Categorize();

  // Initialize optional parameters
  void InitBoolParam(
    const std::string& name,
    bool value,
    const std::string& usage,
    const std::string& fileNotes,
    bool visible
  );
  void InitIntParam(
    const std::string& name,
    int value,
    int min,
    int max,
    const std::string& usage,
    const std::string& fileNotes,
    bool visible
  );
  void InitIntParam(
    const std::string& name,
    int value,
    const std::string& usage,
    const std::string& fileNotes,
    bool visible
  );
  void InitDoubleParam(
    const std::string& name,
    double value,
    double min,
    double max,
    const std::string& usage,
    const std::string& fileNotes,
    bool visible
  );
  void InitDoubleParam(
    const std::string& name,
    double value,
    const std::string& usage,
    const std::string& fileNotes,
    bool visible
  );
  void InitStringParam(
    const std::string& name,
    const std::string& value,
    const std::string& validValues,
    const std::string& usage,
    const std::string& fileNotes,
    bool visible
  );
  void InitStringParam(
    const std::string& name,
    const std::string& value,
    const std::string& usage,
    const std::string& fileNotes,
    bool visible
  );
  void InitArgParam(
    const std::string& name,
    const std::string& usage
  );

  // Get a parameter by name
  // Throws exception if the parameter doesn't exist
  static Param* Require(const std::string& name);

  // Get a parameter by name; returns NULL if the parameter doesn't exist
  Param* Get(const std::string& name);

  // Add a parameter to the container
  // Throws exception if parameter has an invalid value, already exists,
  // or parameters have been finalized
  void Add(Param* param);

  // Add parameter category
  void AddCategory(const std::string& name, const std::set<std::string>& params);

  // Lock parameters and prevent them from being modified
  void FinalizeParams();

  // Throw exception if the parameters have been finalized
  void CanModifyCheck() const;

  std::map<std::string, Param*> params_;
  std::vector<const Param*> paramsOrdered_;
  std::vector<ParamCategory> categories_;
  bool finalized_;
};

class Param {
 public:
  Param(const std::string& name,
        const std::string& usage,
        const std::string& fileNotes,
        bool visible);
  virtual ~Param();

  // Get the name of the parameter
  std::string GetName() const;

  // Get the usage statement for the parameter
  std::string GetUsage() const;

  // Get additional notes about the parameter for the parameter file
  std::string GetFileNotes() const;

  // Get whether the parameter is visible to the user or not
  bool IsVisible() const;

  // Get whether the parameter is an argument or an option
  virtual bool IsArgument() const;

  // Get whether the parameter's value is valid or not
  virtual void ThrowIfInvalid() const;

  // Get the type of the parameter
  virtual std::string GetAcceptedValues() const = 0;

  // Get whether the parameter has been modified from its original value or not
  virtual bool IsDefault() const = 0;

  // Get the value of the parameter
  virtual bool GetBool() const = 0;
  virtual int GetInt() const = 0;
  virtual double GetDouble() const = 0;
  virtual std::string GetString() const = 0;

  // Get the default value of the parameter
  virtual bool GetBoolDefault() const = 0;
  virtual int GetIntDefault() const = 0;
  virtual double GetDoubleDefault() const = 0;
  virtual std::string GetStringDefault() const = 0;

  // Set the value of the parameter
  virtual void Set(bool value);
  virtual void Set(int value);
  virtual void Set(double value);
  virtual void Set(const char* value);
  virtual void Set(const std::string& value);

  // Get the parameter as a string to write to a parameter file
  virtual std::string GetParamFileString(bool defaultValue = false) const;
 protected:
  std::string name_, usage_, fileNotes_;
  bool visible_;
};

class BoolParam : public Param {
 public:
  BoolParam(const std::string& name,
            const std::string& usage,
            const std::string& fileNotes,
            bool visible,
            bool value);
  std::string GetAcceptedValues() const;
  bool IsDefault() const;
  bool GetBool() const;
  int GetInt() const;
  double GetDouble() const;
  std::string GetString() const;
  bool GetBoolDefault() const;
  int GetIntDefault() const;
  double GetDoubleDefault() const;
  std::string GetStringDefault() const;
  void Set(bool value);
  void Set(int value);
  void Set(double value);
  void Set(const std::string& value);

  static bool From(int i);
  static bool From(double d);
  static bool From(std::string s);
 protected:
  bool value_, original_;
};

class IntParam : public Param {
 public:
  IntParam(const std::string& name,
           const std::string& usage,
           const std::string& fileNotes,
           bool visible,
           int value,
           int min = -std::numeric_limits<int>::max(),
           int max = std::numeric_limits<int>::max());
  std::string GetAcceptedValues() const;
  void ThrowIfInvalid() const;
  bool IsDefault() const;
  bool GetBool() const;
  int GetInt() const;
  double GetDouble() const;
  std::string GetString() const;
  bool GetBoolDefault() const;
  int GetIntDefault() const;
  double GetDoubleDefault() const;
  std::string GetStringDefault() const;
  void Set(bool value);
  void Set(int value);
  void Set(double value);
  void Set(const std::string& value);

  static int From(bool b);
  static int From(double d);
  static int From(const std::string& s);
 protected:
  int value_, min_, max_, original_;
};

class DoubleParam : public Param {
 public:
  DoubleParam(const std::string& name,
              const std::string& usage,
              const std::string& fileNotes,
              bool visible,
              double value,
              double min = -std::numeric_limits<double>::max(),
              double max = std::numeric_limits<double>::max());
  std::string GetAcceptedValues() const;
  void ThrowIfInvalid() const;
  bool IsDefault() const;
  bool GetBool() const;
  int GetInt() const;
  double GetDouble() const;
  std::string GetString() const;
  bool GetBoolDefault() const;
  int GetIntDefault() const;
  double GetDoubleDefault() const;
  std::string GetStringDefault() const;
  void Set(bool value);
  void Set(int value);
  void Set(double value);
  void Set(const std::string& value);

  static double From(bool b);
  static double From(int i);
  static double From(const std::string& s);
 protected:
  double value_, min_, max_, original_;
};

class StringParam : public Param {
 public:
  StringParam(const std::string& name,
              const std::string& usage,
              const std::string& fileNotes,
              bool visible,
              const std::string& value,
              const std::vector<std::string>& validvalues = std::vector<std::string>());
  std::string GetAcceptedValues() const;
  void ThrowIfInvalid() const;
  bool IsDefault() const;
  bool GetBool() const;
  int GetInt() const;
  double GetDouble() const;
  std::string GetString() const;
  bool GetBoolDefault() const;
  int GetIntDefault() const;
  double GetDoubleDefault() const;
  std::string GetStringDefault() const;
  void Set(bool value);
  void Set(int value);
  void Set(double value);
  void Set(const std::string& value);

  static std::string From(bool b);
  static std::string From(int i);
  static std::string From(double d);
 protected:
  std::string value_, original_;
  std::vector<std::string> validValues_;
};

class ArgParam : public Param {
 public:
  ArgParam(const std::string& name, const std::string& usage);
  bool IsArgument() const;
  std::string GetAcceptedValues() const;
  bool IsDefault() const;
  bool GetBool() const;
  int GetInt() const;
  double GetDouble() const;
  std::string GetString() const;
  bool GetBoolDefault() const;
  int GetIntDefault() const;
  double GetDoubleDefault() const;
  std::string GetStringDefault() const;

  const std::vector<std::string>& GetStrings() const;
  void Set(bool value);
  void Set(int value);
  void Set(double value);
  void Set(const std::string& value);
  void AddValue(const std::string& value);
 protected:
  std::vector<std::string> values_;
};

#endif

