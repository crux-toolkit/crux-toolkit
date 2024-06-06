#ifndef _CNPXINPUTFILE_H
#define _CNPXINPUTFILE_H

#include <string>
#include <vector>

class CnpxInputFile {
public:

  void write(FILE* f);

  std::string directory;
  std::string name;

private:

};

#endif